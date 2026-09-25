"""Sérialisation des résultats : TSV, GFF3, FASTA protéique et JSON."""

from __future__ import annotations

import io
import json
import math
from collections.abc import Sequence
from dataclasses import asdict, is_dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from orfeval.models import STRAND_SYMBOL, Genome, PredictedGene

PREDICTION_COLUMNS = [
    "gene_id",
    "seqid",
    "start",
    "end",
    "strand",
    "length_nt",
    "start_codon",
    "stop_codon",
    "score_bits",
    "start_score_bits",
    "rbs_motif",
    "rbs_spacer",
    "n_alternative_starts",
]


def gene_ids(genes: Sequence[PredictedGene], prefix: str) -> list[str]:
    """Identifiants stables ``<prefix>_00001`` dans l'ordre des gènes fournis."""
    width = max(5, len(str(len(genes))))
    return [f"{prefix}_{index:0{width}d}" for index in range(1, len(genes) + 1)]


def predictions_dataframe(
    genes: Sequence[PredictedGene], genome: Genome, prefix: str
) -> pd.DataFrame:
    """Tableau des gènes prédits (coordonnées 1-based inclusives, comme GFF3)."""
    rows = [
        {
            "gene_id": gene_id,
            "seqid": genome.seq_id,
            "start": gene.interval.left + 1,
            "end": gene.interval.right,
            "strand": STRAND_SYMBOL[gene.interval.strand],
            "length_nt": gene.length,
            "start_codon": gene.start_codon,
            "stop_codon": gene.stop_codon,
            "score_bits": round(gene.score, 3),
            "start_score_bits": round(gene.start_score, 3),
            "rbs_motif": gene.rbs_motif or "",
            "rbs_spacer": gene.rbs_spacer if gene.rbs_spacer is not None else pd.NA,
            "n_alternative_starts": gene.n_alternative_starts,
        }
        for gene_id, gene in zip(gene_ids(genes, prefix), genes, strict=True)
    ]
    frame = pd.DataFrame(rows, columns=PREDICTION_COLUMNS)
    return frame.astype({"rbs_spacer": "Int64"})


def gff3_text(genes: Sequence[PredictedGene], genome: Genome, prefix: str) -> str:
    """Annotation GFF3 des gènes prédits.

    Pour un génome circulaire, un gène qui passe l'origine a une fin supérieure à la
    longueur de la séquence, conformément à la spécification GFF3 (``Is_circular``).
    """
    lines = ["##gff-version 3", f"##sequence-region {genome.seq_id} 1 {genome.length}"]
    region_attributes = "ID=region_1"
    if genome.is_circular:
        region_attributes += ";Is_circular=true"
    lines.append(
        f"{genome.seq_id}\torfeval\tregion\t1\t{genome.length}\t.\t+\t.\t{region_attributes}"
    )
    for gene_id, gene in zip(gene_ids(genes, prefix), genes, strict=True):
        attributes = [
            f"ID={gene_id}",
            f"start_codon={gene.start_codon}",
            f"start_score={gene.start_score:.2f}",
        ]
        if gene.rbs_motif:
            attributes.append(f"rbs_motif={gene.rbs_motif};rbs_spacer={gene.rbs_spacer}")
        lines.append(
            "\t".join(
                [
                    genome.seq_id,
                    "orfeval",
                    "CDS",
                    str(gene.interval.left + 1),
                    str(gene.interval.right),
                    f"{gene.score:.2f}",
                    STRAND_SYMBOL[gene.interval.strand],
                    "0",
                    ";".join(attributes),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def proteins_fasta_text(
    genes: Sequence[PredictedGene], genome: Genome, prefix: str, table_id: int
) -> str:
    """Protéines traduites avec le code génétique utilisé pour la prédiction.

    La traduction utilise ``cds=True`` : le premier codon est traduit en méthionine
    (y compris GTG/TTG) et Biopython vérifie l'absence de codon stop interne.
    """
    records = []
    for gene_id, gene in zip(gene_ids(genes, prefix), genes, strict=True):
        protein = Seq(gene.sequence).translate(table=table_id, cds=True)
        interval = gene.interval
        description = (
            f"{genome.seq_id}:{interval.left + 1}-{interval.right}"
            f"({STRAND_SYMBOL[interval.strand]}) table={table_id} score={gene.score:.1f}"
        )
        records.append(SeqRecord(protein, id=gene_id, description=description))
    buffer = io.StringIO()
    SeqIO.write(records, buffer, "fasta")
    return buffer.getvalue()


def _to_jsonable(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _to_jsonable(asdict(value))
    if isinstance(value, dict):
        return {str(key): _to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_to_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return None if math.isnan(number) or math.isinf(number) else number
    if isinstance(value, np.ndarray):
        return _to_jsonable(value.tolist())
    if isinstance(value, pd.DataFrame):
        return _to_jsonable(value.to_dict(orient="records"))
    if value is pd.NA:
        return None
    return value


def write_json(data: Any, path: Path) -> None:
    """Écrit un objet en JSON lisible (NaN → null, types numpy convertis)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_to_jsonable(data), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def write_text(text: str, path: Path) -> None:
    """Écrit un texte UTF-8 en créant le dossier parent si nécessaire."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Écrit un tableau TSV sans index."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
