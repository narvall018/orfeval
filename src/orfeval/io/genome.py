"""Lecture et validation des génomes (GenBank ou FASTA, compressés ou non).

Le fichier GenBank fournit à la fois la séquence (entrée du prédicteur) et
l'annotation des CDS (référence pour l'évaluation). Les deux sont séparées ici :
le prédicteur ne reçoit que l'objet :class:`~orfeval.models.Genome`.
"""

from __future__ import annotations

import gzip
import hashlib
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Literal

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    SeqFeature,
    SimpleLocation,
)
from Bio.SeqRecord import SeqRecord

from orfeval.exceptions import InputFormatError
from orfeval.logging_utils import get_logger
from orfeval.models import AnnotatedFeature, Genome, Interval, ReferenceGene, Topology

logger = get_logger(__name__)

FileFormat = Literal["genbank", "fasta"]
RNA_FEATURE_TYPES = frozenset({"tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA"})
GENBANK_SUFFIXES = frozenset({".gb", ".gbk", ".gbff", ".genbank"})
FASTA_SUFFIXES = frozenset({".fa", ".fasta", ".fna", ".fas", ".fsa"})
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
MIN_GENOME_LENGTH = 1_000
SHORT_GENOME_WARNING = 20_000
MAX_AMBIGUOUS_FRACTION = 0.05


@dataclass(frozen=True, slots=True)
class AnnotationSummary:
    """Bilan de l'extraction des CDS de référence.

    Chaque CDS exclue l'est pour une raison explicite, afin que l'évaluation
    reste interprétable (voir ``docs/methodology.md``).
    """

    n_cds: int
    n_pseudo: int
    n_partial: int
    n_unsupported: int
    n_compound: int
    n_reference: int
    declared_tables: dict[int, int] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class LoadedGenome:
    """Génome chargé, avec son annotation éventuelle et sa provenance."""

    genome: Genome
    references: list[ReferenceGene]
    other_features: list[AnnotatedFeature]
    annotation: AnnotationSummary | None
    source: Path
    file_format: FileFormat
    sha256: str
    n_ambiguous: int

    @property
    def has_annotation(self) -> bool:
        """Vrai si des CDS de référence sont disponibles pour l'évaluation."""
        return bool(self.references)


def sha256sum(path: Path) -> str:
    """Somme de contrôle SHA-256 d'un fichier."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text(path: Path) -> IO[str]:
    """Ouvre un fichier texte, compressé en gzip ou non (détection par signature)."""
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def detect_format(path: Path) -> FileFormat:
    """Détermine le format à partir de l'extension, puis du contenu si besoin."""
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]
    if suffixes and suffixes[-1] in GENBANK_SUFFIXES:
        return "genbank"
    if suffixes and suffixes[-1] in FASTA_SUFFIXES:
        return "fasta"
    try:
        with open_text(path) as handle:
            for line in handle:
                if line.strip():
                    if line.startswith("LOCUS"):
                        return "genbank"
                    if line.startswith(">"):
                        return "fasta"
                    break
    except (OSError, UnicodeDecodeError) as exc:
        raise InputFormatError(f"lecture impossible de {path} : {exc}") from exc
    raise InputFormatError(f"format non reconnu pour {path} (GenBank ou FASTA attendu)")


def _read_single_record(path: Path, file_format: FileFormat) -> SeqRecord:
    try:
        with open_text(path) as handle:
            records: list[SeqRecord] = list(SeqIO.parse(handle, file_format))
    except (ValueError, OSError, UnicodeDecodeError) as exc:
        raise InputFormatError(f"{path} n'est pas un fichier {file_format} valide : {exc}") from exc
    if not records:
        raise InputFormatError(f"aucune séquence trouvée dans {path}")
    if len(records) > 1:
        ids = ", ".join(record.id or "?" for record in records[:5])
        raise InputFormatError(
            f"{path} contient {len(records)} séquences ({ids}...). "
            "Cette version analyse une séquence par fichier : séparez les réplicons."
        )
    return records[0]


def _clean_sequence(record: SeqRecord, path: Path) -> tuple[str, int]:
    try:
        sequence = str(record.seq).upper()
    except UndefinedSequenceError as exc:  # enregistrement NCBI de type CON
        raise InputFormatError(
            f"{path} ne contient pas la séquence nucléotidique. Pour un enregistrement NCBI, "
            "téléchargez-le au format 'gbwithparts' (orfeval fetch le fait)."
        ) from exc
    invalid = set(sequence) - IUPAC_DNA
    if invalid:
        shown = "".join(sorted(invalid))[:10]
        raise InputFormatError(f"caractères non ADN dans {path} : {shown!r}")
    length = len(sequence)
    if length < MIN_GENOME_LENGTH:
        raise InputFormatError(f"séquence trop courte ({length} nt < {MIN_GENOME_LENGTH} nt)")
    n_ambiguous = length - sum(sequence.count(base) for base in "ACGT")
    if n_ambiguous / length > MAX_AMBIGUOUS_FRACTION:
        raise InputFormatError(
            f"{n_ambiguous / length:.1%} de bases ambiguës : "
            "séquence trop incomplète pour l'analyse"
        )
    if length < SHORT_GENOME_WARNING:
        logger.warning(
            "Génome court (%d nt) : l'auto-apprentissage dispose de peu de gènes.", length
        )
    return sequence, n_ambiguous


def _is_partial(location: SimpleLocation | CompoundLocation) -> bool:
    return any(
        isinstance(part.start, (BeforePosition, AfterPosition))
        or isinstance(part.end, (BeforePosition, AfterPosition))
        for part in location.parts
    )


def _feature_interval(
    location: SimpleLocation | CompoundLocation, genome_length: int, circular: bool
) -> tuple[Interval, bool] | None:
    """Convertit une localisation Biopython en :class:`Interval`.

    Renvoie ``(intervalle, composé)`` ou ``None`` si la localisation n'est pas
    prise en charge (brins mélangés, par exemple).
    """
    strand = location.strand
    if strand not in (1, -1):
        return None
    if not isinstance(location, CompoundLocation):
        return Interval(int(location.start), int(location.end), strand), False

    parts = location.parts
    if circular and len(parts) == 2:
        ends_at_origin = [part for part in parts if int(part.end) == genome_length]
        starts_at_zero = [part for part in parts if int(part.start) == 0]
        if len(ends_at_origin) == 1 and len(starts_at_zero) == 1:
            left = int(ends_at_origin[0].start)
            right = genome_length + int(starts_at_zero[0].end)
            return Interval(left, right, strand), False
    # Autres CDS composées (décalage de cadre programmé...) : l'extrémité 3'
    # reste définie par la partie la plus en aval, ce qui suffit à l'évaluation.
    left = min(int(part.start) for part in parts)
    right = max(int(part.end) for part in parts)
    return Interval(left, right, strand), True


def _label(feature: SeqFeature, default: str) -> str:
    qualifiers = feature.qualifiers
    for key in ("locus_tag", "gene", "product"):
        if key in qualifiers:
            return str(qualifiers[key][0])
    return default


def extract_reference_genes(
    record: SeqRecord, genome_length: int, *, circular: bool
) -> tuple[list[ReferenceGene], list[AnnotatedFeature], AnnotationSummary]:
    """Extrait les CDS de référence utilisables pour l'évaluation.

    Sont exclues : les pseudogènes (``/pseudo``, ``/pseudogene``), les CDS partielles
    (positions floues ``<`` ou ``>``) et les localisations non prises en charge.
    Les pseudogènes, CDS partielles et gènes d'ARN sont renvoyés à part, comme
    éléments de contexte pour interpréter les prédictions non appariées.
    """
    references: list[ReferenceGene] = []
    others: list[AnnotatedFeature] = []
    tables: Counter[int] = Counter()
    n_cds = n_pseudo = n_partial = n_unsupported = n_compound = 0

    for index, feature in enumerate(record.features):
        location = feature.location
        if feature.type in RNA_FEATURE_TYPES and location is not None:
            converted_rna = _feature_interval(location, genome_length, circular)
            if converted_rna is not None:
                others.append(
                    AnnotatedFeature(converted_rna[0], feature.type, _label(feature, "ARN"))
                )
            continue
        if feature.type != "CDS":
            continue
        n_cds += 1
        qualifiers = feature.qualifiers
        for value in qualifiers.get("transl_table", []):
            try:
                tables[int(value)] += 1
            except ValueError:
                logger.warning("transl_table illisible : %r", value)
        if location is None:
            n_unsupported += 1
            continue
        converted = _feature_interval(location, genome_length, circular)
        if "pseudo" in qualifiers or "pseudogene" in qualifiers:
            n_pseudo += 1
            if converted is not None:
                others.append(AnnotatedFeature(converted[0], "pseudogène", _label(feature, "")))
            continue
        if _is_partial(location):
            n_partial += 1
            if converted is not None:
                others.append(AnnotatedFeature(converted[0], "CDS partielle", _label(feature, "")))
            continue
        if converted is None:
            n_unsupported += 1
            continue
        interval, compound = converted
        n_compound += int(compound)
        locus_tag = qualifiers.get("locus_tag", [f"CDS_{index}"])[0]
        references.append(
            ReferenceGene(
                interval=interval,
                locus_tag=locus_tag,
                product=qualifiers.get("product", [""])[0],
                gene_name=qualifiers.get("gene", [None])[0],
                compound=compound,
            )
        )

    summary = AnnotationSummary(
        n_cds=n_cds,
        n_pseudo=n_pseudo,
        n_partial=n_partial,
        n_unsupported=n_unsupported,
        n_compound=n_compound,
        n_reference=len(references),
        declared_tables=dict(tables),
    )
    return references, others, summary


def load_genome(path: Path, *, topology: Topology | None = None) -> LoadedGenome:
    """Charge un génome GenBank ou FASTA (éventuellement compressé en gzip).

    Parameters
    ----------
    path
        Fichier d'entrée.
    topology
        Force la topologie. Par défaut : celle du fichier GenBank, sinon linéaire.

    Raises
    ------
    InputFormatError
        Fichier absent, mal formé, multi-séquences ou séquence inexploitable.
    """
    path = Path(path)
    if not path.is_file():
        raise InputFormatError(f"fichier introuvable : {path}")
    file_format = detect_format(path)
    record = _read_single_record(path, file_format)
    sequence, n_ambiguous = _clean_sequence(record, path)

    declared_topology = record.annotations.get("topology")
    resolved: Topology
    if topology is not None:
        resolved = topology
    elif declared_topology == "circular":
        resolved = "circular"
    else:
        resolved = "linear"
        if declared_topology != "linear":
            logger.info("Topologie non renseignée pour %s : génome supposé linéaire.", record.id)
    circular = resolved == "circular"

    references: list[ReferenceGene] = []
    other_features: list[AnnotatedFeature] = []
    annotation: AnnotationSummary | None = None
    declared_table: int | None = None
    if file_format == "genbank":
        references, other_features, annotation = extract_reference_genes(
            record, len(sequence), circular=circular
        )
        if annotation.declared_tables:
            declared_table = max(annotation.declared_tables.items(), key=lambda kv: kv[1])[0]
            if len(annotation.declared_tables) > 1:
                logger.warning(
                    "Plusieurs codes génétiques déclarés %s : %d retenu (majoritaire).",
                    annotation.declared_tables,
                    declared_table,
                )
        if annotation.n_cds == 0:
            annotation = None

    organism = record.annotations.get("organism")
    genome = Genome(
        seq_id=record.id or path.name,
        description=record.description or "",
        sequence=sequence,
        topology=resolved,
        organism=organism if isinstance(organism, str) else None,
        declared_table=declared_table,
    )
    logger.info(
        "Génome chargé : %s (%d nt, %s, %d CDS de référence)",
        genome.seq_id,
        genome.length,
        genome.topology,
        len(references),
    )
    return LoadedGenome(
        genome=genome,
        references=references,
        other_features=other_features,
        annotation=annotation,
        source=path,
        file_format=file_format,
        sha256=sha256sum(path),
        n_ambiguous=n_ambiguous,
    )
