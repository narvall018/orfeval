"""Usage des codons : comptages, RSCU et GC en troisième position.

Le **RSCU** (*Relative Synonymous Codon Usage*, Sharp et al. 1986) rapporte l'usage d'un
codon à l'usage moyen des codons synonymes : 1 = pas de préférence, > 1 = codon
préféré. Comparer le RSCU des gènes prédits et des gènes annotés vérifie que
l'ensemble prédit a la même « signature » que l'annotation.

Le **GC3** (GC en 3e position des codons, souvent synonyme) reflète surtout le biais
mutationnel du génome.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
from Bio.Data import CodonTable

from orfeval.orfs.genetic_code import CODONS, encode_codons


def codon_counts(sequences: Iterable[str]) -> np.ndarray:
    """Compte les codons internes (hors codon start et codon stop) de séquences codantes.

    Le premier codon est exclu car un start alternatif (GTG, TTG) code quand même une
    méthionine ; le dernier car c'est le stop.
    """
    counts = np.zeros(64, dtype=np.int64)
    for sequence in sequences:
        inner = encode_codons(sequence)[1:-1]
        counts += np.bincount(inner[inner < 64], minlength=64)
    return counts


def rscu_table(counts: np.ndarray, table_id: int) -> pd.DataFrame:
    """RSCU de chaque codon sens selon le code génétique ``table_id``.

    Returns
    -------
    pandas.DataFrame
        Colonnes ``codon``, ``amino_acid``, ``count``, ``rscu`` (NaN si l'acide aminé
        n'est jamais observé).
    """
    forward = CodonTable.unambiguous_dna_by_id[table_id].forward_table
    families: dict[str, list[str]] = defaultdict(list)
    for codon, amino_acid in forward.items():
        families[amino_acid].append(codon)
    index = {codon: position for position, codon in enumerate(CODONS)}
    rows = []
    for amino_acid, codons in sorted(families.items()):
        family_counts = np.array([counts[index[codon]] for codon in codons], dtype=float)
        mean = family_counts.mean()
        for codon, count in zip(codons, family_counts, strict=True):
            rows.append(
                {
                    "codon": codon,
                    "amino_acid": amino_acid,
                    "count": int(count),
                    "rscu": count / mean if mean > 0 else np.nan,
                }
            )
    return pd.DataFrame(rows)


def gc3_values(sequences: Sequence[str]) -> np.ndarray:
    """GC en 3e position des codons internes, pour chaque séquence."""
    values = []
    for sequence in sequences:
        third = sequence[3:-3][2::3]
        values.append((third.count("G") + third.count("C")) / len(third) if third else np.nan)
    return np.array(values, dtype=float)


def compare_rscu(predicted: Sequence[str], reference: Sequence[str], table_id: int) -> pd.DataFrame:
    """RSCU des gènes prédits et des gènes de référence, côte à côte."""
    predicted_table = rscu_table(codon_counts(predicted), table_id)
    reference_table = rscu_table(codon_counts(reference), table_id)
    merged = predicted_table.merge(
        reference_table, on=["codon", "amino_acid"], suffixes=("_predicted", "_reference")
    )
    return merged
