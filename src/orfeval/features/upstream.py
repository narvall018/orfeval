"""Profil nucléotidique en amont des codons start (``Bio.motifs``).

Les séquences situées juste en amont des starts sont alignées sur le codon start ; on
mesure, position par position, l'écart de composition à la composition du génome
(entropie relative, en bits). Un pic autour de -12 à -6 riche en A/G signale un motif
de Shine-Dalgarno ; son absence suggère que ce signal est peu utilisé par le génome,
ou que les starts retenus sont souvent faux.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio import motifs
from Bio.Seq import Seq

from orfeval.orfs.null_model import strand_symmetric_composition


@dataclass(frozen=True, slots=True)
class UpstreamProfile:
    """Fréquences par position et contenu en information (bits)."""

    frequencies: pd.DataFrame
    information: pd.Series
    n_sequences: int
    consensus: str


def upstream_profile(sequences: Sequence[str], genome_sequence: str) -> UpstreamProfile | None:
    """Construit le profil de séquences amont de même longueur.

    Les séquences contenant des bases ambiguës sont ignorées. Renvoie ``None`` s'il
    reste moins de 10 séquences.

    Parameters
    ----------
    sequences
        Séquences amont (5'→3'), toutes de même longueur, adjacentes au codon start.
    genome_sequence
        Génome, pour la composition de fond.
    """
    clean = [seq for seq in sequences if seq and set(seq) <= set("ACGT")]
    if len(clean) < 10:
        return None
    length = len(clean[0])
    if any(len(seq) != length for seq in clean):
        raise ValueError("les séquences amont doivent avoir la même longueur")
    motif = motifs.create([Seq(seq) for seq in clean], alphabet="ACGT")
    counts = pd.DataFrame({base: np.asarray(motif.counts[base], dtype=float) for base in "ACGT"})
    frequencies = counts.div(counts.sum(axis=1), axis=0)
    frequencies.index = pd.Index(np.arange(-length, 0), name="position")

    background = strand_symmetric_composition(genome_sequence)
    information = pd.Series(0.0, index=frequencies.index, name="bits")
    for base in "ACGT":
        p = frequencies[base].to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            term = np.where(p > 0, p * np.log2(p / background[base]), 0.0)
        information += term
    return UpstreamProfile(
        frequencies=frequencies,
        information=information,
        n_sequences=len(clean),
        consensus=str(motif.consensus),
    )
