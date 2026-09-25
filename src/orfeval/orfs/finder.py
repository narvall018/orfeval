"""Recherche des cadres ouverts de lecture (ORF) sur les six cadres.

Un **ORF** est défini ici comme la région comprise entre deux codons stop consécutifs
d'un même cadre, tronquée au codon start le plus en amont. Tous les codons start
internes sont conservés comme débuts alternatifs : le choix du « bon » start est
fait plus tard (:mod:`orfeval.orfs.selection`).

Génomes circulaires
-------------------
La séquence de chaque brin est concaténée à elle-même. Seuls les codons stop situés
dans la seconde copie sont retenus : chaque stop est ainsi traité une seule fois,
avec tout son contexte amont, y compris lorsque l'ORF traverse l'origine.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from Bio.Seq import reverse_complement

from orfeval.models import Interval
from orfeval.orfs.genetic_code import GeneticCode, encode_codons


@dataclass(frozen=True, slots=True)
class Orf:
    """ORF maximal et ses codons start alternatifs.

    Attributes
    ----------
    strand
        Brin (+1 ou -1).
    left, right
        Coordonnées (brin direct) de l'ORF maximal, du start le plus en amont au stop inclus.
    genome_length
        Longueur du génome (nécessaire pour les coordonnées circulaires).
    sequence
        Séquence 5'→3' sur le brin de l'ORF, du premier start au stop inclus.
    start_offsets
        Décalages (nt, multiples de 3, croissants, le premier vaut 0) des starts possibles.
    upstream
        Séquence située juste en amont du premier start (5'→3').
    """

    strand: int
    left: int
    right: int
    genome_length: int
    sequence: str
    start_offsets: np.ndarray
    upstream: str

    @property
    def length(self) -> int:
        """Longueur maximale (nt, stop inclus)."""
        return self.right - self.left

    def interval_for(self, offset: int) -> Interval:
        """Intervalle du gène obtenu en démarrant au start situé à ``offset``."""
        new_length = self.length - offset
        if self.strand == 1:
            new_left = (self.left + offset) % self.genome_length
            return Interval(new_left, new_left + new_length, 1)
        return Interval(self.left, self.left + new_length, -1)

    def upstream_of(self, offset: int, window: int) -> str:
        """Les ``window`` nucléotides en amont du start situé à ``offset``."""
        context = self.upstream + self.sequence
        position = len(self.upstream) + offset
        return context[max(0, position - window) : position]


@dataclass(frozen=True, slots=True)
class OrfSearch:
    """Résultat de la recherche d'ORF.

    Attributes
    ----------
    orfs
        ORF d'au moins ``min_length`` nucléotides.
    region_lengths
        Nombre de codons sens entre deux stops consécutifs, pour **toutes** les régions
        (y compris très courtes) : sert à comparer au modèle nul géométrique.
    """

    orfs: list[Orf]
    region_lengths: np.ndarray
    min_length: int


def _strand_string(sequence: str, strand: int, circular: bool) -> str:
    local = sequence if strand == 1 else reverse_complement(sequence)
    return local + local + local[:3] if circular else local


def find_orfs(
    sequence: str,
    code: GeneticCode,
    *,
    min_length: int,
    circular: bool,
    upstream_length: int = 30,
) -> OrfSearch:
    """Recherche tous les ORF d'au moins ``min_length`` nt (stop inclus).

    Parameters
    ----------
    sequence
        Génome en majuscules.
    code
        Code génétique (codons stop et codons start retenus).
    min_length
        Longueur minimale en nucléotides, codon stop inclus.
    circular
        Traite le génome comme circulaire.
    upstream_length
        Longueur de contexte amont conservée pour l'analyse du RBS.
    """
    n = len(sequence)
    min_codons = -(-min_length // 3)  # arrondi supérieur
    # Un ORF ne peut pas être plus long que le génome (cas pathologique d'un cadre
    # circulaire sans aucun codon stop en amont).
    max_codons = n // 3
    orfs: list[Orf] = []
    regions: list[np.ndarray] = []

    for strand in (1, -1):
        text = _strand_string(sequence, strand, circular)
        for frame in range(3):
            ids = encode_codons(text, frame)
            stop_idx = np.flatnonzero(code.stop_mask[ids])
            if stop_idx.size == 0:
                continue
            start_idx = np.flatnonzero(code.start_mask[ids])
            previous = np.concatenate(([-1], stop_idx[:-1]))
            stop_nt = frame + 3 * stop_idx

            # Chaque stop n'est retenu qu'une fois (seconde copie si circulaire).
            if circular:
                keep = (stop_nt >= n) & (stop_nt < 2 * n)
            else:
                keep = np.ones(stop_idx.size, dtype=bool)
            regions.append((stop_idx - previous - 1)[keep & (previous >= 0)])

            if start_idx.size == 0:
                continue
            lower_bound = np.maximum(previous + 1, stop_idx + 1 - max_codons)
            lo = np.searchsorted(start_idx, lower_bound, side="left")
            hi = np.searchsorted(start_idx, stop_idx, side="left")
            first_start = start_idx[np.minimum(lo, start_idx.size - 1)]
            n_codons = stop_idx - first_start + 1
            selected = np.flatnonzero(keep & (hi > lo) & (n_codons >= min_codons))

            for k in selected:
                starts = start_idx[lo[k] : hi[k]]
                starts = starts[stop_idx[k] - starts + 1 >= min_codons]
                begin = frame + 3 * int(starts[0])
                end = frame + 3 * int(stop_idx[k]) + 3
                length = end - begin
                left = begin % n if strand == 1 else (n - end) % n
                orfs.append(
                    Orf(
                        strand=strand,
                        left=left,
                        right=left + length,
                        genome_length=n,
                        sequence=text[begin:end],
                        start_offsets=(starts - starts[0]) * 3,
                        upstream=text[max(0, begin - upstream_length) : begin],
                    )
                )

    region_lengths = np.concatenate(regions) if regions else np.zeros(0, dtype=np.int64)
    return OrfSearch(orfs=orfs, region_lengths=region_lengths, min_length=min_length)
