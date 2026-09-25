"""Structures de données partagées et utilitaires de coordonnées.

Convention de coordonnées utilisée dans tout le package
-------------------------------------------------------
Les intervalles sont **0-based et demi-ouverts** ``[left, right)`` sur le brin direct,
comme en Python et dans Biopython. ``left`` est toujours dans ``[0, N)`` ;
``right`` peut dépasser ``N`` lorsqu'un gène chevauche l'origine d'un génome
circulaire (la partie excédentaire se lit alors au début de la séquence).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from Bio.Seq import reverse_complement

Topology = Literal["linear", "circular"]
STRAND_SYMBOL: dict[int, str] = {1: "+", -1: "-"}


@dataclass(frozen=True, slots=True)
class Genome:
    """Séquence génomique nettoyée, prête pour l'analyse.

    Attributes
    ----------
    seq_id
        Identifiant de la séquence (ex. ``NC_001416.1``).
    description
        Description libre issue du fichier d'entrée.
    sequence
        Séquence en majuscules (alphabet IUPAC ADN).
    topology
        ``"linear"`` ou ``"circular"``.
    organism
        Nom de l'organisme s'il est connu.
    declared_table
        Code génétique déclaré par l'annotation (qualifiant ``/transl_table``).
    """

    seq_id: str
    description: str
    sequence: str
    topology: Topology
    organism: str | None = None
    declared_table: int | None = None

    @property
    def length(self) -> int:
        """Longueur de la séquence en nucléotides."""
        return len(self.sequence)

    @property
    def is_circular(self) -> bool:
        """Vrai si le génome est circulaire."""
        return self.topology == "circular"


@dataclass(frozen=True, slots=True)
class Interval:
    """Intervalle orienté sur le génome (voir la convention du module)."""

    left: int
    right: int
    strand: int

    def __post_init__(self) -> None:
        if self.strand not in (1, -1):
            raise ValueError(f"brin invalide : {self.strand!r} (attendu 1 ou -1)")
        if self.left < 0 or self.right <= self.left:
            raise ValueError(f"intervalle invalide : [{self.left}, {self.right})")

    @property
    def length(self) -> int:
        """Longueur en nucléotides."""
        return self.right - self.left

    def three_prime(self, genome_length: int) -> int:
        """Position 0-based (modulo N) de la dernière base du codon stop."""
        position = self.right - 1 if self.strand == 1 else self.left
        return position % genome_length

    def five_prime(self, genome_length: int) -> int:
        """Position 0-based (modulo N) de la première base du codon start."""
        position = self.left if self.strand == 1 else self.right - 1
        return position % genome_length

    def positions(self, genome_length: int) -> np.ndarray:
        """Positions génomiques couvertes (modulo N, gère le passage de l'origine)."""
        return np.arange(self.left, self.right, dtype=np.int64) % genome_length


@dataclass(frozen=True, slots=True)
class ReferenceGene:
    """Gène codant (CDS) issu de l'annotation de référence."""

    interval: Interval
    locus_tag: str
    product: str = ""
    gene_name: str | None = None
    compound: bool = False


@dataclass(frozen=True, slots=True)
class AnnotatedFeature:
    """Autre élément annoté (pseudogène, CDS partielle, gène d'ARN...).

    Ces éléments ne font pas partie de la référence mais aident à interpréter les
    prédictions qui ne correspondent à aucun gène de référence.
    """

    interval: Interval
    kind: str
    label: str


@dataclass(frozen=True, slots=True)
class PredictedGene:
    """Gène prédit par orfeval.

    Attributes
    ----------
    interval
        Coordonnées du gène (codon start et codon stop inclus).
    sequence
        Séquence nucléotidique 5'→3' sur le brin du gène.
    score
        Score de **détection** (bits) : meilleur log-rapport de vraisemblance codant /
        fond parmi tous les starts possibles de l'ORF. Sert à classer et filtrer.
    start_score
        Log-rapport de vraisemblance calculé depuis le start retenu (bits).
    rbs_motif, rbs_spacer
        Motif de Shine-Dalgarno détecté en amont du start retenu et distance au start.
    n_alternative_starts
        Nombre de codons start possibles dans l'ORF.
    start_adjusted
        Vrai si le start a été décalé vers l'aval pour limiter un chevauchement.
    """

    interval: Interval
    sequence: str
    score: float
    start_score: float
    rbs_motif: str | None
    rbs_spacer: int | None
    n_alternative_starts: int
    start_adjusted: bool = False

    @property
    def length(self) -> int:
        """Longueur en nucléotides (codon stop inclus)."""
        return self.interval.length

    @property
    def start_codon(self) -> str:
        """Codon d'initiation."""
        return self.sequence[:3]

    @property
    def stop_codon(self) -> str:
        """Codon de terminaison."""
        return self.sequence[-3:]


def extract_sequence(sequence: str, interval: Interval) -> str:
    """Extrait la séquence d'un intervalle, orientée 5'→3' sur son brin.

    Gère les intervalles qui dépassent la fin de la séquence (génome circulaire).
    """
    n = len(sequence)
    if interval.right <= n:
        segment = sequence[interval.left : interval.right]
    else:
        segment = sequence[interval.left :] + sequence[: interval.right - n]
    return segment if interval.strand == 1 else reverse_complement(segment)


def upstream_sequence(
    sequence: str, interval: Interval, length: int, *, circular: bool
) -> str | None:
    """Renvoie les ``length`` nucléotides situés en amont du codon start.

    La séquence est orientée 5'→3' sur le brin du gène ; son dernier nucléotide est
    adjacent au codon start. Renvoie ``None`` si la région sort d'un génome linéaire.
    """
    n = len(sequence)
    if interval.strand == 1:
        start, end = interval.left - length, interval.left
    else:
        start, end = interval.right, interval.right + length
    if not circular and (start < 0 or end > n):
        return None
    indices = np.arange(start, end) % n
    segment = "".join(sequence[i] for i in indices)
    return segment if interval.strand == 1 else reverse_complement(segment)
