"""Modèle nul de la longueur des ORF.

Dans une séquence aléatoire (nucléotides indépendants, composition du génome), un codon
est un stop avec une probabilité ``p``. Le nombre de codons sens entre deux stops suit
alors une loi géométrique : ``P(L >= k) = (1 - p)^k``.

Ce modèle sert à deux choses :

1. choisir *sans utiliser l'annotation* une longueur au-delà de laquelle les ORF sont
   très probablement de vrais gènes (ensemble d'auto-apprentissage) ;
2. visualiser l'excès d'ORF longs dans le génome réel par rapport au hasard.

Limite : un génome réel n'est pas une suite de nucléotides indépendants (biais de
dinucléotides, répétitions) ; le modèle est une approximation de référence.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np


def strand_symmetric_composition(sequence: str) -> dict[str, float]:
    """Fréquences nucléotidiques moyennées sur les deux brins.

    Les ORF étant recherchés sur les deux brins, on utilise ``p(A) = p(T)`` et
    ``p(C) = p(G)``. Les bases ambiguës sont ignorées.
    """
    counts = {base: sequence.count(base) for base in "ACGT"}
    total = sum(counts.values())
    if total == 0:
        raise ValueError("séquence sans base A, C, G ou T")
    at = (counts["A"] + counts["T"]) / (2 * total)
    gc = (counts["C"] + counts["G"]) / (2 * total)
    return {"A": at, "T": at, "C": gc, "G": gc}


def codon_probability(codon: str, composition: dict[str, float]) -> float:
    """Probabilité d'un codon sous l'hypothèse de nucléotides indépendants."""
    return math.prod(composition[base] for base in codon)


@dataclass(frozen=True, slots=True)
class NullModel:
    """Modèle géométrique de la longueur des régions sans stop.

    Attributes
    ----------
    genome_length
        Longueur du génome (nt).
    p_stop
        Probabilité qu'un codon aléatoire soit un codon stop.
    composition
        Composition nucléotidique symétrique utilisée.
    """

    genome_length: int
    p_stop: float
    composition: dict[str, float]

    @classmethod
    def from_sequence(cls, sequence: str, stop_codons: Iterable[str]) -> NullModel:
        """Estime le modèle à partir de la composition d'un génome."""
        composition = strand_symmetric_composition(sequence)
        p_stop = sum(codon_probability(codon, composition) for codon in stop_codons)
        if not 0.0 < p_stop < 1.0:
            raise ValueError(f"probabilité de stop dégénérée : {p_stop}")
        return cls(genome_length=len(sequence), p_stop=p_stop, composition=composition)

    @property
    def expected_regions(self) -> float:
        """Nombre attendu de régions entre stops (6 cadres × N/3 codons × p)."""
        return 2.0 * self.genome_length * self.p_stop

    def survival(self, n_codons: np.ndarray | float) -> np.ndarray:
        """``P(L >= n_codons)`` pour une région aléatoire."""
        return np.asarray(np.power(1.0 - self.p_stop, np.asarray(n_codons, dtype=float)))

    def expected_count(self, n_codons: np.ndarray | float) -> np.ndarray:
        """Nombre attendu de régions aléatoires d'au moins ``n_codons`` codons sens."""
        return self.expected_regions * self.survival(n_codons)

    def significant_length(self, expected: float = 1.0) -> int:
        """Longueur (nt, stop inclus) au-delà de laquelle un ORF est rare par hasard.

        Renvoie la plus petite longueur pour laquelle moins de ``expected`` ORF
        aléatoires sont attendus dans tout le génome.
        """
        if expected <= 0:
            raise ValueError("expected doit être > 0")
        ratio = expected / self.expected_regions
        n_codons = 0 if ratio >= 1 else math.ceil(math.log(ratio) / math.log(1.0 - self.p_stop))
        return 3 * (n_codons + 1)
