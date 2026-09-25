"""Modèle statistique du potentiel codant.

Le score d'un gène candidat est un **log-rapport de vraisemblance** (en bits) :

    S = somme sur les codons i de log2( P_codant(c_i | c_i-1) / P_fond(c_i) )

- ``P_codant`` est appris sur le génome lui-même (auto-apprentissage) : fréquences des
  codons (ordre 0, modèle « codon ») ou des codons conditionnellement au précédent
  (ordre 1, modèle « dicodon »), avec pseudo-comptes.
- ``P_fond`` suppose des nucléotides indépendants de même composition que le génome,
  **conditionnée à l'absence de codon stop** : un ORF ne contient par définition aucun
  stop interne, et sans ce conditionnement tout ORF long paraîtrait codant.

Le codon start et le codon stop ne sont pas comptés : leur distribution est
contrainte par la définition même d'un ORF.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from orfeval.orfs.genetic_code import (
    CODONS,
    INVALID_CODON,
    N_CODON_STATES,
    GeneticCode,
    encode_codons,
)
from orfeval.orfs.null_model import codon_probability

ModelKind = Literal["codon", "dicodon"]


def background_log_probabilities(code: GeneticCode, composition: dict[str, float]) -> np.ndarray:
    """``log2 P_fond(codon)`` sur les 64 codons, renormalisé sur les codons sens.

    Les codons stop reçoivent 0 (ils ne sont jamais évalués).
    """
    sense = code.sense_mask[:64]
    probabilities = np.array([codon_probability(codon, composition) for codon in CODONS])
    probabilities[~sense] = 0.0
    probabilities /= probabilities.sum()
    log_bg = np.zeros(64)
    log_bg[sense] = np.log2(probabilities[sense])
    return log_bg


@dataclass(frozen=True, slots=True)
class CodingModel:
    """Modèle codant entraîné.

    Attributes
    ----------
    kind
        ``"codon"`` (ordre 0) ou ``"dicodon"`` (ordre 1).
    log_ratio
        Scores par codon en bits : forme (65,) pour l'ordre 0, (65, 65) ``[précédent,
        courant]`` pour l'ordre 1. L'indice 64 (codon ambigu) vaut 0.
    n_genes, n_codons
        Taille de l'ensemble d'entraînement.
    """

    kind: ModelKind
    log_ratio: np.ndarray
    n_genes: int
    n_codons: int

    @classmethod
    def train(
        cls,
        sequences: Sequence[str],
        code: GeneticCode,
        composition: dict[str, float],
        *,
        kind: ModelKind = "codon",
        pseudocount: float = 1.0,
    ) -> CodingModel:
        """Entraîne le modèle sur des séquences codantes (start → stop inclus).

        Raises
        ------
        ValueError
            Si aucune séquence d'entraînement n'est fournie.
        """
        if not sequences:
            raise ValueError("ensemble d'entraînement vide")
        sense = code.sense_mask[:64]
        log_bg = background_log_probabilities(code, composition)

        single = np.zeros(64)
        pairs = np.zeros((64, 64))
        n_codons = 0
        for sequence in sequences:
            ids = encode_codons(sequence)[:-1]  # retire le stop
            inner = ids[1:]  # retire le start
            valid_inner = inner[inner < 64]
            single += np.bincount(valid_inner, minlength=64)
            n_codons += valid_inner.size
            if kind == "dicodon":
                previous, current = ids[:-1], ids[1:]
                ok = (previous < 64) & (current < 64)
                np.add.at(pairs, (previous[ok], current[ok]), 1.0)

        single = single + pseudocount * sense
        single_log = np.zeros(64)
        single_log[sense] = np.log2(single[sense] / single.sum()) - log_bg[sense]

        if kind == "codon":
            log_ratio = np.zeros(N_CODON_STATES)
            log_ratio[:64] = single_log
        else:
            pairs = pairs + pseudocount * sense[np.newaxis, :]
            conditional = pairs / pairs.sum(axis=1, keepdims=True)
            log_ratio = np.zeros((N_CODON_STATES, N_CODON_STATES))
            block = np.zeros((64, 64))
            block[:, sense] = np.log2(conditional[:, sense]) - log_bg[np.newaxis, sense]
            log_ratio[:64, :64] = block
            log_ratio[INVALID_CODON, :64] = single_log  # contexte inconnu : ordre 0
        return cls(kind=kind, log_ratio=log_ratio, n_genes=len(sequences), n_codons=n_codons)

    def codon_scores(self, ids: np.ndarray) -> np.ndarray:
        """Score (bits) de chaque codon d'une séquence encodée, dans son contexte."""
        if self.kind == "codon":
            return np.asarray(self.log_ratio[ids], dtype=float)
        previous = np.concatenate(([INVALID_CODON], ids[:-1]))
        return np.asarray(self.log_ratio[previous, ids], dtype=float)

    def start_scores(self, ids: np.ndarray, start_codon_indices: np.ndarray) -> np.ndarray:
        """Score codant pour chaque codon start possible d'un ORF.

        Parameters
        ----------
        ids
            Codons encodés de l'ORF, du premier start jusqu'au stop inclus.
        start_codon_indices
            Indices (en codons) des starts possibles dans ``ids``.

        Returns
        -------
        numpy.ndarray
            Pour chaque start ``j`` : somme des scores des codons ``j+1`` à l'avant-dernier.
        """
        scores = self.codon_scores(ids)[:-1]  # le stop n'est pas évalué
        suffix = np.concatenate((np.cumsum(scores[::-1])[::-1], [0.0]))
        return np.asarray(suffix[start_codon_indices + 1], dtype=float)

    def per_codon_table(self) -> dict[str, float]:
        """Contribution de chaque codon (ordre 0) : utile pour inspecter le modèle."""
        table = self.log_ratio if self.kind == "codon" else self.log_ratio[INVALID_CODON]
        return {codon: float(table[index]) for index, codon in enumerate(CODONS)}
