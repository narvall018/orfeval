"""Appariement des gènes prédits avec les gènes de référence.

Critère standard pour les gènes procaryotes (utilisé notamment pour évaluer Prodigal) :
une prédiction est **correcte** si elle partage le brin et la position du codon stop
(extrémité 3') d'un gène de référence. Le codon start est évalué séparément, car sa
position est beaucoup plus difficile à prédire.

Cas particulier : deux gènes emboîtés dans le même cadre (starts différents, même stop)
partagent la même extrémité 3'. Une prédiction unique les « retrouve » alors tous deux ;
seul l'un peut avoir le bon start.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from orfeval.models import Interval, ReferenceGene

StopKey = tuple[int, int]


def stop_key(interval: Interval, genome_length: int) -> StopKey:
    """Clé d'appariement : (brin, position de la dernière base du codon stop)."""
    return interval.strand, interval.three_prime(genome_length)


def index_by_stop(
    references: Sequence[ReferenceGene], genome_length: int
) -> dict[StopKey, list[int]]:
    """Indices des gènes de référence regroupés par extrémité 3'."""
    index: dict[StopKey, list[int]] = defaultdict(list)
    for position, reference in enumerate(references):
        index[stop_key(reference.interval, genome_length)].append(position)
    return dict(index)


@dataclass(frozen=True, slots=True)
class GeneMatching:
    """Résultat de l'appariement prédictions ↔ référence.

    Attributes
    ----------
    ref_found, ref_start_correct
        Pour chaque gène de référence : retrouvé (même 3') et start exact.
    pred_correct, pred_start_correct
        Pour chaque prédiction : appariée à la référence et start exact.
    pred_reference_tags
        Pour chaque prédiction, les ``locus_tag`` appariés (séparés par ``;``).
    """

    ref_found: np.ndarray
    ref_start_correct: np.ndarray
    pred_correct: np.ndarray
    pred_start_correct: np.ndarray
    pred_reference_tags: list[str]


def match_genes(
    predictions: Sequence[Interval], references: Sequence[ReferenceGene], genome_length: int
) -> GeneMatching:
    """Apparie des intervalles prédits aux gènes de référence par extrémité 3'."""
    by_stop = index_by_stop(references, genome_length)
    ref_found = np.zeros(len(references), dtype=bool)
    ref_start = np.zeros(len(references), dtype=bool)
    pred_correct = np.zeros(len(predictions), dtype=bool)
    pred_start = np.zeros(len(predictions), dtype=bool)
    tags: list[str] = []

    for position, interval in enumerate(predictions):
        matched = by_stop.get(stop_key(interval, genome_length), [])
        tags.append(";".join(references[i].locus_tag for i in matched))
        if not matched:
            continue
        pred_correct[position] = True
        five_prime = interval.five_prime(genome_length)
        for i in matched:
            ref_found[i] = True
            if references[i].interval.five_prime(genome_length) == five_prime:
                ref_start[i] = True
                pred_start[position] = True

    return GeneMatching(
        ref_found=ref_found,
        ref_start_correct=ref_start,
        pred_correct=pred_correct,
        pred_start_correct=pred_start,
        pred_reference_tags=tags,
    )
