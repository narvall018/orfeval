"""Choix du codon start et résolution des chevauchements entre gènes candidats.

Deux décisions distinctes sont prises pour chaque ORF :

1. **Détection** : l'ORF contient-il un gène ? On utilise le meilleur score codant
   parmi tous ses starts possibles (voir :mod:`orfeval.orfs.predictor`).
2. **Placement du start** : parmi les starts possibles, lequel retenir ? Une stratégie
   fournit un start préféré ; en cas de chevauchement excessif avec un gène déjà
   accepté, le start est décalé vers l'aval (comme dans Glimmer) tant que le gène
   garde la longueur minimale.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from orfeval.models import Interval
from orfeval.orfs.rbs import RbsHit

StartStrategy = Literal["longest", "rbs", "score"]


def rbs_bonuses(hits: Sequence[RbsHit | None], weight: float) -> np.ndarray:
    """Bonus (bits) associé à chaque start : ``weight`` × longueur du motif SD."""
    return np.array([weight * hit.length if hit is not None else 0.0 for hit in hits], dtype=float)


def preferred_start(
    coding_scores: np.ndarray,
    hits: Sequence[RbsHit | None],
    strategy: StartStrategy,
    rbs_weight: float = 1.0,
) -> int:
    """Indice du start préféré parmi les starts possibles (ordonnés de l'amont vers l'aval).

    - ``"longest"`` : le start le plus en amont (ORF le plus long) ;
    - ``"rbs"`` : le start le plus en amont précédé d'un motif de Shine-Dalgarno,
      à défaut le plus en amont ;
    - ``"score"`` : le start qui maximise score codant + bonus RBS.
    """
    if strategy == "longest":
        return 0
    if strategy == "rbs":
        return next((index for index, hit in enumerate(hits) if hit is not None), 0)
    return int(np.argmax(coding_scores + rbs_bonuses(hits, rbs_weight)))


@dataclass(frozen=True, slots=True)
class SelectionCandidate:
    """ORF candidat pour la sélection gloutonne.

    Attributes
    ----------
    intervals
        Intervalle obtenu pour chaque start possible (amont → aval).
    preferred
        Indice du start préféré.
    """

    intervals: Sequence[Interval]
    preferred: int


def resolve_overlaps(
    candidates: Sequence[SelectionCandidate],
    priority: Sequence[int],
    genome_length: int,
    max_overlap: int,
) -> list[tuple[int, int]]:
    """Sélection gloutonne de gènes compatibles.

    Les candidats sont examinés dans l'ordre ``priority`` (meilleur d'abord). Pour chacun,
    on essaie le start préféré puis les starts situés en aval ; le premier dont le
    nombre de positions déjà couvertes par des gènes acceptés (tous brins confondus)
    ne dépasse pas ``max_overlap`` est retenu. Sinon, le candidat est rejeté.

    Propriété utile : l'acceptation d'un candidat ne dépend que des candidats mieux
    classés. Appliquer un seuil de score *après* la sélection donne donc le même résultat
    que l'appliquer avant, ce qui permet de tracer une courbe précision-rappel en une passe.

    Returns
    -------
    list[tuple[int, int]]
        Couples ``(indice du candidat, indice du start retenu)`` dans l'ordre de priorité.
    """
    occupied = np.zeros(genome_length, dtype=bool)
    accepted: list[tuple[int, int]] = []
    for index in priority:
        candidate = candidates[index]
        for start in range(candidate.preferred, len(candidate.intervals)):
            positions = candidate.intervals[start].positions(genome_length)
            if np.count_nonzero(occupied[positions]) <= max_overlap:
                occupied[positions] = True
                accepted.append((index, start))
                break
    return accepted
