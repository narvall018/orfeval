"""Figures matplotlib. Ce module ne fait aucun calcul métier : il reçoit des résultats."""

from orfeval.plotting.composition import orf_length_null_figure, rscu_figure, upstream_figure
from orfeval.plotting.evaluation import (
    comparison_figure,
    error_breakdown_figure,
    precision_recall_figure,
    recall_by_length_figure,
    score_length_figure,
)
from orfeval.plotting.genome import EVALUATION_LEGEND, gc_skew_figure, genome_map
from orfeval.plotting.style import save_figure

__all__ = [
    "EVALUATION_LEGEND",
    "comparison_figure",
    "error_breakdown_figure",
    "gc_skew_figure",
    "genome_map",
    "orf_length_null_figure",
    "precision_recall_figure",
    "recall_by_length_figure",
    "rscu_figure",
    "save_figure",
    "score_length_figure",
    "upstream_figure",
]
