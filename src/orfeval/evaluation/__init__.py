"""Évaluation des prédictions face à une annotation de référence."""

from orfeval.evaluation.matching import GeneMatching, match_genes, stop_key
from orfeval.evaluation.metrics import EvaluationResult, Proportion, evaluate, wilson_interval

__all__ = [
    "EvaluationResult",
    "GeneMatching",
    "Proportion",
    "evaluate",
    "match_genes",
    "stop_key",
    "wilson_interval",
]
