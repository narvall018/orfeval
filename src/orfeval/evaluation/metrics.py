"""Métriques d'évaluation des prédictions face à l'annotation de référence.

Toutes les proportions sont accompagnées d'un intervalle de confiance à 95 % de
Wilson, mieux adapté que l'approximation normale pour les proportions proches de 0
ou 1 et les petits effectifs (classes de longueur, petits génomes).

Ces intervalles décrivent l'incertitude d'échantillonnage *sur ce génome*. Ils ne
disent rien de la performance sur d'autres génomes, ni des erreurs de l'annotation.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any

import numpy as np
import pandas as pd

from orfeval.evaluation.matching import GeneMatching, index_by_stop, match_genes, stop_key
from orfeval.exceptions import AnalysisError
from orfeval.io.genome import RNA_FEATURE_TYPES, LoadedGenome
from orfeval.models import (
    STRAND_SYMBOL,
    AnnotatedFeature,
    Interval,
    ReferenceGene,
    extract_sequence,
)
from orfeval.orfs.genetic_code import GeneticCode
from orfeval.orfs.predictor import PredictionResult

Z_95 = 1.959963984540054

STATUS_EXACT = "correct (start exact)"
STATUS_OTHER_START = "correct (start différent)"
FP_CATEGORIES = (
    "pseudogène ou CDS partielle",
    "gène d'ARN annoté",
    "même brin qu'un gène annoté",
    "antisens d'un gène annoté",
    "hors gènes annotés",
)
FN_CATEGORIES = (
    "stop non canonique (référence)",
    "score sous le seuil",
    "éliminé par chevauchement",
    "plus court que la longueur minimale",
    "aucun ORF candidat",
)
FOUND_EXACT = "retrouvé (start exact)"
FOUND_OTHER_START = "retrouvé (start différent)"


def wilson_interval(successes: int, total: int, z: float = Z_95) -> tuple[float, float]:
    """Intervalle de confiance de Wilson pour une proportion binomiale."""
    if total <= 0:
        return math.nan, math.nan
    if not 0 <= successes <= total:
        raise ValueError(f"effectifs incohérents : {successes}/{total}")
    p = successes / total
    denominator = 1.0 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    half_width = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    # Les bornes encadrent toujours p (évite les erreurs d'arrondi quand p vaut 0 ou 1).
    return max(0.0, min(p, centre - half_width)), min(1.0, max(p, centre + half_width))


@dataclass(frozen=True, slots=True)
class Proportion:
    """Proportion ``successes / total`` avec son intervalle de Wilson à 95 %."""

    successes: int
    total: int

    @property
    def value(self) -> float:
        """Valeur de la proportion (NaN si l'effectif est nul)."""
        return self.successes / self.total if self.total else math.nan

    @property
    def interval(self) -> tuple[float, float]:
        """Intervalle de confiance à 95 %."""
        return wilson_interval(self.successes, self.total)

    def to_dict(self) -> dict[str, Any]:
        """Représentation sérialisable."""
        low, high = self.interval
        return {
            "value": self.value,
            "successes": self.successes,
            "total": self.total,
            "ci95_low": low,
            "ci95_high": high,
        }

    def __str__(self) -> str:
        if not self.total:
            return "n/a (0/0)"
        low, high = self.interval
        return f"{self.value:.3f} [{low:.3f}–{high:.3f}] ({self.successes}/{self.total})"


def f1_score(precision: float, recall: float) -> float:
    """Moyenne harmonique de la précision et du rappel."""
    if math.isnan(precision) or math.isnan(recall) or precision + recall == 0:
        return math.nan if math.isnan(precision) or math.isnan(recall) else 0.0
    return 2 * precision * recall / (precision + recall)


def coverage_mask(intervals: Iterable[Interval], genome_length: int) -> np.ndarray:
    """Masque booléen (2, N) des positions couvertes : ligne 0 brin +, ligne 1 brin -."""
    mask = np.zeros((2, genome_length), dtype=bool)
    for interval in intervals:
        mask[0 if interval.strand == 1 else 1, interval.positions(genome_length)] = True
    return mask


def ranking_curve(
    ranked: Sequence[Interval],
    references: Sequence[ReferenceGene],
    genome_length: int,
    scores: Sequence[float] | None = None,
) -> pd.DataFrame:
    """Précision et rappel des ``k`` premières prédictions d'un classement, pour tout ``k``."""
    by_stop = index_by_stop(references, genome_length)
    gains = np.array(
        [len(by_stop.get(stop_key(interval, genome_length), [])) for interval in ranked],
        dtype=float,
    )
    ranks = np.arange(1, len(ranked) + 1)
    n_reference = max(len(references), 1)
    return pd.DataFrame(
        {
            "rank": ranks,
            "score": np.asarray(scores, dtype=float) if scores is not None else np.nan,
            "precision": np.cumsum(gains > 0) / ranks if len(ranked) else np.zeros(0),
            "recall": np.cumsum(gains) / n_reference,
        }
    )


def average_precision(curve: pd.DataFrame) -> float:
    """Aire sous la courbe précision-rappel (somme des gains de rappel × précision)."""
    if curve.empty:
        return math.nan
    recall = curve["recall"].to_numpy()
    precision = curve["precision"].to_numpy()
    gains = np.diff(np.concatenate(([0.0], recall)))
    return float(np.sum(gains * precision))


def length_class_labels(bins: Sequence[int]) -> list[str]:
    """Libellés des classes de longueur (nt)."""
    labels = [f"< {bins[0]}"]
    labels += [f"{low}–{high - 1}" for low, high in pairwise(bins)]
    labels.append(f"≥ {bins[-1]}")
    return labels


def recall_by_length(lengths: np.ndarray, found: np.ndarray, bins: Sequence[int]) -> pd.DataFrame:
    """Rappel par classe de longueur des gènes de référence, avec IC de Wilson."""
    labels = length_class_labels(bins)
    classes = np.digitize(lengths, bins)
    rows = []
    for index, label in enumerate(labels):
        in_class = classes == index
        proportion = Proportion(int(found[in_class].sum()), int(in_class.sum()))
        low, high = proportion.interval
        rows.append(
            {
                "classe": label,
                "n_reference": proportion.total,
                "retrouves": proportion.successes,
                "rappel": proportion.value,
                "ic95_bas": low,
                "ic95_haut": high,
            }
        )
    return pd.DataFrame(rows)


def classify_false_positives(
    intervals: Sequence[Interval],
    references: Sequence[ReferenceGene],
    other_features: Sequence[AnnotatedFeature],
    genome_length: int,
    min_fraction: float = 0.5,
) -> list[str]:
    """Catégorise les prédictions non appariées selon ce qu'elles recouvrent.

    Une prédiction « fausse » au regard de la référence n'est pas forcément une erreur :
    elle peut tomber sur un pseudogène (exclu de la référence), ou dans une région
    sans annotation (gène éventuellement manqué par l'annotation, ou faux positif).
    """
    excluded = coverage_mask(
        (f.interval for f in other_features if f.kind not in RNA_FEATURE_TYPES), genome_length
    ).any(axis=0)
    rna = coverage_mask(
        (f.interval for f in other_features if f.kind in RNA_FEATURE_TYPES), genome_length
    ).any(axis=0)
    reference_mask = coverage_mask((ref.interval for ref in references), genome_length)

    categories = []
    for interval in intervals:
        positions = interval.positions(genome_length)
        same = 0 if interval.strand == 1 else 1
        if excluded[positions].mean() >= min_fraction:
            categories.append(FP_CATEGORIES[0])
        elif rna[positions].mean() >= min_fraction:
            categories.append(FP_CATEGORIES[1])
        elif reference_mask[same, positions].mean() >= min_fraction:
            categories.append(FP_CATEGORIES[2])
        elif reference_mask[1 - same, positions].mean() >= min_fraction:
            categories.append(FP_CATEGORIES[3])
        else:
            categories.append(FP_CATEGORIES[4])
    return categories


@dataclass(frozen=True, slots=True)
class EvaluationResult:
    """Ensemble des métriques d'une prédiction."""

    n_reference: int
    n_predicted: int
    sensitivity: Proportion
    precision: Proportion
    f1: float
    start_accuracy: Proportion
    sensitivity_reachable: Proportion
    start_accuracy_reachable: Proportion
    nucleotide_sensitivity: float
    nucleotide_precision: float
    average_precision_model: float
    average_precision_length: float
    baseline_same_k: dict[str, float]
    by_length: pd.DataFrame
    pr_model: pd.DataFrame
    pr_length: pd.DataFrame
    reference_table: pd.DataFrame
    prediction_status: list[str]
    fp_categories: dict[str, int] = field(default_factory=dict)
    fn_categories: dict[str, int] = field(default_factory=dict)
    matching: GeneMatching | None = None

    def summary(self) -> dict[str, Any]:
        """Métriques principales sous forme sérialisable."""
        return {
            "n_reference": self.n_reference,
            "n_predicted": self.n_predicted,
            "sensitivity": self.sensitivity.to_dict(),
            "precision": self.precision.to_dict(),
            "f1": self.f1,
            "start_accuracy": self.start_accuracy.to_dict(),
            "sensitivity_reachable": self.sensitivity_reachable.to_dict(),
            "start_accuracy_reachable": self.start_accuracy_reachable.to_dict(),
            "nucleotide_sensitivity": self.nucleotide_sensitivity,
            "nucleotide_precision": self.nucleotide_precision,
            "average_precision_model": self.average_precision_model,
            "average_precision_length": self.average_precision_length,
            "baseline_same_k": self.baseline_same_k,
            "fp_categories": self.fp_categories,
            "fn_categories": self.fn_categories,
            "by_length": self.by_length,
        }


def _ordered_counts(values: Iterable[str], order: Sequence[str]) -> dict[str, int]:
    counts = Counter(values)
    return {category: counts.get(category, 0) for category in order}


def evaluate(
    prediction: PredictionResult,
    loaded: LoadedGenome,
    length_bins: Sequence[int] = (150, 300, 600, 1200),
) -> EvaluationResult:
    """Évalue une prédiction face aux CDS de référence du génome.

    Raises
    ------
    AnalysisError
        Si le génome ne comporte aucune CDS de référence.
    """
    references = loaded.references
    if not references:
        raise AnalysisError("aucune CDS de référence : évaluation impossible (fichier FASTA ?)")
    genome = loaded.genome
    n = genome.length
    genes = prediction.genes
    intervals = [gene.interval for gene in genes]
    matching = match_genes(intervals, references, n)

    code = GeneticCode.from_table(prediction.table_id, prediction.config.start_codons)
    min_length = prediction.config.min_length
    ref_sequences = [
        None if ref.compound else extract_sequence(genome.sequence, ref.interval)
        for ref in references
    ]
    canonical_stop = np.array(
        [seq is None or seq[-3:] in code.stop_codons for seq in ref_sequences], dtype=bool
    )
    allowed_start = np.array(
        [seq is not None and seq[:3] in code.start_codons for seq in ref_sequences], dtype=bool
    )
    lengths = np.array([ref.interval.length for ref in references])
    reachable = canonical_stop & (lengths >= min_length)

    found = matching.ref_found
    start_ok = matching.ref_start_correct
    sensitivity = Proportion(int(found.sum()), len(references))
    precision = Proportion(int(matching.pred_correct.sum()), len(genes))

    # Catégories d'erreurs --------------------------------------------------------
    ranked_keys = {stop_key(gene.interval, n) for gene in prediction.ranked}
    fn_status: list[str] = []
    for index, ref in enumerate(references):
        if found[index]:
            fn_status.append(FOUND_EXACT if start_ok[index] else FOUND_OTHER_START)
            continue
        key = stop_key(ref.interval, n)
        if not canonical_stop[index]:
            fn_status.append(FN_CATEGORIES[0])
        elif key in ranked_keys:
            fn_status.append(FN_CATEGORIES[1])
        elif key in prediction.orf_stop_keys:
            fn_status.append(FN_CATEGORIES[2])
        elif ref.interval.length < min_length:
            fn_status.append(FN_CATEGORIES[3])
        else:
            fn_status.append(FN_CATEGORIES[4])

    false_positive_idx = np.flatnonzero(~matching.pred_correct)
    fp_labels = classify_false_positives(
        [intervals[i] for i in false_positive_idx], references, loaded.other_features, n
    )
    prediction_status = [
        STATUS_EXACT if matching.pred_start_correct[i] else STATUS_OTHER_START
        for i in range(len(genes))
    ]
    for index, label in zip(false_positive_idx.tolist(), fp_labels, strict=True):
        prediction_status[index] = label

    # Niveau nucléotidique ---------------------------------------------------------
    ref_mask = coverage_mask((ref.interval for ref in references), n)
    pred_mask = coverage_mask(intervals, n)
    overlap = int(np.count_nonzero(ref_mask & pred_mask))
    n_ref_bases = int(np.count_nonzero(ref_mask))
    n_pred_bases = int(np.count_nonzero(pred_mask))
    nt_sensitivity = overlap / max(n_ref_bases, 1)
    nt_precision = overlap / n_pred_bases if n_pred_bases else math.nan

    # Courbes précision-rappel : modèle codant vs longueur seule -------------------
    pr_model = ranking_curve(
        [gene.interval for gene in prediction.ranked],
        references,
        n,
        scores=[gene.score for gene in prediction.ranked],
    )
    pr_length = ranking_curve(
        prediction.length_ranked,
        references,
        n,
        scores=[interval.length for interval in prediction.length_ranked],
    )
    k = len(genes)
    baseline = {"n_predictions": k, "precision": math.nan, "recall": math.nan, "f1": math.nan}
    if 0 < k <= len(pr_length):
        row = pr_length.iloc[k - 1]
        baseline = {
            "n_predictions": k,
            "precision": float(row["precision"]),
            "recall": float(row["recall"]),
            "f1": f1_score(float(row["precision"]), float(row["recall"])),
        }

    reference_table = pd.DataFrame(
        {
            "locus_tag": [ref.locus_tag for ref in references],
            "gene": [ref.gene_name or "" for ref in references],
            "product": [ref.product for ref in references],
            "start": [ref.interval.left + 1 for ref in references],
            "end": [ref.interval.right for ref in references],
            "strand": [STRAND_SYMBOL[ref.interval.strand] for ref in references],
            "length_nt": lengths,
            "start_codon": [seq[:3] if seq else "" for seq in ref_sequences],
            "stop_codon": [seq[-3:] if seq else "" for seq in ref_sequences],
            "found": found,
            "start_correct": start_ok,
            "status": fn_status,
        }
    )

    return EvaluationResult(
        n_reference=len(references),
        n_predicted=len(genes),
        sensitivity=sensitivity,
        precision=precision,
        f1=f1_score(precision.value, sensitivity.value),
        start_accuracy=Proportion(int(start_ok.sum()), int(found.sum())),
        sensitivity_reachable=Proportion(int(found[reachable].sum()), int(reachable.sum())),
        start_accuracy_reachable=Proportion(
            int(start_ok[found & allowed_start].sum()), int((found & allowed_start).sum())
        ),
        nucleotide_sensitivity=float(nt_sensitivity),
        nucleotide_precision=float(nt_precision),
        average_precision_model=average_precision(pr_model),
        average_precision_length=average_precision(pr_length),
        baseline_same_k=baseline,
        by_length=recall_by_length(lengths, found, length_bins),
        pr_model=pr_model,
        pr_length=pr_length,
        reference_table=reference_table,
        prediction_status=prediction_status,
        fp_categories=_ordered_counts(fp_labels, FP_CATEGORIES),
        fn_categories=_ordered_counts(
            (status for status in fn_status if status in FN_CATEGORIES), FN_CATEGORIES
        ),
        matching=matching,
    )
