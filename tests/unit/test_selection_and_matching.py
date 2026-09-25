"""Choix du start, résolution des chevauchements et appariement à la référence."""

import numpy as np

from orfeval.evaluation.matching import match_genes
from orfeval.models import Interval, ReferenceGene
from orfeval.orfs.rbs import find_rbs
from orfeval.orfs.selection import SelectionCandidate, preferred_start, resolve_overlaps


# --------------------------------------------------------------------------- RBS
def test_rbs_detects_full_motif_with_spacer():
    upstream = "TTTTTTAGGAGGTTTTTTT"  # motif suivi de 7 nt avant le start
    hit = find_rbs(upstream)
    assert hit is not None
    assert (hit.motif, hit.spacer) == ("AGGAGG", 7)


def test_rbs_rejects_motif_too_close_to_start():
    # toute sous-chaîne du motif (AGG, GGA...) finit à moins de 4 nt du start
    assert find_rbs("TTTTTTTTTTTTTTTAGGAGG", spacer_min=4) is None


def test_rbs_prefers_longest_then_optimal_spacer():
    hit = find_rbs("GGAGTTTTTTTTTTTAGGAGTTTTTTTT", spacer_min=3, spacer_max=25)
    assert hit is not None
    assert hit.motif == "AGGAG"


# --------------------------------------------------------------------------- start
def test_preferred_start_strategies():
    scores = np.array([5.0, 9.0, 2.0])
    hits = [None, None, find_rbs("TTTTAGGAGGTTTTTTT")]
    assert preferred_start(scores, hits, "longest") == 0
    assert preferred_start(scores, hits, "rbs") == 2
    assert preferred_start(scores, hits, "score", rbs_weight=0.5) == 1
    assert preferred_start(scores, [None] * 3, "rbs") == 0


# --------------------------------------------------------------------------- chevauchements
def _single(left: int, right: int, strand: int = 1) -> SelectionCandidate:
    return SelectionCandidate(intervals=[Interval(left, right, strand)], preferred=0)


def test_overlap_resolution_keeps_best_candidate():
    candidates = [_single(0, 300), _single(200, 600, -1)]
    assert resolve_overlaps(candidates, [1, 0], 1000, max_overlap=60) == [(1, 0)]
    assert resolve_overlaps(candidates, [1, 0], 1000, max_overlap=100) == [(1, 0), (0, 0)]


def test_start_is_shifted_downstream_to_fit():
    upstream_gene = _single(0, 300)
    shiftable = SelectionCandidate(
        intervals=[Interval(250, 900, 1), Interval(290, 900, 1)], preferred=0
    )
    accepted = resolve_overlaps([upstream_gene, shiftable], [0, 1], 1000, max_overlap=20)
    assert accepted == [(0, 0), (1, 1)]


def test_threshold_after_selection_equals_threshold_before():
    """Propriété utilisée pour la courbe précision-rappel (voir selection.resolve_overlaps)."""
    rng = np.random.default_rng(5)
    candidates, scores = [], []
    for _ in range(200):
        left = int(rng.integers(0, 9000))
        candidates.append(
            _single(left, left + int(rng.integers(90, 900)), int(rng.choice([1, -1])))
        )
        scores.append(float(rng.normal()))
    order = sorted(range(len(candidates)), key=lambda i: -scores[i])
    for threshold in (-1.0, 0.0, 0.8):
        after = [
            i for i, _ in resolve_overlaps(candidates, order, 10_000, 60) if scores[i] > threshold
        ]
        before_order = [i for i in order if scores[i] > threshold]
        before = [i for i, _ in resolve_overlaps(candidates, before_order, 10_000, 60)]
        assert after == before


# --------------------------------------------------------------------------- appariement
def _ref(tag: str, left: int, right: int, strand: int) -> ReferenceGene:
    return ReferenceGene(Interval(left, right, strand), tag)


def test_matching_uses_strand_and_stop_position():
    references = [_ref("A", 100, 400, 1), _ref("B", 500, 800, -1)]
    predictions = [
        Interval(130, 400, 1),  # même stop que A, start différent
        Interval(500, 800, -1),  # identique à B
        Interval(100, 400, -1),  # mauvais brin
    ]
    result = match_genes(predictions, references, 1000)
    assert list(result.pred_correct) == [True, True, False]
    assert list(result.pred_start_correct) == [False, True, False]
    assert list(result.ref_found) == [True, True]
    assert list(result.ref_start_correct) == [False, True]
    assert result.pred_reference_tags == ["A", "B", ""]


def test_nested_genes_sharing_a_stop_are_both_found():
    references = [_ref("long", 100, 700, 1), _ref("court", 400, 700, 1)]
    result = match_genes([Interval(400, 700, 1)], references, 1000)
    assert list(result.ref_found) == [True, True]
    assert list(result.ref_start_correct) == [False, True]
    assert result.pred_reference_tags == ["long;court"]
