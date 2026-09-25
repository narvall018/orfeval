"""Modèle nul, modèle codant et métriques statistiques."""

import math

import numpy as np
import pandas as pd
import pytest

from orfeval.evaluation.metrics import (
    Proportion,
    average_precision,
    f1_score,
    recall_by_length,
    wilson_interval,
)
from orfeval.orfs.coding_model import CodingModel
from orfeval.orfs.genetic_code import CODON_INDEX, GeneticCode, encode_codons
from orfeval.orfs.null_model import NullModel, strand_symmetric_composition

CODE_11 = GeneticCode.from_table(11, ("ATG",))
UNIFORM = dict.fromkeys("ACGT", 0.25)


# --------------------------------------------------------------------------- modèle nul
def test_composition_is_strand_symmetric():
    composition = strand_symmetric_composition("AAAAGG")
    assert composition["A"] == composition["T"] == pytest.approx(4 / 12)
    assert composition["C"] == composition["G"] == pytest.approx(2 / 12)


def test_stop_probability_uniform_composition():
    model = NullModel.from_sequence("ACGT" * 100, CODE_11.stop_codons)
    assert model.p_stop == pytest.approx(3 / 64)


def test_significant_length_is_the_first_length_below_target():
    model = NullModel.from_sequence("ACGT" * 25_000, CODE_11.stop_codons)
    length = model.significant_length(expected=1.0)
    n_codons = length // 3 - 1
    assert model.expected_count(n_codons) <= 1.0
    assert model.expected_count(n_codons - 1) > 1.0


def test_null_model_matches_simulated_random_sequence():
    """Sur une séquence réellement aléatoire, la loi géométrique doit être respectée."""
    from orfeval.orfs.finder import find_orfs

    rng = np.random.default_rng(1)
    sequence = "".join(rng.choice(list("ACGT"), size=60_000))
    model = NullModel.from_sequence(sequence, CODE_11.stop_codons)
    regions = find_orfs(sequence, CODE_11, min_length=3000, circular=True).region_lengths
    for n_codons in (10, 30):
        observed = np.mean(regions >= n_codons)
        assert observed == pytest.approx(model.survival(n_codons), rel=0.1)


# --------------------------------------------------------------------------- modèle codant
def test_coding_model_learns_codon_preference():
    training = ["ATG" + "AAA" * 50 + "TAA"] * 5
    model = CodingModel.train(training, CODE_11, UNIFORM, kind="codon", pseudocount=0.5)
    assert model.log_ratio[CODON_INDEX["AAA"]] > 0
    assert model.log_ratio[CODON_INDEX["CCC"]] < 0
    assert model.n_codons == 250


def test_start_scores_are_suffix_sums():
    training = ["ATG" + "AAA" * 50 + "TAA"] * 5
    model = CodingModel.train(training, CODE_11, UNIFORM, kind="codon")
    ids = encode_codons("ATG" + "AAA" + "ATG" + "AAA" * 2 + "TAA")
    per_codon = model.codon_scores(ids)
    scores = model.start_scores(ids, np.array([0, 2]))
    assert scores[0] == pytest.approx(per_codon[1:5].sum())
    assert scores[1] == pytest.approx(per_codon[3:5].sum())


def test_dicodon_model_trains_and_scores():
    training = ["ATG" + "AAACCC" * 30 + "TAA"] * 4
    model = CodingModel.train(training, CODE_11, UNIFORM, kind="dicodon")
    ids = encode_codons("ATG" + "AAACCC" * 5 + "TAA")
    assert model.log_ratio.shape == (65, 65)
    assert model.start_scores(ids, np.array([0]))[0] > 0


def test_coding_model_requires_training_data():
    with pytest.raises(ValueError):
        CodingModel.train([], CODE_11, UNIFORM)


# --------------------------------------------------------------------------- métriques
def test_wilson_interval_reference_values():
    low, high = wilson_interval(8, 10)
    assert low == pytest.approx(0.4902, abs=1e-4)
    assert high == pytest.approx(0.9433, abs=1e-4)


@pytest.mark.parametrize(("k", "n"), [(0, 7), (7, 7), (504, 504)])
def test_wilson_bounds_contain_the_estimate(k, n):
    low, high = wilson_interval(k, n)
    assert 0.0 <= low <= k / n <= high <= 1.0


def test_wilson_edge_cases():
    assert all(math.isnan(bound) for bound in wilson_interval(0, 0))
    with pytest.raises(ValueError):
        wilson_interval(5, 3)


def test_proportion_and_f1():
    proportion = Proportion(3, 4)
    assert proportion.value == 0.75
    assert "3/4" in str(proportion)
    assert f1_score(0.5, 1.0) == pytest.approx(2 / 3)
    assert f1_score(0.0, 0.0) == 0.0
    assert math.isnan(Proportion(0, 0).value)


def test_average_precision_of_a_perfect_ranking():
    curve = pd.DataFrame({"precision": [1.0, 1.0, 1.0], "recall": [1 / 3, 2 / 3, 1.0]})
    assert average_precision(curve) == pytest.approx(1.0)


def test_recall_by_length_bins():
    lengths = np.array([100, 200, 250, 700, 2000])
    found = np.array([False, True, False, True, True])
    table = recall_by_length(lengths, found, (150, 300, 600, 1200))
    assert list(table["n_reference"]) == [1, 2, 0, 1, 1]
    assert list(table["retrouves"]) == [0, 1, 0, 1, 1]
    assert math.isnan(table.loc[2, "rappel"])
