import numpy as np
import pytest

from orfeval.features.codon_usage import codon_counts, gc3_values, rscu_table
from orfeval.features.composition import analyse_composition, auto_window, gc_skew_profile
from orfeval.features.upstream import upstream_profile
from orfeval.models import Genome, Interval, ReferenceGene
from orfeval.orfs.genetic_code import CODON_INDEX


def test_gc_skew_profile_values():
    sequence = "G" * 100 + "C" * 100 + "GGCC" * 25
    profile = gc_skew_profile(sequence, 100)
    assert list(profile["skew"]) == pytest.approx([1.0, -1.0, 0.0])
    assert list(profile["cumulative_skew"]) == pytest.approx([1.0, 0.0, 0.0])


def test_short_last_window_is_dropped():
    profile = gc_skew_profile("G" * 100 + "C" * 20, 100)
    assert len(profile) == 1


def test_auto_window():
    assert auto_window(48_502) == 200
    assert auto_window(5_000) == 100


def test_replication_landmarks_only_for_circular_genomes():
    # skew négatif puis positif : minimum du cumul au milieu
    sequence = "C" * 5000 + "G" * 5000
    reference = ReferenceGene(Interval(4900, 5200, 1), "D1", gene_name="dnaA")
    circular = analyse_composition(Genome("c", "", sequence, "circular"), 500, [reference])
    assert circular.landmarks is not None
    assert circular.landmarks.origin_candidate == 4750
    assert circular.landmarks.dnaa_distance == 150
    linear = analyse_composition(Genome("l", "", sequence, "linear"), 500)
    assert linear.landmarks is None


def test_codon_counts_exclude_start_and_stop():
    counts = codon_counts(["ATG" + "AAA" * 3 + "GAA" + "TAA"])
    assert counts.sum() == 4
    assert counts[CODON_INDEX["AAA"]] == 3
    assert counts[CODON_INDEX["ATG"]] == 0


def test_rscu_sums_to_family_size():
    counts = codon_counts(["ATG" + "AAA" * 9 + "AAG" + "TGG" * 2 + "TAA"])
    table = rscu_table(counts, 11)
    lysine = table[table["amino_acid"] == "K"]
    assert lysine["rscu"].sum() == pytest.approx(2.0)
    assert lysine.set_index("codon").loc["AAA", "rscu"] == pytest.approx(1.8)
    # code 4 : TGA code le tryptophane, famille de 2 codons
    assert set(rscu_table(counts, 4).query("amino_acid == 'W'")["codon"]) == {"TGG", "TGA"}


def test_gc3():
    values = gc3_values(["ATG" + "AAG" + "AAA" + "TAA"])
    assert values[0] == pytest.approx(0.5)


def test_upstream_profile_detects_planted_motif():
    rng = np.random.default_rng(0)
    sequences = []
    for _ in range(200):
        background = "".join(rng.choice(list("ACGT"), size=20))
        sequences.append(background[:6] + "AGGAGG" + background[12:])
    profile = upstream_profile(sequences, "ACGT" * 1000)
    assert profile is not None
    assert profile.consensus[6:12] == "AGGAGG"
    peak = profile.information.idxmax()
    assert -14 <= peak <= -9
    assert upstream_profile(sequences[:5], "ACGT" * 100) is None
