import pytest

from orfeval.models import Interval, extract_sequence, upstream_sequence


def test_interval_rejects_invalid_values():
    with pytest.raises(ValueError):
        Interval(10, 10, 1)
    with pytest.raises(ValueError):
        Interval(-1, 5, 1)
    with pytest.raises(ValueError):
        Interval(0, 5, 0)


@pytest.mark.parametrize(
    ("interval", "five", "three"),
    [
        (Interval(10, 40, 1), 10, 39),
        (Interval(10, 40, -1), 39, 10),
        # gène à cheval sur l'origine d'un génome de 100 nt
        (Interval(90, 130, 1), 90, 29),
        (Interval(90, 130, -1), 29, 90),
    ],
)
def test_gene_ends(interval, five, three):
    assert interval.five_prime(100) == five
    assert interval.three_prime(100) == three


def test_positions_wrap_around_origin():
    assert list(Interval(98, 102, 1).positions(100)) == [98, 99, 0, 1]


def test_extract_sequence_handles_strand_and_origin():
    sequence = "AAAACCCCGGGGTTTT"
    assert extract_sequence(sequence, Interval(4, 8, 1)) == "CCCC"
    assert extract_sequence(sequence, Interval(4, 8, -1)) == "GGGG"
    assert extract_sequence(sequence, Interval(14, 18, 1)) == "TTAA"


def test_upstream_sequence():
    sequence = "ACGTACGTAC" + "ATG" + "AAACCC"
    gene = Interval(10, 19, 1)
    assert upstream_sequence(sequence, gene, 4, circular=False) == "GTAC"
    assert upstream_sequence(sequence, Interval(2, 11, 1), 4, circular=False) is None
    assert upstream_sequence(sequence, Interval(2, 11, 1), 4, circular=True) == "CCAC"
    # brin − : l'amont se lit à droite du gène, en complément inverse
    assert upstream_sequence(sequence, Interval(0, 6, -1), 3, circular=False) == "TAC"


def test_minus_strand_across_origin():
    """Gène du brin − qui chevauche l'origine d'un génome circulaire de 10 nt."""
    sequence = "AACCGGTTAC"
    # brin + : "AC" (positions 8-9) + "AA" (0-1) ; brin − : complément inverse
    assert extract_sequence(sequence, Interval(8, 12, -1)) == "TTGT"
    gene = Interval(5, 10, -1)
    # l'amont d'un gène du brin − se lit à droite, ici après l'origine
    assert upstream_sequence(sequence, gene, 3, circular=True) == "GTT"
    assert upstream_sequence(sequence, gene, 3, circular=False) is None
