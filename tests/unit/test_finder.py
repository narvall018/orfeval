import numpy as np
import pytest
from Bio.Seq import reverse_complement

from orfeval.orfs.finder import find_orfs
from orfeval.orfs.genetic_code import CODON_INDEX, INVALID_CODON, GeneticCode, encode_codons

CODE_11 = GeneticCode.from_table(11, ("ATG", "GTG", "TTG"))
CODE_4 = GeneticCode.from_table(4, ("ATG", "GTG", "TTG"))
BODY = "GCA" * 30  # 30 codons sans stop ni start


def _pad(core: str) -> str:
    return "CC" + core + "CC"


def test_encode_codons_values_and_frames():
    assert list(encode_codons("ATGAAA")) == [CODON_INDEX["ATG"], CODON_INDEX["AAA"]]
    assert list(encode_codons("CATGA", frame=1)) == [CODON_INDEX["ATG"]]
    assert encode_codons("ANG")[0] == INVALID_CODON
    with pytest.raises(ValueError):
        encode_codons("ATG", frame=3)


def test_genetic_codes_differ_on_tga():
    tga = CODON_INDEX["TGA"]
    assert CODE_11.stop_mask[tga]
    assert not CODE_4.stop_mask[tga]


def test_start_codon_must_be_allowed_by_table():
    with pytest.raises(ValueError, match="non autorisé"):
        GeneticCode.from_table(11, ("AAA",))


def test_forward_orf_coordinates():
    gene = "ATG" + BODY + "TAA"
    search = find_orfs(_pad(gene), CODE_11, min_length=60, circular=False)
    orf = next(o for o in search.orfs if o.strand == 1 and o.length == len(gene))
    assert (orf.left, orf.right) == (2, 2 + len(gene))
    assert orf.sequence == gene


def test_reverse_orf_coordinates():
    gene = "ATG" + BODY + "TAA"
    sequence = _pad(reverse_complement(gene))
    search = find_orfs(sequence, CODE_11, min_length=60, circular=False)
    orf = next(o for o in search.orfs if o.strand == -1 and o.length == len(gene))
    assert (orf.left, orf.right) == (2, 2 + len(gene))
    assert orf.sequence == gene


def test_alternative_starts_are_listed():
    gene = "ATG" + "GCA" + "GTG" + BODY + "TAG"
    search = find_orfs(_pad(gene), CODE_11, min_length=60, circular=False)
    orf = next(o for o in search.orfs if o.strand == 1 and o.length == len(gene))
    assert list(orf.start_offsets) == [0, 6]
    shifted = orf.interval_for(6)
    assert (shifted.left, shifted.right) == (orf.left + 6, orf.right)


def test_min_length_filter():
    gene = "ATG" + "GCA" * 10 + "TAA"
    assert not [
        o
        for o in find_orfs(_pad(gene), CODE_11, min_length=90, circular=False).orfs
        if o.strand == 1
    ]


def test_genetic_code_4_reads_through_tga():
    gene = "ATG" + BODY + "TGA" + BODY + "TAA"
    long_11 = find_orfs(_pad(gene), CODE_11, min_length=150, circular=False).orfs
    long_4 = find_orfs(_pad(gene), CODE_4, min_length=150, circular=False).orfs
    assert not [o for o in long_11 if o.strand == 1]
    assert any(o.strand == 1 and o.length == len(gene) for o in long_4)


def test_circular_orf_across_origin():
    gene = "ATG" + BODY + "TAA"
    filler = "C" * 40
    # le gène commence 20 nt avant la fin et se termine au début de la séquence
    sequence = gene[20:] + filler + gene[:20]
    n = len(sequence)
    linear = find_orfs(sequence, CODE_11, min_length=60, circular=False)
    circular = find_orfs(sequence, CODE_11, min_length=60, circular=True)
    assert not any(o.strand == 1 and o.length == len(gene) for o in linear.orfs)
    orf = next(o for o in circular.orfs if o.strand == 1 and o.length == len(gene))
    assert orf.left == n - 20
    assert orf.right == n - 20 + len(gene)
    assert orf.sequence == gene


def test_circular_search_reports_each_stop_once():
    rng = np.random.default_rng(0)
    sequence = "".join(rng.choice(list("ACGT"), size=3001))
    search = find_orfs(sequence, CODE_11, min_length=60, circular=True)
    keys = [(o.strand, o.interval_for(0).three_prime(len(sequence))) for o in search.orfs]
    assert len(keys) == len(set(keys))
    # une région par codon stop : nombre total = nombre de stops sur les 6 cadres
    assert search.region_lengths.size > 0
