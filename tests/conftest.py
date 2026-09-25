"""Fixtures partagées : chemins de démonstration, génome synthétique, GenBank minimal."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from orfeval.models import Genome, Interval
from orfeval.orfs.genetic_code import CODONS

ROOT = Path(__file__).resolve().parents[1]
DEMO_DIR = ROOT / "data" / "demo"
LAMBDA_PATH = DEMO_DIR / "NC_001416.1.gb.gz"
MGEN_PATH = DEMO_DIR / "NC_000908.2.gb.gz"


@pytest.fixture(scope="session")
def lambda_path() -> Path:
    return LAMBDA_PATH


@pytest.fixture(scope="session")
def mgen_path() -> Path:
    return MGEN_PATH


# --------------------------------------------------------------------------- synthétique
@dataclass(frozen=True)
class SyntheticGenome:
    genome: Genome
    genes: list[Interval]


def _codon_weights(rng: np.random.Generator) -> np.ndarray:
    """Usage des codons fortement biaisé (sans stop ni start interne ATG/GTG/TTG)."""
    weights = rng.gamma(0.4, 1.0, size=64)
    for codon in ("TAA", "TAG", "TGA", "ATG", "GTG", "TTG"):
        weights[CODONS.index(codon)] = 0.0
    return weights / weights.sum()


def make_synthetic_genome(seed: int = 7, n_genes: int = 80) -> SyntheticGenome:
    """Génome linéaire artificiel : gènes plantés séparés par de l'intergénique aléatoire.

    Les gènes suivent un usage des codons biaisé, l'intergénique est tiré uniformément :
    c'est la situation idéale où le modèle codant doit fonctionner.
    """
    rng = np.random.default_rng(seed)
    weights = _codon_weights(rng)
    parts: list[str] = []
    genes: list[Interval] = []
    position = 0
    for _ in range(n_genes):
        intergenic = "".join(rng.choice(list("ACGT"), size=int(rng.integers(60, 250))))
        parts.append(intergenic)
        position += len(intergenic)
        n_codons = int(rng.integers(100, 450))
        body = "".join(rng.choice(CODONS, size=n_codons, p=weights))
        gene = "ATG" + body + "TAA"
        strand = 1 if rng.random() < 0.6 else -1
        parts.append(gene if strand == 1 else reverse_complement(gene))
        genes.append(Interval(position, position + len(gene), strand))
        position += len(gene)
    parts.append("".join(rng.choice(list("ACGT"), size=200)))
    sequence = "".join(parts)
    genome = Genome("SYNTH_1", "génome synthétique", sequence, "linear", declared_table=11)
    return SyntheticGenome(genome=genome, genes=genes)


@pytest.fixture(scope="session")
def synthetic() -> SyntheticGenome:
    return make_synthetic_genome()


# --------------------------------------------------------------------------- GenBank minimal
def write_toy_genbank(path: Path, *, circular: bool = True) -> Path:
    """Écrit un petit GenBank couvrant les cas particuliers d'annotation.

    - CDS normal (+) 100..400, CDS (−) 500..800 ;
    - pseudogène 900..1200 ; CDS partielle <1300..1500 ;
    - ARNt 1600..1675 ;
    - CDS à cheval sur l'origine join(1901..2000,1..60) si circulaire.
    """
    rng = np.random.default_rng(3)
    sequence = "".join(rng.choice(list("ACGT"), size=2000))
    record = SeqRecord(Seq(sequence), id="TOY_1.1", name="TOY_1", description="toy genome")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular" if circular else "linear"
    record.annotations["organism"] = "Testus minimus"
    qualifiers = {"transl_table": ["11"]}
    record.features = [
        SeqFeature(
            SimpleLocation(99, 400, 1), "CDS", qualifiers={**qualifiers, "locus_tag": ["T1"]}
        ),
        SeqFeature(
            SimpleLocation(499, 800, -1),
            "CDS",
            qualifiers={**qualifiers, "locus_tag": ["T2"], "gene": ["dnaA"]},
        ),
        SeqFeature(
            SimpleLocation(899, 1200, 1),
            "CDS",
            qualifiers={**qualifiers, "locus_tag": ["T3"], "pseudo": [""]},
        ),
        SeqFeature(
            SimpleLocation(BeforePosition(1299), 1500, 1),
            "CDS",
            qualifiers={**qualifiers, "locus_tag": ["T4"]},
        ),
        SeqFeature(SimpleLocation(1599, 1675, 1), "tRNA", qualifiers={"locus_tag": ["R1"]}),
    ]
    if circular:
        location = CompoundLocation([SimpleLocation(1900, 2000, 1), SimpleLocation(0, 60, 1)])
        record.features.append(
            SeqFeature(location, "CDS", qualifiers={**qualifiers, "locus_tag": ["T5"]})
        )
    SeqIO.write(record, path, "genbank")
    return path


@pytest.fixture
def toy_genbank(tmp_path: Path) -> Path:
    return write_toy_genbank(tmp_path / "toy.gb")
