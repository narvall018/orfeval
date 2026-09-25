"""Lecture des génomes, écriture des résultats et configuration."""

import gzip
import json
import math
from pathlib import Path

import pytest
import yaml

from orfeval.config import PredictorConfig, load_config
from orfeval.exceptions import ConfigError, InputFormatError
from orfeval.io.genome import detect_format, load_genome
from orfeval.io.ncbi import validate_accession, write_gzip_deterministic
from orfeval.io.writers import gff3_text, proteins_fasta_text, write_json
from orfeval.models import Genome, Interval, PredictedGene

ROOT = Path(__file__).resolve().parents[2]


# --------------------------------------------------------------------------- GenBank
def test_toy_genbank_reference_extraction(toy_genbank):
    loaded = load_genome(toy_genbank)
    tags = {ref.locus_tag: ref for ref in loaded.references}
    assert set(tags) == {"T1", "T2", "T5"}
    assert loaded.annotation.n_pseudo == 1
    assert loaded.annotation.n_partial == 1
    assert loaded.genome.declared_table == 11
    assert loaded.genome.topology == "circular"
    # CDS à cheval sur l'origine : un seul intervalle qui dépasse N
    assert (tags["T5"].interval.left, tags["T5"].interval.right) == (1900, 2060)
    kinds = sorted(feature.kind for feature in loaded.other_features)
    assert kinds == ["CDS partielle", "pseudogène", "tRNA"]


def test_gzip_and_fasta_inputs(tmp_path, toy_genbank):
    compressed = tmp_path / "toy.gb.gz"
    compressed.write_bytes(gzip.compress(toy_genbank.read_bytes()))
    assert len(load_genome(compressed).references) == 3

    fasta = tmp_path / "genome.fna"
    fasta.write_text(">seq1 test\n" + "ACGT" * 400 + "\n")
    loaded = load_genome(fasta, topology="circular")
    assert loaded.file_format == "fasta"
    assert not loaded.has_annotation
    assert loaded.genome.topology == "circular"


def test_format_detection_by_content(tmp_path):
    unknown = tmp_path / "genome.txt"
    unknown.write_text(">x\nACGT\n")
    assert detect_format(unknown) == "fasta"
    unknown.write_text("pas un génome\n")
    with pytest.raises(InputFormatError):
        detect_format(unknown)


@pytest.mark.parametrize(
    ("content", "message"),
    [
        (">a\n" + "ACGT" * 300 + "\n>b\n" + "ACGT" * 300 + "\n", "2 séquences"),
        (">a\n" + "ACGU" * 300 + "\n", "non ADN"),
        (">a\n" + "ACGT" * 10 + "\n", "trop courte"),
        (">a\n" + "N" * 1500 + "\n", "ambiguës"),
    ],
)
def test_invalid_genomes_are_rejected(tmp_path, content, message):
    path = tmp_path / "bad.fasta"
    path.write_text(content)
    with pytest.raises(InputFormatError, match=message):
        load_genome(path)


def test_missing_file():
    with pytest.raises(InputFormatError, match="introuvable"):
        load_genome(Path("nope.gb"))


# --------------------------------------------------------------------------- NCBI (hors ligne)
def test_accession_validation():
    assert validate_accession(" nc_001416.1 ") == "NC_001416.1"
    with pytest.raises(Exception, match="invalide"):
        validate_accession("pas-un-numéro")


def test_deterministic_gzip(tmp_path):
    first, second = tmp_path / "a.gz", tmp_path / "b.gz"
    write_gzip_deterministic("LOCUS test\n", first)
    write_gzip_deterministic("LOCUS test\n", second)
    assert first.read_bytes() == second.read_bytes()


# --------------------------------------------------------------------------- écriture
def _gene(sequence: str, interval: Interval) -> PredictedGene:
    return PredictedGene(interval, sequence, 12.5, 10.0, "AGGAGG", 7, 2)


def test_gff3_and_protein_outputs():
    genome = Genome("G1", "test", "A" * 1000, "circular")
    genes = [
        _gene("GTG" + "GCA" * 29 + "TAA", Interval(10, 103, 1)),
        _gene("ATG" + "AAA" * 29 + "TGA", Interval(980, 1073, -1)),
    ]
    gff = gff3_text(genes, genome, "run").splitlines()
    assert gff[0] == "##gff-version 3"
    assert "Is_circular=true" in gff[2]
    fields = gff[3].split("\t")
    assert fields[2:5] == ["CDS", "11", "103"]
    assert gff[4].split("\t")[4] == "1073"  # fin > N autorisée pour un génome circulaire

    proteins = proteins_fasta_text(genes, genome, "run", table_id=11)
    records = proteins.split(">")[1:]
    assert records[0].splitlines()[1].startswith("MAAA")  # GTG initiateur traduit en M
    assert records[1].splitlines()[1] == "M" + "K" * 29


def test_write_json_converts_special_values(tmp_path):
    import numpy as np

    path = tmp_path / "out.json"
    write_json({"nan": math.nan, "int": np.int64(3), "path": Path("x")}, path)
    assert json.loads(path.read_text()) == {"nan": None, "int": 3, "path": "x"}


# --------------------------------------------------------------------------- configuration
def test_demo_config_loads_and_resolves_paths():
    config = load_config(ROOT / "config" / "demo.yaml")
    assert [run.name for run in config.runs][:2] == ["lambda", "mgenitalium"]
    assert all(run.genome.is_absolute() and run.genome.is_file() for run in config.runs)
    rbs_run = next(run for run in config.runs if run.name == "mgenitalium_rbs")
    assert rbs_run.resolve_predictor(config.predictor).start_strategy == "rbs"


def test_default_config_is_valid():
    config = load_config(ROOT / "config" / "default.yaml", check_files=False)
    assert config.predictor == PredictorConfig()


def _write_config(tmp_path: Path, data: dict) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


@pytest.mark.parametrize(
    "data",
    [
        {"runs": []},
        {"runs": [{"name": "a", "genome": "x.gb"}, {"name": "a", "genome": "x.gb"}]},
        {"runs": [{"name": "a", "genome": "x.gb", "translation_table": 99}]},
        {"runs": [{"name": "a", "genome": "x.gb"}], "predictor": {"start_codons": ["AXG"]}},
        {"runs": [{"name": "a", "genome": "x.gb"}], "predictor": {"unknown": 1}},
        {"runs": [{"name": "a", "genome": "x.gb", "predictor": {"max_overlap": -5}}]},
    ],
)
def test_invalid_configurations(tmp_path, data):
    with pytest.raises(ConfigError):
        load_config(_write_config(tmp_path, data), check_files=False)


def test_missing_genome_is_reported(tmp_path):
    path = _write_config(tmp_path, {"runs": [{"name": "a", "genome": "absent.gb"}]})
    with pytest.raises(ConfigError, match="introuvable"):
        load_config(path)
