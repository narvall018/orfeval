"""Téléchargement NCBI (simulé, sans réseau) et cas limites du prédicteur."""

import io
import logging

import pytest

from orfeval.config import PredictorConfig
from orfeval.exceptions import AnalysisError, DataFetchError
from orfeval.io import ncbi
from orfeval.io.genome import load_genome
from orfeval.orfs.predictor import GenePredictor


class _FakeHandle(io.StringIO):
    """Imite le handle renvoyé par Bio.Entrez.efetch."""


def test_fetch_writes_a_reproducible_gzip(monkeypatch, tmp_path, toy_genbank):
    text = toy_genbank.read_text()
    calls = []

    def fake_efetch(**kwargs):
        calls.append(kwargs)
        return _FakeHandle(text)

    monkeypatch.setattr(ncbi.Entrez, "efetch", fake_efetch)
    path = ncbi.fetch_genbank("NC_999999.1", tmp_path, email="test@example.org")
    assert path.name == "NC_999999.1.gb.gz"
    assert calls[0]["rettype"] == "gbwithparts"
    assert len(load_genome(path).references) == 3

    first_bytes = path.read_bytes()
    assert ncbi.fetch_genbank("NC_999999.1", tmp_path) == path  # déjà présent : pas d'appel
    assert len(calls) == 1
    ncbi.fetch_genbank("NC_999999.1", tmp_path, overwrite=True)
    assert path.read_bytes() == first_bytes


def test_fetch_rejects_unexpected_answers(monkeypatch, tmp_path):
    monkeypatch.setattr(ncbi.Entrez, "efetch", lambda **_: _FakeHandle("Error: bad id"))
    with pytest.raises(DataFetchError, match="inattendue"):
        ncbi.fetch_genbank("NC_000001.1", tmp_path)

    def failing(**_):
        raise OSError("réseau indisponible")

    monkeypatch.setattr(ncbi.Entrez, "efetch", failing)
    with pytest.raises(DataFetchError, match="échec"):
        ncbi.fetch_genbank("NC_000001.1", tmp_path)


def test_training_falls_back_to_longest_orfs(synthetic, caplog):
    config = PredictorConfig(train_min_length=20_000)
    with caplog.at_level(logging.WARNING, logger="orfeval"):
        result = GenePredictor(config, 11).predict(synthetic.genome)
    assert result.training_rounds[0].n_training_genes == config.min_training_genes
    assert result.train_min_length < 20_000
    assert "plus longs" in caplog.text


def test_too_few_orfs_is_an_explicit_error(synthetic):
    config = PredictorConfig(min_length=6000)
    with pytest.raises(AnalysisError, match="ORF"):
        GenePredictor(config, 11).predict(synthetic.genome)


def test_invalid_start_codon_for_table_is_reported():
    with pytest.raises(AnalysisError, match="non autorisé"):
        GenePredictor(PredictorConfig(start_codons=("CTG", "AAA")), 11)
