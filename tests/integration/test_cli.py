"""Interface en ligne de commande, exécutée comme un utilisateur le ferait."""

import json

import pytest
from typer.testing import CliRunner

from orfeval import __version__
from orfeval.cli import app

pytestmark = pytest.mark.integration
runner = CliRunner()


def test_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_validate_demo_genome(mgen_path):
    result = runner.invoke(app, ["validate", str(mgen_path)])
    assert result.exit_code == 0, result.output
    assert "580" in result.output
    assert "Fichier valide" in result.output


def test_validate_reports_errors_cleanly(tmp_path):
    bad = tmp_path / "bad.fasta"
    bad.write_text(">x\nACGU\n")
    result = runner.invoke(app, ["validate", str(bad)])
    assert result.exit_code == 1
    assert "Traceback" not in result.output


def test_predict_with_forced_genetic_code(lambda_path, tmp_path):
    outdir = tmp_path / "pred"
    result = runner.invoke(
        app,
        ["-q", "predict", str(lambda_path), "-o", str(outdir), "--table", "11", "--name", "lam"],
    )
    assert result.exit_code == 0, result.output
    summary = json.loads((outdir / "run_summary.json").read_text())
    assert summary["genetic_code"]["source"] == "configuration"
    assert (outdir / "predictions.gff3").read_text().count("lam_00001") == 1


def test_predict_rejects_invalid_parameter(lambda_path, tmp_path):
    result = runner.invoke(
        app, ["predict", str(lambda_path), "-o", str(tmp_path), "--start-strategy", "magie"]
    )
    assert result.exit_code == 1


def test_demo_runs_end_to_end(tmp_path):
    """La démonstration complète, telle que documentée dans le README."""
    result = runner.invoke(app, ["-q", "demo", "--output", str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "report.html").is_file()
    assert (tmp_path / "comparison.png").is_file()
    summary = (tmp_path / "summary.tsv").read_text().splitlines()
    assert len(summary) == 1 + 4
    rebuilt = runner.invoke(app, ["report", str(tmp_path), "-o", str(tmp_path / "copie.html")])
    assert rebuilt.exit_code == 0
    assert (tmp_path / "copie.html").is_file()
