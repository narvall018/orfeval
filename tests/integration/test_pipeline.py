"""Pipeline complet sur le phage lambda (données de démonstration réelles)."""

import json
import math
from pathlib import Path

import pandas as pd
import pytest
import yaml

from orfeval.config import load_config
from orfeval.pipeline import run_project
from orfeval.report.builder import build_report

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def lambda_results(tmp_path_factory, lambda_path):
    workdir = tmp_path_factory.mktemp("lambda")
    config_path = workdir / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "output_dir": "results",
                "runs": [{"name": "lambda", "genome": str(lambda_path)}],
            }
        )
    )
    outcome = run_project(load_config(config_path), config_path=config_path)
    return outcome


def test_expected_files_are_written(lambda_results):
    run_dir = lambda_results.output_dir / "lambda"
    for name in (
        "predictions.tsv",
        "predictions.gff3",
        "proteins.faa",
        "reference_comparison.tsv",
        "recall_by_length.tsv",
        "run_summary.json",
    ):
        assert (run_dir / name).stat().st_size > 0, name
    figures = {path.stem for path in (run_dir / "figures").glob("*.png")}
    assert figures >= {"genome_map", "precision_recall", "recall_by_length", "gc_skew", "upstream"}
    assert (lambda_results.output_dir / "summary.tsv").is_file()
    assert lambda_results.report_path is not None and lambda_results.report_path.is_file()


def test_outputs_are_consistent(lambda_results):
    run_dir = lambda_results.output_dir / "lambda"
    predictions = pd.read_csv(run_dir / "predictions.tsv", sep="\t")
    summary = json.loads((run_dir / "run_summary.json").read_text())
    proteins = (run_dir / "proteins.faa").read_text().count(">")
    gff_cds = sum(
        1 for line in (run_dir / "predictions.gff3").read_text().splitlines() if "\tCDS\t" in line
    )
    n_genes = summary["prediction"]["n_genes"]
    assert len(predictions) == proteins == gff_cds == n_genes
    assert (predictions["length_nt"] % 3 == 0).all()
    assert set(predictions["strand"]) <= {"+", "-"}
    assert summary["genetic_code"] == {
        "table": 11,
        "source": "annotation (/transl_table)",
        "declared": 11,
    }


def test_metrics_are_valid_proportions(lambda_results):
    summary = json.loads((lambda_results.output_dir / "lambda" / "run_summary.json").read_text())
    evaluation = summary["evaluation"]
    for key in ("sensitivity", "precision", "start_accuracy"):
        metric = evaluation[key]
        assert 0.0 <= metric["ci95_low"] <= metric["value"] <= metric["ci95_high"] <= 1.0
    assert evaluation["n_reference"] == 73
    assert (
        sum(evaluation["fp_categories"].values())
        == evaluation["n_predicted"] - evaluation["precision"]["successes"]
    )


def test_regression_floor_on_lambda(lambda_results):
    """Garde-fou de non-régression (valeurs mesurées en v0.1.0 : 0,78 / 0,81).

    Ce n'est pas une affirmation scientifique : seulement une alerte si une
    modification du code dégrade nettement les résultats sur la démonstration.
    """
    row = lambda_results.summary_table.iloc[0]
    assert row["sensitivity"] >= 0.70
    assert row["precision"] >= 0.70
    assert not math.isnan(row["f1"])


def test_report_can_be_rebuilt_from_files(lambda_results, tmp_path):
    output = build_report(lambda_results.output_dir, tmp_path / "rapport.html")
    html = output.read_text(encoding="utf-8")
    for section in ("Données", "Paramètres", "Méthodes", "Limites", "Reproductibilité"):
        assert section in html
    assert html.count("data:image/png;base64,") >= 8
    assert "{{" not in html and "}}" not in html


def test_report_requires_results(tmp_path):
    from orfeval.exceptions import OrfevalError

    with pytest.raises(OrfevalError, match="introuvable"):
        build_report(Path(tmp_path))
