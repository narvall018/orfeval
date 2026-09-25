"""Orchestration des analyses.

Deux niveaux, pour séparer calcul et effets de bord :

- :func:`analyze_genome` : calcul pur (aucune écriture), utilisé par le dashboard ;
- :func:`write_run` / :func:`run_project` : écriture des tableaux, figures, JSON et du
  rapport, utilisées par la CLI.
"""

from __future__ import annotations

import platform
import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from orfeval import __version__
from orfeval.config import AnalysisConfig, PredictorConfig, ProjectConfig, RunConfig
from orfeval.evaluation.matching import match_genes
from orfeval.evaluation.metrics import EvaluationResult, evaluate
from orfeval.features.codon_usage import compare_rscu, gc3_values
from orfeval.features.composition import CompositionResult, analyse_composition
from orfeval.features.upstream import UpstreamProfile, upstream_profile
from orfeval.io.genome import LoadedGenome, load_genome
from orfeval.io.writers import (
    gff3_text,
    predictions_dataframe,
    proteins_fasta_text,
    write_json,
    write_text,
    write_tsv,
)
from orfeval.logging_utils import get_logger
from orfeval.models import extract_sequence, upstream_sequence
from orfeval.orfs.predictor import GenePredictor, PredictionResult
from orfeval.plotting import (
    EVALUATION_LEGEND,
    comparison_figure,
    error_breakdown_figure,
    gc_skew_figure,
    genome_map,
    orf_length_null_figure,
    precision_recall_figure,
    recall_by_length_figure,
    rscu_figure,
    save_figure,
    score_length_figure,
    upstream_figure,
)
from orfeval.plotting.style import COLOR_CORRECT, COLOR_ERROR, COLOR_OTHER_START, COLOR_REFERENCE

logger = get_logger(__name__)

DEFAULT_TABLE = 11
TRACKED_PACKAGES = ("biopython", "numpy", "pandas", "matplotlib", "pydantic", "jinja2", "typer")


@dataclass(frozen=True, slots=True)
class RunResult:
    """Résultats complets (en mémoire) d'une analyse."""

    run: RunConfig
    predictor_config: PredictorConfig
    loaded: LoadedGenome
    table_id: int
    table_source: str
    prediction: PredictionResult
    evaluation: EvaluationResult | None
    composition: CompositionResult
    rscu: pd.DataFrame | None
    gc3_predicted: np.ndarray
    upstream_predicted: UpstreamProfile | None
    upstream_reference: UpstreamProfile | None
    elapsed_seconds: float


@dataclass(frozen=True, slots=True)
class ProjectOutcome:
    """Bilan d'une exécution de projet."""

    results: list[RunResult]
    output_dir: Path
    summary_table: pd.DataFrame
    report_path: Path | None


def environment_versions() -> dict[str, str]:
    """Versions de Python, du système et des dépendances principales."""
    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "orfeval": __version__,
    }
    for package in TRACKED_PACKAGES:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "absent"
    return versions


def resolve_table(run: RunConfig, loaded: LoadedGenome) -> tuple[int, str]:
    """Code génétique à utiliser et sa provenance (configuration > annotation > défaut)."""
    if run.translation_table is not None:
        return run.translation_table, "configuration"
    if loaded.genome.declared_table is not None:
        return loaded.genome.declared_table, "annotation (/transl_table)"
    logger.warning("Code génétique non déclaré : code %d utilisé par défaut.", DEFAULT_TABLE)
    return DEFAULT_TABLE, f"défaut ({DEFAULT_TABLE})"


def _upstream_sequences(loaded: LoadedGenome, intervals: list[Any], length: int) -> list[str]:
    genome = loaded.genome
    sequences = (
        upstream_sequence(genome.sequence, interval, length, circular=genome.is_circular)
        for interval in intervals
    )
    return [sequence for sequence in sequences if sequence is not None]


def analyze_genome(
    loaded: LoadedGenome,
    run: RunConfig,
    predictor_config: PredictorConfig,
    analysis_config: AnalysisConfig,
) -> RunResult:
    """Prédit, évalue et décrit un génome déjà chargé (aucune écriture sur disque)."""
    started = time.perf_counter()
    table_id, table_source = resolve_table(run, loaded)
    prediction = GenePredictor(predictor_config, table_id).predict(loaded.genome)
    evaluation = (
        evaluate(prediction, loaded, analysis_config.length_bins) if loaded.has_annotation else None
    )
    window = None if analysis_config.gc_window == "auto" else int(analysis_config.gc_window)
    composition = analyse_composition(loaded.genome, window, loaded.references)

    predicted_sequences = [gene.sequence for gene in prediction.genes]
    rscu = None
    upstream_reference = None
    if loaded.has_annotation:
        reference_sequences = [
            extract_sequence(loaded.genome.sequence, ref.interval)
            for ref in loaded.references
            if not ref.compound
        ]
        rscu = compare_rscu(predicted_sequences, reference_sequences, table_id)
        upstream_reference = upstream_profile(
            _upstream_sequences(
                loaded, [ref.interval for ref in loaded.references], analysis_config.upstream_length
            ),
            loaded.genome.sequence,
        )
    upstream_predicted = upstream_profile(
        _upstream_sequences(
            loaded, [gene.interval for gene in prediction.genes], analysis_config.upstream_length
        ),
        loaded.genome.sequence,
    )
    return RunResult(
        run=run,
        predictor_config=predictor_config,
        loaded=loaded,
        table_id=table_id,
        table_source=table_source,
        prediction=prediction,
        evaluation=evaluation,
        composition=composition,
        rscu=rscu,
        gc3_predicted=gc3_values(predicted_sequences),
        upstream_predicted=upstream_predicted,
        upstream_reference=upstream_reference,
        elapsed_seconds=time.perf_counter() - started,
    )


def prediction_table(result: RunResult) -> pd.DataFrame:
    """Tableau des gènes prédits, enrichi du statut d'évaluation s'il existe."""
    frame = predictions_dataframe(result.prediction.genes, result.loaded.genome, result.run.name)
    frame["start_adjusted"] = [gene.start_adjusted for gene in result.prediction.genes]
    frame["gc3"] = np.round(result.gc3_predicted, 4)
    if result.evaluation is not None and result.evaluation.matching is not None:
        frame["status"] = result.evaluation.prediction_status
        frame["reference_locus_tags"] = result.evaluation.matching.pred_reference_tags
    return frame


# ---------------------------------------------------------------------------- figures
def build_figures(result: RunResult) -> dict[str, Figure]:
    """Construit toutes les figures d'une analyse (clé → figure)."""
    genome = result.loaded.genome
    prediction = result.prediction
    evaluation = result.evaluation
    name = result.run.name
    figures: dict[str, Figure] = {}

    intervals = [gene.interval for gene in prediction.genes]
    if evaluation is not None and evaluation.matching is not None:
        matching = evaluation.matching
        prediction_colors = [
            COLOR_CORRECT if exact else COLOR_OTHER_START if correct else COLOR_ERROR
            for correct, exact in zip(
                matching.pred_correct, matching.pred_start_correct, strict=True
            )
        ]
        figures["genome_map"] = genome_map(
            genome.length,
            intervals,
            prediction_colors,
            [ref.interval for ref in result.loaded.references],
            [COLOR_REFERENCE if found else COLOR_ERROR for found in matching.ref_found],
            legend=EVALUATION_LEGEND,
            title=f"{name} — gènes prédits et gènes de référence",
        )
    else:
        figures["genome_map"] = genome_map(
            genome.length,
            intervals,
            [COLOR_CORRECT] * len(intervals),
            title=f"{name} — gènes prédits",
        )

    figures["gc_skew"] = gc_skew_figure(
        result.composition.profile,
        result.composition.landmarks,
        title=f"{name} — composition en GC et GC skew",
    )
    figures["orf_length_null"] = orf_length_null_figure(
        prediction.region_lengths,
        prediction.null_model,
        prediction.train_min_length,
        title=f"{name} — longueur des ORF : génome vs hasard",
    )
    figures["upstream"] = upstream_figure(
        result.upstream_predicted,
        result.upstream_reference,
        title=f"{name} — signal nucléotidique en amont des starts",
    )

    if evaluation is not None:
        figures["precision_recall"] = precision_recall_figure(
            evaluation.pr_model,
            evaluation.pr_length,
            evaluation.n_predicted,
            title=f"{name} — précision-rappel",
        )
        figures["recall_by_length"] = recall_by_length_figure(
            evaluation.by_length, title=f"{name} — rappel selon la longueur des gènes"
        )
        ranked_match = match_genes(
            [gene.interval for gene in prediction.ranked],
            result.loaded.references,
            genome.length,
        )
        figures["score_length"] = score_length_figure(
            [gene.length for gene in prediction.ranked],
            [gene.score for gene in prediction.ranked],
            list(ranked_match.pred_correct),
            prediction.config.score_threshold,
            title=f"{name} — score de détection des candidats",
        )
        figures["errors"] = error_breakdown_figure(
            evaluation.fp_categories,
            evaluation.fn_categories,
            title=f"{name} — origine des écarts avec l'annotation",
        )
    if result.rscu is not None:
        figures["codon_usage"] = rscu_figure(result.rscu, title=f"{name} — usage des codons (RSCU)")
    return figures


# ---------------------------------------------------------------------------- écriture
def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        return path.name


def run_summary(result: RunResult, figure_files: dict[str, str]) -> dict[str, Any]:
    """Résumé sérialisable d'une analyse (entrée unique du rapport)."""
    loaded = result.loaded
    genome = loaded.genome
    prediction = result.prediction
    genes = prediction.genes
    coding_positions = sum(gene.length for gene in genes)
    rscu_r = None
    if result.rscu is not None:
        data = result.rscu.dropna(subset=["rscu_predicted", "rscu_reference"])
        if len(data) > 2:
            rscu_r = float(np.corrcoef(data["rscu_reference"], data["rscu_predicted"])[0, 1])

    def upstream_info(profile: UpstreamProfile | None) -> dict[str, Any] | None:
        if profile is None:
            return None
        return {
            "n_sequences": profile.n_sequences,
            "max_bits": float(profile.information.max()),
            "max_position": int(profile.information.index[int(profile.information.argmax())]),
            "consensus": profile.consensus,
        }

    landmarks = result.composition.landmarks
    return {
        "name": result.run.name,
        "description": result.run.description,
        "input": {
            "path": _display_path(loaded.source),
            "sha256": loaded.sha256,
            "format": loaded.file_format,
            "seq_id": genome.seq_id,
            "description": genome.description,
            "organism": genome.organism,
            "length": genome.length,
            "topology": genome.topology,
            "n_ambiguous": loaded.n_ambiguous,
        },
        "annotation": asdict(loaded.annotation) if loaded.annotation else None,
        "genetic_code": {
            "table": result.table_id,
            "source": result.table_source,
            "declared": genome.declared_table,
        },
        "parameters": result.predictor_config.model_dump(mode="json"),
        "prediction": {
            "n_orfs": prediction.n_orfs,
            "n_genes": len(genes),
            "p_stop": prediction.null_model.p_stop,
            "train_min_length": prediction.train_min_length,
            "training_rounds": [asdict(item) for item in prediction.training_rounds],
            "coding_fraction": coding_positions / genome.length,
            "mean_length": float(np.mean([gene.length for gene in genes])) if genes else None,
            "start_codons": dict(Counter(gene.start_codon for gene in genes).most_common()),
            "n_start_adjusted": sum(gene.start_adjusted for gene in genes),
            "n_with_rbs_motif": sum(gene.rbs_motif is not None for gene in genes),
        },
        "evaluation": result.evaluation.summary() if result.evaluation else None,
        "composition": {
            "gc_content": result.composition.gc_content,
            "window": result.composition.window,
            "landmarks": asdict(landmarks) if landmarks else None,
            "gc3_mean_predicted": float(np.nanmean(result.gc3_predicted))
            if len(result.gc3_predicted)
            else None,
        },
        "upstream": {
            "predicted": upstream_info(result.upstream_predicted),
            "reference": upstream_info(result.upstream_reference),
        },
        "rscu_correlation": rscu_r,
        "figures": figure_files,
        "files": {
            "predictions": "predictions.tsv",
            "gff3": "predictions.gff3",
            "proteins": "proteins.faa",
            "reference_comparison": "reference_comparison.tsv" if result.evaluation else None,
            "recall_by_length": "recall_by_length.tsv" if result.evaluation else None,
        },
        "elapsed_seconds": result.elapsed_seconds,
    }


def write_run(
    result: RunResult, run_dir: Path, *, dpi: int = 130, figures: bool = True
) -> dict[str, Any]:
    """Écrit tous les fichiers d'une analyse et renvoie son résumé.

    Avec ``figures=False``, seuls les tableaux et le résumé sont écrits (plus rapide pour
    les grands génomes).
    """
    genome = result.loaded.genome
    genes = result.prediction.genes
    name = result.run.name
    run_dir.mkdir(parents=True, exist_ok=True)

    write_tsv(prediction_table(result), run_dir / "predictions.tsv")
    write_text(gff3_text(genes, genome, name), run_dir / "predictions.gff3")
    write_text(proteins_fasta_text(genes, genome, name, result.table_id), run_dir / "proteins.faa")
    if result.evaluation is not None:
        write_tsv(result.evaluation.reference_table, run_dir / "reference_comparison.tsv")
        write_tsv(result.evaluation.by_length, run_dir / "recall_by_length.tsv")
        write_tsv(result.evaluation.pr_model, run_dir / "precision_recall_model.tsv")
        write_tsv(result.evaluation.pr_length, run_dir / "precision_recall_length.tsv")
    write_tsv(result.composition.profile, run_dir / "gc_profile.tsv")
    if result.rscu is not None:
        write_tsv(result.rscu, run_dir / "codon_usage.tsv")

    figure_files: dict[str, str] = {}
    for key, figure in (build_figures(result) if figures else {}).items():
        relative = f"figures/{key}.png"
        save_figure(figure, run_dir / relative, dpi=dpi)
        figure_files[key] = relative

    summary = run_summary(result, figure_files)
    write_json(summary, run_dir / "run_summary.json")
    return summary


SUMMARY_COLUMNS = (
    "run",
    "genome",
    "table",
    "n_genes",
    "sensitivity",
    "precision",
    "f1",
    "start_accuracy",
    "ap_model",
    "ap_length",
)


def summary_row(summary: dict[str, Any]) -> dict[str, Any]:
    """Ligne du tableau récapitulatif à partir d'un résumé d'analyse."""
    evaluation = summary.get("evaluation") or {}

    def metric(key: str) -> float | None:
        value = evaluation.get(key)
        return value.get("value") if isinstance(value, dict) else value

    return {
        "run": summary["name"],
        "genome": summary["input"]["seq_id"],
        "table": summary["genetic_code"]["table"],
        "n_genes": summary["prediction"]["n_genes"],
        "sensitivity": metric("sensitivity"),
        "precision": metric("precision"),
        "f1": metric("f1"),
        "start_accuracy": metric("start_accuracy"),
        "ap_model": metric("average_precision_model"),
        "ap_length": metric("average_precision_length"),
    }


def run_project(
    config: ProjectConfig,
    *,
    output_dir: Path | None = None,
    config_path: Path | None = None,
    make_report: bool = True,
) -> ProjectOutcome:
    """Exécute toutes les analyses d'une configuration et écrit les résultats."""
    output = (output_dir or config.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = datetime.now(UTC)
    results: list[RunResult] = []
    summaries: list[dict[str, Any]] = []
    for run in config.runs:
        logger.info("=== Analyse « %s » ===", run.name)
        loaded = load_genome(run.genome, topology=run.topology)
        predictor_config = run.resolve_predictor(config.predictor)
        result = analyze_genome(loaded, run, predictor_config, config.analysis)
        summaries.append(write_run(result, output / run.name, dpi=config.analysis.figure_dpi))
        results.append(result)
        if result.evaluation is not None:
            logger.info(
                "%s : %d gènes prédits | sensibilité %s | précision %s",
                run.name,
                len(result.prediction.genes),
                result.evaluation.sensitivity,
                result.evaluation.precision,
            )

    table = pd.DataFrame([summary_row(summary) for summary in summaries], columns=SUMMARY_COLUMNS)
    write_tsv(table.round(4), output / "summary.tsv")
    evaluated = [summary for summary in summaries if summary.get("evaluation")]
    if evaluated:
        figure = comparison_figure(
            [{"name": s["name"], **s["evaluation"]} for s in evaluated],
            title="Comparaison des analyses (IC de Wilson à 95 %)",
        )
        save_figure(figure, output / "comparison.png", dpi=config.analysis.figure_dpi)

    write_json(
        {
            "created_utc": started,
            "finished_utc": datetime.now(UTC),
            "command": " ".join(
                Path(sys.argv[0]).name if i == 0 else arg for i, arg in enumerate(sys.argv)
            ),
            "config_file": _display_path(config_path) if config_path else None,
            "runs": [run.name for run in config.runs],
            "predictor_defaults": config.predictor.model_dump(mode="json"),
            "analysis": config.analysis.model_dump(mode="json"),
            "versions": environment_versions(),
            "comparison_figure": "comparison.png" if evaluated else None,
        },
        output / "project_summary.json",
    )

    report_path = None
    if make_report:
        from orfeval.report.builder import build_report

        report_path = build_report(output)
    return ProjectOutcome(results, output, table, report_path)
