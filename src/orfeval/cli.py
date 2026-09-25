"""Interface en ligne de commande (Typer).

La CLI ne contient aucune logique scientifique : elle lit les arguments, appelle le
package et affiche les résultats. Les erreurs attendues (:class:`OrfevalError`) sont
affichées proprement avec un code de sortie non nul.
"""

from __future__ import annotations

import logging
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Annotated

import typer
from rich.console import Console
from rich.table import Table

from orfeval import __version__
from orfeval.config import PredictorConfig, ProjectConfig, RunConfig, load_config
from orfeval.exceptions import OrfevalError
from orfeval.logging_utils import setup_logging

if TYPE_CHECKING:
    import pandas as pd

app = typer.Typer(
    name="orfeval",
    help="Prédiction auto-entraînée de gènes procaryotes et évaluation face à une annotation.",
    no_args_is_help=True,
    add_completion=False,
    rich_markup_mode="rich",
)
console = Console()
error_console = Console(stderr=True)

DEMO_CONFIG = Path("config") / "demo.yaml"


@contextmanager
def _handle_errors() -> Iterator[None]:
    """Convertit les erreurs attendues en message lisible et code de sortie 1."""
    try:
        yield
    except OrfevalError as exc:
        error_console.print(f"[bold red]Erreur :[/] {exc}")
        raise typer.Exit(code=1) from exc


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"orfeval {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Affiche les messages de débogage.")
    ] = False,
    quiet: Annotated[
        bool, typer.Option("--quiet", "-q", help="N'affiche que les avertissements.")
    ] = False,
    log_file: Annotated[Path | None, typer.Option(help="Copie le journal dans ce fichier.")] = None,
    version: Annotated[
        bool,
        typer.Option(
            "--version", callback=_version_callback, is_eager=True, help="Affiche la version."
        ),
    ] = False,
) -> None:
    """Trouve les gènes d'un génome procaryote et mesure la qualité du résultat."""
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    setup_logging(level, log_file)


# ---------------------------------------------------------------------------- fetch
@app.command()
def fetch(
    accessions: Annotated[list[str], typer.Argument(help="Numéros d'accession NCBI.")],
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Dossier de sortie.")] = Path(
        "data/raw"
    ),
    email: Annotated[
        str | None,
        typer.Option(envvar="NCBI_EMAIL", help="E-mail transmis au NCBI (recommandé)."),
    ] = None,
    api_key: Annotated[
        str | None, typer.Option(envvar="NCBI_API_KEY", help="Clé d'API NCBI (optionnelle).")
    ] = None,
    force: Annotated[bool, typer.Option(help="Retélécharge les fichiers existants.")] = False,
) -> None:
    """Télécharge des génomes GenBank complets depuis le NCBI (Bio.Entrez)."""
    from orfeval.io.genome import sha256sum
    from orfeval.io.ncbi import fetch_genbank

    with _handle_errors():
        for accession in accessions:
            path = fetch_genbank(accession, outdir, email=email, api_key=api_key, overwrite=force)
            console.print(f"[green]✓[/] {path}  sha256={sha256sum(path)}")


# ---------------------------------------------------------------------------- validate
@app.command()
def validate(
    genome: Annotated[Path, typer.Argument(help="Fichier GenBank ou FASTA (.gz accepté).")],
) -> None:
    """Vérifie qu'un génome est lisible et résume son contenu."""
    from orfeval.io.genome import load_genome

    with _handle_errors():
        loaded = load_genome(genome)
    data = loaded.genome
    table = Table(title=f"Validation de {genome.name}", show_header=False)
    table.add_row("Format", loaded.file_format)
    table.add_row("Séquence", f"{data.seq_id} — {data.description}")
    table.add_row("Organisme", data.organism or "non renseigné")
    table.add_row("Longueur", f"{data.length:,} nt".replace(",", " "))
    table.add_row("Topologie", data.topology)
    table.add_row("Bases ambiguës", str(loaded.n_ambiguous))
    table.add_row(
        "Code génétique déclaré", str(data.declared_table) if data.declared_table else "aucun"
    )
    if loaded.annotation is not None:
        annotation = loaded.annotation
        table.add_row("CDS annotées", str(annotation.n_cds))
        table.add_row(
            "CDS utilisables comme référence",
            f"{annotation.n_reference} (exclues : {annotation.n_pseudo} pseudogènes, "
            f"{annotation.n_partial} partielles, {annotation.n_unsupported} non prises en charge)",
        )
    else:
        table.add_row("Annotation", "aucune : prédiction possible, évaluation impossible")
    table.add_row("SHA-256", loaded.sha256)
    console.print(table)
    console.print("[green]✓ Fichier valide.[/]")


# ---------------------------------------------------------------------------- predict
@app.command()
def predict(
    genome: Annotated[Path, typer.Argument(help="Fichier GenBank ou FASTA (.gz accepté).")],
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Dossier de sortie.")] = Path(
        "results/predict"
    ),
    table: Annotated[
        int | None, typer.Option("--table", "-t", help="Code génétique NCBI (ex. 11, 4).")
    ] = None,
    min_length: Annotated[
        int | None, typer.Option(help="Longueur minimale des gènes (nt, stop inclus).")
    ] = None,
    start_strategy: Annotated[
        str | None, typer.Option(help="Choix du start : longest, rbs ou score.")
    ] = None,
    topology: Annotated[
        str | None, typer.Option(help="Force la topologie : linear ou circular.")
    ] = None,
    name: Annotated[str | None, typer.Option(help="Nom de l'analyse (préfixe des gènes).")] = None,
) -> None:
    """Prédit les gènes d'un génome et écrit TSV, GFF3, protéines et figures."""
    from orfeval.io.genome import load_genome
    from orfeval.pipeline import analyze_genome, write_run

    overrides = {
        key: value
        for key, value in {"min_length": min_length, "start_strategy": start_strategy}.items()
        if value is not None
    }
    with _handle_errors():
        try:
            predictor = PredictorConfig.model_validate(overrides)
            run = RunConfig(
                name=name or genome.name.split(".")[0],
                genome=genome,
                translation_table=table,
                topology=topology,
            )
        except ValueError as exc:
            raise OrfevalError(f"paramètre invalide :\n{exc}") from exc
        loaded = load_genome(genome, topology=run.topology)
        project = ProjectConfig(runs=[run])
        result = analyze_genome(loaded, run, predictor, project.analysis)
        write_run(result, outdir, dpi=project.analysis.figure_dpi)

    console.print(
        f"[green]✓[/] {len(result.prediction.genes)} gènes prédits "
        f"(code génétique {result.table_id}, {result.table_source}) → {outdir}"
    )
    if result.evaluation is not None:
        evaluation = result.evaluation
        console.print(f"  sensibilité : {evaluation.sensitivity}")
        console.print(f"  précision   : {evaluation.precision}")
        console.print(f"  start exact : {evaluation.start_accuracy}")


# ---------------------------------------------------------------------------- analyze
def _print_summary(summary: pd.DataFrame) -> None:
    from pandas import isna

    table = Table(title="Récapitulatif")
    for column in ("run", "table", "n_genes", "sensitivity", "precision", "f1", "start_accuracy"):
        table.add_column(column, justify="right" if column != "run" else "left")
    for _, row in summary.iterrows():
        cells = [str(row["run"]), str(row["table"]), str(row["n_genes"])]
        for column in ("sensitivity", "precision", "f1", "start_accuracy"):
            value = row[column]
            cells.append("n/a" if isna(value) else f"{value:.3f}")
        table.add_row(*cells)
    console.print(table)


def _run_config(config_path: Path, output: Path | None, report: bool) -> None:
    from orfeval.pipeline import run_project

    with _handle_errors():
        config = load_config(config_path)
        outcome = run_project(
            config, output_dir=output, config_path=config_path, make_report=report
        )
    _print_summary(outcome.summary_table)
    console.print(f"[green]✓[/] Résultats : {outcome.output_dir}")
    if outcome.report_path is not None:
        console.print(f"[green]✓[/] Rapport : {outcome.report_path}")


@app.command()
def analyze(
    config: Annotated[Path, typer.Argument(help="Fichier de configuration YAML.")],
    output: Annotated[
        Path | None, typer.Option("--output", "-o", help="Remplace output_dir de la configuration.")
    ] = None,
    report: Annotated[bool, typer.Option(help="Génère le rapport HTML.")] = True,
) -> None:
    """Exécute toutes les analyses décrites dans une configuration YAML."""
    _run_config(config, output, report)


@app.command()
def report(
    results: Annotated[Path, typer.Argument(help="Dossier de résultats d'une analyse.")],
    output: Annotated[Path | None, typer.Option("--output", "-o", help="Fichier HTML.")] = None,
) -> None:
    """Régénère le rapport HTML à partir d'un dossier de résultats (sans recalcul)."""
    from orfeval.report.builder import build_report

    with _handle_errors():
        path = build_report(results, output)
    console.print(f"[green]✓[/] Rapport : {path}")


def find_demo_config(start: Path | None = None) -> Path:
    """Cherche ``config/demo.yaml`` depuis le dossier courant, puis ses parents."""
    here = (start or Path.cwd()).resolve()
    for folder in (here, *here.parents):
        candidate = folder / DEMO_CONFIG
        if candidate.is_file():
            return candidate
    source_tree = Path(__file__).resolve().parents[2] / DEMO_CONFIG
    if source_tree.is_file():
        return source_tree
    raise OrfevalError(
        "config/demo.yaml introuvable : lancez la commande depuis le dépôt orfeval "
        "(ou utilisez 'orfeval analyze <config.yaml>')."
    )


@app.command()
def demo(
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="Dossier de sortie (défaut : results/demo)."),
    ] = None,
) -> None:
    """Lance la démonstration complète sur les génomes fournis dans data/demo."""
    with _handle_errors():
        config_path = find_demo_config()
    console.print(f"Configuration de démonstration : {config_path}")
    _run_config(config_path, output, report=True)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(app())
