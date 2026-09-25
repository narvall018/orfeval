"""Construction du rapport HTML autonome.

Le rapport est construit **uniquement à partir des fichiers de résultats**
(``project_summary.json``, ``<analyse>/run_summary.json`` et figures PNG) : il peut
être régénéré sans relancer les calculs, et il ne contient que des valeurs réellement
calculées. Les images sont intégrées en base64 : le fichier HTML se suffit à lui-même.
"""

from __future__ import annotations

import base64
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from jinja2 import Environment, PackageLoader, StrictUndefined, select_autoescape

from orfeval.exceptions import OrfevalError
from orfeval.logging_utils import get_logger

logger = get_logger(__name__)

FIGURE_INFO: dict[str, tuple[str, str]] = {
    "genome_map": (
        "Carte des gènes",
        "Chaque ligne couvre une portion du génome ; pistes : référence (+), prédiction (+), "
        "prédiction (−), référence (−). Les erreurs (orange) montrent où se concentrent les "
        "écarts : régions denses en gènes chevauchants, petits gènes, pseudogènes.",
    ),
    "precision_recall": (
        "Courbes précision-rappel",
        "Chaque point correspond à un seuil de score. La courbe bleue classe les candidats "
        "par score codant, l'orange par longueur seule (référence naïve). Si la courbe bleue "
        "domine, le modèle codant apporte une information que la longueur ne donne pas.",
    ),
    "recall_by_length": (
        "Rappel selon la longueur",
        "Proportion de gènes annotés retrouvés par classe de longueur, avec IC de Wilson à "
        "95 %. Les gènes courts sont les plus difficiles : peu de codons, donc peu de signal "
        "statistique. Les classes à faible effectif ont des intervalles larges.",
    ),
    "score_length": (
        "Score de détection des candidats",
        "Tous les candidats retenus après résolution des chevauchements, y compris sous le "
        "seuil. Le score (log-rapport de vraisemblance, en bits) croît avec la longueur des "
        "vrais gènes ; les candidats sans correspondance se concentrent aux faibles scores.",
    ),
    "errors": (
        "Origine des écarts avec l'annotation",
        "Gauche : ce que recouvrent les prédictions sans correspondance (une prédiction sur "
        "un pseudogène n'est pas forcément absurde). Droite : pourquoi des gènes annotés "
        "n'ont pas été retrouvés.",
    ),
    "orf_length_null": (
        "Longueur des ORF : génome réel vs hasard",
        "Nombre de régions sans codon stop d'au moins L nucléotides (échelle log). L'écart "
        "entre le génome et le modèle aléatoire de même composition traduit la présence de "
        "gènes. Le trait vertical marque le seuil d'auto-apprentissage choisi sans annotation.",
    ),
    "gc_skew": (
        "Composition et GC skew",
        "GC skew par fenêtre et cumulé. Pour un chromosome bactérien circulaire, le minimum "
        "du skew cumulé est un indicateur classique de l'origine de réplication. "
        "Indicateur sans objet pour un génome linéaire de phage.",
    ),
    "codon_usage": (
        "Usage des codons (RSCU)",
        "Un point par codon sens : RSCU des gènes annotés (x) et des gènes prédits (y). Des "
        "points proches de la diagonale indiquent que l'ensemble prédit a la même signature "
        "de codons que l'annotation.",
    ),
    "upstream": (
        "Signal en amont des codons start",
        "Écart de composition à celle du génome (entropie relative, bits) pour chaque "
        "position avant le start. Un pic vers −12 à −6 signale un motif de Shine-Dalgarno ; "
        "un profil plat indique un signal faible ou des starts incorrects.",
    ),
}
FIGURE_ORDER = tuple(FIGURE_INFO)


def _image_data_uri(path: Path) -> str | None:
    if not path.is_file():
        logger.warning("Figure absente : %s", path)
        return None
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _load_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise OrfevalError(f"lecture impossible de {path} : {exc}") from exc
    if not isinstance(data, dict):
        raise OrfevalError(f"{path} : dictionnaire JSON attendu")
    return data


def _percent(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "n/a"
    return f"{100 * value:.{digits}f} %".replace(".", ",")


def _metric(value: dict[str, Any] | None) -> str:
    if not value or value.get("value") is None:
        return "n/a"
    low, high = value.get("ci95_low"), value.get("ci95_high")
    interval = f" [{_percent(low)} – {_percent(high)}]" if low is not None else ""
    return f"{_percent(value['value'])}{interval} ({value['successes']}/{value['total']})"


def _number(value: float | int | None, digits: int = 3) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int):
        return f"{value:,}".replace(",", " ")
    return f"{value:.{digits}f}".replace(".", ",")


def _parameter_overrides(run: dict[str, Any], defaults: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in run["parameters"].items() if defaults.get(key) != value}


def _facts(run: dict[str, Any]) -> list[str]:
    """Constats statistiques factuels (aucune interprétation biologique)."""
    facts = []
    evaluation = run.get("evaluation")
    if evaluation:
        facts.append(
            "Classement par score codant : précision moyenne (AP) de "
            f"{_number(evaluation['average_precision_model'])} contre "
            f"{_number(evaluation['average_precision_length'])} pour la longueur seule."
        )
        baseline = evaluation.get("baseline_same_k") or {}
        if baseline.get("f1") is not None:
            facts.append(
                f"À nombre de prédictions égal ({baseline['n_predictions']}), le classement par "
                f"longueur seule obtient un F1 de {_number(baseline['f1'])} contre "
                f"{_number(evaluation['f1'])} pour orfeval."
            )
        facts.append(
            "Au niveau nucléotidique : "
            f"{_percent(evaluation['nucleotide_sensitivity'])} des bases codantes annotées "
            f"sont couvertes par une prédiction sur le bon brin, et "
            f"{_percent(evaluation['nucleotide_precision'])} des bases prédites sont annotées "
            "comme codantes."
        )
    upstream = run["upstream"].get("predicted")
    if upstream:
        facts.append(
            f"Signal amont maximal (starts prédits) : {_number(upstream['max_bits'], 2)} bit "
            f"à la position {upstream['max_position']}."
        )
    landmarks = (run.get("composition") or {}).get("landmarks")
    if landmarks:
        text = (
            "Minimum du GC skew cumulé vers "
            f"{_number(landmarks['origin_candidate'] // 1000)} kb, maximum vers "
            f"{_number(landmarks['terminus_candidate'] // 1000)} kb."
        )
        if landmarks.get("dnaa_distance") is not None:
            text += (
                f" Le gène dnaA annoté se trouve à {_number(landmarks['dnaa_distance'])} nt "
                "du minimum."
            )
        facts.append(text)
    if run.get("rscu_correlation") is not None:
        facts.append(
            "Corrélation (Pearson) entre RSCU des gènes prédits et annotés : "
            f"{_number(run['rscu_correlation'])}."
        )
    return facts


def load_results(results_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Charge le résumé de projet et les résumés d'analyses d'un dossier de résultats."""
    project_file = results_dir / "project_summary.json"
    if not project_file.is_file():
        raise OrfevalError(
            f"{project_file} introuvable : lancez d'abord 'orfeval analyze' ou 'orfeval demo'."
        )
    project = _load_json(project_file)
    runs = []
    for name in project.get("runs", []):
        summary_file = results_dir / name / "run_summary.json"
        if not summary_file.is_file():
            raise OrfevalError(f"résumé d'analyse manquant : {summary_file}")
        runs.append(_load_json(summary_file))
    if not runs:
        raise OrfevalError(f"aucune analyse trouvée dans {results_dir}")
    return project, runs


def build_report(results_dir: Path, output: Path | None = None) -> Path:
    """Génère ``report.html`` à partir d'un dossier de résultats.

    Returns
    -------
    Path
        Chemin du rapport écrit.
    """
    results_dir = results_dir.resolve()
    project, runs = load_results(results_dir)
    defaults = project.get("predictor_defaults", {})
    for run in runs:
        base = results_dir / run["name"]
        run["figure_blocks"] = [
            {
                "key": key,
                "title": FIGURE_INFO[key][0],
                "caption": FIGURE_INFO[key][1],
                "src": _image_data_uri(base / run["figures"][key]),
            }
            for key in FIGURE_ORDER
            if key in run.get("figures", {})
        ]
        run["overrides"] = _parameter_overrides(run, defaults)
        run["facts"] = _facts(run)

    comparison = project.get("comparison_figure")
    environment = Environment(
        loader=PackageLoader("orfeval.report", "templates"),
        autoescape=select_autoescape(["html", "j2"]),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    environment.filters.update(percent=_percent, metric=_metric, number=_number)
    template = environment.get_template("report.html.j2")
    html = template.render(
        project=project,
        runs=runs,
        defaults=defaults,
        comparison_src=_image_data_uri(results_dir / comparison) if comparison else None,
        generated=datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"),
    )
    output = output or results_dir / "report.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    logger.info("Rapport écrit : %s", output)
    return output
