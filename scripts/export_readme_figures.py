"""Copie une sélection de figures de la démo vers docs/images pour le README.

Les figures du README sont ainsi toujours issues d'une exécution réelle de
``orfeval demo`` (aucune figure n'est produite à la main).

Usage : ``python scripts/export_readme_figures.py [dossier_de_résultats]``
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELECTION = {
    "comparison.png": "comparison.png",
    "mgenitalium/figures/genome_map.png": "mgenitalium_genome_map.png",
    "lambda/figures/precision_recall.png": "lambda_precision_recall.png",
    "mgenitalium/figures/orf_length_null.png": "mgenitalium_orf_length_null.png",
    "lambda/figures/upstream.png": "lambda_upstream.png",
    "mgenitalium/figures/upstream.png": "mgenitalium_upstream.png",
    "mgenitalium_code11/figures/errors.png": "mgenitalium_code11_errors.png",
    "mgenitalium/figures/gc_skew.png": "mgenitalium_gc_skew.png",
    "lambda/figures/recall_by_length.png": "lambda_recall_by_length.png",
}


def main(results: Path) -> int:
    destination = ROOT / "docs" / "images"
    destination.mkdir(parents=True, exist_ok=True)
    missing = [source for source in SELECTION if not (results / source).is_file()]
    if missing:
        print(f"Figures absentes dans {results} : {missing}. Lancez d'abord 'orfeval demo'.")
        return 1
    for source, target in SELECTION.items():
        shutil.copyfile(results / source, destination / target)
        print(f"{source} -> docs/images/{target}")
    return 0


if __name__ == "__main__":
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "demo"
    sys.exit(main(folder))
