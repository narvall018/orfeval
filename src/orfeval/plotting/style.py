"""Style graphique commun (palette validée pour les daltonismes, axes discrets).

Les figures sont construites avec :class:`matplotlib.figure.Figure` (sans ``pyplot``) :
aucun état global, compatible avec un serveur (Streamlit) et les exécutions parallèles.

Palette : slots catégoriels ordonnés (bleu, orange, aqua...) dont les trois premiers
sont distinguables deux à deux, y compris en cas de deutéranopie ou de tritanopie.
Le rouge/vert n'est jamais utilisé seul pour opposer « correct » et « erroné ».
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.layout_engine import ConstrainedLayoutEngine

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
BLUE, ORANGE, AQUA = SERIES[0], SERIES[1], SERIES[2]

# Rôles sémantiques de l'évaluation (identité portée aussi par la légende et la forme).
COLOR_CORRECT = BLUE
COLOR_OTHER_START = AQUA
COLOR_ERROR = ORANGE
COLOR_REFERENCE = MUTED

LINE_WIDTH = 1.6
FONT_SIZE = 9.5


def new_figure(width: float, height: float) -> Figure:
    """Figure vide avec fond et mise en page homogènes."""
    figure = Figure(figsize=(width, height), facecolor=SURFACE, layout="constrained")
    engine = figure.get_layout_engine()
    if isinstance(engine, ConstrainedLayoutEngine):
        engine.set(w_pad=0.08, h_pad=0.08)
    return figure


def style_axes(ax: Axes, grid: Literal["x", "y", "both"] | None = "y") -> None:
    """Applique le style discret : axes fins, grille pleine très claire, texte gris."""
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=FONT_SIZE - 1, length=3)
    ax.xaxis.label.set_color(INK_SECONDARY)
    ax.yaxis.label.set_color(INK_SECONDARY)
    ax.xaxis.label.set_fontsize(FONT_SIZE)
    ax.yaxis.label.set_fontsize(FONT_SIZE)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.7, linestyle="-")
        ax.set_axisbelow(True)


def set_title(ax: Axes, title: str) -> None:
    """Titre aligné à gauche, en encre principale."""
    ax.set_title(title, loc="left", fontsize=FONT_SIZE + 2, fontweight="bold", color=INK)


def set_suptitle(figure: Figure, title: str) -> None:
    """Titre de figure multi-panneaux, aligné à gauche."""
    figure.suptitle(
        title, x=0.01, ha="left", fontsize=FONT_SIZE + 2.5, fontweight="bold", color=INK
    )


def style_legend(ax: Axes, **kwargs: Any) -> None:
    """Légende sans cadre, texte en encre secondaire."""
    legend = ax.legend(frameon=False, fontsize=FONT_SIZE - 0.5, **kwargs)
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)


def save_figure(figure: Figure, path: Path, dpi: int = 130) -> Path:
    """Enregistre la figure en PNG (dossier créé si besoin)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, facecolor=SURFACE)
    return path
