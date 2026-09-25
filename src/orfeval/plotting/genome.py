"""Figures à l'échelle du génome : carte des gènes et GC skew."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

from orfeval.features.composition import ReplicationLandmarks
from orfeval.models import Interval
from orfeval.plotting.style import (
    BLUE,
    COLOR_CORRECT,
    COLOR_ERROR,
    COLOR_OTHER_START,
    COLOR_REFERENCE,
    INK_SECONDARY,
    LINE_WIDTH,
    MUTED,
    ORANGE,
    SURFACE,
    new_figure,
    set_suptitle,
    set_title,
    style_axes,
)

KB = FuncFormatter(lambda value, _: f"{value / 1000:,.0f}".replace(",", " "))

TRACKS = {"ref+": 3, "pred+": 2, "pred-": 1, "ref-": 0}
TRACK_LABELS = ["Réf. (−)", "Préd. (−)", "Préd. (+)", "Réf. (+)"]


def _segments(interval: Interval, genome_length: int) -> list[tuple[int, int]]:
    """Découpe un intervalle qui passe l'origine en segments dans [0, N)."""
    if interval.right <= genome_length:
        return [(interval.left, interval.right)]
    return [(interval.left, genome_length), (0, interval.right - genome_length)]


def genome_map(
    genome_length: int,
    predictions: Sequence[Interval],
    prediction_colors: Sequence[str],
    references: Sequence[Interval] | None = None,
    reference_colors: Sequence[str] | None = None,
    *,
    legend: Sequence[tuple[str, str]] = (),
    title: str = "Carte des gènes",
    max_row_length: int = 150_000,
) -> Figure:
    """Carte linéaire des gènes, en plusieurs lignes pour les grands génomes.

    Quatre pistes par ligne : référence (+), prédiction (+), prédiction (−),
    référence (−). Les couleurs sont fournies par l'appelant (statut d'évaluation).
    """
    n_rows = max(1, math.ceil(genome_length / max_row_length))
    row_length = math.ceil(genome_length / n_rows)
    figure = new_figure(12, 0.9 + 1.15 * n_rows)
    axes = figure.subplots(n_rows, 1, squeeze=False)[:, 0]

    boxes: dict[tuple[int, int, str], list[tuple[int, int]]] = defaultdict(list)

    def add(interval: Interval, track: int, color: str) -> None:
        for start, end in _segments(interval, genome_length):
            position = start
            while position < end:
                row = position // row_length
                row_end = min(end, (row + 1) * row_length)
                boxes[(row, track, color)].append((position, row_end - position))
                position = row_end

    for interval, color in zip(predictions, prediction_colors, strict=True):
        add(interval, TRACKS["pred+" if interval.strand == 1 else "pred-"], color)
    if references is not None and reference_colors is not None:
        for interval, color in zip(references, reference_colors, strict=True):
            add(interval, TRACKS["ref+" if interval.strand == 1 else "ref-"], color)

    for (row, track, color), spans in boxes.items():
        # Un liseré couleur de fond sépare visuellement deux gènes adjacents.
        axes[row].broken_barh(
            spans, (track - 0.36, 0.72), facecolors=color, edgecolors=SURFACE, linewidth=0.6
        )

    for row, ax in enumerate(axes):
        style_axes(ax, grid=None)
        start = row * row_length
        ax.set_xlim(start, min(genome_length, start + row_length))
        ax.set_ylim(-0.6, 3.6)
        ax.set_yticks(range(4), TRACK_LABELS)
        ax.xaxis.set_major_formatter(KB)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[-1].set_xlabel("Position (kb)")
    set_suptitle(figure, title)
    if legend:
        handles = [Patch(facecolor=color, label=label) for label, color in legend]
        figure.legend(
            handles=handles,
            loc="outside lower center",
            ncols=len(handles),
            frameon=False,
            fontsize=8.5,
            labelcolor=INK_SECONDARY,
            handlelength=1.2,
        )
    return figure


EVALUATION_LEGEND: tuple[tuple[str, str], ...] = (
    ("Prédiction correcte, start exact", COLOR_CORRECT),
    ("Prédiction correcte, start différent", COLOR_OTHER_START),
    ("Prédiction sans correspondance / gène manqué", COLOR_ERROR),
    ("Gène de référence retrouvé", COLOR_REFERENCE),
)


def gc_skew_figure(
    profile: pd.DataFrame,
    landmarks: ReplicationLandmarks | None,
    *,
    title: str = "Composition : GC skew",
) -> Figure:
    """GC skew par fenêtre et skew cumulé, en deux panneaux (une seule échelle chacun)."""
    figure = new_figure(10, 5.2)
    top, bottom = figure.subplots(2, 1, sharex=True)
    midpoints = ((profile["start"] + profile["end"]) / 2).to_numpy()

    skew = profile["skew"].to_numpy()
    top.axhline(0, color=MUTED, linewidth=0.8)
    top.fill_between(midpoints, skew, 0, where=skew >= 0, color=BLUE, alpha=0.25, linewidth=0)
    top.fill_between(midpoints, skew, 0, where=skew < 0, color=ORANGE, alpha=0.25, linewidth=0)
    top.plot(midpoints, skew, color=INK_SECONDARY, linewidth=0.8)
    top.set_ylabel("(G − C) / (G + C)")
    style_axes(top)
    set_title(top, "GC skew par fenêtre")

    bottom.plot(midpoints, profile["cumulative_skew"], color=BLUE, linewidth=LINE_WIDTH)
    bottom.set_ylabel("Skew cumulé")
    bottom.set_xlabel("Position (kb)")
    bottom.xaxis.set_major_formatter(KB)
    style_axes(bottom)
    set_title(bottom, "GC skew cumulé")
    if landmarks is not None:
        for position, label in (
            (landmarks.origin_candidate, "min. : origine candidate"),
            (landmarks.terminus_candidate, "max. : terminus candidat"),
        ):
            bottom.axvline(position, color=MUTED, linewidth=0.9)
            bottom.annotate(
                label,
                xy=(position, 1.0),
                xycoords=("data", "axes fraction"),
                xytext=(4, -12),
                textcoords="offset points",
                fontsize=8.5,
                color=INK_SECONDARY,
            )
        if landmarks.dnaa_position is not None:
            value = float(np.interp(landmarks.dnaa_position, midpoints, profile["cumulative_skew"]))
            bottom.plot(
                [landmarks.dnaa_position],
                [value],
                marker="D",
                markersize=7,
                color=ORANGE,
                markeredgecolor="white",
                markeredgewidth=1.5,
                linestyle="none",
                label="gène dnaA annoté",
            )
            bottom.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
    set_suptitle(figure, title)
    return figure
