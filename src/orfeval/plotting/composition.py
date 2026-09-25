"""Figures de composition : modèle nul des ORF, usage des codons, profil amont."""

from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from orfeval.features.upstream import UpstreamProfile
from orfeval.orfs.null_model import NullModel
from orfeval.plotting.style import (
    BLUE,
    INK_SECONDARY,
    LINE_WIDTH,
    MUTED,
    ORANGE,
    SURFACE,
    new_figure,
    set_title,
    style_axes,
    style_legend,
)


def orf_length_null_figure(
    region_lengths: np.ndarray, null_model: NullModel, train_min_length: int, *, title: str
) -> Figure:
    """Nombre de régions sans stop d'au moins L codons : génome réel vs modèle nul."""
    figure = new_figure(7.5, 4.6)
    ax = figure.subplots()
    lengths = np.sort(np.asarray(region_lengths, dtype=int))
    grid = np.arange(0, max(int(lengths.max()) if lengths.size else 1, 1) + 1)
    observed = lengths.size - np.searchsorted(lengths, grid, side="left")
    expected = null_model.expected_count(grid)
    nucleotides = (grid + 1) * 3
    visible = observed > 0
    ax.plot(
        nucleotides[visible],
        observed[visible],
        color=BLUE,
        linewidth=LINE_WIDTH,
        label="Génome réel",
    )
    shown = expected >= 0.01
    ax.plot(
        nucleotides[shown],
        expected[shown],
        color=ORANGE,
        linewidth=LINE_WIDTH,
        label="Séquence aléatoire de même composition",
    )
    ax.axvline(train_min_length, color=MUTED, linewidth=0.9)
    ax.annotate(
        f"seuil d'auto-apprentissage\n{train_min_length} nt",
        xy=(train_min_length, 1.0),
        xycoords=("data", "axes fraction"),
        xytext=(5, -24),
        textcoords="offset points",
        fontsize=8.5,
        color=INK_SECONDARY,
    )
    ax.set_yscale("log")
    ax.set_xlabel("Longueur L entre deux codons stop (nt)")
    ax.set_ylabel("Nombre de régions de longueur ≥ L")
    style_axes(ax, grid="both")
    style_legend(ax, loc="upper right")
    set_title(ax, title)
    return figure


def rscu_figure(rscu: pd.DataFrame, *, title: str) -> Figure:
    """RSCU des gènes prédits vs gènes annotés (un point par codon sens)."""
    figure = new_figure(5.6, 5.4)
    ax = figure.subplots()
    data = rscu.dropna(subset=["rscu_predicted", "rscu_reference"])
    limit = float(max(data["rscu_predicted"].max(), data["rscu_reference"].max(), 1.0)) * 1.08
    ax.plot([0, limit], [0, limit], color=MUTED, linewidth=0.9)
    ax.scatter(
        data["rscu_reference"],
        data["rscu_predicted"],
        s=30,
        color=BLUE,
        edgecolors=SURFACE,
        linewidths=1.2,
        zorder=3,
    )
    if len(data) > 2:
        r = float(np.corrcoef(data["rscu_reference"], data["rscu_predicted"])[0, 1])
        ax.annotate(
            f"r de Pearson = {r:.3f}\n{len(data)} codons sens",
            xy=(0.03, 0.97),
            xycoords="axes fraction",
            va="top",
            fontsize=9,
            color=INK_SECONDARY,
        )
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_aspect("equal")
    ax.set_xlabel("RSCU — gènes annotés")
    ax.set_ylabel("RSCU — gènes prédits")
    style_axes(ax, grid="both")
    set_title(ax, title)
    return figure


def upstream_figure(
    predicted: UpstreamProfile | None, reference: UpstreamProfile | None, *, title: str
) -> Figure:
    """Information (bits) par position en amont des codons start."""
    figure = new_figure(7.5, 4.2)
    ax = figure.subplots()
    for profile, color, marker, label in (
        (predicted, BLUE, "o", "Starts prédits"),
        (reference, ORANGE, "s", "Starts annotés"),
    ):
        if profile is None:
            continue
        ax.plot(
            profile.information.index,
            profile.information.to_numpy(),
            color=color,
            linewidth=LINE_WIDTH,
            marker=marker,
            markersize=5,
            markeredgecolor=SURFACE,
            label=f"{label} (n = {profile.n_sequences})",
        )
    ax.set_xlabel("Position par rapport au codon start (nt)")
    ax.set_ylabel("Écart à la composition du génome (bits)")
    ax.set_ylim(bottom=0)
    style_axes(ax, grid="y")
    style_legend(ax, loc="upper left")
    set_title(ax, title)
    return figure
