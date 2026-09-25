"""Figures d'évaluation : rappel par longueur, précision-rappel, erreurs, comparaison."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

from orfeval.plotting.style import (
    BLUE,
    COLOR_CORRECT,
    COLOR_ERROR,
    INK_SECONDARY,
    LINE_WIDTH,
    MUTED,
    ORANGE,
    SURFACE,
    new_figure,
    set_suptitle,
    set_title,
    style_axes,
    style_legend,
)


def recall_by_length_figure(by_length: pd.DataFrame, *, title: str) -> Figure:
    """Rappel par classe de longueur des gènes de référence, avec IC de Wilson à 95 %."""
    figure = new_figure(7.5, 4.2)
    ax = figure.subplots()
    x = np.arange(len(by_length))
    recall = by_length["rappel"].to_numpy(dtype=float)
    lower = np.clip(recall - by_length["ic95_bas"].to_numpy(dtype=float), 0, None)
    upper = np.clip(by_length["ic95_haut"].to_numpy(dtype=float) - recall, 0, None)
    ax.bar(x, np.nan_to_num(recall), width=0.4, color=BLUE, zorder=2)
    ax.errorbar(
        x, recall, yerr=[lower, upper], fmt="none", ecolor=INK_SECONDARY, elinewidth=1, capsize=3
    )
    for position, total in zip(x, by_length["n_reference"], strict=True):
        ax.annotate(
            f"n = {total}",
            xy=(position, 0),
            xytext=(0, -26),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            color=MUTED,
            annotation_clip=False,
        )
    ax.set_xticks(x, by_length["classe"])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Rappel (gènes retrouvés / gènes annotés)")
    ax.set_xlabel("Longueur du gène de référence (nt)", labelpad=16)
    style_axes(ax)
    set_title(ax, title)
    return figure


def precision_recall_figure(
    pr_model: pd.DataFrame, pr_length: pd.DataFrame, n_predicted: int, *, title: str
) -> Figure:
    """Courbes précision-rappel : score du modèle codant vs longueur seule."""
    figure = new_figure(6.8, 4.8)
    ax = figure.subplots()
    ax.plot(
        pr_model["recall"],
        pr_model["precision"],
        color=BLUE,
        linewidth=LINE_WIDTH,
        label="Classement par score codant",
    )
    ax.plot(
        pr_length["recall"],
        pr_length["precision"],
        color=ORANGE,
        linewidth=LINE_WIDTH,
        label="Classement par longueur seule",
    )
    if 0 < n_predicted <= len(pr_model):
        point = pr_model.iloc[n_predicted - 1]
        ax.plot(
            [point["recall"]],
            [point["precision"]],
            marker="o",
            markersize=8,
            color=BLUE,
            markeredgecolor=SURFACE,
            markeredgewidth=2,
            linestyle="none",
            label="Seuil retenu",
        )
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Rappel (sensibilité)")
    ax.set_ylabel("Précision")
    style_axes(ax, grid="both")
    style_legend(ax, loc="lower left")
    set_title(ax, title)
    return figure


def score_length_figure(
    lengths: Sequence[int],
    scores: Sequence[float],
    correct: Sequence[bool],
    threshold: float,
    *,
    title: str,
) -> Figure:
    """Score de détection en fonction de la longueur, par statut d'appariement."""
    figure = new_figure(7.5, 4.8)
    ax = figure.subplots()
    lengths_array = np.asarray(lengths, dtype=float)
    scores_array = np.asarray(scores, dtype=float)
    correct_array = np.asarray(correct, dtype=bool)
    for mask, color, marker, label in (
        (correct_array, COLOR_CORRECT, "o", "Apparié à un gène annoté"),
        (~correct_array, COLOR_ERROR, "^", "Sans correspondance"),
    ):
        ax.scatter(
            lengths_array[mask],
            scores_array[mask],
            s=18,
            color=color,
            marker=marker,
            alpha=0.75,
            linewidths=0.4,
            edgecolors=SURFACE,
            label=f"{label} (n = {int(mask.sum())})",
        )
    ax.axhline(threshold, color=MUTED, linewidth=0.9)
    ax.annotate(
        f"seuil = {threshold:g} bits",
        xy=(1.0, threshold),
        xycoords=("axes fraction", "data"),
        xytext=(-4, 4),
        textcoords="offset points",
        ha="right",
        fontsize=8.5,
        color=INK_SECONDARY,
    )
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xlabel("Longueur du candidat (nt, échelle log)")
    ax.set_ylabel("Score de détection (bits, échelle symlog)")
    style_axes(ax, grid="both")
    style_legend(ax, loc="upper left")
    set_title(ax, title)
    return figure


def error_breakdown_figure(
    false_positives: Mapping[str, int], false_negatives: Mapping[str, int], *, title: str
) -> Figure:
    """Répartition des prédictions non appariées et des gènes manqués par cause."""
    figure = new_figure(10, 3.6)
    left, right = figure.subplots(1, 2)
    for ax, counts, subtitle in (
        (left, false_positives, "Prédictions sans correspondance"),
        (right, false_negatives, "Gènes de référence manqués"),
    ):
        labels = list(counts)[::-1]
        values = [counts[label] for label in labels]
        y = np.arange(len(labels))
        ax.barh(y, values, height=0.55, color=COLOR_ERROR, zorder=2)
        for position, value in zip(y, values, strict=True):
            ax.annotate(
                str(value),
                xy=(value, position),
                xytext=(4, 0),
                textcoords="offset points",
                va="center",
                fontsize=8.5,
                color=INK_SECONDARY,
            )
        ax.set_yticks(y, labels)
        ax.set_xlim(0, max([*values, 1]) * 1.18)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel("Nombre")
        style_axes(ax, grid="x")
        set_title(ax, subtitle)
    set_suptitle(figure, title)
    return figure


COMPARISON_METRICS = (
    ("sensitivity", "Sensibilité"),
    ("precision", "Précision"),
    ("start_accuracy", "Start exact (parmi retrouvés)"),
)


def _field(rows: Sequence[Any], key: str) -> np.ndarray:
    """Extrait un champ numérique de dictionnaires de métriques (NaN si absent)."""
    values = [
        row[key] if isinstance(row, Mapping) and row.get(key) is not None else np.nan
        for row in rows
    ]
    return np.asarray(values, dtype=float)


def comparison_figure(summaries: Sequence[Mapping[str, Any]], *, title: str) -> Figure:
    """Comparaison des analyses : chaque métrique avec son IC de Wilson à 95 %.

    ``summaries`` contient, pour chaque analyse, ``name`` et un dictionnaire par
    métrique (``value``, ``ci95_low``, ``ci95_high``).
    """
    figure = new_figure(11, 1.4 + 0.45 * len(summaries))
    axes = figure.subplots(1, len(COMPARISON_METRICS), sharey=True)
    names = [str(summary["name"]) for summary in summaries][::-1]
    y = np.arange(len(names))
    for ax, (key, label) in zip(axes, COMPARISON_METRICS, strict=True):
        rows = [summary.get(key) for summary in summaries][::-1]
        values = _field(rows, "value")
        lows = _field(rows, "ci95_low")
        highs = _field(rows, "ci95_high")
        ax.hlines(y, lows, highs, color=INK_SECONDARY, linewidth=1.2)
        ax.plot(
            values,
            y,
            marker="o",
            markersize=8,
            color=BLUE,
            markeredgecolor=SURFACE,
            markeredgewidth=2,
            linestyle="none",
        )
        for position, value in zip(y, values, strict=True):
            if not np.isnan(value):
                ax.annotate(
                    f"{value:.2f}",
                    xy=(value, position),
                    xytext=(0, 7),
                    textcoords="offset points",
                    ha="center",
                    fontsize=8,
                    color=INK_SECONDARY,
                )
        ax.set_xlim(0, 1.05)
        ax.set_ylim(-0.6, len(names) - 0.3)
        ax.set_yticks(y, names)
        style_axes(ax, grid="x")
        set_title(ax, label)
    set_suptitle(figure, title)
    return figure
