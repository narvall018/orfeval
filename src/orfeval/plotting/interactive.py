"""Carte interactive des gènes (Plotly), utilisée par le dashboard.

Plotly est une dépendance optionnelle (extra ``app``) : il n'est importé qu'à l'appel.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pandas as pd

from orfeval.plotting.style import (
    COLOR_CORRECT,
    COLOR_ERROR,
    COLOR_OTHER_START,
    COLOR_REFERENCE,
    GRID,
    INK_SECONDARY,
    SURFACE,
)

if TYPE_CHECKING:
    import plotly.graph_objects as go

TRACK_ORDER = ["Réf. (−)", "Préd. (−)", "Préd. (+)", "Réf. (+)"]
STATUS_COLORS = {
    "correct (start exact)": COLOR_CORRECT,
    "correct (start différent)": COLOR_OTHER_START,
    "retrouvé (start exact)": COLOR_REFERENCE,
    "retrouvé (start différent)": COLOR_REFERENCE,
}


def _segments(frame: pd.DataFrame, genome_length: int) -> pd.DataFrame:
    """Découpe les gènes qui passent l'origine (fin > N) en deux segments affichables."""
    rows = []
    for record in frame.to_dict(orient="records"):
        left, right = int(record["start"]) - 1, int(record["end"])
        if right <= genome_length:
            rows.append({**record, "left": left, "right": right})
        else:
            rows.append({**record, "left": left, "right": genome_length})
            rows.append({**record, "left": 0, "right": right - genome_length})
    return pd.DataFrame(rows)


def _color(status: str) -> str:
    return STATUS_COLORS.get(status, COLOR_ERROR)


def genome_browser(
    genome_length: int,
    predictions: pd.DataFrame,
    references: pd.DataFrame | None = None,
    region: tuple[int, int] | None = None,
) -> go.Figure:
    """Carte des gènes prédits (et de référence) avec survol détaillé.

    Parameters
    ----------
    genome_length
        Longueur du génome.
    predictions
        Tableau des prédictions (colonnes ``gene_id``, ``start``, ``end``, ``strand``,
        ``score_bits`` et, si disponible, ``status``) en coordonnées 1-based.
    references
        Tableau de comparaison à la référence (``locus_tag``, ``start``, ``end``,
        ``strand``, ``product``, ``status``).
    region
        Fenêtre affichée ``(début, fin)`` en nt ; tout le génome par défaut.
    """
    import plotly.graph_objects as go

    figure = go.Figure()
    layers: list[tuple[pd.DataFrame, str, str]] = []
    prediction_data = predictions.copy()
    if "status" not in prediction_data:
        prediction_data["status"] = "prédiction"
    prediction_data["track"] = prediction_data["strand"].map({"+": "Préd. (+)", "-": "Préd. (−)"})
    prediction_data["label"] = prediction_data["gene_id"]
    prediction_data["detail"] = prediction_data["score_bits"].map(lambda s: f"score {s:.1f} bits")
    layers.append((prediction_data, "Prédictions", "status"))
    if references is not None:
        reference_data = references.copy()
        reference_data["track"] = reference_data["strand"].map({"+": "Réf. (+)", "-": "Réf. (−)"})
        reference_data["label"] = reference_data["locus_tag"]
        reference_data["detail"] = reference_data["product"].fillna("")
        layers.append((reference_data, "Référence", "status"))

    for data, group, status_column in layers:
        segments = _segments(data, genome_length)
        if segments.empty:
            continue
        for status, subset in segments.groupby(status_column, sort=False):
            custom: Any = subset[["label", "start", "end", "strand", "detail"]].to_numpy()
            figure.add_trace(
                go.Bar(
                    x=subset["right"] - subset["left"],
                    base=subset["left"],
                    y=subset["track"],
                    orientation="h",
                    marker={"color": _color(str(status)), "line": {"color": SURFACE, "width": 1}},
                    name=f"{group} — {status}",
                    legendgroup=group,
                    customdata=custom,
                    hovertemplate=(
                        "<b>%{customdata[0]}</b><br>%{customdata[1]:,}–%{customdata[2]:,} "
                        "(%{customdata[3]})<br>%{customdata[4]}<br>" + f"{status}<extra></extra>"
                    ),
                )
            )

    start, end = region or (0, genome_length)
    figure.update_layout(
        barmode="overlay",
        bargap=0.35,
        height=360,
        margin={"l": 10, "r": 10, "t": 30, "b": 10},
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font={"color": INK_SECONDARY},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "x": 0},
        hoverlabel={"bgcolor": "white"},
    )
    figure.update_xaxes(
        range=[start, end], title="Position (nt)", gridcolor=GRID, tickformat=",d", zeroline=False
    )
    figure.update_yaxes(categoryorder="array", categoryarray=TRACK_ORDER, gridcolor=GRID)
    return figure
