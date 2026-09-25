"""Composition nucléotidique et biais GC (GC skew).

Le **GC skew** ``(G - C) / (G + C)`` change de signe au niveau de l'origine et du
terminus de réplication de nombreux chromosomes bactériens circulaires (asymétrie de
mutation entre brin précoce et brin tardif). Le minimum du skew **cumulé** est un
indicateur classique de la position de l'origine.

Ce n'est qu'un indicateur : il est sans objet pour un génome linéaire de phage, et il
peut être brouillé par des réarrangements ou un transfert horizontal récent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio.SeqUtils import GC_skew, gc_fraction

from orfeval.models import Genome, ReferenceGene


def auto_window(genome_length: int) -> int:
    """Taille de fenêtre par défaut : ~200 fenêtres, arrondie à 100 nt (minimum 100)."""
    return max(100, round(genome_length / 200 / 100) * 100)


@dataclass(frozen=True, slots=True)
class ReplicationLandmarks:
    """Positions extrêmes du GC skew cumulé (génomes circulaires uniquement)."""

    origin_candidate: int
    terminus_candidate: int
    dnaa_position: int | None
    dnaa_distance: int | None


@dataclass(frozen=True, slots=True)
class CompositionResult:
    """Composition globale et profil de GC skew."""

    gc_content: float
    window: int
    profile: pd.DataFrame
    landmarks: ReplicationLandmarks | None


def gc_skew_profile(sequence: str, window: int) -> pd.DataFrame:
    """Profil de GC et de GC skew par fenêtres consécutives (``Bio.SeqUtils``).

    La dernière fenêtre est ignorée si elle fait moins de la moitié de la taille
    demandée : calculé sur quelques nucléotides, son skew serait surtout du bruit.

    Returns
    -------
    pandas.DataFrame
        Colonnes ``start`` (0-based), ``end``, ``gc``, ``skew``, ``cumulative_skew``.
    """
    if window <= 0:
        raise ValueError("la fenêtre doit être > 0")
    skews = np.asarray(GC_skew(sequence, window=window), dtype=float)
    starts = np.arange(0, len(sequence), window)
    ends = np.minimum(starts + window, len(sequence))
    if len(starts) > 1 and ends[-1] - starts[-1] < window / 2:
        skews, starts, ends = skews[:-1], starts[:-1], ends[:-1]
    gc = [gc_fraction(sequence[start:end]) for start, end in zip(starts, ends, strict=True)]
    return pd.DataFrame(
        {
            "start": starts,
            "end": ends,
            "gc": gc,
            "skew": skews,
            "cumulative_skew": np.cumsum(skews),
        }
    )


def _circular_distance(a: int, b: int, length: int) -> int:
    difference = abs(a - b) % length
    return min(difference, length - difference)


def replication_landmarks(
    profile: pd.DataFrame, genome_length: int, references: list[ReferenceGene]
) -> ReplicationLandmarks:
    """Candidats origine/terminus d'après le minimum/maximum du skew cumulé.

    Si un gène ``dnaA`` est annoté, sa distance au candidat origine est rapportée à titre
    de comparaison (``dnaA`` est souvent, mais pas toujours, proche de l'origine).
    """
    midpoints = ((profile["start"] + profile["end"]) // 2).to_numpy()
    cumulative = profile["cumulative_skew"].to_numpy()
    origin = int(midpoints[int(np.argmin(cumulative))])
    terminus = int(midpoints[int(np.argmax(cumulative))])
    dnaa = next((ref for ref in references if (ref.gene_name or "").lower() == "dnaa"), None)
    dnaa_position = dnaa.interval.five_prime(genome_length) if dnaa else None
    distance = (
        _circular_distance(origin, dnaa_position, genome_length)
        if dnaa_position is not None
        else None
    )
    return ReplicationLandmarks(origin, terminus, dnaa_position, distance)


def analyse_composition(
    genome: Genome, window: int | None = None, references: list[ReferenceGene] | None = None
) -> CompositionResult:
    """Calcule la composition du génome et, s'il est circulaire, les repères de réplication."""
    window = window or auto_window(genome.length)
    profile = gc_skew_profile(genome.sequence, window)
    landmarks = (
        replication_landmarks(profile, genome.length, references or [])
        if genome.is_circular
        else None
    )
    return CompositionResult(
        gc_content=float(gc_fraction(genome.sequence)),
        window=window,
        profile=profile,
        landmarks=landmarks,
    )
