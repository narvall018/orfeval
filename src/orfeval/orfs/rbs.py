"""Détection simplifiée du site de fixation du ribosome (Shine-Dalgarno).

Chez de nombreuses bactéries, l'ARNr 16S s'apparie à un motif riche en purines situé
quelques nucléotides en amont du codon start (consensus ``AGGAGG``). On recherche ici
la plus longue sous-chaîne du consensus (au moins 3 nt) placée à une distance
(« spacer ») plausible du codon start.

Limite : certaines lignées (dont les mycoplasmes) utilisent peu ce mécanisme ;
le motif est alors peu informatif, ce que l'analyse du profil amont permet de vérifier.
"""

from __future__ import annotations

from dataclasses import dataclass

SD_CONSENSUS = "AGGAGG"
OPTIMAL_SPACER = 7


def _consensus_submotifs(consensus: str, min_length: int) -> list[str]:
    motifs = {
        consensus[i:j]
        for i in range(len(consensus))
        for j in range(i + min_length, len(consensus) + 1)
    }
    return sorted(motifs, key=lambda motif: (-len(motif), motif))


_ALL_MOTIFS = _consensus_submotifs(SD_CONSENSUS, 3)


@dataclass(frozen=True, slots=True)
class RbsHit:
    """Motif détecté : sous-chaîne du consensus et distance au codon start."""

    motif: str
    spacer: int

    @property
    def length(self) -> int:
        """Longueur du motif (3 à 6 nt)."""
        return len(self.motif)


def find_rbs(
    upstream: str,
    *,
    min_motif: int = 3,
    spacer_min: int = 3,
    spacer_max: int = 15,
) -> RbsHit | None:
    """Cherche le meilleur motif de Shine-Dalgarno dans une séquence amont.

    Parameters
    ----------
    upstream
        Séquence 5'→3' dont le dernier nucléotide est adjacent au codon start.
    min_motif
        Longueur minimale du motif retenu.
    spacer_min, spacer_max
        Bornes de la distance entre la fin du motif et le codon start.

    Returns
    -------
    RbsHit | None
        Le motif le plus long ; à longueur égale, celui dont le spacer est le plus
        proche de 7 nt. ``None`` si aucun motif ne convient.
    """
    total = len(upstream)
    for size in range(len(SD_CONSENSUS), min_motif - 1, -1):
        best: RbsHit | None = None
        for motif in _ALL_MOTIFS:
            if len(motif) != size:
                continue
            position = upstream.find(motif)
            while position != -1:
                spacer = total - (position + size)
                if spacer_min <= spacer <= spacer_max:
                    hit = RbsHit(motif=motif, spacer=spacer)
                    if best is None or abs(spacer - OPTIMAL_SPACER) < abs(
                        best.spacer - OPTIMAL_SPACER
                    ):
                        best = hit
                position = upstream.find(motif, position + 1)
        if best is not None:
            return best
    return None
