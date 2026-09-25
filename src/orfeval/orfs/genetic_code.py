"""Encodage numérique des codons et codes génétiques.

Chaque codon est représenté par un entier ``16*b0 + 4*b1 + b2`` (A=0, C=1, G=2, T=3),
soit 0..63 ; la valeur 64 désigne un codon contenant une base ambiguë. Ce codage
permet de traiter des génomes entiers avec des opérations numpy vectorisées.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from itertools import product

import numpy as np
from Bio.Data import CodonTable

from orfeval.config import check_table

NUCLEOTIDES = "ACGT"
CODONS: tuple[str, ...] = tuple("".join(bases) for bases in product(NUCLEOTIDES, repeat=3))
CODON_INDEX: dict[str, int] = {codon: index for index, codon in enumerate(CODONS)}
INVALID_CODON = 64
N_CODON_STATES = 65

_ENCODING = np.full(256, 4, dtype=np.uint8)
for _index, _base in enumerate(NUCLEOTIDES.encode("ascii")):
    _ENCODING[_base] = _index


def encode_codons(sequence: str, frame: int = 0) -> np.ndarray:
    """Encode les codons successifs d'une séquence à partir du décalage ``frame``.

    Parameters
    ----------
    sequence
        Séquence en majuscules.
    frame
        Décalage de lecture (0, 1 ou 2).

    Returns
    -------
    numpy.ndarray
        Tableau d'entiers (0..64) de longueur ``(len(sequence) - frame) // 3``.
    """
    if frame not in (0, 1, 2):
        raise ValueError(f"cadre de lecture invalide : {frame}")
    bases = _ENCODING[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    n_codons = max(0, (len(bases) - frame) // 3)
    triplets = bases[frame : frame + 3 * n_codons].reshape(n_codons, 3).astype(np.int16)
    ids = triplets[:, 0] * 16 + triplets[:, 1] * 4 + triplets[:, 2]
    ids[(triplets >= 4).any(axis=1)] = INVALID_CODON
    return ids.astype(np.int64)


@dataclass(frozen=True, slots=True)
class GeneticCode:
    """Code génétique restreint aux codons start retenus par l'analyse.

    Attributes
    ----------
    table_id
        Identifiant NCBI du code (11 = bactéries, 4 = mycoplasmes...).
    stop_codons, start_codons
        Codons stop du code et codons start retenus.
    stop_mask, start_mask, sense_mask
        Masques booléens indexés par l'identifiant de codon (taille 65).
    """

    table_id: int
    stop_codons: tuple[str, ...]
    start_codons: tuple[str, ...]
    stop_mask: np.ndarray
    start_mask: np.ndarray
    sense_mask: np.ndarray

    @classmethod
    def from_table(cls, table_id: int, start_codons: Iterable[str]) -> GeneticCode:
        """Construit le code à partir des tables NCBI de Biopython.

        Raises
        ------
        ValueError
            Si un codon start demandé n'est pas un codon d'initiation de ce code.
        """
        table = CodonTable.unambiguous_dna_by_id[check_table(table_id)]
        starts = tuple(start_codons)
        not_allowed = [codon for codon in starts if codon not in table.start_codons]
        if not_allowed:
            raise ValueError(
                f"codon(s) {not_allowed} non autorisé(s) comme start par le code {table_id} "
                f"(autorisés : {sorted(table.start_codons)})"
            )
        stop_mask = np.zeros(N_CODON_STATES, dtype=bool)
        start_mask = np.zeros(N_CODON_STATES, dtype=bool)
        for codon in table.stop_codons:
            stop_mask[CODON_INDEX[codon]] = True
        for codon in starts:
            start_mask[CODON_INDEX[codon]] = True
        sense_mask = ~stop_mask
        sense_mask[INVALID_CODON] = False
        return cls(
            table_id=table_id,
            stop_codons=tuple(sorted(table.stop_codons)),
            start_codons=starts,
            stop_mask=stop_mask,
            start_mask=start_mask,
            sense_mask=sense_mask,
        )

    def forward_table(self) -> dict[str, str]:
        """Correspondance codon → acide aminé (codons sens uniquement)."""
        return dict(CodonTable.unambiguous_dna_by_id[self.table_id].forward_table)
