"""Prédicteur de gènes auto-entraîné.

Étapes (voir ``docs/methodology.md``) :

1. recherche des ORF sur les six cadres ;
2. modèle nul → longueur au-delà de laquelle un ORF est très improbable par hasard ;
3. entraînement du modèle codant sur ces longs ORF ;
4. score de détection de chaque ORF et start préféré ;
5. sélection gloutonne des gènes compatibles (décalage du start si besoin), seuil ;
6. ré-entraînement sur les gènes prédits (itérations).

Le prédicteur ne reçoit que la séquence : il n'a **aucun accès** à l'annotation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orfeval.config import PredictorConfig
from orfeval.exceptions import AnalysisError
from orfeval.logging_utils import get_logger
from orfeval.models import Genome, Interval, PredictedGene
from orfeval.orfs.coding_model import CodingModel
from orfeval.orfs.finder import Orf, OrfSearch, find_orfs
from orfeval.orfs.genetic_code import GeneticCode, encode_codons
from orfeval.orfs.null_model import NullModel
from orfeval.orfs.rbs import RbsHit, find_rbs
from orfeval.orfs.selection import SelectionCandidate, preferred_start, resolve_overlaps

logger = get_logger(__name__)


@dataclass(frozen=True, slots=True)
class TrainingRound:
    """Bilan d'une itération d'auto-apprentissage."""

    iteration: int
    n_training_genes: int
    n_training_codons: int
    n_predicted: int


@dataclass(frozen=True, slots=True)
class PredictionResult:
    """Sortie complète du prédicteur.

    Attributes
    ----------
    genes
        Gènes prédits (score > seuil), triés par position.
    ranked
        Tous les candidats retenus par la sélection gloutonne, triés par score
        décroissant (y compris sous le seuil) : base de la courbe précision-rappel.
    length_ranked
        Référence naïve sans modèle : ORF classés par longueur, même règle de
        chevauchement. Sert à mesurer l'apport du modèle codant.
    """

    genes: list[PredictedGene]
    ranked: list[PredictedGene]
    length_ranked: list[Interval]
    table_id: int
    config: PredictorConfig
    null_model: NullModel
    train_min_length: int
    training_rounds: list[TrainingRound]
    model: CodingModel
    region_lengths: np.ndarray
    n_orfs: int
    orf_stop_keys: frozenset[tuple[int, int]]


class GenePredictor:
    """Prédit les gènes codants d'un génome procaryote à partir de sa seule séquence.

    Parameters
    ----------
    config
        Paramètres du prédicteur.
    table_id
        Code génétique NCBI à utiliser.
    """

    def __init__(self, config: PredictorConfig, table_id: int) -> None:
        self.config = config
        self.table_id = table_id
        try:
            self.code = GeneticCode.from_table(table_id, config.start_codons)
        except ValueError as exc:
            raise AnalysisError(str(exc)) from exc

    # ------------------------------------------------------------------ étapes
    def _training_threshold(self, null_model: NullModel) -> int:
        if self.config.train_min_length == "auto":
            return max(null_model.significant_length(expected=1.0), self.config.min_length)
        return int(self.config.train_min_length)

    def _initial_training_set(self, orfs: list[Orf], threshold: int) -> tuple[list[str], int]:
        long_orfs = [orf for orf in orfs if orf.length >= threshold]
        needed = self.config.min_training_genes
        if len(long_orfs) >= needed:
            return [orf.sequence for orf in long_orfs], threshold
        if len(orfs) < needed:
            raise AnalysisError(
                f"seulement {len(orfs)} ORF trouvés : impossible d'entraîner le modèle "
                f"(minimum {needed}). Génome trop court ou paramètres trop stricts ?"
            )
        longest = sorted(orfs, key=lambda orf: orf.length, reverse=True)[:needed]
        fallback = min(orf.length for orf in longest)
        logger.warning(
            "Seulement %d ORF >= %d nt : entraînement sur les %d plus longs (>= %d nt).",
            len(long_orfs),
            threshold,
            needed,
            fallback,
        )
        return [orf.sequence for orf in longest], fallback

    def _rbs_hits(self, orfs: list[Orf]) -> list[list[RbsHit | None]]:
        cfg = self.config
        return [
            [
                find_rbs(
                    orf.upstream_of(int(offset), cfg.rbs_window),
                    min_motif=cfg.rbs_min_motif,
                    spacer_min=cfg.rbs_spacer_min,
                    spacer_max=cfg.rbs_spacer_max,
                )
                for offset in orf.start_offsets
            ]
            for orf in orfs
        ]

    def _select(
        self,
        orfs: list[Orf],
        intervals: list[list[Interval]],
        hits: list[list[RbsHit | None]],
        start_scores: list[np.ndarray],
        genome_length: int,
    ) -> list[PredictedGene]:
        """Sélection gloutonne par score de détection décroissant."""
        cfg = self.config
        detection = np.array([scores.max() for scores in start_scores])
        candidates = [
            SelectionCandidate(
                intervals=orf_intervals,
                preferred=preferred_start(scores, orf_hits, cfg.start_strategy, cfg.rbs_weight),
            )
            for orf_intervals, scores, orf_hits in zip(intervals, start_scores, hits, strict=True)
        ]
        priority = sorted(
            range(len(orfs)),
            key=lambda i: (-detection[i], -orfs[i].length, orfs[i].left, orfs[i].strand),
        )
        accepted = resolve_overlaps(candidates, priority, genome_length, cfg.max_overlap)
        genes = []
        for index, start in accepted:
            orf = orfs[index]
            offset = int(orf.start_offsets[start])
            hit = hits[index][start]
            genes.append(
                PredictedGene(
                    interval=intervals[index][start],
                    sequence=orf.sequence[offset:],
                    score=float(detection[index]),
                    start_score=float(start_scores[index][start]),
                    rbs_motif=hit.motif if hit else None,
                    rbs_spacer=hit.spacer if hit else None,
                    n_alternative_starts=len(orf.start_offsets),
                    start_adjusted=start != candidates[index].preferred,
                )
            )
        return genes

    def _length_baseline(
        self, orfs: list[Orf], intervals: list[list[Interval]], genome_length: int
    ) -> list[Interval]:
        candidates = [
            SelectionCandidate(intervals=orf_intervals, preferred=0) for orf_intervals in intervals
        ]
        priority = sorted(
            range(len(orfs)), key=lambda i: (-orfs[i].length, orfs[i].left, orfs[i].strand)
        )
        accepted = resolve_overlaps(candidates, priority, genome_length, self.config.max_overlap)
        return [intervals[index][start] for index, start in accepted]

    # ------------------------------------------------------------------ API
    def find(self, genome: Genome) -> OrfSearch:
        """Recherche des ORF avec les paramètres du prédicteur."""
        return find_orfs(
            genome.sequence,
            self.code,
            min_length=self.config.min_length,
            circular=genome.is_circular,
            upstream_length=self.config.rbs_window,
        )

    def predict(self, genome: Genome) -> PredictionResult:
        """Prédit les gènes du génome.

        Raises
        ------
        AnalysisError
            Si le génome ne fournit pas assez d'ORF pour l'auto-apprentissage.
        """
        cfg = self.config
        search = self.find(genome)
        orfs = search.orfs
        logger.info(
            "%d ORF >= %d nt (code génétique %d, génome %s)",
            len(orfs),
            cfg.min_length,
            self.table_id,
            genome.topology,
        )
        null_model = NullModel.from_sequence(genome.sequence, self.code.stop_codons)
        threshold = self._training_threshold(null_model)
        training, threshold = self._initial_training_set(orfs, threshold)
        logger.info(
            "Modèle nul : p(stop) = %.4f ; seuil d'auto-apprentissage = %d nt (%d ORF)",
            null_model.p_stop,
            threshold,
            len(training),
        )

        encoded = [encode_codons(orf.sequence) for orf in orfs]
        start_indices = [orf.start_offsets // 3 for orf in orfs]
        intervals = [
            [orf.interval_for(int(offset)) for offset in orf.start_offsets] for orf in orfs
        ]
        hits = self._rbs_hits(orfs)

        rounds: list[TrainingRound] = []
        model: CodingModel | None = None
        ranked: list[PredictedGene] = []
        for iteration in range(1, cfg.self_training_iterations + 1):
            model = CodingModel.train(
                training,
                self.code,
                null_model.composition,
                kind=cfg.coding_model,
                pseudocount=cfg.pseudocount,
            )
            start_scores = [
                model.start_scores(ids, indices)
                for ids, indices in zip(encoded, start_indices, strict=True)
            ]
            ranked = self._select(orfs, intervals, hits, start_scores, genome.length)
            predicted = [gene for gene in ranked if gene.score > cfg.score_threshold]
            rounds.append(
                TrainingRound(
                    iteration=iteration,
                    n_training_genes=model.n_genes,
                    n_training_codons=model.n_codons,
                    n_predicted=len(predicted),
                )
            )
            logger.info(
                "Itération %d : entraînement sur %d séquences → %d gènes prédits",
                iteration,
                model.n_genes,
                len(predicted),
            )
            retrain = [gene.sequence for gene in predicted if gene.length >= threshold]
            if len(retrain) >= cfg.min_training_genes:
                training = retrain

        if model is None:  # pragma: no cover - au moins une itération (validé par pydantic)
            raise AnalysisError("aucune itération d'entraînement effectuée")
        genes = sorted(
            (gene for gene in ranked if gene.score > cfg.score_threshold),
            key=lambda gene: (gene.interval.left, gene.interval.strand),
        )
        return PredictionResult(
            genes=genes,
            ranked=ranked,
            length_ranked=self._length_baseline(orfs, intervals, genome.length),
            table_id=self.table_id,
            config=cfg,
            null_model=null_model,
            train_min_length=threshold,
            training_rounds=rounds,
            model=model,
            region_lengths=search.region_lengths,
            n_orfs=len(orfs),
            orf_stop_keys=frozenset(
                (orf.strand, orf_intervals[0].three_prime(genome.length))
                for orf, orf_intervals in zip(orfs, intervals, strict=True)
            ),
        )
