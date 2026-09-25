"""Le prédicteur sur un génome synthétique dont on connaît tous les gènes.

C'est le seul cadre où la vérité est parfaitement connue : il vérifie que la méthode
fonctionne quand ses hypothèses sont respectées (usage des codons biaisé dans les
gènes, intergénique aléatoire). Ce n'est pas une mesure de performance sur du réel.
"""

from orfeval.config import PredictorConfig
from orfeval.evaluation.matching import match_genes
from orfeval.models import ReferenceGene
from orfeval.orfs.predictor import GenePredictor


def test_planted_genes_are_recovered(synthetic):
    genome = synthetic.genome
    result = GenePredictor(PredictorConfig(), 11).predict(genome)
    references = [ReferenceGene(gene, f"g{i}") for i, gene in enumerate(synthetic.genes)]
    matching = match_genes([g.interval for g in result.genes], references, genome.length)
    sensitivity = matching.ref_found.mean()
    precision = matching.pred_correct.mean()
    assert sensitivity >= 0.95
    assert precision >= 0.85
    # l'auto-apprentissage a bien démarré sur les longs ORF puis itéré
    assert len(result.training_rounds) == 2
    assert result.train_min_length > PredictorConfig().min_length


def test_prediction_is_deterministic(synthetic):
    first = GenePredictor(PredictorConfig(), 11).predict(synthetic.genome)
    second = GenePredictor(PredictorConfig(), 11).predict(synthetic.genome)
    assert [(g.interval, g.score) for g in first.genes] == [
        (g.interval, g.score) for g in second.genes
    ]


def test_predictor_never_reads_the_annotation(synthetic):
    """Garantie structurelle : l'API de prédiction ne reçoit qu'un objet Genome."""
    import inspect

    parameters = inspect.signature(GenePredictor.predict).parameters
    assert list(parameters) == ["self", "genome"]


def test_every_strategy_runs(synthetic):
    for strategy in ("longest", "rbs", "score"):
        config = PredictorConfig(start_strategy=strategy, coding_model="dicodon")
        result = GenePredictor(config, 11).predict(synthetic.genome)
        assert result.genes
        assert all(gene.score > config.score_threshold for gene in result.genes)
