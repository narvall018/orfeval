"""orfeval — prédiction auto-entraînée de gènes procaryotes et évaluation.

Le package est organisé en couches indépendantes :

- :mod:`orfeval.io` : lecture des génomes (GenBank/FASTA), écriture des résultats ;
- :mod:`orfeval.orfs` : recherche d'ORF, modèles statistiques et prédicteur ;
- :mod:`orfeval.evaluation` : comparaison avec l'annotation de référence ;
- :mod:`orfeval.features` : descripteurs du génome (composition, codons, amont) ;
- :mod:`orfeval.plotting` et :mod:`orfeval.report` : figures et rapport HTML ;
- :mod:`orfeval.pipeline` : orchestration, utilisée par la CLI et le dashboard.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("orfeval")
except PackageNotFoundError:  # pragma: no cover - package non installé
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
