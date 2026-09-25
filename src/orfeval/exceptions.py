"""Exceptions propres à orfeval.

Toutes dérivent de :class:`OrfevalError`, ce qui permet à la CLI de distinguer
les erreurs attendues (entrée invalide, configuration incorrecte) des bogues.
"""


class OrfevalError(Exception):
    """Erreur de base du package."""


class InputFormatError(OrfevalError):
    """Fichier d'entrée illisible, mal formé ou non pris en charge."""


class ConfigError(OrfevalError):
    """Configuration invalide ou incohérente."""


class AnalysisError(OrfevalError):
    """Analyse impossible avec les données et paramètres fournis."""


class DataFetchError(OrfevalError):
    """Échec du téléchargement de données publiques."""
