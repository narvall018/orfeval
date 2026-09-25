"""Configuration centralisée du logging."""

from __future__ import annotations

import logging
from pathlib import Path

LOGGER_NAME = "orfeval"
_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
_DATE_FORMAT = "%H:%M:%S"


def setup_logging(level: int | str = logging.INFO, log_file: Path | None = None) -> logging.Logger:
    """Configure le logger racine du package.

    Parameters
    ----------
    level
        Niveau de journalisation (``"DEBUG"``, ``"INFO"``, ...).
    log_file
        Fichier optionnel qui reçoit une copie des messages.

    Returns
    -------
    logging.Logger
        Le logger ``orfeval`` configuré.
    """
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    formatter = logging.Formatter(_FORMAT, datefmt=_DATE_FORMAT)
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.propagate = False
    return logger


def get_logger(name: str) -> logging.Logger:
    """Renvoie un logger enfant de ``orfeval`` (``name`` est typiquement ``__name__``)."""
    if name == LOGGER_NAME or name.startswith(f"{LOGGER_NAME}."):
        return logging.getLogger(name)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")
