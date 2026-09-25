"""Téléchargement de génomes depuis le NCBI (E-utilities via ``Bio.Entrez``).

Les fichiers sont écrits compressés avec un en-tête gzip déterministe (date nulle),
si bien qu'un même contenu produit toujours la même somme SHA-256.
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError

from Bio import Entrez, SeqIO

from orfeval.exceptions import DataFetchError
from orfeval.io.genome import sha256sum
from orfeval.logging_utils import get_logger

logger = get_logger(__name__)

ACCESSION_PATTERN = re.compile(r"^[A-Z]{1,6}_?[A-Z]{0,2}\d{5,12}(\.\d+)?$")


def validate_accession(accession: str) -> str:
    """Vérifie la forme d'un numéro d'accession NCBI (ex. ``NC_001416.1``)."""
    accession = accession.strip().upper()
    if not ACCESSION_PATTERN.match(accession):
        raise DataFetchError(f"numéro d'accession invalide : {accession!r}")
    return accession


def download_genbank_text(accession: str, *, email: str | None, api_key: str | None) -> str:
    """Récupère un enregistrement GenBank complet (``gbwithparts``) sous forme de texte."""
    accession = validate_accession(accession)
    entrez: Any = Entrez  # attributs de module non typés par Biopython
    entrez.tool = "orfeval"
    entrez.email = email
    if api_key:
        entrez.api_key = api_key
    if not email:
        logger.warning("Aucun e-mail fourni : le NCBI recommande d'en indiquer un (--email).")
    try:
        with Entrez.efetch(
            db="nuccore", id=accession, rettype="gbwithparts", retmode="text"
        ) as handle:
            content = handle.read()
    except (HTTPError, URLError, OSError) as exc:
        raise DataFetchError(f"échec du téléchargement de {accession} : {exc}") from exc
    text = content.decode("utf-8") if isinstance(content, bytes) else str(content)
    if not text.startswith("LOCUS"):
        raise DataFetchError(f"réponse inattendue du NCBI pour {accession} : {text[:120]!r}")
    try:
        record = SeqIO.read(io.StringIO(text), "genbank")
        _ = str(record.seq)[:10]
    except Exception as exc:
        raise DataFetchError(f"enregistrement {accession} illisible : {exc}") from exc
    return text


def write_gzip_deterministic(text: str, path: Path) -> None:
    """Écrit un texte compressé en gzip sans horodatage (sortie reproductible)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as gz:
        gz.write(text.encode("utf-8"))


def fetch_genbank(
    accession: str,
    outdir: Path,
    *,
    email: str | None = None,
    api_key: str | None = None,
    overwrite: bool = False,
) -> Path:
    """Télécharge un génome GenBank et l'écrit sous ``<outdir>/<accession>.gb.gz``.

    Returns
    -------
    Path
        Chemin du fichier écrit (ou existant si ``overwrite`` est faux).
    """
    accession = validate_accession(accession)
    destination = outdir / f"{accession}.gb.gz"
    if destination.exists() and not overwrite:
        logger.info("%s existe déjà : téléchargement ignoré.", destination)
        return destination
    text = download_genbank_text(accession, email=email, api_key=api_key)
    write_gzip_deterministic(text, destination)
    logger.info("Écrit : %s (sha256 %s)", destination, sha256sum(destination)[:12])
    return destination
