"""Configuration du projet : modèles pydantic et chargement des fichiers YAML.

Les chemins relatifs d'un fichier de configuration sont résolus par rapport au
dossier qui contient ce fichier, ce qui rend les configurations portables.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Literal

import yaml
from Bio.Data import CodonTable
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

from orfeval.exceptions import ConfigError
from orfeval.models import Topology

_RUN_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def available_tables() -> list[int]:
    """Identifiants des codes génétiques connus de Biopython."""
    return sorted(CodonTable.unambiguous_dna_by_id)


def check_table(table_id: int) -> int:
    """Vérifie qu'un code génétique NCBI existe et le renvoie."""
    if table_id not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f"code génétique inconnu : {table_id} (connus : {available_tables()})")
    return table_id


class _StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class PredictorConfig(_StrictModel):
    """Paramètres du prédicteur de gènes.

    Les valeurs par défaut sont fixées *a priori* à partir de la littérature
    (voir ``docs/methodology.md``), et non optimisées sur les génomes de démonstration.
    """

    min_length: int = Field(90, ge=30, description="Longueur minimale d'un gène (nt, stop inclus).")
    start_codons: tuple[str, ...] = Field(("ATG", "GTG", "TTG"), min_length=1)
    start_strategy: Literal["longest", "rbs", "score"] = "longest"
    coding_model: Literal["codon", "dicodon"] = "codon"
    pseudocount: float = Field(1.0, gt=0)
    train_min_length: int | Literal["auto"] = "auto"
    min_training_genes: int = Field(20, ge=5)
    self_training_iterations: int = Field(2, ge=1, le=10)
    max_overlap: int = Field(60, ge=0)
    score_threshold: float = 0.0
    rbs_weight: float = Field(
        1.0, ge=0, description="Stratégie 'score' : bonus (bits) par nucléotide du motif SD."
    )
    rbs_min_motif: int = Field(3, ge=3, le=6)
    rbs_spacer_min: int = Field(3, ge=0)
    rbs_spacer_max: int = Field(15, ge=1)
    rbs_window: int = Field(21, ge=6, le=60)

    @field_validator("start_codons", mode="before")
    @classmethod
    def _normalise_codons(cls, value: Any) -> tuple[str, ...]:
        if isinstance(value, str):
            value = [value]
        codons = tuple(str(codon).upper() for codon in value)
        for codon in codons:
            if len(codon) != 3 or set(codon) - set("ACGT"):
                raise ValueError(f"codon start invalide : {codon!r}")
        if len(set(codons)) != len(codons):
            raise ValueError("codons start en double")
        return codons

    @field_validator("train_min_length")
    @classmethod
    def _check_train_length(cls, value: int | str) -> int | str:
        if isinstance(value, int) and value < 90:
            raise ValueError("train_min_length doit être >= 90 nt ou 'auto'")
        return value

    @model_validator(mode="after")
    def _check_spacer(self) -> PredictorConfig:
        if self.rbs_spacer_min > self.rbs_spacer_max:
            raise ValueError("rbs_spacer_min doit être <= rbs_spacer_max")
        if self.rbs_spacer_max + self.rbs_min_motif > self.rbs_window:
            raise ValueError("rbs_window trop court pour rbs_spacer_max + rbs_min_motif")
        return self


class AnalysisConfig(_StrictModel):
    """Paramètres des analyses descriptives et de l'évaluation."""

    gc_window: int | Literal["auto"] = "auto"
    upstream_length: int = Field(25, ge=6, le=100)
    length_bins: tuple[int, ...] = (150, 300, 600, 1200)
    figure_dpi: int = Field(130, ge=50, le=600)

    @field_validator("length_bins")
    @classmethod
    def _sorted_bins(cls, value: tuple[int, ...]) -> tuple[int, ...]:
        if not value or list(value) != sorted(set(value)) or value[0] <= 0:
            raise ValueError("length_bins doit être une liste strictement croissante d'entiers > 0")
        return value

    @field_validator("gc_window")
    @classmethod
    def _check_window(cls, value: int | str) -> int | str:
        if isinstance(value, int) and value < 50:
            raise ValueError("gc_window doit être >= 50 nt ou 'auto'")
        return value


class RunConfig(_StrictModel):
    """Une analyse : un génome et d'éventuelles surcharges de paramètres."""

    name: str
    genome: Path
    description: str = ""
    translation_table: int | None = None
    topology: Topology | None = None
    predictor: dict[str, Any] = Field(default_factory=dict)

    @field_validator("name")
    @classmethod
    def _check_name(cls, value: str) -> str:
        if not _RUN_NAME.match(value):
            raise ValueError(
                f"nom d'analyse invalide : {value!r} (lettres, chiffres, '_', '-', '.')"
            )
        return value

    @field_validator("translation_table")
    @classmethod
    def _check_table(cls, value: int | None) -> int | None:
        return None if value is None else check_table(value)

    def resolve_predictor(self, base: PredictorConfig) -> PredictorConfig:
        """Fusionne les surcharges de l'analyse avec la configuration globale."""
        if not self.predictor:
            return base
        try:
            return PredictorConfig.model_validate({**base.model_dump(), **self.predictor})
        except ValidationError as exc:
            raise ConfigError(
                f"paramètres invalides pour l'analyse '{self.name}' :\n{exc}"
            ) from exc


class ProjectConfig(_StrictModel):
    """Configuration complète d'un projet d'analyse."""

    output_dir: Path = Path("results")
    predictor: PredictorConfig = Field(default_factory=PredictorConfig)
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    runs: list[RunConfig] = Field(min_length=1)

    @model_validator(mode="after")
    def _unique_names(self) -> ProjectConfig:
        names = [run.name for run in self.runs]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ValueError(f"noms d'analyse en double : {duplicates}")
        return self


def _resolve(path: Path, base_dir: Path) -> Path:
    return path if path.is_absolute() else (base_dir / path).resolve()


def load_config(path: Path, *, check_files: bool = True) -> ProjectConfig:
    """Charge et valide un fichier de configuration YAML.

    Parameters
    ----------
    path
        Chemin du fichier YAML.
    check_files
        Si vrai, vérifie que chaque génome référencé existe.

    Raises
    ------
    ConfigError
        Si le fichier est absent, illisible ou invalide.
    """
    if not path.is_file():
        raise ConfigError(f"fichier de configuration introuvable : {path}")
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"YAML invalide dans {path} : {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"{path} doit contenir un dictionnaire YAML")

    try:
        config = ProjectConfig.model_validate(raw)
    except ValidationError as exc:
        raise ConfigError(f"configuration invalide ({path}) :\n{exc}") from exc

    base_dir = path.resolve().parent
    runs = [
        run.model_copy(update={"genome": _resolve(run.genome, base_dir)}) for run in config.runs
    ]
    config = config.model_copy(
        update={"runs": runs, "output_dir": _resolve(config.output_dir, base_dir)}
    )
    if check_files:
        missing = [str(run.genome) for run in config.runs if not run.genome.is_file()]
        if missing:
            raise ConfigError("génome(s) introuvable(s) :\n  " + "\n  ".join(missing))
    # Valide dès maintenant les surcharges propres à chaque analyse.
    for run in config.runs:
        run.resolve_predictor(config.predictor)
    return config
