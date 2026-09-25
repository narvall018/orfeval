# Journal des modifications

Format inspiré de [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) ; versions
selon [Semantic Versioning](https://semver.org/lang/fr/).

## [Non publié]

### Ajouté

- `orfeval validate` affiche le taux de GC, et accepte `--json` pour les scripts.
- `orfeval predict --no-figures` : tableaux, GFF3 et protéines seulement (environ un tiers
  plus rapide sur *E. coli*).
- `make check` : lint, types et tests en une commande.
- Marqueur `py.typed` (PEP 561).
- `.editorconfig`, `.gitattributes`.
- Dependabot pour les GitHub Actions.
- Tests : gène du brin − à cheval sur l'origine ; recherche de `config/demo.yaml`.

### Modifié

- CI : durée maximale par job.
- README : lien vers le projet msannot.

## [0.1.0] — 2026-09-25

Première version publique.

### Ajouté

- **Prédicteur auto-entraîné** :
  - recherche d'ORF vectorisée sur six cadres, pour tout code génétique NCBI ;
  - génomes circulaires ;
  - modèle nul géométrique pour le seuil d'auto-apprentissage ;
  - modèle codant (codons ou dicodons) en log-rapport de vraisemblance ;
  - stratégies de choix du start (`longest`, `rbs`, `score`) ;
  - sélection gloutonne avec décalage du start ;
  - ré-entraînement itératif.
- **Évaluation** :
  - appariement par extrémité 3′ ;
  - exactitude du start ;
  - intervalles de confiance de Wilson ;
  - courbes précision-rappel face à la longueur seule ;
  - métriques nucléotidiques ;
  - catégorisation des faux positifs et faux négatifs.
- **Descripteurs** : GC skew et repères de réplication, RSCU, GC3, profil amont des starts.
- **Entrées/sorties** :
  - lecture GenBank/FASTA (gzip) ;
  - téléchargement NCBI ;
  - sorties TSV, GFF3, FASTA protéique, JSON ;
  - rapport HTML autonome.
- **CLI** `orfeval` : `fetch`, `validate`, `predict`, `analyze`, `report`, `demo`.
- **Dashboard** Streamlit avec carte interactive Plotly.
- **Données** : phage lambda (NC_001416.1) et *M. genitalium* G37 (NC_000908.2), avec leurs
  sommes SHA-256.
- **Qualité** : 95 tests (unitaires, intégration, génome synthétique, CLI, dashboard), ruff,
  mypy, pre-commit, CI GitHub Actions (Python 3.11–3.13, démo, Docker).
- **Documentation** : README, 9 documents dans `docs/`, analyse de sensibilité
  (`config/sensitivity.yaml`).

### Choix documentés

- Stratégie de start par défaut : `longest`, retenue après comparaison sur les génomes de
  démonstration (voir `docs/methodology.md`, section 11).

[Non publié]: https://github.com/narvall018/orfeval/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/narvall018/orfeval/releases/tag/v0.1.0
