# Reproductibilité

L'objectif : que n'importe qui obtienne **exactement** les chiffres de la documentation, et
sache comment ils ont été produits.

## Ce qui est garanti

| Élément | Mécanisme |
|---|---|
| Données d'entrée identiques | génomes versionnés dans `data/demo/`, sommes SHA-256, gzip déterministe |
| Calcul identique | prédicteur **déterministe** (aucun aléa, ordre de tri total) ; test dédié |
| Paramètres traçables | configuration YAML validée ; paramètres complets écrits dans `run_summary.json` et le rapport |
| Environnement traçable | versions de Python, du système et des dépendances dans `project_summary.json` et le rapport |
| Environnement reconstructible | `requirements.txt` (versions exactes), `Dockerfile` (Python 3.12-slim) |
| Provenance des entrées | chemin et SHA-256 de chaque fichier dans le rapport |
| Figures du README | copiées d'une exécution réelle (`make figures`), jamais produites à la main |
| Vérification continue | GitHub Actions : lint, typage, tests sur Python 3.11/3.12/3.13, démo, image Docker |

## Reproduire les résultats documentés

```bash
make install-locked         # versions exactes de requirements.txt
source .venv/bin/activate
cd data/demo && sha256sum -c SHA256SUMS && cd ../..
orfeval demo                              # docs/results.md
orfeval analyze config/sensitivity.yaml   # tableau de docs/methodology.md, section 11
```

Ou, sans installer Python :

```bash
docker build -t orfeval .
docker run --rm -v "$PWD/results/docker:/app/results" orfeval
```

## Ce qui peut varier

- **Versions des dépendances.** Avec `pip install -e .` (bornes minimales plutôt que versions
  exactes), une autre version de numpy ou de pandas peut en théorie modifier un arrondi. Les
  tests ne dépendent pas des derniers chiffres significatifs.
- **Données NCBI en ligne.** Retélécharger un génome peut donner une annotation plus récente ;
  les résultats documentés utilisent les fichiers du dépôt.
- **Rendu des figures.** Selon la version de matplotlib et les polices installées, les PNG
  peuvent différer au pixel près ; les données tracées restent identiques (voir les TSV).
- **Temps de calcul.** Ils dépendent de la machine ; ils sont donnés à titre indicatif.

## Environnement de référence

Les résultats documentés ont été produits avec :

```
Python 3.12.3 — Linux 6.8 (glibc 2.39), x86_64
orfeval 0.1.0 · biopython 1.88 · numpy 2.5.3 · pandas 3.0.6 · matplotlib 3.11.2
pydantic 2.13.5 · jinja2 3.1.6 · typer 0.27.2
```

La suite de tests passe également sous Python 3.11.16 et 3.13.15 (vérifié localement et en CI).
