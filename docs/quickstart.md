# Démarrage rapide

En trois commandes, depuis un clone du dépôt :

```bash
make install && source .venv/bin/activate
orfeval demo                      # environ 15 s
xdg-open results/demo/report.html # macOS : open ; Windows : start
```

## Ce que fait `orfeval demo`

La démo lit [`config/demo.yaml`](../config/demo.yaml) et lance quatre analyses sur deux
génomes publics versionnés dans `data/demo/` :

| Analyse | Génome | But |
|---|---|---|
| `lambda` | phage lambda, 48,5 kb | petit génome dense |
| `mgenitalium` | *M. genitalium*, 580 kb, code 4 | génome bactérien minimal |
| `mgenitalium_code11` | le même, code 11 forcé | contrôle : effet d'un mauvais code génétique |
| `mgenitalium_rbs` | le même, start guidé par le RBS | analyse de sensibilité |

## Ce qui est produit

```
results/demo/
├── report.html              # rapport autonome (données, paramètres, résultats, limites, versions)
├── summary.tsv              # une ligne par analyse
├── comparison.png           # métriques de toutes les analyses, avec IC à 95 %
├── project_summary.json     # commande, versions, configuration
└── lambda/                  # un dossier par analyse
    ├── predictions.tsv      # gènes prédits + statut d'évaluation
    ├── predictions.gff3     # annotation lisible par IGV, Artemis, JBrowse…
    ├── proteins.faa         # protéines traduites avec le bon code génétique
    ├── reference_comparison.tsv  # chaque CDS annotée : retrouvée ? start exact ? cause ?
    ├── recall_by_length.tsv, precision_recall_*.tsv, gc_profile.tsv, codon_usage.tsv
    ├── run_summary.json     # toutes les métriques et paramètres
    └── figures/*.png        # 9 figures
```

## Et ensuite ?

- Explorer interactivement : `make run`, puis ouvrir <http://localhost:8501>.
- Analyser un autre génome : `orfeval fetch NC_000913.3` (*E. coli* K-12), puis
  `orfeval predict data/raw/NC_000913.3.gb.gz -o results/ecoli`.
- Comprendre la méthode : [methodology.md](methodology.md).
- Lire les résultats commentés : [results.md](results.md).
