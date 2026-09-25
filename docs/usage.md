# Utilisation

## Ligne de commande

```text
orfeval [--verbose | --quiet] [--log-file FICHIER] COMMANDE [OPTIONS]
```

| Commande | Rôle |
|---|---|
| `fetch ACCESSION…` | télécharge des génomes GenBank complets depuis le NCBI (`Bio.Entrez`) |
| `validate GÉNOME` | vérifie un fichier et résume son contenu (taille, GC, topologie, code génétique, CDS) ; `--json` pour les scripts |
| `predict GÉNOME` | prédit (et évalue si le fichier est annoté) un génome, sans fichier de configuration ; `--no-figures` pour les tableaux seuls |
| `analyze CONFIG.yaml` | exécute toutes les analyses d'une configuration, puis génère le rapport |
| `report DOSSIER` | régénère le rapport HTML à partir d'un dossier de résultats, sans recalcul |
| `demo` | lance la démonstration (`config/demo.yaml`) |

`orfeval COMMANDE --help` détaille les options. Les erreurs attendues (fichier illisible,
paramètre invalide…) produisent un message en français et un code de sortie 1, sans trace
Python. `--verbose` affiche les messages de débogage.

### Exemples

```bash
# Télécharger un génome (l'e-mail est recommandé par le NCBI)
orfeval fetch NC_000913.3 --email vous@exemple.org --outdir data/raw

# Vérifier un fichier (tableau lisible, ou JSON pour un script)
orfeval validate data/raw/NC_000913.3.gb.gz
orfeval -q validate data/raw/NC_000913.3.gb.gz --json | jq .gc

# Prédire avec des paramètres explicites
orfeval predict data/raw/NC_000913.3.gb.gz -o results/ecoli --min-length 120 --start-strategy longest

# Tableaux, GFF3 et protéines seulement (environ un tiers plus rapide sur E. coli)
orfeval predict data/raw/NC_000913.3.gb.gz -o results/ecoli --no-figures

# FASTA sans annotation : prédiction seule, code génétique et topologie à préciser
orfeval predict mon_genome.fna -o results/mon_genome --table 11 --topology circular
```

## Fichier de configuration

[`config/default.yaml`](../config/default.yaml) documente **toutes** les options avec leur
valeur par défaut. Structure :

```yaml
output_dir: ../results/mon_projet      # relatif au fichier de configuration
predictor:                             # paramètres communs à toutes les analyses
  min_length: 90
  start_strategy: longest
analysis:
  gc_window: auto
runs:
  - name: ecoli
    genome: ../data/raw/NC_000913.3.gb.gz
    description: "Escherichia coli K-12 MG1655"
  - name: ecoli_rbs
    genome: ../data/raw/NC_000913.3.gb.gz
    predictor: {start_strategy: rbs}   # surcharge propre à cette analyse
    translation_table: 11              # force le code génétique
```

La configuration est validée par pydantic avant tout calcul : une clé inconnue, un codon
invalide, un code génétique inexistant, des noms d'analyse en double ou un fichier absent sont
signalés avec un message explicite.

Priorité du code génétique : `translation_table` de l'analyse > `/transl_table` de
l'annotation > 11 par défaut (avec avertissement). La source retenue est indiquée dans le
rapport.

## Formats de sortie

| Fichier | Format | Contenu |
|---|---|---|
| `predictions.tsv` | TSV | coordonnées 1-based inclusives, brin, codons start/stop, scores (bits), motif RBS, statut |
| `predictions.gff3` | GFF3 | type `CDS` ; `Is_circular=true` et fin > longueur pour un gène qui passe l'origine |
| `proteins.faa` | FASTA | traduction `cds=True` (start alternatif traduit en M), code génétique de l'analyse |
| `reference_comparison.tsv` | TSV | chaque CDS annotée : retrouvée, start exact, cause si manquée |
| `run_summary.json` | JSON | métriques (avec IC), paramètres, provenance (SHA-256), catégories d'écarts |

Colonnes `status` de `predictions.tsv` : `correct (start exact)`, `correct (start différent)`,
ou la catégorie de la prédiction non appariée (`pseudogène ou CDS partielle`,
`antisens d'un gène annoté`, `hors gènes annotés`…).

## Dashboard

```bash
make run          # ou : streamlit run app/streamlit_app.py
```

- Barre latérale : génome de démonstration ou fichier importé (GenBank ou FASTA, gzip
  accepté), et paramètres du prédicteur.
- Onglets :
  - **Résumé** : métriques et intervalles de confiance ;
  - **Gènes prédits** : tableau filtrable, téléchargements TSV, GFF3 et FASTA ;
  - **Carte interactive** : navigateur de gènes Plotly avec survol ;
  - **Évaluation** : courbes et gènes manqués, avec leur cause ;
  - **Composition** ;
  - **Méthode**.

Le dashboard ne contient aucune logique scientifique : il appelle `orfeval.pipeline.analyze_genome`
et les fonctions de figures du package, et met les résultats en cache.

## Utilisation comme bibliothèque

```python
from pathlib import Path
from orfeval.config import PredictorConfig
from orfeval.evaluation import evaluate
from orfeval.io import load_genome
from orfeval.orfs import GenePredictor

loaded = load_genome(Path("data/demo/NC_000908.2.gb.gz"))
genome = loaded.genome  # le prédicteur ne reçoit que la séquence
predictor = GenePredictor(PredictorConfig(), table_id=genome.declared_table)
prediction = predictor.predict(genome)
evaluation = evaluate(prediction, loaded)  # l'annotation ne sert qu'ici
print(len(prediction.genes), evaluation.sensitivity, evaluation.precision)
```
