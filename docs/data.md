# Données

## Données de démonstration (`data/demo/`)

| Fichier | Organisme | Taille | Topologie | Code génétique | CDS (référence) |
|---|---|---|---|---|---|
| `NC_001416.1.gb.gz` | *Enterobacteria phage lambda* | 48 502 nt | linéaire | 11 | 73 (73) |
| `NC_000908.2.gb.gz` | *Mycoplasmoides genitalium* G37 | 580 076 nt | circulaire | 4 | 524 (504 hors pseudogènes) |

- **Source** : NCBI RefSeq (base `nuccore`), format `gbwithparts`, séquence et annotation
  complètes.
- **Téléchargement** : 25/09/2026, avec `orfeval fetch NC_001416.1 NC_000908.2 --outdir data/demo`.
- **Intégrité** : `data/demo/SHA256SUMS` (vérifier avec
  `cd data/demo && sha256sum -c SHA256SUMS`).
- **Taille totale** : environ 450 Ko compressés. Les fichiers gzip sont écrits sans
  horodatage, si bien qu'un même contenu donne toujours la même somme.
- **Pourquoi ces génomes ?**
  - Lambda est petit, dense et historiquement très étudié (Sanger et al., 1982).
  - *M. genitalium* est un génome bactérien minimal (Fraser et al., 1995), riche en AT, qui
    utilise le code génétique 4. Il permet l'expérience de contrôle sur le code génétique.
- **Conditions d'utilisation** : le NCBI n'impose pas de restriction à l'utilisation ni à la
  redistribution des données GenBank, mais rappelle que certains déposants peuvent revendiquer
  des droits sur leurs données
  ([politique du NCBI](https://www.ncbi.nlm.nih.gov/home/about/policies/)). Les deux
  génomes sont des séquences de référence publiques, largement utilisées en enseignement.
- **Mises à jour** : le NCBI peut modifier l'annotation d'un enregistrement sans changer son
  numéro de version. L'enregistrement de *M. genitalium* a été mis à jour en décembre 2025.
  Les résultats documentés correspondent aux fichiers versionnés ici ;
  `scripts/refresh_demo_data.sh` signale toute différence avec la version en ligne.

## Utiliser ses propres données

| Entrée | Prédiction | Évaluation |
|---|---|---|
| GenBank annoté (`.gb`, `.gbk`, `.gbff`, `.genbank`, avec ou sans `.gz`) | oui | oui (CDS du fichier) |
| FASTA (`.fa`, `.fasta`, `.fna`, `.fas`, avec ou sans `.gz`) | oui | non |

Contraintes de la version 0.1 :

- **Une séquence par fichier**. Séparez chromosome et plasmides ; un assemblage fragmenté
  n'est pas pris en charge.
- **Alphabet ADN IUPAC**. Plus de 5 % de bases ambiguës est refusé, car la séquence serait
  trop incomplète.
- **Au moins 1 000 nt** ; un avertissement est émis sous 20 kb, où l'auto-apprentissage
  dispose de peu de gènes.
- Pour un FASTA, précisez `--table` (code génétique) et `--topology` si besoin : par défaut,
  le code 11 et une topologie linéaire sont utilisés.

## Données générées

Tout ce qui est produit va dans `results/`, ignoré par Git sauf le dossier lui-même. Les
figures du README sont copiées depuis une exécution réelle de la démo par
`scripts/export_readme_figures.py` (`make figures`).

Les données téléchargées hors démonstration vont dans `data/raw/`, également ignoré par Git.
