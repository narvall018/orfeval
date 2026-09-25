# Données

| Dossier | Contenu | Versionné |
|---|---|---|
| `demo/` | génomes de démonstration (GenBank compressé) et `SHA256SUMS` | oui |
| `raw/` | génomes téléchargés avec `orfeval fetch` | non (`.gitignore`) |

## Génomes de démonstration

| Fichier | Organisme | Taille | Code génétique | SHA-256 (début) |
|---|---|---|---|---|
| `demo/NC_001416.1.gb.gz` | *Enterobacteria phage lambda* | 48 502 nt, linéaire | 11 | `fbb33cbe…` |
| `demo/NC_000908.2.gb.gz` | *Mycoplasmoides genitalium* G37 | 580 076 nt, circulaire | 4 | `b9bd2042…` |

- **Source** : NCBI RefSeq (`nuccore`, format `gbwithparts`), téléchargés le 25/09/2026 avec
  `orfeval fetch`.
- **Vérifier** : `cd data/demo && sha256sum -c SHA256SUMS`.
- **Retélécharger et comparer** : `make data`.

Provenance détaillée, conditions d'utilisation et formats acceptés :
[../docs/data.md](../docs/data.md).
