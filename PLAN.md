# PLAN — orfeval

Feuille de route du projet : ce document fixe les objectifs, les choix de conception
et les étapes de construction. Il est mis à jour au fil du développement.

## 1. Question scientifique

> Peut-on retrouver les gènes codants d'un génome procaryote **à partir de sa seule
> séquence**, et comment mesurer honnêtement la qualité de cette prédiction ?

- Le prédicteur ne voit **jamais** l'annotation : il s'entraîne sur le génome lui-même
  (auto-apprentissage, principe de GeneMark-S / Prodigal, en version simplifiée).
- L'annotation GenBank/RefSeq sert **uniquement** à l'évaluation, après la prédiction.
- L'annotation de référence est elle-même en partie issue de prédictions : c'est une
  référence imparfaite, et cette limite est documentée.

## 2. Données

| Génome | Accession | Taille | Topologie | Code génétique | Rôle |
|---|---|---|---|---|---|
| Phage lambda | NC_001416.1 | 48,5 kb | linéaire | 11 | petit génome dense, gènes chevauchants |
| *Mycoplasmoides* (*Mycoplasma*) *genitalium* G37 | NC_000908.2 | 580 kb | circulaire | **4** (TGA = Trp) | génome bactérien minimal, riche en AT |

- Fichiers GenBank compressés (< 1 Mo au total) versionnés dans `data/demo/`
  avec provenance, date de téléchargement et sommes SHA-256.
- Commande `orfeval fetch` (Biopython `Entrez`) pour retélécharger ou ajouter des génomes.

## 3. Méthode (résumé)

1. **ORF** : recherche vectorisée (numpy) sur les 6 cadres, code génétique paramétrable,
   gestion des génomes circulaires (ORF à cheval sur l'origine).
2. **Modèle nul** : longueur des ORF dans une séquence aléatoire de même composition
   (loi géométrique) → seuil de longueur « significatif » pour l'auto-apprentissage.
3. **Modèle codant** : usage des codons (ou dicodons) appris sur les longs ORF,
   comparé à un modèle de fond conditionné à l'absence de codon stop → score en bits.
4. **Codon start** : ORF le plus long, ou score combinant potentiel codant et
   motif de Shine-Dalgarno.
5. **Chevauchements** : sélection gloutonne par score avec chevauchement maximal toléré.
6. **Auto-apprentissage itératif** : ré-entraînement sur les gènes prédits.
7. **Évaluation** : appariement par extrémité 3′ (codon stop), sensibilité, précision,
   F1, exactitude du codon start, intervalles de Wilson, rappel par classe de longueur,
   courbe précision-rappel (score vs longueur seule), niveau nucléotidique.
8. **Descripteurs du génome** : GC, GC skew et biais cumulé, usage des codons (RSCU),
   GC3, profil nucléotidique en amont des codons start (`Bio.motifs`).

Expérience de contrôle : *M. genitalium* analysé avec le code 11 au lieu du code 4,
pour mesurer l'effet d'un mauvais code génétique.

## 4. Architecture

```
src/orfeval/
├── cli.py              # Typer : fetch | validate | predict | analyze | report | demo
├── config.py           # modèles pydantic + chargement YAML
├── io/                 # GenBank/FASTA, écriture TSV/GFF3/FASTA/JSON, NCBI
├── orfs/               # finder, null_model, coding_model, rbs, selection, predictor
├── evaluation/         # appariement, métriques
├── features/           # composition, usage des codons, profil amont
├── plotting/           # figures matplotlib (aucun calcul métier)
├── report/             # rapport HTML autonome (Jinja2)
└── pipeline.py         # orchestration d'une analyse
app/                    # dashboard Streamlit (appelle uniquement l'API du package)
```

## 5. Étapes

- [x] Inspection de l'environnement (Python 3.12, réseau, GitHub)
- [x] Choix du sujet et des données de démonstration
- [ ] Squelette : pyproject, configuration qualité, arborescence
- [ ] Entrées/sorties et validation des génomes
- [ ] Recherche d'ORF + modèle nul
- [ ] Modèle codant, RBS, sélection, prédicteur auto-entraîné
- [ ] Évaluation
- [ ] Descripteurs du génome
- [ ] Figures, pipeline, écriture des résultats
- [ ] Rapport HTML
- [ ] CLI
- [ ] Tests unitaires et d'intégration (exécutés)
- [ ] Dashboard Streamlit
- [ ] Docker, Makefile, CI GitHub Actions, pre-commit
- [ ] Documentation (`docs/`) et README avec les résultats réels de la démo
- [ ] Vérification finale (lint, typage, tests, démo) puis publication GitHub

## 6. Règles scientifiques du projet

- Aucun résultat n'est écrit à la main : les chiffres de la documentation proviennent
  d'exécutions réelles de `orfeval demo`, avec la version et la date.
- Trois niveaux distingués partout : **technique** (le code fonctionne),
  **statistique** (ce que mesurent les métriques), **biologique** (interprétation prudente).
- Les paramètres par défaut sont fixés *a priori* ; toute comparaison faite sur les
  génomes de démonstration est présentée comme telle (risque de sur-ajustement).
