# Méthodologie

Ce document décrit ce que fait orfeval, pourquoi, et avec quelles hypothèses. Les chiffres
cités proviennent de `orfeval demo` et `orfeval analyze config/sensitivity.yaml`
(orfeval 0.1.0, 25/09/2026) ; ils sont reproductibles avec les fichiers du dépôt.

## 1. Question et principe

> Retrouver les gènes codant des protéines d'un génome procaryote **à partir de sa seule
> séquence**, puis mesurer honnêtement l'accord avec une annotation de référence.

Deux règles structurent le projet :

1. **Le prédicteur ne voit jamais l'annotation.** `GenePredictor.predict()` ne reçoit qu'un
   objet `Genome` (séquence, topologie) : c'est une garantie de conception, vérifiée par un
   test. Le modèle est appris sur le génome lui-même (*auto-apprentissage*), comme dans
   GeneMarkS ou Prodigal, en version volontairement simplifiée.
2. **L'annotation sert uniquement à l'évaluation**, après coup. Elle est traitée comme une
   référence imparfaite, et non comme une vérité absolue.

## 2. Vocabulaire

| Terme | Définition utilisée ici |
|---|---|
| Codon stop | Triplet qui termine la traduction ; dépend du **code génétique** (TGA est un stop dans le code 11, mais code le tryptophane dans le code 4). |
| ORF | *Open Reading Frame* : région d'un cadre de lecture entre deux codons stop, tronquée au codon start le plus en amont. |
| CDS | Séquence codante annotée dans le fichier GenBank. |
| Extrémité 3′ | Position du codon stop : c'est la partie d'un gène la plus facile à identifier. |
| Code génétique | Table NCBI (`transl_table`) ; 11 pour la plupart des bactéries, 4 pour les mycoplasmes. |

## 3. Recherche des ORF (`orfs/finder.py`)

- Les six cadres (3 par brin) sont parcourus en encodant chaque codon par un entier (0–63, et
  64 pour un codon ambigu). Les stops et les starts sont repérés par des masques numpy, puis
  chaque région entre deux stops est associée à ses starts par `searchsorted` :
  le coût est linéaire dans la taille du génome (moins d'une seconde pour 580 kb).
- Les starts retenus par défaut sont **ATG, GTG, TTG**. orfeval vérifie qu'ils sont des codons
  d'initiation valides pour le code génétique choisi.
- **Génomes circulaires** : chaque brin est concaténé à lui-même, et seuls les stops de la
  seconde copie sont retenus. Chaque stop est ainsi traité une fois, avec tout son contexte
  amont, y compris quand un gène traverse l'origine. Un ORF ne peut pas dépasser la longueur
  du génome.
- Longueur minimale : **90 nt**, codon stop inclus.

## 4. Modèle nul de la longueur des ORF (`orfs/null_model.py`)

Dans une séquence aléatoire (nucléotides indépendants, composition du génome), un codon est un
stop avec une probabilité `p`. On utilise une composition **symétrique** entre les brins
(`p(A) = p(T)`, `p(C) = p(G)`), puisque les deux brins sont analysés. Le nombre `L` de codons
sens entre deux stops suit alors une loi géométrique :

```
P(L ≥ k) = (1 − p)^k        nombre attendu de régions ≥ k  ≈  2·N·p · (1 − p)^k
```

Le **seuil d'auto-apprentissage** est la plus petite longueur pour laquelle on attend moins
d'**un** ORF aléatoire dans tout le génome. Les ORF plus longs sont très probablement codants :
ils forment l'ensemble d'entraînement initial, choisi **sans l'annotation**.

| Génome | p(stop) | Seuil obtenu | ORF d'entraînement |
|---|---|---|---|
| Lambda (code 11, GC 49,9 %) | 0,0471 | 528 nt | 43 |
| *M. genitalium* (code 4, GC 31,7 %) | 0,0583 | 561 nt | 405 |
| *M. genitalium* (code 11, contrôle) | 0,0768 | 432 nt | 376 |

S'il y a moins de 20 ORF au-delà du seuil (petit génome), les 20 ORF les plus longs sont
utilisés, et un avertissement est émis.

*Hypothèse :* un génome réel n'est pas une suite de nucléotides indépendants (biais de
dinucléotides, répétitions). Le modèle nul est une référence, pas une description exacte.
La figure `orf_length_null` compare les deux distributions.

## 5. Modèle codant (`orfs/coding_model.py`)

Le score d'un candidat est un **log-rapport de vraisemblance** en bits :

```
S = Σ_i log2( P_codant(c_i | c_i−1) / P_fond(c_i) )
```

- `P_codant` : fréquences des codons (modèle `codon`, ordre 0) ou des codons conditionnées au
  précédent (modèle `dicodon`, ordre 1), avec pseudo-comptes (1 par défaut).
- `P_fond` : nucléotides indépendants de même composition, **renormalisés sur les codons
  sens**. Un ORF ne contient par définition aucun stop interne : sans ce conditionnement, tout
  ORF long gagnerait environ −log2(1 − p) bit par codon, et paraîtrait codant à tort.
- Le codon start et le codon stop ne sont pas comptés, car leur distribution est contrainte
  par la définition même d'un ORF.
- Pour un ORF, les scores de tous les starts possibles s'obtiennent en une passe, par sommes
  cumulées inversées.

## 6. Détection et placement du codon start (`orfs/selection.py`)

Deux décisions distinctes sont prises pour chaque ORF :

1. **Détection** : le score de l'ORF est le **maximum** du score codant sur ses starts
   possibles. Il mesure la présence d'un signal codant, indépendamment du start.
2. **Placement du start**, selon trois stratégies :
   - `longest` (défaut) : le start le plus en amont ;
   - `rbs` : le start le plus en amont précédé d'un motif de Shine-Dalgarno, c'est-à-dire une
     sous-chaîne d'au moins 3 nt de `AGGAGG`, située à 3–15 nt du start ;
   - `score` : le start qui maximise score codant + bonus RBS (1 bit par nucléotide du motif).

Séparer ces deux décisions a été nécessaire. Dans une première version, le start et la
détection partageaient le même score : la stratégie `longest` perdait alors des gènes, parce
qu'un ORF prolongé vers l'amont chevauche davantage ses voisins.

## 7. Sélection gloutonne et chevauchements

Les ORF sont examinés par score de détection décroissant. Un candidat est accepté si ses
positions déjà couvertes par des gènes acceptés, tous brins confondus, ne dépassent pas
`max_overlap` (**60 nt**). Sinon, le start est décalé vers l'aval, comme dans Glimmer,
tant que le gène garde la longueur minimale.

Cette sélection a une propriété utile : l'acceptation d'un candidat ne dépend que des
candidats mieux classés. Appliquer le seuil de score (0 bit) **après** la sélection donne donc
le même ensemble que l'appliquer avant. C'est ce qui permet de tracer toute la courbe
précision-rappel en une seule passe. La propriété est vérifiée par un test sur 200 candidats
aléatoires.

## 8. Auto-apprentissage itératif

1. Entraînement sur les longs ORF (section 4).
2. Prédiction.
3. Ré-entraînement sur les gènes prédits d'au moins la même longueur, ce qui élimine en
   particulier les ORF « ombres » situés en antisens de vrais gènes.
4. Nouvelle prédiction. Deux itérations sont faites par défaut.

Le prédicteur est **déterministe** : aucun tirage aléatoire, et un ordre de tri total en cas
d'égalité.

## 9. Évaluation (`evaluation/`)

**Référence.** Les CDS du GenBank, en excluant :

- les pseudogènes (`/pseudo`, `/pseudogene`) : 20 chez *M. genitalium* ;
- les CDS partielles (positions `<` ou `>`) ;
- les localisations non prises en charge.

Les CDS qui traversent l'origine (`join(…,1..x)`) sont reconstituées. Les éléments exclus et
les gènes d'ARN sont conservés comme contexte d'interprétation.

**Appariement.** Une prédiction est correcte si elle partage le brin et la position du codon
stop d'une CDS annotée. C'est le critère courant pour les gènes procaryotes (voir par exemple
Hyatt et al., 2010). Le codon start est évalué séparément.

**Métriques.**

| Métrique | Définition |
|---|---|
| Sensibilité | CDS de référence retrouvées / CDS de référence |
| Précision | prédictions appariées / prédictions |
| F1 | moyenne harmonique des deux |
| Start exact | CDS retrouvées avec le bon start / CDS retrouvées |
| Niveau nucléotidique | bases codantes (brin compris) partagées / annotées, et / prédites |
| AP | aire sous la courbe précision-rappel (somme des gains de rappel × précision) |

- Toutes les proportions ont un **intervalle de confiance à 95 % de Wilson** (Wilson, 1927),
  mieux adapté que l'approximation normale près de 0 ou 1 et sur de petits effectifs.
- **Référence naïve** : les mêmes ORF classés par longueur seule, avec la même règle de
  chevauchement. L'écart entre les deux courbes mesure l'apport du modèle codant.
- **Analyse des écarts.** Chaque prédiction non appariée est classée selon ce qu'elle recouvre
  à au moins 50 % : pseudogène ou CDS partielle, gène d'ARN, gène annoté sur le même brin, sur
  le brin opposé (antisens), ou région sans annotation. Chaque gène manqué est classé selon la
  cause : stop non canonique, score sous le seuil, élimination par chevauchement, gène trop
  court, ou absence d'ORF candidat.
- **Sous-ensembles atteignables.** Les gènes plus courts que la longueur minimale, ou dont le
  stop n'est pas canonique, ne peuvent pas être retrouvés : une sensibilité restreinte aux
  gènes atteignables est aussi rapportée.

## 10. Descripteurs du génome (`features/`)

- **GC skew** `(G − C)/(G + C)` par fenêtres (`Bio.SeqUtils.GC_skew`), et son **cumul**. Chez de
  nombreux chromosomes bactériens circulaires, le minimum du cumul est proche de l'origine de
  réplication, et le maximum proche du terminus (Lobry, 1996 ; Grigoriev, 1998). La distance au
  gène *dnaA* annoté est rapportée à titre de comparaison. Cet indicateur n'est pas calculé
  pour un génome linéaire.
- **RSCU** (Sharp et al., 1986), selon le code génétique réellement utilisé. Avec le code 4, le
  tryptophane a deux codons, TGG et TGA. Comparer les gènes prédits aux gènes annotés vérifie
  que l'ensemble prédit a la même signature de codons.
- **Profil amont des starts** (`Bio.motifs`) : entropie relative de chaque position, par
  rapport à la composition du génome. Un pic vers −12 à −6 signale un motif de Shine-Dalgarno
  (Shine & Dalgarno, 1974).

## 11. Choix des paramètres et analyse de sensibilité

Les paramètres ont été fixés *a priori* (longueur minimale de 90 nt, chevauchement de 60 nt,
seuil de 0 bit, modèle `codon`), à une exception près, documentée ci-dessous.

**Stratégie de choix du start.** La valeur par défaut prévue était `rbs`. Elle a été remplacée
par `longest` **après** comparaison sur les génomes de démonstration :

| Analyse | Gènes | Sensibilité | Précision | F1 | Start exact | AP |
|---|---|---|---|---|---|---|
| lambda — longest | 68 | 0,781 | 0,809 | 0,795 | 0,737 | 0,821 |
| lambda — rbs | 70 | 0,795 | 0,800 | 0,797 | 0,741 | 0,831 |
| lambda — score | 70 | 0,781 | 0,786 | 0,783 | 0,667 | 0,826 |
| lambda — dicodon | 72 | 0,781 | 0,778 | 0,779 | 0,737 | 0,815 |
| *M. genitalium* — longest | 564 | 0,970 | 0,867 | 0,916 | 0,824 | 0,963 |
| *M. genitalium* — rbs | 791 | 0,982 | 0,626 | 0,765 | 0,289 | 0,969 |
| *M. genitalium* — score | 575 | 0,980 | 0,859 | 0,916 | 0,739 | 0,973 |
| *M. genitalium* — dicodon | 567 | 0,976 | 0,868 | 0,919 | 0,825 | 0,959 |

*Reproduire :* `orfeval analyze config/sensitivity.yaml`.

**Ce que montre ce tableau.**

- *Technique.* Les trois stratégies et les deux modèles fonctionnent sur les deux génomes.
- *Statistique.* Sur lambda, les écarts entre variantes sont faibles au regard des
  intervalles de confiance (73 gènes). Sur *M. genitalium*, la stratégie `rbs` fait chuter
  l'exactitude du start de 0,82 à 0,29 et la précision de 0,87 à 0,63 : les starts sont
  placés trop en aval, et les gènes raccourcis laissent passer des ORF antisens
  (211 prédictions antisens, contre 16 avec `longest`).
- *Biologique (prudent).* Ce résultat est **cohérent** avec le profil amont mesuré sur les
  starts annotés : signal maximal de 0,39 bit à −10 chez lambda, contre 0,10 bit chez
  *M. genitalium*. Le motif de Shine-Dalgarno semble donc peu informatif pour ce génome. Cette
  observation porte sur un seul génome et ne suffit pas à généraliser.

`longest` a été retenu parce qu'il est le plus **robuste** sur les deux génomes. Les
performances de la démonstration peuvent donc être légèrement optimistes : la valeur par
défaut a été choisie en voyant ces données.

**Vérification sur un génome indépendant.** Sur *E. coli* K-12, qui n'a servi à aucun choix,
`rbs` améliore au contraire l'exactitude du start : 0,752 contre 0,649, avec des intervalles
disjoints (voir `docs/results.md`). Le signal amont y est fort (0,45 bit). Le choix optimal
dépend donc du génome : `longest` est un compromis robuste, pas un optimum. Une sélection
automatique de la stratégie, d'après la force du signal amont mesurée sur l'ensemble
d'auto-apprentissage, est inscrite dans la roadmap.

## 12. Limites

- La référence (GenBank/RefSeq) est en partie issue de prédictions automatiques. Une
  prédiction « sans correspondance » peut être un vrai gène non annoté, et inversement.
- Le modèle codant se limite aux codons ou aux dicodons, sans modèle de start entraîné. Les
  outils de référence combinent des modèles plus riches : modèles de Markov interpolés de
  Glimmer, modèles d'ordre élevé de GeneMarkS, statistiques d'hexamères et modèle de site
  d'initiation appris de Prodigal. orfeval est un projet pédagogique et méthodologique : il
  ne prétend pas les remplacer.
- Les gènes courts (< 300 nt) apportent peu de codons, donc peu de signal : sur lambda, le
  rappel est de 0,55 pour 150–299 nt (IC 95 % : 0,34–0,74).
- Les chevauchements supérieurs à 60 nt sont impossibles : c'est la première cause de gènes
  manqués sur les deux génomes (10 sur lambda, 13 chez *M. genitalium*).
- Les décalages de cadre programmés, les sélénoprotéines et les introns ne sont pas modélisés.
- Une seule séquence est analysée par fichier : pas de plasmides, pas d'assemblages
  fragmentés.
- Deux génomes de démonstration ne suffisent pas à estimer une performance générale.

## Références

- Besemer J., Lomsadze A., Borodovsky M. (2001). GeneMarkS: a self-training method for
  prediction of gene starts in microbial genomes. *Nucleic Acids Research* 29(12):2607–2618.
- Cock P.J.A. et al. (2009). Biopython: freely available Python tools for computational
  molecular biology and bioinformatics. *Bioinformatics* 25(11):1422–1423.
- Delcher A.L., Bratke K.A., Powers E.C., Salzberg S.L. (2007). Identifying bacterial genes
  and endosymbiont DNA with Glimmer. *Bioinformatics* 23(6):673–679.
- Fraser C.M. et al. (1995). The minimal gene complement of *Mycoplasma genitalium*.
  *Science* 270(5235):397–403.
- Grigoriev A. (1998). Analyzing genomes with cumulative skew diagrams. *Nucleic Acids
  Research* 26(10):2286–2290.
- Hyatt D. et al. (2010). Prodigal: prokaryotic gene recognition and translation initiation
  site identification. *BMC Bioinformatics* 11:119.
- Lobry J.R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria.
  *Molecular Biology and Evolution* 13(5):660–665.
- Sanger F. et al. (1982). Nucleotide sequence of bacteriophage λ DNA. *Journal of Molecular
  Biology* 162(4):729–773.
- Sharp P.M., Tuohy T.M.F., Mosurski K.R. (1986). Codon usage in yeast: cluster analysis
  clearly differentiates highly and lowly expressed genes. *Nucleic Acids Research*
  14(13):5125–5143.
- Shine J., Dalgarno L. (1974). The 3′-terminal sequence of *Escherichia coli* 16S ribosomal
  RNA: complementarity to nonsense triplets and ribosome binding sites. *PNAS*
  71(4):1342–1346.
- Wilson E.B. (1927). Probable inference, the law of succession, and statistical inference.
  *Journal of the American Statistical Association* 22(158):209–212.
- Yamao F. et al. (1985). UGA is read as tryptophan in *Mycoplasma capricolum*. *PNAS*
  82(8):2306–2309.
