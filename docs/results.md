# Résultats de la démonstration

Les résultats ci-dessous sont ceux de `orfeval demo`, avec orfeval 0.1.0, Python 3.12.3,
Biopython 1.88, numpy 2.5.3 et pandas 3.0.6, exécuté le 25/09/2026. Ils sont
**déterministes** : la même commande sur les mêmes fichiers donne les mêmes nombres. Le rapport
HTML complet est produit dans `results/demo/report.html`, et publié comme artefact de chaque
exécution de la CI.

Chaque analyse est lue à trois niveaux : **technique** (ce que le code a produit),
**statistique** (ce que mesurent les métriques) et **biologique** (interprétation prudente).

## Vue d'ensemble

| Analyse | Code | Gènes prédits | Sensibilité [IC 95 %] | Précision [IC 95 %] | F1 | Start exact |
|---|---|---|---|---|---|---|
| `lambda` | 11 | 68 | 0,781 [0,673–0,860] | 0,809 [0,700–0,885] | 0,795 | 0,737 |
| `mgenitalium` | 4 | 564 | 0,970 [0,951–0,982] | 0,867 [0,837–0,893] | 0,916 | 0,824 |
| `mgenitalium_code11` | 11 *(forcé)* | 1 195 | 0,687 [0,645–0,725] | 0,290 [0,265–0,316] | 0,407 | 0,347 |
| `mgenitalium_rbs` | 4 | 791 | 0,982 [0,966–0,991] | 0,626 [0,592–0,659] | 0,764 | 0,289 |

![Comparaison des analyses](images/comparison.png)

## Phage lambda (NC_001416.1)

**Technique.** 502 ORF d'au moins 90 nt ont été trouvés. Le seuil d'auto-apprentissage, tiré
du modèle nul, vaut 528 nt : 43 ORF d'entraînement, puis 29 après ré-entraînement. 68 gènes
sont prédits en 0,1 s.

**Statistique.**

- **Sensibilité et précision.** 57 des 73 CDS annotées sont retrouvées, et 55 des
  68 prédictions sont appariées. Les deux nombres diffèrent parce que des gènes emboîtés
  partagent un même codon stop.
- **Apport du modèle codant.** Son classement atteint une aire sous la courbe
  précision-rappel (AP) de 0,821, contre 0,682 pour la longueur seule. À nombre de
  prédictions égal (68), la longueur seule obtient un F1 de 0,709 contre 0,795.
- **Gènes manqués** (16) : 10 sont éliminés par la règle de chevauchement, 5 ont un score
  sous le seuil, 1 est plus court que 90 nt.
- **Prédictions sans correspondance** (13) : 10 tombent hors des gènes annotés, 3 sont en
  antisens d'un gène annoté.
- **Effet de la longueur.** Le rappel est de 0,55 pour les gènes de 150 à 299 nt
  (IC : 0,34–0,74), contre 0,86 à 0,90 au-delà de 300 nt.

![Précision-rappel, lambda](images/lambda_precision_recall.png)

**Biologique (prudent).** Le génome de lambda est très dense : les prédictions couvrent
82,6 % de sa longueur, et de nombreux gènes se chevauchent. La règle des 60 nt est la
principale cause de gènes manqués : c'est une limite de la méthode, pas une propriété de
lambda. En amont des starts annotés, le profil nucléotidique montre un pic de 0,39 bit à la
position −10, avec un consensus riche en G (`…GAAGGGGG…`). Ce signal est compatible avec un
motif de Shine-Dalgarno.

![Profil amont, lambda](images/lambda_upstream.png)

## *Mycoplasmoides* (*Mycoplasma*) *genitalium* G37 (NC_000908.2)

**Technique.** L'annotation déclare le code génétique 4 pour ses 524 CDS, et orfeval l'utilise
automatiquement. 20 pseudogènes sont exclus de la référence, qui compte donc 504 CDS. On
obtient 4 664 ORF, un seuil d'auto-apprentissage de 561 nt, 405 puis 374 séquences
d'entraînement, et 564 gènes prédits en 0,8 s.

**Statistique.**

- **Sensibilité** : 0,970, soit 489 gènes sur 504. Tous les gènes d'au moins 1 200 nt sont
  retrouvés (152 sur 152).
- **Précision** : 0,867.
- **Start exact** : 0,824 des gènes retrouvés.
- **Modèle codant contre longueur seule** : AP de 0,963 contre 0,918.
- **Gènes manqués** (15) : 13 par chevauchement, 2 sous le seuil.
- **Prédictions sans correspondance** (75), réparties ainsi :
  - 22 recouvrent un pseudogène exclu de la référence ;
  - 7 recouvrent un gène d'ARN ;
  - 16 sont en antisens d'un gène annoté ;
  - 1 est sur le même brin qu'un gène annoté ;
  - 29 tombent hors des gènes annotés.
- **Usage des codons** : le RSCU des gènes prédits et celui des gènes annotés sont presque
  identiques (r = 0,9999).

![Carte des gènes, M. genitalium](images/mgenitalium_genome_map.png)

![Longueur des ORF face au modèle nul](images/mgenitalium_orf_length_null.png)

**Biologique (prudent).**

- **Faux positifs apparents.** Près d'un tiers des prédictions « fausses » (22 sur 75)
  recouvrent des pseudogènes. Ces régions ont gardé une composition codante, ce qui est
  plausible pour des gènes inactivés (décalage de cadre, troncature). On ne peut pas conclure,
  à partir de ces seules données, que les 29 prédictions hors annotation soient de vrais gènes.
- **Réplication.** Le minimum du GC skew cumulé se situe à environ 1,5 kb, à 2 946 nt du gène
  *dnaA* annoté, compte tenu de la circularité. Le maximum est vers 294 kb, soit à peu près à
  l'opposé sur ce génome de 580 kb. C'est cohérent avec le schéma attendu pour un chromosome
  circulaire à réplication bidirectionnelle, mais cela ne localise pas expérimentalement
  l'origine.
- **Signal amont.** Il est faible : 0,10 bit au maximum, contre 0,39 bit chez lambda. Cela
  explique en partie l'échec de la stratégie `rbs` sur ce génome (voir plus bas).

![GC skew, M. genitalium](images/mgenitalium_gc_skew.png)

![Profil amont, M. genitalium](images/mgenitalium_upstream.png)

## Expérience de contrôle : mauvais code génétique (`mgenitalium_code11`)

Le même génome est analysé en forçant le code 11, où TGA est un codon stop.

**Technique.** La probabilité de stop passe de 0,058 à 0,077. Les ORF sont plus nombreux et
plus courts : la longueur moyenne des gènes prédits est de 396 nt, contre 986 nt avec le
code 4.

**Statistique.**

- La précision chute de 0,867 à 0,290.
- La sensibilité passe de 0,970 à 0,687.
- 649 prédictions sans correspondance sont sur le même brin qu'un gène annoté : ce sont des
  fragments de gènes coupés par TGA.
- 144 gènes n'ont plus aucun ORF candidat qui se termine à leur codon stop.

![Écarts avec le mauvais code génétique](images/mgenitalium_code11_errors.png)

**Biologique.** Ce résultat illustre, sur des données réelles, la réassignation du codon UGA
en tryptophane chez les mycoplasmes (Yamao et al., 1985). Sur ce génome, la précision de la
prédiction dépend directement du bon code génétique. orfeval lit ce code dans l'annotation
quand il y est déclaré, et le rapport en indique la provenance.

## Analyse de sensibilité : codon start guidé par le RBS (`mgenitalium_rbs`)

La sensibilité ne bouge presque pas (0,982). En revanche, l'exactitude du start tombe à
0,289 et la précision à 0,626, avec 211 prédictions antisens. Les starts sont placés trop en
aval, et les gènes raccourcis laissent de la place à des ORF antisens. Sur lambda, la même
stratégie est comparable à `longest` (voir `docs/methodology.md`, section 11). La valeur par
défaut `longest` a été choisie à la suite de cette comparaison.

## Au-delà de la démonstration : *E. coli* K-12 MG1655 (NC_000913.3)

Ce génome **n'a servi à choisir aucun paramètre**. Il fait 4,6 Mb, soit huit fois la taille de
*M. genitalium*, et n'est pas versionné dans le dépôt (3 Mo compressés). Pour reproduire :

```bash
orfeval fetch NC_000913.3 --outdir data/raw
orfeval predict data/raw/NC_000913.3.gb.gz -o results/ecoli                          # longest
orfeval predict data/raw/NC_000913.3.gb.gz -o results/ecoli_rbs --start-strategy rbs
```

| Stratégie | Gènes | Sensibilité [IC 95 %] | Précision [IC 95 %] | Start exact [IC 95 %] |
|---|---|---|---|---|
| `longest` (défaut) | 4 850 | 0,914 [0,906–0,922] | 0,809 [0,798–0,820] | 0,649 [0,634–0,664] |
| `rbs` | 4 976 | 0,919 [0,911–0,927] | 0,792 [0,781–0,803] | 0,752 [0,738–0,765] |

**Technique.** L'analyse complète (prédiction, évaluation et figures) prend 15 s et 370 Mo de
mémoire sur la machine de développement.

**Statistique.**

- **Référence.** 4 300 CDS sont utilisées, sur 4 318 annotées (18 pseudogènes exclus). Le
  modèle codant atteint une AP de 0,920, contre 0,861 pour la longueur seule.
- **Rappel selon la longueur.** Il est de 0,30 sous 150 nt (47 gènes sur 157), 0,78 entre 150
  et 299 nt, et au moins 0,95 au-delà de 600 nt.
- **Prédictions sans correspondance.** Sur 927, 591 tombent hors des gènes annotés.
- **Choix du start.** Ici, la stratégie `rbs` améliore nettement l'exactitude du start
  (+0,10, intervalles disjoints), au prix d'une précision légèrement plus faible.

**Biologique (prudent).**

- **Signal amont.** En amont des starts annotés, il atteint 0,45 bit à −9, contre 0,10 bit
  chez *M. genitalium*. Le motif de Shine-Dalgarno est donc informatif pour *E. coli* : c'est
  cohérent avec le gain de la stratégie `rbs` sur ce génome, et avec sa perte sur
  *M. genitalium*. La meilleure stratégie de start **dépend du génome**. Une piste
  d'amélioration, inscrite dans la roadmap, est de choisir la stratégie d'après la force du
  signal amont mesuré sur l'ensemble d'auto-apprentissage.
- **Réplication.** Le minimum du GC skew cumulé se trouve à 25 kb du gène *dnaA* annoté.

## Ce que ces résultats ne montrent pas

- Ils ne montrent pas une performance générale : trois génomes, dont un phage, ne
  représentent pas la diversité procaryote.
- Ils ne mettent pas orfeval en concurrence avec Prodigal ou GeneMarkS : aucune comparaison
  n'a été faite dans cette version (c'est dans la roadmap).
- Ils ne valident pas les prédictions hors annotation : il faudrait des données
  indépendantes (homologie, protéomique, transcriptomique).
