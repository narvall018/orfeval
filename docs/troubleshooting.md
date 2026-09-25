# Dépannage

### `make install` : « Python >= 3.11 requis »

Votre `python3` est trop ancien, par exemple Anaconda en 3.9. Indiquez un interpréteur
explicitement :

```bash
make install PYTHON=/usr/bin/python3.12
```

Sous Ubuntu 24.04, `python3.12` est déjà installé. Sinon, `pyenv`, `uv python install 3.12` ou
un environnement conda (`conda create -n orfeval python=3.12`) conviennent.

### `orfeval: command not found`

L'environnement virtuel n'est pas activé : `source .venv/bin/activate`, ou appelez
`.venv/bin/orfeval`.

### `orfeval demo` : « config/demo.yaml introuvable »

La démo utilise les fichiers du dépôt : lancez la commande depuis le dossier du dépôt ou l'un
de ses sous-dossiers. Avec une installation hors dépôt, utilisez
`orfeval analyze <votre_config.yaml>`.

### « ne contient pas la séquence nucléotidique »

Le fichier GenBank est un enregistrement de type `CON`, qui liste des contigs sans leur
séquence. Retéléchargez-le avec `orfeval fetch <accession>`, qui demande le format
`gbwithparts`.

### « contient N séquences »

La version 0.1 analyse une séquence par fichier. Séparez les réplicons, par exemple avec
Biopython :

```python
from Bio import SeqIO

for record in SeqIO.parse("assemblage.gbff", "genbank"):
    SeqIO.write(record, f"{record.id}.gb", "genbank")
```

### Résultats très mauvais sur un génome

Vérifiez d'abord le **code génétique** : `orfeval validate` affiche le code déclaré. Un
mycoplasme analysé avec le code 11 perd plus de la moitié de sa précision (voir
`docs/results.md`). Pour un FASTA, précisez `--table`. Vérifiez aussi la **topologie** : un
génome circulaire traité comme linéaire perd les gènes qui traversent l'origine.

### Docker : `permission denied while trying to connect to the docker API`

Votre utilisateur n'appartient pas au groupe `docker`. Deux solutions :
`sudo docker …`, ou `sudo usermod -aG docker $USER` puis déconnexion et reconnexion. Cette
seconde option donne des droits équivalents à root : à réserver à une machine personnelle.
Sans Docker, l'installation Python suffit : la CI construit et teste l'image à chaque push.

### Le dashboard ne s'ouvre pas

- Port déjà utilisé : `make run PORT=8502`.
- Sur une machine distante : `streamlit run app/streamlit_app.py --server.address 0.0.0.0`,
  puis ouvrir le port. Le dashboard n'a pas d'authentification : ne l'exposez pas sur
  Internet.

### `findfont` ou avertissements de police matplotlib

Ils sont sans conséquence sur les résultats. Les figures utilisent DejaVu Sans, fournie avec
matplotlib.

### `pre-commit` : `mypy: command not found`

Le hook mypy utilise l'environnement du projet : lancez `make pre-commit`, ou activez `.venv`
avant `pre-commit run`.
