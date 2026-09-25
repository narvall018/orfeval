# Installation

## Prérequis

| Outil | Version | Remarque |
|---|---|---|
| Python | **3.11, 3.12 ou 3.13** | testé sur les trois en CI |
| Git | toute version récente | pour cloner le dépôt |
| make | optionnel | raccourcis `make install`, `make demo`… |
| Docker | optionnel | exécution sans installer Python |

Aucun outil bioinformatique externe n'est nécessaire : tout repose sur des paquets Python
installables avec `pip`, dont Biopython.

> **Attention aux Python trop anciens.** Une distribution Anaconda ancienne peut fournir
> `python3` en version 3.9, ce qui ne suffit pas. Le Makefile cherche automatiquement
> `python3.12`, puis `python3.13`, `python3.11` et enfin `python3`, et s'arrête avec un message
> clair si la version est insuffisante. Vérifiez avec `python3 --version`.

## Installation recommandée (Makefile)

```bash
git clone https://github.com/narvall018/orfeval.git
cd orfeval
make install                 # crée .venv et installe orfeval + dashboard + outils de dev
source .venv/bin/activate
orfeval --version
```

## Installation manuelle (pip)

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e ".[app,dev]"   # ou simplement "pip install -e ." pour la CLI seule
```

Extras disponibles :

| Extra | Contenu |
|---|---|
| *(aucun)* | CLI, pipeline, figures, rapport HTML |
| `app` | dashboard Streamlit et carte interactive Plotly |
| `dev` | pytest, ruff, mypy, pre-commit, stubs de typage |

## Versions exactes (reproductibilité)

`requirements.txt` fige les versions testées (Python 3.12, Linux) :

```bash
make install-locked
# ou : pip install -r requirements.txt && pip install --no-deps -e .
```

## Docker

```bash
docker build -t orfeval .
docker run --rm -v "$PWD/results/docker:/app/results" orfeval     # lance la démo
```

Sous Linux, si `docker` renvoie `permission denied … docker.sock`, votre utilisateur n'a pas
accès au démon Docker. Voir [troubleshooting.md](troubleshooting.md).

## Vérifier l'installation

```bash
orfeval validate data/demo/NC_001416.1.gb.gz   # lecture d'un génome
make test                                       # suite de tests (~30 s)
```
