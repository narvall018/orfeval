# Commandes courantes d'orfeval. `make help` liste les cibles.
# Python >= 3.11 requis : on privilégie python3.12, puis 3.13, 3.11, puis python3.

PYTHON ?= $(shell command -v python3.12 || command -v python3.13 || command -v python3.11 || command -v python3)
VENV   ?= .venv
BIN    := $(VENV)/bin
PORT   ?= 8501
IMAGE  ?= orfeval:latest

.DEFAULT_GOAL := help
.PHONY: help venv install install-locked demo test test-fast coverage lint format typecheck check \
        pre-commit run report figures data docker-build docker-demo docker-app clean

help: ## Affiche cette aide
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

$(BIN)/python:
	@$(PYTHON) -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else "Python >= 3.11 requis (trouvé : " + sys.version.split()[0] + ")")'
	$(PYTHON) -m venv $(VENV)
	$(BIN)/python -m pip install --quiet --upgrade pip

venv: $(BIN)/python ## Crée l'environnement virtuel (.venv)

install: venv ## Installe orfeval + dashboard + outils de développement
	$(BIN)/pip install --quiet -e ".[app,dev]"
	@echo "OK : activez l'environnement avec 'source $(VENV)/bin/activate'"

install-locked: venv ## Installe les versions exactes testées (requirements.txt)
	$(BIN)/pip install --quiet -r requirements.txt
	$(BIN)/pip install --quiet --no-deps -e .
	$(BIN)/pip install --quiet -e ".[dev]"

demo: ## Lance la démonstration complète (résultats dans results/demo)
	$(BIN)/orfeval demo

test: ## Lance tous les tests
	$(BIN)/pytest

test-fast: ## Lance uniquement les tests unitaires
	$(BIN)/pytest -m "not integration"

coverage: ## Tests avec rapport de couverture
	$(BIN)/pytest --cov=orfeval --cov-report=term-missing

lint: ## Vérifie le style (ruff) et les types (mypy)
	$(BIN)/ruff check src tests app scripts
	$(BIN)/ruff format --check src tests app scripts
	$(BIN)/mypy

check: lint test ## Lint, types et tests : à lancer avant un commit ou une PR

format: ## Reformate le code
	$(BIN)/ruff format src tests app scripts
	$(BIN)/ruff check --fix src tests app scripts

typecheck: ## Vérifie uniquement les types
	$(BIN)/mypy

pre-commit: ## Exécute tous les hooks pre-commit
	PATH="$(abspath $(BIN)):$$PATH" $(BIN)/pre-commit run --all-files

run: ## Lance le dashboard Streamlit (http://localhost:8501)
	$(BIN)/streamlit run app/streamlit_app.py --server.port $(PORT)

report: ## Régénère le rapport HTML de la démo sans recalcul
	$(BIN)/orfeval report results/demo

figures: demo ## Met à jour les figures du README (docs/images)
	$(BIN)/python scripts/export_readme_figures.py

data: ## Retélécharge les génomes de démonstration depuis le NCBI
	PATH="$(abspath $(BIN)):$$PATH" bash scripts/refresh_demo_data.sh

docker-build: ## Construit l'image Docker
	docker build -t $(IMAGE) .

docker-demo: ## Lance la démo dans Docker (résultats dans results/docker)
	mkdir -p results/docker
	docker run --rm -v "$(CURDIR)/results/docker:/app/results" $(IMAGE)

docker-app: ## Lance le dashboard dans Docker (http://localhost:8501)
	docker run --rm -p $(PORT):8501 --entrypoint streamlit $(IMAGE) \
		run app/streamlit_app.py --server.address 0.0.0.0 --server.headless true

clean: ## Supprime caches, builds et résultats générés
	rm -rf build dist .pytest_cache .mypy_cache .ruff_cache .coverage coverage.xml htmlcov
	find . -path ./$(VENV) -prune -o -name "__pycache__" -type d -print -exec rm -rf {} +
	find results -mindepth 1 ! -name .gitkeep -exec rm -rf {} +
