# Contribuer à orfeval

Merci de votre intérêt ! Ce guide décrit comment proposer une modification.

## Environnement de développement

```bash
git clone https://github.com/narvall018/orfeval.git
cd orfeval
make install              # .venv + orfeval en mode éditable + outils de dev
source .venv/bin/activate
pre-commit install        # vérifications automatiques à chaque commit
make check                # lint, types et tests : tout doit passer avant de commencer
```

## Déroulement

1. Ouvrez une *issue* pour discuter d'un changement important avant de coder.
2. Créez une branche depuis `main` : `git switch -c feat/strategie-start-adaptative`.
3. Faites des commits petits et explicites, en français à l'impératif ou au nominatif
   (« Ajout du modèle de start appris »).
4. Vérifiez localement : `make check`, et `orfeval demo` si le comportement
   scientifique change.
5. Ouvrez une *pull request* en remplissant le modèle. La CI doit être verte.

## Conventions

| Sujet | Règle |
|---|---|
| Langue | documentation, docstrings, messages et commentaires en **français** ; identifiants en anglais |
| Style | `ruff format` (100 colonnes), `ruff check` sans erreur |
| Typage | annotations partout dans `src/` ; `mypy` sans erreur |
| Docstrings | style NumPy (`Parameters`, `Returns`, `Raises`) pour les fonctions publiques |
| Erreurs | lever une sous-classe de `OrfevalError` pour une erreur attendue (entrée, configuration) |
| Journal | `logger = get_logger(__name__)`, jamais de `print` dans le package (sauf la CLI) |
| Chemins | `pathlib.Path` exclusivement |
| Architecture | aucune logique scientifique dans `cli.py` ni `app/` ; les figures ne calculent rien |

## Tests

- Toute correction de bug s'accompagne d'un test qui échouait avant.
- Les tests unitaires utilisent des séquences construites à la main ; les tests
  d'intégration sont marqués `@pytest.mark.integration` (`make test-fast` les exclut).
- Pas de réseau dans les tests : simulez `Bio.Entrez` avec `monkeypatch`, comme dans
  `tests/unit/test_ncbi_and_fallbacks.py`.
- Si une modification change les résultats de la démonstration, mettez à jour
  `docs/results.md`, les figures (`make figures`) et le `CHANGELOG.md`, puis expliquez
  pourquoi dans la PR.

## Exemple : ajouter une stratégie de choix du start

1. `config.py` : ajouter la valeur au `Literal` de `PredictorConfig.start_strategy`.
2. `orfs/selection.py` : l'ajouter à `StartStrategy` et à `preferred_start()`.
3. `tests/unit/test_selection_and_matching.py` : tester le choix sur un cas construit.
4. `config/sensitivity.yaml` : ajouter une analyse pour la comparer aux autres.
5. `docs/methodology.md` : décrire la stratégie et ses hypothèses.

L'évaluation, le rapport et le dashboard la prennent en charge sans autre modification.

## Règles scientifiques

- Ne jamais présenter un résultat qui n'a pas été produit par le code du dépôt.
- Ne pas optimiser un paramètre sur les génomes de démonstration sans le signaler.
- Distinguer ce qui est mesuré (statistique) de ce qui est interprété (biologique).
