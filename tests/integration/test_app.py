"""Le dashboard Streamlit s'exécute sans erreur sur les génomes de démonstration."""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
pytest.importorskip("plotly")
from streamlit.testing.v1 import AppTest

pytestmark = pytest.mark.integration
APP = Path(__file__).resolve().parents[2] / "app" / "streamlit_app.py"


def test_dashboard_on_lambda():
    app = AppTest.from_file(str(APP), default_timeout=180)
    app.run()
    assert not app.exception, app.exception
    assert len(app.tabs) == 6
    labels = [metric.label for metric in app.metric]
    assert labels[:3] == ["Gènes prédits", "Sensibilité", "Précision"]


def test_dashboard_switches_genome_and_parameters():
    app = AppTest.from_file(str(APP), default_timeout=180)
    app.run()
    app.sidebar.selectbox[0].select_index(1).run()  # M. genitalium
    assert not app.exception, app.exception
    app.sidebar.selectbox[2].set_value("rbs").run()  # stratégie de choix du start
    assert not app.exception, app.exception
    assert any("NC_000908" in caption.value for caption in app.caption)
