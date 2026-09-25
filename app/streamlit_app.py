"""Dashboard Streamlit d'orfeval.

Couche d'interface uniquement : chargement, paramètres, affichage et téléchargements.
Tous les calculs proviennent du package ``orfeval`` (``analyze_genome``) ; aucune
logique scientifique n'est dupliquée ici.

Lancement : ``streamlit run app/streamlit_app.py`` (ou ``make run``).
"""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from orfeval import __version__
from orfeval.config import AnalysisConfig, PredictorConfig, RunConfig, available_tables
from orfeval.exceptions import OrfevalError
from orfeval.io.genome import LoadedGenome, load_genome
from orfeval.io.writers import gff3_text, proteins_fasta_text
from orfeval.pipeline import RunResult, analyze_genome, build_figures, prediction_table
from orfeval.plotting.interactive import genome_browser

ROOT = Path(__file__).resolve().parents[1]
DEMO_GENOMES = {
    "Phage lambda — NC_001416.1 (48,5 kb)": ROOT / "data" / "demo" / "NC_001416.1.gb.gz",
    "M. genitalium G37 — NC_000908.2 (580 kb)": ROOT / "data" / "demo" / "NC_000908.2.gb.gz",
}
UPLOAD_TYPES = ["gb", "gbk", "gbff", "genbank", "fa", "fasta", "fna", "fas", "gz"]
AUTO_TABLE = "Déclaré dans l'annotation"

st.set_page_config(page_title="orfeval", page_icon="🧬", layout="wide")


# ---------------------------------------------------------------------------- calculs (cache)
@st.cache_resource(max_entries=6, show_spinner="Lecture du génome…")
def load_cached(content_hash: str, data: bytes, suffix: str) -> LoadedGenome:
    """Écrit le fichier dans un dossier temporaire puis le charge (clé : empreinte)."""
    folder = Path(tempfile.gettempdir()) / "orfeval-app"
    folder.mkdir(exist_ok=True)
    path = folder / f"{content_hash}{suffix}"
    path.write_bytes(data)
    return load_genome(path)


@st.cache_resource(max_entries=12, show_spinner="Prédiction et évaluation…")
def analyze_cached(
    content_hash: str, _loaded: LoadedGenome, table: int | None, predictor_json: str
) -> RunResult:
    """Analyse mise en cache par (fichier, code génétique, paramètres)."""
    predictor = PredictorConfig.model_validate_json(predictor_json)
    run = RunConfig(name="app", genome=_loaded.source, translation_table=table)
    return analyze_genome(_loaded, run, predictor, AnalysisConfig())


@st.cache_resource(max_entries=12, show_spinner="Figures…")
def figures_cached(content_hash: str, predictor_json: str, table: int | None, _result: RunResult):
    """Figures matplotlib du package, calculées une fois par analyse."""
    return build_figures(_result)


# ---------------------------------------------------------------------------- barre latérale
def sidebar() -> tuple[bytes, str, str, int | None, PredictorConfig] | None:
    st.sidebar.title("🧬 orfeval")
    st.sidebar.caption(f"v{__version__} — prédiction de gènes procaryotes")
    source = st.sidebar.radio("Génome", ["Démonstration", "Importer un fichier"], horizontal=True)
    if source == "Démonstration":
        choice = st.sidebar.selectbox("Génome de démonstration", list(DEMO_GENOMES))
        path = DEMO_GENOMES[choice]
        data, name = path.read_bytes(), path.name
    else:
        uploaded = st.sidebar.file_uploader(
            "GenBank (évaluation possible) ou FASTA", type=UPLOAD_TYPES
        )
        if uploaded is None:
            st.info("Importez un fichier GenBank ou FASTA pour commencer, ou choisissez la démo.")
            return None
        data, name = uploaded.getvalue(), uploaded.name
    suffix = "".join(Path(name).suffixes[-2:]) or ".gb"

    st.sidebar.subheader("Paramètres")
    table_choice = st.sidebar.selectbox(
        "Code génétique",
        [AUTO_TABLE, *available_tables()],
        index=0,
        help="Par défaut, le code déclaré par l'annotation (/transl_table), sinon 11.",
    )
    table = None if table_choice == AUTO_TABLE else int(table_choice)
    min_length = st.sidebar.slider("Longueur minimale (nt)", 60, 300, 90, step=3)
    strategy = st.sidebar.selectbox(
        "Choix du codon start",
        ["longest", "rbs", "score"],
        help="longest : ORF le plus long · rbs : motif de Shine-Dalgarno · score : score + RBS",
    )
    model = st.sidebar.selectbox("Modèle codant", ["codon", "dicodon"])
    max_overlap = st.sidebar.slider("Chevauchement maximal (nt)", 0, 200, 60, step=5)
    threshold = st.sidebar.number_input("Seuil de score (bits)", value=0.0, step=1.0)
    try:
        predictor = PredictorConfig(
            min_length=min_length,
            start_strategy=strategy,
            coding_model=model,
            max_overlap=max_overlap,
            score_threshold=threshold,
        )
    except ValueError as exc:
        st.sidebar.error(str(exc))
        return None
    return data, name, suffix, table, predictor


# ---------------------------------------------------------------------------- onglets
def percent(value: float | None) -> str:
    return "n/a" if value is None or pd.isna(value) else f"{100 * value:.1f} %"


def tab_summary(result: RunResult) -> None:
    genome = result.loaded.genome
    evaluation = result.evaluation
    columns = st.columns(5)
    columns[0].metric("Gènes prédits", f"{len(result.prediction.genes):,}".replace(",", " "))
    if evaluation is not None:
        columns[1].metric("Sensibilité", percent(evaluation.sensitivity.value))
        columns[2].metric("Précision", percent(evaluation.precision.value))
        columns[3].metric("F1", f"{evaluation.f1:.3f}")
        columns[4].metric("Start exact", percent(evaluation.start_accuracy.value))
        st.caption(
            "Intervalles de confiance à 95 % (Wilson) : "
            f"sensibilité {evaluation.sensitivity} · précision {evaluation.precision}"
        )
    else:
        st.warning("Pas d'annotation de référence : prédiction sans évaluation (fichier FASTA).")

    info = {
        "Séquence": f"{genome.seq_id} — {genome.description}",
        "Organisme": genome.organism or "non renseigné",
        "Longueur": f"{genome.length:,} nt".replace(",", " "),
        "Topologie": genome.topology,
        "GC": percent(result.composition.gc_content),
        "Code génétique": f"{result.table_id} ({result.table_source})",
        "ORF candidats": str(result.prediction.n_orfs),
        "Seuil d'auto-apprentissage": f"{result.prediction.train_min_length} nt",
        "Temps de calcul": f"{result.elapsed_seconds:.1f} s",
    }
    st.table(pd.DataFrame({"Valeur": info}))
    st.info(
        "**Lecture prudente.** Les métriques mesurent l'accord avec l'annotation de "
        "référence, elle-même en partie prédite. Une prédiction « sans correspondance » "
        "n'est pas forcément fausse, et un accord n'est pas une validation expérimentale."
    )


def tab_genes(result: RunResult) -> None:
    table = prediction_table(result)
    left, middle, right = st.columns(3)
    strands = left.multiselect("Brin", ["+", "-"], default=["+", "-"])
    statuses = sorted(table["status"].unique()) if "status" in table else []
    chosen = middle.multiselect("Statut", statuses, default=statuses) if statuses else []
    min_score = right.number_input("Score minimal (bits)", value=float(table["score_bits"].min()))
    mask = table["strand"].isin(strands) & (table["score_bits"] >= min_score)
    if statuses:
        mask &= table["status"].isin(chosen)
    st.dataframe(table[mask], hide_index=True)
    st.caption(f"{int(mask.sum())} gènes affichés sur {len(table)}.")

    genome = result.loaded.genome
    genes = result.prediction.genes
    downloads = st.columns(3)
    downloads[0].download_button(
        "Tableau (TSV)",
        table.to_csv(sep="\t", index=False),
        "predictions.tsv",
        "text/tab-separated-values",
    )
    downloads[1].download_button(
        "Annotation (GFF3)", gff3_text(genes, genome, "orfeval"), "predictions.gff3", "text/plain"
    )
    downloads[2].download_button(
        "Protéines (FASTA)",
        proteins_fasta_text(genes, genome, "orfeval", result.table_id),
        "proteins.faa",
        "text/plain",
    )


def tab_map(result: RunResult) -> None:
    length = result.loaded.genome.length
    default_end = min(length, 60_000)
    start, end = st.slider(
        "Région affichée (nt)", 0, length, (0, default_end), step=max(1, length // 1000)
    )
    references = result.evaluation.reference_table if result.evaluation is not None else None
    st.plotly_chart(genome_browser(length, prediction_table(result), references, (start, end)))
    st.caption("Survolez un gène pour afficher ses coordonnées, son score ou son produit annoté.")


def tab_evaluation(result: RunResult, figures: dict) -> None:
    evaluation = result.evaluation
    if evaluation is None:
        st.warning("Évaluation impossible sans annotation de référence.")
        return
    left, right = st.columns(2)
    left.pyplot(figures["precision_recall"])
    right.pyplot(figures["recall_by_length"])
    st.pyplot(figures["errors"])
    st.subheader("Gènes de référence")
    references = evaluation.reference_table
    missed_only = st.toggle("Uniquement les gènes manqués", value=True)
    shown = references[~references["found"]] if missed_only else references
    st.dataframe(shown, hide_index=True)
    st.download_button(
        "Comparaison à la référence (TSV)",
        references.to_csv(sep="\t", index=False),
        "reference_comparison.tsv",
        "text/tab-separated-values",
    )


def tab_composition(figures: dict) -> None:
    st.pyplot(figures["gc_skew"])
    left, right = st.columns(2)
    left.pyplot(figures["orf_length_null"])
    right.pyplot(figures["upstream"])
    if "codon_usage" in figures:
        st.pyplot(figures["codon_usage"], width="content")


def tab_method() -> None:
    st.markdown(
        """
        1. **ORF** recherchés sur les six cadres avec le code génétique choisi.
        2. **Modèle nul** : longueur au-delà de laquelle moins d'un ORF est attendu par
           hasard → ensemble d'auto-apprentissage, construit **sans l'annotation**.
        3. **Modèle codant** (usage des codons) → score en bits pour chaque candidat.
        4. **Sélection** gloutonne limitant les chevauchements, puis ré-entraînement.
        5. **Évaluation** : même brin et même codon stop qu'une CDS annotée.

        **Limites principales** : référence imparfaite, gènes courts difficiles, modèle
        plus simple que Prodigal ou GeneMarkS, une séquence par fichier.
        Détails : `docs/methodology.md` du dépôt.
        """
    )


# ---------------------------------------------------------------------------- page
def main() -> None:
    selection = sidebar()
    st.title("Prédiction et évaluation de gènes procaryotes")
    if selection is None:
        return
    data, name, suffix, table, predictor = selection
    content_hash = hashlib.sha256(data).hexdigest()[:16]
    predictor_json = predictor.model_dump_json()
    try:
        loaded = load_cached(content_hash, data, suffix)
        result = analyze_cached(content_hash, loaded, table, predictor_json)
    except OrfevalError as exc:
        st.error(f"Analyse impossible : {exc}")
        return
    st.caption(f"Fichier : {name}")
    figures = figures_cached(content_hash, predictor_json, table, result)

    tabs = st.tabs(
        ["Résumé", "Gènes prédits", "Carte interactive", "Évaluation", "Composition", "Méthode"]
    )
    with tabs[0]:
        tab_summary(result)
    with tabs[1]:
        tab_genes(result)
    with tabs[2]:
        tab_map(result)
    with tabs[3]:
        tab_evaluation(result, figures)
    with tabs[4]:
        tab_composition(figures)
    with tabs[5]:
        tab_method()


main()
