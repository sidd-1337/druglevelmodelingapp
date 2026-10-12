"""Headless smoke test of the Streamlit app using Streamlit's AppTest."""
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

ROOT_APP = str(Path(__file__).resolve().parent.parent / "app.py")


def _chart_traces(at: AppTest, index: int = 0) -> list[str]:
    """Trace names of a plotly chart (AppTest exposes only the serialised proto)."""
    return [tr.get("name") for tr in json.loads(at.get("plotly_chart")[index].proto.spec)["data"]]


def _run() -> AppTest:
    at = AppTest.from_file(ROOT_APP, default_timeout=60)
    at.run()
    return at


def test_app_runs_without_exception():
    at = _run()
    assert not at.exception
    assert len(at.metric) >= 6


def test_app_all_drugs_both_languages():
    at = _run()
    drugs = list(at.selectbox[0].options)
    for lang in ("English", "Čeština"):
        at.radio[0].set_value(lang).run()
        for d in drugs:
            at.selectbox[0].set_value(d).run()
            assert not at.exception, (lang, d, at.exception)


def test_app_missed_dose_and_semilog():
    at = _run()
    next(c for c in at.checkbox if "missed" in c.label.lower()).check().run()
    assert not at.exception
    at.radio[2].set_value(True).run()  # semilog axis
    assert not at.exception


def test_app_renal_function_changes_lithium_half_life():
    at = _run()
    at.selectbox[0].set_value("lithium").run()
    t_half = lambda: next(m.value for m in at.metric if m.label == "t½")
    normal = float(t_half().split()[0])
    next(c for c in at.checkbox if "renal" in c.label.lower()).check().run()
    at.number_input(key="scr").set_value(200.0).run()
    assert not at.exception
    assert float(t_half().split()[0]) > 1.5 * normal


def test_app_measured_level_fit():
    at = _run()
    at.selectbox[0].set_value("lithium").run()
    next(c for c in at.checkbox if "measured" in c.label.lower()).check().run()
    assert not at.exception
    at.number_input(key="fit_c0_lithium").set_value(1.2).run()
    assert not at.exception
    fitted = next(m for m in at.metric if m.label.startswith("Clearance (fitted)"))
    assert fitted.delta.startswith("-")


def test_app_sukl_search_integration():
    at = _run()
    sukl_input = at.sidebar.text_input[0]
    sukl_input.set_value("Livizux").run()
    assert not at.exception


def test_app_sukl_search_livizux_both_packages():
    at = _run()
    sukl_input = at.sidebar.text_input[0]
    sukl_input.set_value("Livizux 70mg").run()
    assert not at.exception
    options = list(at.sidebar.selectbox[0].options)
    assert len(options) >= 2
    joined = " ".join(options)
    assert "0265887" in joined and "0265888" in joined


def test_app_sukl_search_ritalin_matches_preset():
    at = _run()
    sukl_input = at.sidebar.text_input[0]
    sukl_input.set_value("Ritalin").run()
    assert not at.exception
    # Drug selectbox automatically switches to matched preset 'ritalin'
    assert at.sidebar.selectbox[1].value == "ritalin"


def test_app_sukl_search_ritalin_alias_code_matches_preset():
    at = _run()
    sukl_input = at.sidebar.text_input[0]
    sukl_input.set_value("0015642").run()
    assert not at.exception
    assert at.sidebar.selectbox[1].value == "ritalin"



def test_app_glossary_tab_and_tooltips():
    at = _run()
    # Check that KPI metrics have help tooltips
    metrics_with_help = [m for m in at.metric if m.help]
    assert len(metrics_with_help) >= 5



def test_app_model_tab_shows_validation_table():
    at = _run()
    assert not at.exception
    assert any("within tolerance" in m.label for m in at.metric)


def test_app_fit_with_edited_dose_log_runs():
    at = _run()
    at.selectbox[0].set_value("lithium").run()
    next(c for c in at.checkbox if "measured" in c.label.lower()).check().run()
    assert not at.exception
    assert any("Clearance (fitted)" in m.label for m in at.metric)
    assert at.tabs  # later tabs still render after the fit tab


def test_app_co_medication_tab_applies_interaction():
    at = _run()
    at.selectbox[0].set_value("lithium").run()
    at.multiselect[0].set_value(["ibuprofen"]).run()
    assert not at.exception
    assert any("1.15x" in m.value for m in at.metric)


def test_app_co_medication_shows_a_derived_estimate():
    at = _run()
    at.selectbox[0].set_value("mirtazapine").run()
    at.multiselect[1].set_value(["Clarithromycin"]).run()  # no study for this pair
    assert not at.exception
    cells = [str(v).lower() for df in at.dataframe for v in df.value.to_numpy().ravel()]
    assert any(c == "estimated" for c in cells)
    assert any("carries 35 % of this drug's clearance" in c for c in cells)


def test_app_shows_the_metabolite_curve_for_diazepam():
    at = _run()
    at.selectbox[0].set_value("diazepam").run()
    assert not at.exception
    traces = _chart_traces(at, 0)
    assert any("Nordazepam" in str(name) for name in traces)
    assert "Concentration" in traces  # the reported quantity is still diazepam itself
