"""Headless smoke test of the Streamlit app using Streamlit's AppTest."""
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

ROOT_APP = str(Path(__file__).resolve().parent.parent / "app.py")


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

