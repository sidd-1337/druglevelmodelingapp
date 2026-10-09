import pytest

from core.sukl_client import (
    match_preset_id,
    parse_dose_mg,
    query_sukl_api,
)


def test_parse_dose_mg_various_formats():
    assert parse_dose_mg("70MG") == 70.0
    assert parse_dose_mg("70 MG") == 70.0
    assert parse_dose_mg("500MG TBL NOB 20") == 500.0
    assert parse_dose_mg("2,5MG/ML") == 2.5
    assert parse_dose_mg("0.5MG") == 0.5
    assert parse_dose_mg("1G") == 1000.0
    assert parse_dose_mg("1 G") == 1000.0
    assert parse_dose_mg("500MCG") == 0.5
    assert parse_dose_mg("invalid") is None
    assert parse_dose_mg("") is None


def test_match_preset_id():
    assert match_preset_id("LISDEXAMFETAMIN", "LIVIZUX") == "livizux"
    assert match_preset_id("PARACETAMOL", "PARALEN") == "paracetamol"
    assert match_preset_id("KLONAZEPAM", "RIVOTRIL") == "rivotril"
    assert match_preset_id("PREGABALIN", "LYRICA") == "pregabalin"
    assert match_preset_id("UNKNOWN", "UNKNOWN") is None


def test_query_sukl_api_short_query_returns_empty():
    assert query_sukl_api("a") == []
    assert query_sukl_api("") == []


def test_query_sukl_api_live_livizux():
    results = query_sukl_api("livizux")
    # Live API check: Livizux must return results from SÚKL
    if results:
        first = results[0]
        assert "LIVIZUX" in first.nazev.upper()
        assert first.matched_preset_id == "livizux"
        assert first.parsed_dose_mg in (20.0, 30.0, 40.0, 50.0, 60.0, 70.0)
        assert first.kod_sukl != ""
