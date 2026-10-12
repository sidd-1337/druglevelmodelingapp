import pytest

from core.sukl_client import (
    PRESET_DEFAULT_SUKL_CODES,
    fetch_sukl_detail,
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
    assert match_preset_id("METHYLFENIDÁT", "RITALIN") == "ritalin"
    assert match_preset_id("METHYLPHENIDATE", "RITALIN") == "ritalin"
    assert match_preset_id("METHYLFENIDAT", "CONCERTA") == "ritalin"
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


def test_query_sukl_api_live_livizux_expands_all_packages():
    results = query_sukl_api("Livizux 70mg")
    if results:
        codes = [r.kod_sukl for r in results]
        # Must return both 30-cap (0265887) and 100-cap (0265888) packages
        assert "0265887" in codes
        assert "0265888" in codes
        for r in results:
            assert r.parsed_dose_mg == 70.0
            assert r.matched_preset_id == "livizux"


def test_fetch_sukl_detail_invalid():
    assert fetch_sukl_detail("xyz") is None
    assert fetch_sukl_detail("") is None
    assert fetch_sukl_detail("999999999") is None


def test_preset_default_sukl_codes():
    assert "livizux" in PRESET_DEFAULT_SUKL_CODES
    assert "levetiracetam" in PRESET_DEFAULT_SUKL_CODES
    assert "ritalin" in PRESET_DEFAULT_SUKL_CODES
    assert PRESET_DEFAULT_SUKL_CODES["livizux"] == "0265888"
    assert PRESET_DEFAULT_SUKL_CODES["ritalin"] == "0268766"
    assert PRESET_DEFAULT_SUKL_CODES["ritalin"] in ("0015642", "0015643")
    assert PRESET_DEFAULT_SUKL_CODES["ritalin"] == "0015642"


def test_no_two_presets_share_a_default_sukl_code():
    """One code is one product, so two presets sharing it means one points at the wrong drug.

    paracetamol was mapped to IBALGIN 400 (ibuprofen) this way, which showed ibuprofen's
    price, dispensing status and shortage state under paracetamol in the app.
    """
    seen: dict[str, str] = {}
    for preset_id, code in PRESET_DEFAULT_SUKL_CODES.items():
        assert str(code) not in seen, f"{preset_id} and {seen[str(code)]} share SÚKL code {code}"
        seen[str(code)] = preset_id


def test_every_default_code_belongs_to_a_real_preset():
    from core.drug_database import load_drugs

    known = {d.id for d in load_drugs()}
    assert set(PRESET_DEFAULT_SUKL_CODES) <= known


def test_fetch_sukl_detail_live_livizux():
    detail = fetch_sukl_detail("0265888")
    if detail:
        assert detail.kod_sukl == "0265888"
        assert "LIVIZUX" in detail.nazev.upper()
        assert detail.zpusob_vydeje_kod in ("rx", "blue_stripe")
        assert detail.uhrada_rezim == "none"  # Livizux 100 caps is non-reimbursed in CZ


def test_fetch_sukl_detail_live_livizux_shortage_and_reimbursement():
    detail = fetch_sukl_detail("0265887")
    if detail:
        assert detail.kod_sukl == "0265887"
        assert "LIVIZUX" in detail.nazev.upper()
        # 30-cap package is reimbursed and has interrupted supply (dostupnost: 'B' / mrData)
        assert detail.dostupnost_status == "interrupted"
        assert "přerušení" in detail.dostupnost_text.lower() or "výpadek" in detail.dostupnost_text.lower()
        assert detail.uhrada_rezim == "partial"
        assert detail.uhrada_zp is not None and detail.uhrada_zp > 1000.0
        assert detail.max_doplatek is not None and detail.max_doplatek > 0.0
        assert detail.indikacni_omezeni is not None


def test_fetch_sukl_detail_live_ritalin():
    detail = fetch_sukl_detail("0268766")
    if detail:
        assert detail.kod_sukl == "0268766"
        assert "RITALIN" in detail.nazev.upper()
        assert detail.dostupnost_status == "active"
        assert detail.uhrada_rezim == "partial"
        assert detail.uhrada_zp is not None and detail.uhrada_zp > 0.0
        assert detail.preskripcni_omezeni is not None


def test_fetch_sukl_detail_live_ritalin_alias_code():
    detail = fetch_sukl_detail("0015642")
    if detail:
        assert detail.kod_sukl in ("0015642", "0268766")
        assert "RITALIN" in detail.nazev.upper()
        assert detail.uhrada_rezim == "partial"


def test_query_sukl_api_live_unpadded_and_alias():
    res_unpadded = query_sukl_api("265887")
    if res_unpadded:
        codes = [r.kod_sukl for r in res_unpadded]
        assert "0265887" in codes

    res_alias = query_sukl_api("0015642")
    if res_alias:
        assert any("RITALIN" in r.nazev.upper() for r in res_alias)
        assert any(r.matched_preset_id == "ritalin" for r in res_alias)


def test_fetch_sukl_detail_live_levelanz_reimbursement():
    detail = fetch_sukl_detail("0267073")
    if detail:
        assert detail.kod_sukl == "0267073"
        assert detail.uhrada_zp is not None and detail.uhrada_zp > 0
        assert detail.uhrada_rezim == "partial"
        assert detail.preskripcni_omezeni is not None
        assert "neurologie" in detail.preskripcni_omezeni.lower()


def test_fetch_sukl_detail_parses_drzitel_and_substances_for_livizux():
    detail = fetch_sukl_detail("0265887")
    if detail:
        assert "Egis Pharmaceuticals" in detail.drzitel
        assert detail.atc_kod == "N06BA12"
        assert "LISDEXAMFETAMIN" in detail.atc_nazev.upper()
        assert "LISDEXAMFETAMIN" in detail.leciva_latka.upper()
        assert detail.sila != ""


def test_match_preset_id_livizux_trade_name():
    assert match_preset_id("", "LIVIZUX") == "livizux"
    assert match_preset_id("N06BA12", "LIVIZUX 70MG") == "livizux"
    assert match_preset_id("LISDEXAMFETAMIN", "LIVIZUX") == "livizux"


