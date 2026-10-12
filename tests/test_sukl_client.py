import time
import urllib.error
import urllib.request

import pytest

from core import sukl_client
from core.sukl_client import (
    PRESET_ATC_CODES,
    PRESET_DEFAULT_SUKL_CODES,
    USER_AGENT,
    fetch_sukl_detail,
    match_preset_id,
    parse_dose_mg,
    query_sukl_api,
)


def _detail(code: str):
    """fetch_sukl_detail, retried: SÚKL sometimes resets a connection, and None then means transport, not data."""
    for attempt in range(4):
        detail = fetch_sukl_detail(code)
        if detail is not None:
            return detail
        time.sleep(1 + attempt)
    return None


def _query(text: str):
    """query_sukl_api, retried for the same reason as _detail: an empty list may be a reset connection."""
    for attempt in range(4):
        results = query_sukl_api(text)
        if results:
            return results
        time.sleep(1 + attempt)
    return []


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


@pytest.mark.network
def test_query_sukl_api_live_livizux(requires_sukl):
    results = _query("livizux")
    assert results
    # Live API check: Livizux must return results from SÚKL
    first = results[0]
    assert "LIVIZUX" in first.nazev.upper()
    assert first.matched_preset_id == "livizux"
    assert first.parsed_dose_mg in (20.0, 30.0, 40.0, 50.0, 60.0, 70.0)
    assert first.kod_sukl != ""


@pytest.mark.network
def test_query_sukl_api_live_livizux_expands_all_packages(requires_sukl):
    results = _query("Livizux 70mg")
    assert results
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


def test_every_default_code_declares_the_atc_it_must_resolve_to():
    assert set(PRESET_ATC_CODES) == set(PRESET_DEFAULT_SUKL_CODES)


@pytest.mark.network
@pytest.mark.parametrize("preset_id", sorted(PRESET_DEFAULT_SUKL_CODES))
def test_default_code_is_the_presets_substance(requires_sukl, preset_id):
    """diazepam's code was DITHIADEN (bisulepin) and amoxicillin's was co-amoxiclav; the ATC catches both."""
    code = PRESET_DEFAULT_SUKL_CODES[preset_id]
    detail = _detail(code)
    assert detail is not None, f"SÚKL returned nothing for {preset_id} code {code}"
    assert detail.atc_kod == PRESET_ATC_CODES[preset_id], (
        f"{preset_id} code {code} is {detail.nazev} {detail.doplnek}, ATC {detail.atc_kod}"
    )
    assert detail.spc_url, f"{preset_id} code {code} has no SPC on file"


@pytest.mark.network
def test_spc_link_serves_the_pdf(requires_sukl):
    """The old sukl.cz download.php link answered 200 with an HTML notice instead of the PDF."""
    detail = _detail(PRESET_DEFAULT_SUKL_CODES["lithium"])
    assert detail is not None and detail.spc_url
    for attempt in range(4):
        try:
            req = urllib.request.Request(detail.spc_url, headers={"User-Agent": USER_AGENT, "Range": "bytes=0-1023"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                head = resp.read(8)
            break
        except OSError:
            if attempt == 3:
                raise
            time.sleep(1 + attempt)
    assert head.startswith(b"%PDF"), head


@pytest.mark.network
def test_fetch_sukl_detail_live_livizux(requires_sukl):
    detail = _detail("0265888")
    assert detail
    assert detail.kod_sukl == "0265888"
    assert "LIVIZUX" in detail.nazev.upper()
    assert detail.zpusob_vydeje_kod in ("rx", "blue_stripe")
    assert detail.uhrada_rezim == "none"  # Livizux 100 caps is non-reimbursed in CZ


@pytest.mark.network
def test_fetch_sukl_detail_live_livizux_shortage_and_reimbursement(requires_sukl):
    detail = _detail("0265887")
    assert detail
    assert detail.kod_sukl == "0265887"
    assert "LIVIZUX" in detail.nazev.upper()
    # 30-cap package is reimbursed and has interrupted supply (dostupnost: 'B' / mrData)
    assert detail.dostupnost_status == "interrupted"
    assert "přerušení" in detail.dostupnost_text.lower() or "výpadek" in detail.dostupnost_text.lower()
    assert detail.uhrada_rezim == "partial"
    assert detail.uhrada_zp is not None and detail.uhrada_zp > 1000.0
    assert detail.max_doplatek is not None and detail.max_doplatek > 0.0
    assert detail.indikacni_omezeni is not None


@pytest.mark.network
def test_fetch_sukl_detail_live_ritalin(requires_sukl):
    detail = _detail("0268766")
    assert detail
    assert detail.kod_sukl == "0268766"
    assert "RITALIN" in detail.nazev.upper()
    assert detail.dostupnost_status == "active"
    assert detail.uhrada_rezim == "partial"
    assert detail.uhrada_zp is not None and detail.uhrada_zp > 0.0
    assert detail.preskripcni_omezeni is not None


@pytest.mark.network
def test_fetch_sukl_detail_live_ritalin_alias_code(requires_sukl):
    detail = _detail("0015642")
    assert detail
    assert detail.kod_sukl in ("0015642", "0268766")
    assert "RITALIN" in detail.nazev.upper()
    assert detail.uhrada_rezim == "partial"


@pytest.mark.network
def test_query_sukl_api_live_unpadded_and_alias(requires_sukl):
    res_unpadded = _query("265887")
    assert "0265887" in [r.kod_sukl for r in res_unpadded]

    res_alias = _query("0015642")
    assert any("RITALIN" in r.nazev.upper() for r in res_alias)
    assert any(r.matched_preset_id == "ritalin" for r in res_alias)


@pytest.mark.network
def test_fetch_sukl_detail_live_levelanz_reimbursement(requires_sukl):
    detail = _detail("0267073")
    assert detail
    assert detail.kod_sukl == "0267073"
    assert detail.uhrada_zp is not None and detail.uhrada_zp > 0
    assert detail.uhrada_rezim == "partial"
    assert detail.preskripcni_omezeni is not None
    assert "neurologie" in detail.preskripcni_omezeni.lower()


@pytest.mark.network
def test_fetch_sukl_detail_parses_drzitel_and_substances_for_livizux(requires_sukl):
    detail = _detail("0265887")
    assert detail
    assert "Egis Pharmaceuticals" in detail.drzitel
    assert detail.atc_kod == "N06BA12"
    assert "LISDEXAMFETAMIN" in detail.atc_nazev.upper()
    assert "LISDEXAMFETAMIN" in detail.leciva_latka.upper()
    assert detail.sila != ""


def test_match_preset_id_livizux_trade_name():
    assert match_preset_id("", "LIVIZUX") == "livizux"
    assert match_preset_id("N06BA12", "LIVIZUX 70MG") == "livizux"
    assert match_preset_id("LISDEXAMFETAMIN", "LIVIZUX") == "livizux"




class _Response:
    def __init__(self, body: bytes):
        self.body = body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self):
        return self.body


@pytest.mark.parametrize("failure", [
    urllib.error.URLError(ConnectionResetError(104, "Connection reset by peer")),
    urllib.error.URLError(TimeoutError("_ssl.c:1015: The handshake operation timed out")),
    TimeoutError("The read operation timed out"),
], ids=["reset", "handshake-timeout", "read-timeout"])
def test_a_dropped_connection_is_retried_once(monkeypatch, failure):
    """One dropped SÚKL connection must not read as 'no packages'."""
    calls = []

    def urlopen(req, timeout):
        calls.append(req)
        if len(calls) == 1:
            raise failure
        return _Response(b'{"kodSUKL": "0265887"}')

    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(time, "sleep", lambda s: None)
    assert sukl_client._open_json(urllib.request.Request("https://example.invalid")) == {"kodSUKL": "0265887"}
    assert len(calls) == 2


def test_a_second_reset_or_an_http_error_is_not_retried(monkeypatch):
    calls = []

    def reset(req, timeout):
        calls.append(req)
        raise ConnectionResetError(104, "Connection reset by peer")

    monkeypatch.setattr(urllib.request, "urlopen", reset)
    monkeypatch.setattr(time, "sleep", lambda s: None)
    with pytest.raises(ConnectionResetError):
        sukl_client._open_json(urllib.request.Request("https://example.invalid"))
    assert len(calls) == 2

    calls.clear()

    def not_found(req, timeout):
        calls.append(req)
        raise urllib.error.HTTPError(req.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(urllib.request, "urlopen", not_found)
    with pytest.raises(urllib.error.HTTPError):
        sukl_client._open_json(urllib.request.Request("https://example.invalid"))
    assert len(calls) == 1
