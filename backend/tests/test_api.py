import pytest
from fastapi.testclient import TestClient

from backend.app.main import app

client = TestClient(app)


def test_health():
    res = client.get("/health")
    assert res.status_code == 200
    data = res.json()
    assert data["status"] == "ok"
    assert "version" in data


def test_root():
    res = client.get("/")
    assert res.status_code == 200
    data = res.json()
    assert "docs" in data


def test_drugs_presets():
    res = client.get("/api/v1/pk/presets")
    assert res.status_code == 200
    presets = res.json()
    assert isinstance(presets, list)
    assert len(presets) >= 12
    ids = {p["id"] for p in presets}
    assert {"ibuprofen", "ritalin", "livizux", "diazepam", "rivotril"} <= ids

    # Check ritalin has route params
    rit = next(p for p in presets if p["id"] == "ritalin")
    assert rit["t_half_h"] == 2.5
    assert rit["routes"] is not None
    assert "insufflated" in rit["routes"]
    assert "smoked" in rit["routes"]


def test_drugs_search():
    res = client.get("/api/v1/drugs/search?q=ibuprofen&limit=5")
    assert res.status_code == 200
    items = res.json()
    assert isinstance(items, list)
    assert len(items) > 0


def test_drugs_search_preset_fallback():
    res = client.get("/api/v1/drugs/search?q=mirtazapine&limit=5")
    assert res.status_code == 200
    items = res.json()
    assert isinstance(items, list)
    assert len(items) > 0


def test_drugs_detail_existing_preset():
    res = client.get("/api/v1/drugs/ritalin")
    assert res.status_code == 200
    data = res.json()
    assert data["nazev"] == "Ritalin (methylfenidát)"
    assert data["matched_preset"] is not None
    assert data["matched_preset"]["id"] == "ritalin"


def test_drugs_detail_not_found():
    res = client.get("/api/v1/drugs/completely_non_existent_code_99999")
    assert res.status_code == 404
    assert "not found" in res.json()["detail"].lower()


def test_drugs_alternatives():
    res = client.get("/api/v1/drugs/ritalin/alternatives?limit=5")
    assert res.status_code == 200
    alts = res.json()
    assert isinstance(alts, list)


def test_pk_simulate_oral():
    payload = {
        "drug_id": "ibuprofen",
        "dose_mg": 400.0,
        "interval_h": 8.0,
        "duration_h": 24.0,
        "body_weight_kg": 70.0,
        "route": "oral",
    }
    res = client.post("/api/v1/pk/simulate", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["drug_id"] == "ibuprofen"
    assert data["route"] == "oral"
    assert len(data["concentrations"]) > 0
    assert len(data["time_grid"]) == len(data["concentrations"])
    assert data["metrics"]["cmax"] > 0
    assert data["metrics"]["tmax_h"] > 0
    assert data["steady_state"]["accumulation_index"] >= 1.0


def test_pk_simulate_iv_bolus():
    payload = {
        "drug_id": "ibuprofen",
        "dose_mg": 400.0,
        "interval_h": 12.0,
        "duration_h": 6.0,
        "body_weight_kg": 70.0,
        "route": "iv_bolus",
    }
    res = client.post("/api/v1/pk/simulate", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["route"] == "iv_bolus"
    # IV bolus peak occurs immediately at t = 0
    assert data["metrics"]["tmax_h"] == pytest.approx(0.0, abs=0.01)
    assert data["steady_state"]["tmax_ss_h"] == pytest.approx(0.0, abs=0.01)
    # Concentration at t=0 equals D / Vd = 400 / (0.15 * 70) = 400 / 10.5 ≈ 38.095
    assert data["concentrations"][0] == pytest.approx(400.0 / 10.5, rel=1e-2)


def test_pk_simulate_routes_comparison():
    """Verify ordering of peak times: Smoked < Insufflated < Oral."""
    base_payload = {
        "drug_id": "ritalin",
        "dose_mg": 10.0,
        "interval_h": 8.0,
        "duration_h": 24.0,
        "body_weight_kg": 70.0,
    }

    res_oral = client.post("/api/v1/pk/simulate", json={**base_payload, "route": "oral"}).json()
    res_nasal = client.post("/api/v1/pk/simulate", json={**base_payload, "route": "insufflated"}).json()
    res_smoked = client.post("/api/v1/pk/simulate", json={**base_payload, "route": "smoked"}).json()

    tmax_oral = res_oral["metrics"]["tmax_h"]
    tmax_nasal = res_nasal["metrics"]["tmax_h"]
    tmax_smoked = res_smoked["metrics"]["tmax_h"]

    assert tmax_smoked < tmax_nasal < tmax_oral


def test_pk_simulate_missed_dose():
    payload = {
        "drug_id": "ibuprofen",
        "dose_mg": 400.0,
        "interval_h": 8.0,
        "duration_h": 48.0,
        "body_weight_kg": 70.0,
        "route": "oral",
        "missed_dose": {
            "dose_index": 2,
            "delay_h": 4.0,
        },
    }
    res = client.post("/api/v1/pk/simulate", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["missed_concentrations"] is not None
    assert len(data["missed_concentrations"]) == len(data["concentrations"])


def test_harm_reduction_safe_interval_and_stacking():
    payload = {
        "drug_id": "ritalin",
        "dose_mg": 10.0,
        "body_weight_kg": 70.0,
        "route": "oral",
        "elapsed_h": 1.5,  # Redosing near peak
    }
    res = client.post("/api/v1/pk/harm-reduction", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["safe_redosing_interval_h"] > 0
    assert "phases" in data
    assert data["phases"]["onset_h"] < data["phases"]["tmax_h"] < data["phases"]["offset_h"]
    assert len(data["route_timelines"]) == 4
    assert len(data["guidelines"]) >= 5
    assert len(data["single_curve"]) > 0
    assert len(data["stacked_curve"]) > 0


def test_harm_reduction_toxic_breach():
    # Large dose of ritalin with early redose pushes over MTC (mtc = 0.03 ug/mL)
    payload = {
        "drug_id": "ritalin",
        "dose_mg": 30.0,
        "body_weight_kg": 70.0,
        "route": "oral",
        "elapsed_h": 1.5,
        "next_dose_mg": 30.0,
    }
    res = client.post("/api/v1/pk/harm-reduction", json=payload)
    assert res.status_code == 200
    data = res.json()
    assessment = data["stacking_assessment"]
    assert assessment["will_breach_mtc"] is True
    assert assessment["safe_to_redose"] is False
    assert assessment["risk_level"] in ("critical_toxic", "high")
    assert assessment["time_above_mtc_h"] > 0.0


def test_harm_reduction_routes_timelines():
    payload = {
        "drug_id": "ketamine",
        "dose_mg": 50.0,
        "body_weight_kg": 70.0,
        "route": "insufflated",
        "elapsed_h": 3.0,
    }
    res = client.post("/api/v1/pk/harm-reduction", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["route"] == "insufflated"
    routes_data = {r["route"]: r for r in data["route_timelines"]}
    assert "oral" in routes_data
    assert "iv_bolus" in routes_data
    assert "insufflated" in routes_data
    assert "smoked" in routes_data

    # IV has instant onset
    assert routes_data["iv_bolus"]["onset_minutes_range"][0] == 0

    # Ketamine has specific routes in drugs.json: oral F=20%, insufflated F=50%, smoked F=30%
    assert routes_data["oral"]["bioavailability_pct"] == 20.0
    assert routes_data["insufflated"]["bioavailability_pct"] == 50.0
    assert routes_data["smoked"]["bioavailability_pct"] == 30.0


def test_pk_simulate_route_aliases():
    base = {"drug_id": "ibuprofen", "dose_mg": 400.0, "body_weight_kg": 70.0}
    for alias, canonical in [("nasal", "insufflated"), ("inhaled", "smoked"), ("iv", "iv_bolus")]:
        res = client.post("/api/v1/pk/simulate", json={**base, "route": alias})
        assert res.status_code == 200
        assert res.json()["route"] == canonical


def test_drugs_detail_by_sukl_code():
    res = client.get("/api/v1/drugs/0268766")
    assert res.status_code == 200
    data = res.json()
    assert "ritalin" in data["nazev"].lower() or (data["matched_preset"] and data["matched_preset"]["id"] == "ritalin")


def test_harm_reduction_monotonic_phases():
    payload = {
        "drug_id": "livizux",
        "dose_mg": 70.0,
        "body_weight_kg": 70.0,
        "route": "oral",
        "elapsed_h": 4.0,
    }
    res = client.post("/api/v1/pk/harm-reduction", json=payload)
    assert res.status_code == 200
    phases = res.json()["phases"]
    assert phases["onset_h"] <= phases["peak_start_h"] <= phases["tmax_h"] <= phases["peak_end_h"] <= phases["offset_h"] <= phases["comedown_end_h"]


def test_simulate_validation_errors():
    # Negative dose
    res = client.post("/api/v1/pk/simulate", json={"drug_id": "ibuprofen", "dose_mg": -100})
    assert res.status_code == 422

    # Invalid route
    res = client.post("/api/v1/pk/simulate", json={"drug_id": "ibuprofen", "dose_mg": 100, "route": "unknown_route"})
    assert res.status_code == 400

    # Nonexistent drug
    res = client.post("/api/v1/pk/simulate", json={"drug_id": "fake_drug", "dose_mg": 100})
    assert res.status_code == 404


def test_drugs_detail_livizux_30_caps_package():
    res = client.get("/api/v1/drugs/0265887")
    assert res.status_code == 200
    data = res.json()
    assert data["kod_sukl"] == "0265887"
    assert "Egis Pharmaceuticals" in data["drzitel"]
    assert data["atc_kod"] == "N06BA12"
    assert "LISDEXAMFETAMIN" in data["atc_nazev"].upper()
    assert "LISDEXAMFETAMIN" in data["leciva_latka"].upper()
    assert data["matched_preset_id"] == "livizux"
    assert data["matched_preset"] is not None
    assert data["matched_preset"]["dose_conversion"] == 0.294

