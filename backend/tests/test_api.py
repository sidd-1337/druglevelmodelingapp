import numpy as np
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
    assert rit["t_half_h"] == 2.8
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
    # Concentration at t=0 equals D / Vd = 400 / (0.13 * 70) = 400 / 9.1 ≈ 43.956
    assert data["concentrations"][0] == pytest.approx(400.0 / 9.1, rel=1e-2)


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
        "dose_mg": 40.0,
        "body_weight_kg": 70.0,
        "route": "oral",
        "elapsed_h": 1.5,
        "next_dose_mg": 40.0,
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



def test_pk_simulate_renal_function_and_variability_band():
    base = {"drug_id": "lithium", "dose_mg": 300.0, "interval_h": 8.0, "duration_h": 240.0, "body_weight_kg": 70.0}
    normal = client.post("/api/v1/pk/simulate", json=base).json()
    impaired = client.post("/api/v1/pk/simulate", json={**base, "patient": {"crcl_ml_min": 50}}).json()
    assert impaired["crcl_ml_min"] == 50
    assert impaired["metrics"]["t_half_h"] > 1.5 * normal["metrics"]["t_half_h"]
    assert impaired["steady_state"]["cavg_ss"] > 1.5 * normal["steady_state"]["cavg_ss"]
    c, lo, hi = (np.array(normal[k]) for k in ("concentrations", "band_low", "band_high"))
    assert np.all(lo <= c + 1e-9) and np.all(c <= hi + 1e-9)
    assert hi[-1] > 1.3 * c[-1]

    cg = client.post("/api/v1/pk/simulate", json={
        **base, "patient": {"age_y": 80, "female": True, "serum_creatinine_umol_l": 120}}).json()
    assert cg["crcl_ml_min"] == pytest.approx((140 - 80) * 70 * 1.23 / 120 * 0.85, abs=0.1)


def test_pk_simulate_custom_override_keeps_nonlinear_elimination():
    payload = {"drug_id": "mdma", "dose_mg": 100.0, "interval_h": 24.0, "duration_h": 24.0,
               "custom_pk": {"ka_per_h": 2.0}}
    res = client.post("/api/v1/pk/simulate", json=payload)
    assert res.status_code == 200
    presets = {p["id"]: p for p in client.get("/api/v1/pk/presets").json()}
    assert presets["mdma"]["inactivation_fraction"] > 0 and presets["lithium"]["fe_renal"] == 0.95
    assert presets["lithium"]["v1_fraction"] < 1 and presets["paracetamol"]["tlag_h"] > 0


def test_pk_individualize_from_measured_trough():
    base = {"drug_id": "lithium", "dose_mg": 400.0, "interval_h": 12.0, "duration_h": 168.0}
    sim = client.post("/api/v1/pk/simulate", json=base).json()
    i = sim["time_grid"].index(156.0)
    measured = 1.5 * sim["concentrations"][i]  # patient clears lithium slower than average
    res = client.post("/api/v1/pk/individualize",
                      json={**base, "observations": [{"time_h": 156.0, "concentration": measured}]})
    assert res.status_code == 200
    data = res.json()
    assert data["clearance_individual_l_per_h"] < data["clearance_population_l_per_h"]
    assert data["t_half_individual_h"] > data["t_half_population_h"]
    obs = data["observations"][0]
    assert obs["population_predicted"] < obs["individual_predicted"] < obs["observed"]
    assert data["steady_state"]["cavg_ss"] > sim["steady_state"]["cavg_ss"]
    assert len(data["band_low"]) == len(data["time_grid"])


def test_pk_individualize_rejects_level_before_first_dose():
    res = client.post("/api/v1/pk/individualize", json={
        "drug_id": "lithium", "dose_mg": 400.0, "interval_h": 12.0, "duration_h": 48.0,
        "observations": [{"time_h": 6.0, "concentration": 0.5}],
        "missed_dose": {"dose_index": 0}})  # first dose skipped: nothing taken before 12 h
    assert res.status_code == 422


def test_validation_report():
    data = client.get("/api/v1/pk/validation").json()
    assert data["summary"]["exposure"]["within"] == data["summary"]["exposure"]["n"]
    assert all(r["source"].startswith("http") for r in data["rows"])


def test_pk_individualize_uses_actual_dose_history():
    base = {"drug_id": "lithium", "dose_mg": 400.0, "interval_h": 12.0, "duration_h": 168.0,
            "observations": [{"time_h": 156.0, "concentration": 0.6}]}
    regular = client.post("/api/v1/pk/individualize", json=base).json()
    # Same level, but the patient skipped the last three doses: clearance must be lower to explain it.
    history = [{"time_h": 12.0 * k, "dose_mg": 400.0} for k in range(10)]
    skipped = client.post("/api/v1/pk/individualize", json={**base, "dose_history": history}).json()
    assert skipped["dose_times"][-1] == 108.0
    assert skipped["clearance_individual_l_per_h"] < regular["clearance_individual_l_per_h"]


def test_multi_drug_lithium_plus_ibuprofen_raises_lithium():
    payload = {"drugs": [{"drug_id": "lithium", "dose_mg": 400, "interval_h": 12},
                         {"drug_id": "ibuprofen", "dose_mg": 400, "interval_h": 8}], "duration_h": 240}
    data = client.post("/api/v1/pk/multi-drug", json=payload).json()
    li = next(s for s in data["series"] if s["drug_id"] == "lithium")
    ibu = next(s for s in data["series"] if s["drug_id"] == "ibuprofen")
    assert li["exposure_factor"] == pytest.approx(1.15, rel=0.01)
    assert li["interactions"][0]["perpetrator"] == "Ibuprofen"
    assert np.all(np.array(li["concentrations"]) >= np.array(li["concentrations_alone"]) - 1e-9)
    assert ibu["exposure_factor"] == 1.0 and ibu["interactions"] == []
    assert li["unit"] == "mmol/L" and ibu["unit"] == "ug/mL"


def test_multi_drug_named_perpetrator_and_unmatched():
    payload = {"drugs": [{"drug_id": "mirtazapine"}, {"perpetrator": "Carbamazepine 400 mg b.i.d."},
                         {"perpetrator": "Nonexistent drug"}]}
    data = client.post("/api/v1/pk/multi-drug", json=payload).json()
    mirt = data["series"][0]
    assert mirt["exposure_factor"] < 0.5
    assert data["unmatched_perpetrators"] == ["Nonexistent drug"]


def test_interactions_endpoint():
    items = client.get("/api/v1/pk/interactions/diazepam").json()
    assert {i["perpetrator"] for i in items} >= {"Cimetidine", "Omeprazole"}
    assert client.get("/api/v1/pk/interactions/nope").status_code == 404
    assert client.get("/api/v1/pk/interactions/levetiracetam").json() == []


def test_simulate_prolonged_release_formulation():
    base = {"drug_id": "propranolol", "dose_mg": 160.0, "interval_h": 24.0, "duration_h": 120.0}
    ir = client.post("/api/v1/pk/simulate", json=base).json()
    pr = client.post("/api/v1/pk/simulate", json={**base, "formulation": "pr"}).json()
    assert pr["steady_state"]["tmax_ss_h"] > 3 * ir["steady_state"]["tmax_ss_h"]
    assert pr["metrics"]["cmax"] < ir["metrics"]["cmax"]
    assert pr["steady_state"]["fluctuation_pct"] < ir["steady_state"]["fluctuation_pct"]
    presets = {p["id"]: p for p in client.get("/api/v1/pk/presets").json()}
    assert "pr" in presets["propranolol"]["formulations"]
    assert client.post("/api/v1/pk/simulate", json={**base, "formulation": "nope"}).status_code == 400


def test_phenotype_changes_clearance_and_unknown_phenotype_widens_the_band():
    base = {"drug_id": "diazepam", "dose_mg": 10.0, "interval_h": 24.0, "duration_h": 240.0}
    unknown = client.post("/api/v1/pk/simulate", json=base).json()
    poor = client.post("/api/v1/pk/simulate",
                       json={**base, "patient": {"phenotypes": {"CYP2C19": "poor"}}}).json()
    normal = client.post("/api/v1/pk/simulate",
                         json={**base, "patient": {"phenotypes": {"CYP2C19": "normal"}}}).json()
    # a poor metaboliser clears diazepam more slowly
    assert poor["metrics"]["t_half_h"] > 1.4 * normal["metrics"]["t_half_h"]
    # stating "normal" matches the default prediction but narrows the band, because an
    # unstated phenotype really could be any of them
    assert normal["metrics"]["t_half_h"] == pytest.approx(unknown["metrics"]["t_half_h"])
    span = lambda d: max(d["band_high"]) / max(d["band_low"])
    assert span(unknown) > span(normal)


def test_unknown_enzyme_phenotype_is_rejected():
    res = client.post("/api/v1/pk/simulate", json={
        "drug_id": "diazepam", "dose_mg": 10.0, "interval_h": 24.0, "duration_h": 48.0,
        "patient": {"phenotypes": {"CYP2C19": "nonsense"}}})
    assert res.status_code == 400
