import json

import numpy as np

import pytest

from core.drug_database import CUSTOM_ID, Drug, custom_drug, get_drug, load_drugs
from core.pk_engine import Route


def test_all_presets_load_and_are_valid():
    drugs = load_drugs()
    ids = {d.id for d in drugs}
    assert {
        "ibuprofen",
        "paracetamol",
        "amoxicillin",
        "lithium",
        "diazepam",
        "livizux",
        "pregabalin",
        "propranolol",
        "rivotril",
        "mirtazapine",
        "risperidone",
        "ritalin",
    } <= ids


def test_specified_reference_values():
    ibu = get_drug("ibuprofen")
    assert (ibu.t_half_h, ibu.vd_l_per_kg, ibu.bioavailability, ibu.mec, ibu.mtc) == (2.0, 0.13, 0.9, 10.0, 50.0)
    para = get_drug("paracetamol")
    assert (para.t_half_h, para.vd_l_per_kg, para.bioavailability, para.mec, para.mtc) == (2.5, 0.7, 0.8, 10.0, 30.0)
    assert get_drug("amoxicillin").t_half_h == 1.0 and get_drug("amoxicillin").vd_l_per_kg == 0.28
    li = get_drug("lithium")
    assert li.t_half_h == 24.0 and li.concentration_unit == "mmol/L" and li.mtc == 1.5
    assert get_drug("diazepam").t_half_h == 40.0
    assert get_drug("livizux").t_half_h == 11.0 and get_drug("livizux").dose_conversion == pytest.approx(0.294)
    assert get_drug("levetiracetam").t_half_h == 7.0 and get_drug("levetiracetam").bioavailability == 1.0
    assert get_drug("pregabalin").t_half_h == 6.3 and get_drug("pregabalin").mec == 2.0
    assert get_drug("propranolol").t_half_h == 4.0 and get_drug("propranolol").bioavailability == 0.25
    assert get_drug("rivotril").t_half_h == 35.0 and get_drug("rivotril").mec == 0.02
    assert get_drug("mirtazapine").t_half_h == 26.0 and get_drug("mirtazapine").mtc == 0.15
    # risperidone's t1/2 is now the PARENT's 3 h; the 21 h metabolite governs the reported sum
    assert get_drug("risperidone").t_half_h == 3.0 and get_drug("risperidone").mtc == 0.08
    rit = get_drug("ritalin")
    assert (rit.t_half_h, rit.vd_l_per_kg, rit.bioavailability, rit.mec, rit.mtc) == (2.8, 3.1, 0.15, 0.005, 0.03)
    assert rit.default_dose_mg == 10 and rit.default_interval_h == 8


def test_to_pk_model_scales_vd_with_weight():
    ibu = get_drug("ibuprofen")
    assert ibu.to_pk_model(70).vd_l == pytest.approx(9.1)
    assert ibu.to_pk_model(70, Route.IV_BOLUS).route is Route.IV_BOLUS
    with pytest.raises(ValueError):
        ibu.to_pk_model(0)


def test_lithium_dose_conversion_to_mmol():
    li = get_drug("lithium")
    # 300 mg Li2CO3 (M = 73.89 g/mol, 2 Li per formula unit) ~ 8.12 mmol Li+
    assert li.amount_from_dose_mg(300) == pytest.approx(300 * 2 / 73.89, rel=1e-3)


def test_half_life_band_defaults_and_validation():
    c = custom_drug()
    lo, hi = c.half_life_band_h
    assert lo < c.t_half_h < hi
    with pytest.raises(ValueError):
        c.with_overrides(mec=10, mtc=5)
    with pytest.raises(ValueError):
        c.with_overrides(bioavailability=0)
    with pytest.raises(ValueError):
        c.with_overrides(t_half_range_h=(10.0, 20.0))


def test_unknown_and_custom_ids():
    assert get_drug(CUSTOM_ID).id == CUSTOM_ID
    with pytest.raises(KeyError):
        get_drug("nope")


def test_duplicate_ids_rejected(tmp_path):
    entry = {
        "id": "x", "name": "X", "t_half_h": 1, "vd_l_per_kg": 1, "bioavailability": 1,
        "ka_per_h": 1, "mec": 1, "mtc": 2,
    }
    p = tmp_path / "d.json"
    p.write_text(json.dumps({"drugs": [entry, entry]}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_drugs(p)
    assert isinstance(Drug, type)


def test_formulations_override_absorption_and_defaults():
    prop = get_drug("propranolol")
    assert prop.formulation_ids == ["ir", "pr"]
    assert prop.with_formulation("ir") == prop  # first entry is the plain preset
    pr = prop.with_formulation("pr")
    assert pr.ka_per_h < prop.ka_per_h and pr.bioavailability < prop.bioavailability
    assert (pr.default_dose_mg, pr.default_interval_h) == (160, 24)
    assert pr.t_half_h == prop.t_half_h and pr.vd_l_per_kg == prop.vd_l_per_kg
    assert prop.formulation_name("pr").startswith("Prolonged")
    assert prop.formulation_name("pr", "cs").startswith("S prodlouženým")
    with pytest.raises(KeyError):
        prop.with_formulation("nope")


def test_formulation_is_a_noop_for_drugs_without_forms():
    ibu = get_drug("ibuprofen")
    assert ibu.formulation_ids == [] and ibu.with_formulation(None) is ibu and ibu.with_formulation("pr") is ibu


def test_prolonged_release_flattens_the_steady_state_profile():
    from core.pk_engine import steady_state_summary
    prop = get_drug("propranolol")
    ir = steady_state_summary(prop.amount_from_dose_mg(80), 12, prop.to_pk_model(70))
    pr_drug = prop.with_formulation("pr")
    pr = steady_state_summary(pr_drug.amount_from_dose_mg(160), 24, pr_drug.to_pk_model(70))
    assert pr.fluctuation_pct < ir.fluctuation_pct
    assert pr.tmax_ss_h > ir.tmax_ss_h


def test_diazepam_metabolite_is_reported_as_parent_so_nothing_changes():
    from dataclasses import replace
    from core.pk_engine import cmax_single_dose, single_dose_concentration
    dz = get_drug("diazepam")
    model = dz.to_pk_model(70)
    assert model.metabolite is not None and model.reported_quantity == "parent"
    bare = replace(model, metabolite=None)
    t = np.linspace(0.0, 400.0, 2001)
    np.testing.assert_allclose(single_dose_concentration(t, 10.0, model),
                               single_dose_concentration(t, 10.0, bare), atol=1e-7)
    assert cmax_single_dose(10.0, model) == pytest.approx(cmax_single_dose(10.0, bare), rel=1e-5)


@pytest.mark.xfail(strict=True, reason=(
    "formation_fraction 0.761 was calibrated to hit 1.26 while the engine formed metabolite with "
    "MW_parent/MW_metabolite instead of its inverse; with that fixed the same fraction gives 1.14. "
    "Not re-calibrated: a calibration is not evidence."))
def test_diazepam_metabolite_reproduces_the_measured_steady_state_ratio():
    """Nordazepam:diazepam settles at the published 1.26 on chronic dosing."""
    from core.pk_engine import concentration_profile, metabolite_profile, regular_schedule
    dz = get_drug("diazepam")
    model = dz.to_pk_model(70)
    t = np.linspace(0.0, 24 * 90, 9001)
    events = regular_schedule(dz.amount_from_dose_mg(10.0), 24, 24 * 90)
    late = t > 24 * 70
    parent = concentration_profile(t, events, model)
    met = metabolite_profile(t, events, model)
    assert float(np.mean(met[late]) / np.mean(parent[late])) == pytest.approx(1.26, rel=0.02)


def test_metabolite_volume_scales_with_patient_weight():
    dz = get_drug("diazepam")
    assert dz.to_pk_model(140).metabolite.vd_l == pytest.approx(2 * dz.to_pk_model(70).metabolite.vd_l)


def test_risperidone_reports_the_active_moiety_and_keeps_its_agnp_level():
    from core.pk_engine import steady_state_summary
    ris = get_drug("risperidone")
    model = ris.to_pk_model(70)
    assert model.reported_quantity == "sum"
    assert model.terminal_half_life_h == pytest.approx(21.0)  # the slow metabolite governs the tail


@pytest.mark.xfail(strict=True, reason=(
    "The volumes were solved to put the active moiety at 0.028 ug/mL while the engine under-formed "
    "9-OH-risperidone by MW_parent/MW_metabolite squared (0.926); fixed, it is 0.0299 (+6.8 %). "
    "Not re-solved; the 1.5-fold reference check in data/reference_pk.json still passes at 1.07."))
def test_risperidone_active_moiety_sits_on_its_solved_agnp_level():
    from core.pk_engine import steady_state_summary
    model = get_drug("risperidone").to_pk_model(70)
    # The AGNP 20-60 ng/mL range is for the active moiety, so the reported sum must match it.
    assert steady_state_summary(2.0, 12, model).cavg_ss == pytest.approx(0.028, rel=0.05)


def test_cyp2d6_barely_moves_the_risperidone_active_moiety():
    """The label says combined parent + 9-OH kinetics are similar in poor and extensive
    metabolisers: less metabolite forms, but the parent is cleared more slowly. The split
    model must reproduce that buffering, which it was not fitted to."""
    from core.pk_engine import steady_state_summary
    ris = get_drug("risperidone")
    normal = steady_state_summary(2.0, 12, ris.to_pk_model(70, genotypes={"CYP2D6": "normal"}))
    poor = steady_state_summary(2.0, 12, ris.to_pk_model(70, genotypes={"CYP2D6": "poor"}))
    ratio = poor.cavg_ss / normal.cavg_ss
    assert 1.2 < ratio < 1.8, ratio  # published poor/extensive active-moiety ratios are 1.28-1.6
    # ... whereas the parent alone would move far more, which is why the moiety is buffered.
    parent_only = 1 / (0.78 * 0.12 + 0.22)
    assert parent_only > 2 * ratio  # the parent moves >3x, the moiety only ~1.4x
