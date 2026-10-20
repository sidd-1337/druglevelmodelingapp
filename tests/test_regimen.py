"""A regimen simulated as a whole: interaction routing, timing, enzymes and clashes."""
from dataclasses import replace

import numpy as np
import pytest

from core.drug_database import get_drug
from core.enzymes import load_perpetrators
from core.pk_engine import Route, regular_schedule, time_grid
from core.regimen import SUBSTANCE, AddedSubstance, RegimenDrug, simulate_regimen

DAYS = 10
T = time_grid(24.0 * DAYS)


def _entry(drug_id, route=Route.ORAL, **model_changes):
    d = get_drug(drug_id)
    m = d.to_pk_model(70.0, route)
    if model_changes:
        m = replace(m, **model_changes)
    events = tuple(regular_schedule(d.amount_from_dose_mg(d.default_dose_mg), d.default_interval_h, 24.0 * DAYS))
    return RegimenDrug(d, m, events, d.default_interval_h, d.enzyme_shares(70.0))


def _drug(result, drug_id):
    return next(r for r in result.drugs if r.entry.key == drug_id)


def test_measured_factor_between_two_presets():
    res = simulate_regimen(T, [_entry("lithium"), _entry("ibuprofen")])
    assert _drug(res, "lithium").exposure_factor == pytest.approx(1.15, rel=1e-3)
    assert _drug(res, "ibuprofen").exposure_factor == pytest.approx(1.0)


def test_added_substance_defaults_to_already_taken():
    res = simulate_regimen(T, [_entry("risperidone")], [AddedSubstance("Paroxetine 20 mg/day")])
    r = _drug(res, "risperidone")
    assert r.exposure_factor == pytest.approx(1.6, rel=1e-3)  # the measured active-moiety ratio
    assert r.drivers == ["Paroxetine 20 mg/day"] and r.applied[0].timed


def test_started_late_it_has_not_acted_yet_early_on():
    res = simulate_regimen(T, [_entry("risperidone")], [AddedSubstance("Paroxetine 20 mg/day", start_h=24.0 * 8)])
    r = _drug(res, "risperidone")
    early = T < 24.0 * 8
    assert np.allclose(r.concentrations[early], r.alone[early])
    assert 1.0 < r.exposure_factor < 1.6


def test_bioavailability_interaction_skips_an_iv_dose():
    res = simulate_regimen(T, [_entry("ketamine", Route.IV_BOLUS)], [AddedSubstance("Clarithromycin")])
    r = _drug(res, "ketamine")
    assert [a.skipped for a in r.applied] == ["iv"]
    assert r.exposure_factor == pytest.approx(1.0)


def test_bioavailability_interaction_is_capped_at_complete_absorption():
    res = simulate_regimen(T, [_entry("ketamine", bioavailability=0.9)], [AddedSubstance("Clarithromycin")])
    (a,) = _drug(res, "ketamine").applied
    assert a.factor == pytest.approx(1.0 / 0.9)  # the measured x2.6 cannot push F past 1


def test_inhibitor_on_a_destroyed_enzyme_goes_through_the_pool():
    added = [AddedSubstance("Paroxetine 20 mg/day")]
    both = simulate_regimen(T, [_entry("risperidone"), _entry("mdma")], added)
    r = _drug(both, "risperidone")
    assert [(a.interaction.perpetrator, a.via_enzyme) for a in r.applied] == [("Paroxetine 20 mg/day", True)]
    assert both.enzyme_drivers["CYP2D6"] == ["MDMA (ecstasy)", "Paroxetine 20 mg/day"]
    # Bounded by losing the whole CYP2D6 share (1 / (1 - 0.8)), unlike stacking x1.6 on the pool.
    mdma_only = _drug(simulate_regimen(T, [_entry("risperidone"), _entry("mdma")]), "risperidone")
    assert r.exposure_factor < 1.0 / (1.0 - 0.8)
    assert r.exposure_factor < 1.6 * mdma_only.exposure_factor


def test_alias_gives_a_measured_pair_its_time_course():
    res = simulate_regimen(T, [_entry("diazepam")], [AddedSubstance("Fluvoxamine 100-150 mg/day", start_h=0.0)])
    (a,) = _drug(res, "diazepam").applied
    assert a.timed and a.factor == pytest.approx(2.86)


def test_three_inhibitors_of_one_enzyme_stay_below_losing_the_whole_enzyme():
    added = [AddedSubstance(n) for n in ("Paroxetine 20 mg/day", "Fluoxetine 20 mg/day", "Terbinafine")]
    r = _drug(simulate_regimen(T, [_entry("mirtazapine")], added), "mirtazapine")
    ceiling = 1.0 / (1.0 - 0.36)  # mirtazapine's CYP2D6 share
    assert 1.5 < r.exposure_factor < ceiling  # multiplying the three factors gave about x3.4
    assert [a.pooled for a in r.applied] == ["CYP2D6"] * 3
    multiplied, combined = r.pooled["CYP2D6"]
    assert multiplied > 3.0 and combined == pytest.approx(r.exposure_factor, rel=0.02)


def test_two_measured_inhibitors_give_more_than_either_and_less_than_their_product():
    added = [AddedSubstance("Paroxetine 20 mg/day"), AddedSubstance("Fluoxetine 20 mg/day")]
    r = _drug(simulate_regimen(T, [_entry("risperidone")], added), "risperidone")
    assert 1.6 <= r.exposure_factor < 1.6 * 1.4


def test_inhibitor_and_inducer_follow_the_net_effect_model():
    keto, carba = "Ketoconazole 200 mg b.i.d.", "Carbamazepine 400 mg b.i.d."
    r = _drug(simulate_regimen(T, [_entry("diazepam")], [AddedSubstance(keto), AddedSubstance(carba)]), "diazepam")
    residuals = {p.name: p.residual_activity for p in load_perpetrators()}
    fm = get_drug("diazepam").pathways["CYP3A4"]
    # Activities multiply, so the strong inhibitor wins on the induced enzyme. The two
    # factors multiplied would have given x0.80.
    assert r.exposure_factor == pytest.approx(1.0 / (fm * residuals[keto] * residuals[carba] + 1.0 - fm), rel=0.01)


def test_one_substance_per_enzyme_is_not_pooled():
    added = [AddedSubstance("Fluconazole"), AddedSubstance("Ketoconazole 200 mg b.i.d.")]  # CYP2C19, CYP3A4
    r = _drug(simulate_regimen(T, [_entry("diazepam")], added), "diazepam")
    assert r.pooled == {} and all(not a.pooled for a in r.applied)
    assert r.exposure_factor == pytest.approx(np.prod([a.factor for a in r.applied]), rel=1e-3)


def test_two_substances_on_one_enzyme_combine_on_the_chart():
    added = [AddedSubstance("Ketoconazole 200 mg b.i.d."), AddedSubstance("Rifampicin")]
    res = simulate_regimen(T, [_entry("diazepam")], added)
    assert res.enzyme_drivers["CYP3A4"] == ["Ketoconazole 200 mg b.i.d.", "Rifampicin"]
    keto = simulate_regimen(T, [_entry("diazepam")], added[:1]).enzyme_activity["CYP3A4"]
    rif = simulate_regimen(T, [_entry("diazepam")], added[1:]).enzyme_activity["CYP3A4"]
    assert np.allclose(res.enzyme_activity["CYP3A4"], keto * rif)


def test_enzyme_not_used_by_any_drug_is_not_charted():
    res = simulate_regimen(T, [_entry("lithium")], [AddedSubstance("Rifampicin")])  # lithium is renal
    assert res.enzyme_activity == {}
    assert res.unmatched == ["Rifampicin"]


def test_added_substances_take_part_in_the_clash_check():
    res = simulate_regimen(T, [_entry("mdma")], [AddedSubstance("Paroxetine 20 mg/day")])
    key = SUBSTANCE + "Paroxetine 20 mg/day"
    assert key in res.labels
    assert key in res.tags or key in res.unassessed  # never silently left out
