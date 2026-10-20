import json

import numpy as np
import pytest

from core.drug_database import get_drug, load_drugs
from core.interactions import (
    BIOAVAILABILITY,
    CLEARANCE,
    Interaction,
    apply_interactions,
    combined_factor,
    find,
    interactions_for,
    load_interactions,
    mechanistic_exposure_factor,
    mechanistic_exposure_response,
)
from core.pk_engine import Route, steady_state_summary


def test_dataset_loads_and_cites_sources():
    items = load_interactions()
    assert len(items) > 10
    assert all(i.source.startswith("http") for i in items)
    known = {d.id for d in load_drugs()}
    assert {i.victim_id for i in items} <= known


def test_quantified_entries_have_a_plausible_factor():
    for i in load_interactions():
        if not i.qualitative:
            assert 0.05 <= i.exposure_factor <= 10.0, i.perpetrator


def test_clearance_interaction_raises_steady_state_proportionally():
    li = get_drug("lithium")
    model = li.to_pk_model(70)
    nsaid = find("lithium", ["Ibuprofen"])
    assert len(nsaid) == 1 and nsaid[0].direction == "increase"
    inhibited = apply_interactions(model, nsaid)
    amount = li.amount_from_dose_mg(300)
    ratio = steady_state_summary(amount, 8, inhibited).cavg_ss / steady_state_summary(amount, 8, model).cavg_ss
    assert ratio == pytest.approx(nsaid[0].exposure_factor, rel=1e-6)
    assert inhibited.t_half_h > model.t_half_h and inhibited.vd_l == model.vd_l


def test_inducer_lowers_exposure():
    mirt = get_drug("mirtazapine")
    model = mirt.to_pk_model(70)
    cbz = find("mirtazapine", ["Carbamazepine 400 mg b.i.d."])
    induced = apply_interactions(model, cbz)
    assert induced.t_half_h < model.t_half_h
    assert steady_state_summary(30, 24, induced).cavg_ss < steady_state_summary(30, 24, model).cavg_ss


def test_bioavailability_interaction_leaves_half_life_alone():
    ket = get_drug("ketamine")
    model = ket.to_pk_model(70)
    rif = find("ketamine", ["Rifampicin"])
    assert rif[0].site == BIOAVAILABILITY
    induced = apply_interactions(model, rif)
    assert induced.t_half_h == model.t_half_h
    assert induced.bioavailability == pytest.approx(model.bioavailability * rif[0].exposure_factor)


def test_bioavailability_factor_is_capped_at_one():
    lev = get_drug("levetiracetam").to_pk_model(70)  # F = 1.0
    boosted = apply_interactions(lev, [Interaction("levetiracetam", "X", "X", BIOAVAILABILITY, "http://e", 3.0)])
    assert boosted.bioavailability == 1.0


def test_effects_multiply_and_qualitative_entries_do_nothing():
    dz = get_drug("diazepam").to_pk_model(70)
    both = find("diazepam", ["Cimetidine", "Omeprazole"])
    assert len(both) == 2
    expected = both[0].exposure_factor * both[1].exposure_factor
    assert apply_interactions(dz, both).t_half_h == pytest.approx(dz.t_half_h * expected)

    qualitative = [i for i in interactions_for("lithium") if i.qualitative]
    assert qualitative and qualitative[0].direction == "unknown"
    assert apply_interactions(dz, qualitative) is dz
    assert combined_factor(qualitative, CLEARANCE) == 1.0


def test_interaction_validation():
    with pytest.raises(ValueError):
        Interaction("lithium", "X", "X", "absorption", "http://e", 1.2)
    with pytest.raises(ValueError):
        Interaction("lithium", "X", "X", CLEARANCE, "http://e", 0.0)
    with pytest.raises(ValueError):
        Interaction("lithium", "X", "X", CLEARANCE, "http://e", 1.2, qualitative=True)


def test_czech_names_present():
    for i in load_interactions():
        assert i.perpetrator_cs
    assert find("lithium", ["Ibuprofen"])[0].name("cs")


def test_every_entry_states_whether_it_was_verified():
    raw = json.loads((__import__("pathlib").Path("data/interactions.json")).read_text(encoding="utf-8"))
    assert all(isinstance(i["verified"], bool) for i in raw["interactions"])
    assert Route.ORAL.value == "oral"


@pytest.mark.parametrize("drug_id, enzyme, tolerance", [
    ("diazepam", "CYP3A4", 1e-12),  # parent only: closed form
    ("risperidone", "CYP2D6", 0.01),  # parent + 9-OH sum: tabulated simulation
])
def test_exposure_response_matches_the_pointwise_factor(drug_id, enzyme, tolerance):
    drug = get_drug(drug_id)
    activities = np.array([0.004, 0.069, 0.3, 1.0, 3.65, 25.0])
    response = mechanistic_exposure_response(drug, enzyme, activities)
    expected = [mechanistic_exposure_factor(drug, enzyme, float(a))[0] for a in activities]
    assert response == pytest.approx(expected, rel=tolerance)


def test_exposure_response_is_flat_off_the_drugs_pathways():
    assert np.all(mechanistic_exposure_response(get_drug("lithium"), "CYP3A4", [0.01, 1.0, 10.0]) == 1.0)
