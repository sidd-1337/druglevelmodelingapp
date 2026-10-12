import math

import pytest

from core.drug_database import get_drug, load_drugs
from core.enzymes import (
    CLASS_BANDS,
    Enzyme,
    Perpetrator,
    derived_exposure_factor,
    find_perpetrators,
    genotype_activity,
    load_enzymes,
    load_perpetrators,
    polymorphic_enzymes,
    residual_activity,
    resolve_enzyme_activity,
)
from core.interactions import CLEARANCE, derived_interactions, find, interactions_for


def _perp(**kw):
    base = dict(name="X", name_cs="X", enzyme="CYP2D6", drug_class="strong_inhibitor",
                index_substrate="desipramine", index_auc_ratio=5.0, index_fm=0.85, source="http://e")
    return Perpetrator(**{**base, **kw})


def test_registry_loads_and_cites_sources():
    enzymes, perpetrators = load_enzymes(), load_perpetrators()
    assert {"CYP2D6", "CYP2C19", "CYP3A4"} <= set(enzymes)
    assert all(e.source.startswith("http") for e in enzymes.values())
    assert all(p.source.startswith("http") for p in perpetrators)
    assert all(isinstance(e.verified, bool) for e in enzymes.values())


def test_every_drug_pathway_resolves_to_a_known_enzyme_or_a_non_enzyme_route():
    # Non-enzyme routes legitimately have no registry entry.
    non_enzyme = {"renal", "UGT", "sulfation", "hydrolysis", "CES1", "CYP2E1"}
    known = set(load_enzymes()) | non_enzyme
    for d in load_drugs():
        assert set(d.pathways or {}) <= known, d.id


def test_residual_is_derived_from_the_index_substrate_not_stored():
    p = _perp(index_auc_ratio=4.8, index_fm=0.85)
    assert p.residual_activity == pytest.approx((1 / 4.8 - 0.15) / 0.85)
    # ignoring the index substrate's own fm would overstate what activity is left
    assert p.residual_activity < 1 / 4.8


def test_every_derived_residual_sits_in_its_fda_class_band():
    for p in load_perpetrators():
        lo, hi = CLASS_BANDS[p.drug_class]
        assert lo - 1e-9 <= p.residual_activity <= hi + 1e-9, p.name


def test_a_residual_outside_its_class_band_is_rejected():
    with pytest.raises(ValueError, match="band"):
        _perp(drug_class="weak_inhibitor", index_auc_ratio=10.0)


def test_inducers_have_residual_above_one():
    inducers = [p for p in load_perpetrators() if p.drug_class == "inducer"]
    assert inducers and all(p.residual_activity > 1 for p in inducers)
    assert derived_exposure_factor(0.5, inducers[0].residual_activity) < 1.0


def test_exposure_factor_relation():
    assert derived_exposure_factor(0.0, 0.1) == 1.0  # no clearance through the enzyme, no effect
    assert derived_exposure_factor(1.0, 0.5) == pytest.approx(2.0)
    assert derived_exposure_factor(0.5, 0.0) == pytest.approx(2.0)
    # a higher fraction metabolised means a bigger effect from the same inhibitor
    assert derived_exposure_factor(0.8, 0.1) > derived_exposure_factor(0.3, 0.1)
    with pytest.raises(ValueError):
        derived_exposure_factor(1.0, 0.0)  # wholly cleared by a fully inhibited enzyme


def test_genotype_activity_defaults_to_normal():
    assert genotype_activity("CYP2D6", None) == 1.0
    assert genotype_activity("CYP2D6", "normal") == 1.0
    assert genotype_activity("CYP2D6", "poor") < 0.2
    assert genotype_activity("CYP2C19", "ultrarapid") > 1.0
    with pytest.raises(KeyError):
        genotype_activity("CYP2D6", "nonsense")
    assert resolve_enzyme_activity({"CYP2D6": "poor"})["CYP2D6"] < 0.2


def test_enzymes_without_a_phenotype_system_are_not_polymorphic():
    # CYP3A4 and CYP1A2 have no validated poor/ultrarapid system, so no genotype selector.
    assert {"CYP2D6", "CYP2C19", "CYP2C9"} <= polymorphic_enzymes()
    assert not ({"CYP3A4", "CYP1A2"} & polymorphic_enzymes())


def test_enzyme_validation():
    with pytest.raises(ValueError, match="normal"):
        Enzyme("E", {"poor": 0.1}, {}, "http://e")
    with pytest.raises(ValueError, match="sum"):
        Enzyme("E", {"normal": 1.0, "poor": 0.1}, {"normal": 0.2, "poor": 0.2}, "http://e")


def test_combined_residual_multiplies():
    perps = find_perpetrators(["Paroxetine 20 mg/day", "Fluoxetine 20 mg/day"])
    assert len(perps) == 2
    assert residual_activity("CYP2D6", perps) == pytest.approx(math.prod(p.residual_activity for p in perps))
    assert residual_activity("CYP3A4", perps) == 1.0  # a different enzyme is untouched


# --------------------------------------------------------- measured vs derived
def test_measured_entries_win_over_derived():
    mirt = get_drug("mirtazapine")
    rows = find("mirtazapine", ["Ketoconazole 200 mg b.i.d."], drug=mirt)
    assert len(rows) == 1 and not rows[0].estimated  # a measured study exists, so no estimate is added


def test_derived_entry_appears_only_where_no_study_exists():
    mirt = get_drug("mirtazapine")
    rows = find("mirtazapine", ["Clarithromycin"], drug=mirt)
    assert len(rows) == 1 and rows[0].estimated
    assert rows[0].site == CLEARANCE and rows[0].exposure_factor > 1
    assert "CYP3A4" in rows[0].basis and "clearance" in rows[0].basis


def test_no_derived_entry_without_a_shared_pathway():
    # lithium is wholly renal: no CYP inhibitor can act on it mechanistically
    assert derived_interactions(get_drug("lithium"), ["Ketoconazole 200 mg b.i.d."]) == []
    # ... and drugs with no pathways encoded never get estimates
    assert derived_interactions(get_drug("ketamine"), ["Clarithromycin"]) == []


def test_find_without_a_drug_is_unchanged():
    assert find("diazepam", ["Cimetidine"]) == [i for i in interactions_for("diazepam")
                                                if i.perpetrator == "Cimetidine"]
