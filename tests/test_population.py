import math

import numpy as np
import pytest

from core.drug_database import custom_drug, get_drug
from core.pk_engine import Route, steady_state_summary
from core.drug_database import load_drugs
from core.population import (
    RENAL,
    Z_90,
    clearance_factor,
    cockcroft_gault,
    exposure_band_models,
    log_sd,
    pathway_clearance_factor,
)


def _auc(amount, model):
    return model.f_eff * amount / model.clearance_l_per_h


def test_cockcroft_gault_si_units():
    # 1 mg/dL creatinine = 88.4 umol/L; classic formula gives (140 - 40) * 70 / 72 = 97.2 mL/min
    assert cockcroft_gault(40, 70, 88.4, female=False) == pytest.approx(97.2, rel=0.01)
    assert cockcroft_gault(40, 70, 88.4, female=True) == pytest.approx(0.85 * 97.2, rel=0.01)
    with pytest.raises(ValueError):
        cockcroft_gault(40, 70, 0, female=False)


def test_clearance_factor():
    assert clearance_factor(70) == pytest.approx(1.0)
    assert clearance_factor(140) == pytest.approx(2 ** 0.75)
    assert clearance_factor(70, crcl_ml_min=50, fe_renal=1.0) == pytest.approx(0.5)
    assert clearance_factor(70, crcl_ml_min=50, fe_renal=0.5) == pytest.approx(0.75)
    assert clearance_factor(70, crcl_ml_min=50, fe_renal=0.0) == pytest.approx(1.0)


def test_reference_patient_keeps_preset_clearance():
    for drug_id in ("lithium", "ibuprofen", "mdma"):
        d = get_drug(drug_id)
        m = d.to_pk_model(70)
        # engine t1/2 carries CL; with enzyme inactivation it is the fresh-enzyme value
        assert m.t_half_h == pytest.approx(d.t_half_h * (1 - d.inactivation_fraction))
    mdma = get_drug("mdma").to_pk_model(70)
    assert mdma.terminal_half_life_h == pytest.approx(get_drug("mdma").t_half_h)


def test_renal_impairment_prolongs_lithium_and_raises_steady_state():
    li = get_drug("lithium")
    normal, impaired = li.to_pk_model(70), li.to_pk_model(70, crcl_ml_min=50)
    assert impaired.t_half_h == pytest.approx(24.0 / (0.05 + 0.95 * 0.5))
    amount = li.amount_from_dose_mg(300)
    ratio = steady_state_summary(amount, 8, impaired).cavg_ss / steady_state_summary(amount, 8, normal).cavg_ss
    assert ratio == pytest.approx(1 / 0.525, rel=1e-6)


def test_renal_function_ignored_for_hepatically_cleared_drug():
    ibu = get_drug("ibuprofen")
    assert ibu.to_pk_model(70, crcl_ml_min=20).t_half_h == pytest.approx(ibu.t_half_h)


def test_allometric_weight_scaling():
    d = get_drug("levetiracetam")
    m70, m140 = d.to_pk_model(70), d.to_pk_model(140)
    assert m140.vd_l == pytest.approx(2 * m70.vd_l)
    assert m140.clearance_l_per_h == pytest.approx(2 ** 0.75 * m70.clearance_l_per_h)
    assert m140.t_half_h == pytest.approx(2 ** 0.25 * m70.t_half_h)


def test_exposure_band_matches_lognormal_quantiles():
    model = get_drug("propranolol").to_pk_model(70)
    cv_cl, cv_f = 0.25, 0.45
    lo, hi = exposure_band_models(model, cv_cl, cv_f)
    w = math.hypot(log_sd(cv_cl), log_sd(cv_f))
    assert _auc(40, hi) / _auc(40, model) == pytest.approx(math.exp(Z_90 * w))
    assert _auc(40, lo) / _auc(40, model) == pytest.approx(math.exp(-Z_90 * w))
    assert lo.vd_l == hi.vd_l == model.vd_l


def test_exposure_band_caps_bioavailability_and_ignores_f_for_iv():
    model = get_drug("levetiracetam").to_pk_model(70)  # F = 1.0
    _, hi = exposure_band_models(model, 0.2, 0.3)
    assert hi.bioavailability == 1.0
    iv = get_drug("propranolol").to_pk_model(70, Route.IV_BOLUS)
    lo, hi = exposure_band_models(iv, 0.25, 0.45)
    assert _auc(40, hi) / _auc(40, iv) == pytest.approx(math.exp(Z_90 * log_sd(0.25)))


def test_drug_band_falls_back_to_half_life_range_scaled_by_covariates():
    d = custom_drug().with_overrides(cv_cl=None)
    m = d.to_pk_model(140)
    lo, hi = d.exposure_band(m)
    band = d.half_life_band_h
    scale = m.t_half_h / d.t_half_h
    assert (lo.t_half_h, hi.t_half_h) == pytest.approx((band[0] * scale, band[1] * scale))
    assert np.isfinite(lo.clearance_l_per_h)


# ------------------------------------------------------------------ pathways
@pytest.mark.parametrize("weight", [30.0, 55.0, 70.0, 95.0, 140.0])
@pytest.mark.parametrize("crcl", [None, 20.0, 50.0, 100.0, 160.0])
@pytest.mark.parametrize("fe", [0.0, 0.2, 0.7, 0.95, 1.0])
def test_renal_pathway_reproduces_the_scalar_form(weight, crcl, fe):
    """The pathway form must be algebraically identical to the renal-only one."""
    expected = clearance_factor(weight, crcl, fe)
    assert pathway_clearance_factor(weight, {RENAL: fe}, crcl) == pytest.approx(expected, abs=1e-12)


def test_reference_patient_is_exactly_one():
    assert pathway_clearance_factor(70.0, {"CYP2D6": 0.3, RENAL: 0.5}, 100.0) == 1.0
    assert pathway_clearance_factor(70.0) == 1.0


def test_renal_clearance_gets_no_allometric_term():
    # CrCl already carries body size, so a wholly renal drug at normal CrCl is unscaled.
    assert pathway_clearance_factor(140.0, {RENAL: 1.0}, 100.0) == pytest.approx(1.0)
    # ... while a wholly metabolic drug does scale.
    assert pathway_clearance_factor(140.0, {"CYP3A4": 1.0}, 100.0) == pytest.approx(2 ** 0.75)


def test_unassigned_remainder_scales_with_size():
    # 40 % CYP2D6 abolished, 60 % unattributed and still size-scaled.
    assert pathway_clearance_factor(70.0, {"CYP2D6": 0.4}, None, {"CYP2D6": 0.0}) == pytest.approx(0.6)


def test_enzyme_activity_scales_only_its_own_pathway():
    pathways = {"CYP2D6": 0.3, RENAL: 0.5}
    half = pathway_clearance_factor(70.0, pathways, 100.0, {"CYP2D6": 0.5})
    assert half == pytest.approx(0.3 * 0.5 + 0.5 + 0.2)


@pytest.mark.parametrize("bad", [{"a": 0.6, "b": 0.6}, {"a": 1.2}, {"a": -0.1}])
def test_invalid_pathways_rejected(bad):
    with pytest.raises(ValueError):
        pathway_clearance_factor(70.0, bad)


def test_every_preset_keeps_fe_renal_and_pathways_consistent():
    for d in load_drugs():
        if d.pathways is not None:
            assert d.pathways.get(RENAL, 0.0) == pytest.approx(d.fe_renal)
            assert sum(d.pathways.values()) <= 1 + 1e-9


def test_enzyme_activity_does_not_touch_a_renally_cleared_drug():
    li = get_drug("lithium")  # wholly renal: abolishing a CYP must change nothing
    assert li.to_pk_model(70, enzyme_activity={"CYP2D6": 0.0}).t_half_h == pytest.approx(li.t_half_h)


def test_halving_an_enzyme_doubles_the_half_life_of_a_drug_cleared_by_it():
    d = get_drug("lithium").with_overrides(pathways={"CYP2D6": 1.0}, fe_renal=0.0)
    assert d.to_pk_model(70, enzyme_activity={"CYP2D6": 0.5}).t_half_h == pytest.approx(2 * d.t_half_h)
