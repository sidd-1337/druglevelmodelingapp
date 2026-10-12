import pytest
import numpy as np

from core.drug_database import get_drug
from core.harm_reduction import (
    assess_dose_stacking,
    compute_safe_redosing_interval_h,
    get_harm_reduction_guidelines,
    get_route_timeline,
)
from core.metrics import TherapeuticWindow
from core.pk_engine import PKModel, Route, single_dose_concentration, tmax_single_dose


def test_routes_enum():
    assert {Route.ORAL, Route.IV_BOLUS, Route.INSUFFLATED, Route.SMOKED}.issubset(set(Route))
    assert Route.INSUFFLATED.value == "insufflated"
    assert Route.SMOKED.value == "smoked"


def test_route_timelines_ordering():
    drug = get_drug("ibuprofen")
    model = drug.to_pk_model(70.0)
    amount = drug.amount_from_dose_mg(400.0)

    oral_tl = get_route_timeline(Route.ORAL, amount, model)
    iv_tl = get_route_timeline(Route.IV_BOLUS, amount, model)
    nasal_tl = get_route_timeline(Route.INSUFFLATED, amount, model)
    smoked_tl = get_route_timeline(Route.SMOKED, amount, model)

    assert iv_tl.onset_minutes_range[0] == 0
    assert iv_tl.tmax_h == 0.0
    assert iv_tl.bioavailability_pct == 100.0

    # Inhalation/smoked onset is faster than oral
    assert smoked_tl.onset_minutes_range[0] < oral_tl.onset_minutes_range[0]
    # Insufflated onset is faster than oral
    assert nasal_tl.onset_minutes_range[0] < oral_tl.onset_minutes_range[0]
    assert nasal_tl.tmax_h < oral_tl.tmax_h


def test_safe_redosing_interval():
    model = PKModel(t_half_h=4.0, vd_l=40.0, bioavailability=0.8, ka_per_h=1.5)
    window = TherapeuticWindow(mec=10.0, mtc=30.0)
    amount = 500.0

    interval_h = compute_safe_redosing_interval_h(amount, model, window)
    assert interval_h > tmax_single_dose(model)
    assert interval_h > 0

    # At interval_h, verify residual + new peak is below MTC
    residual = float(single_dose_concentration(interval_h, amount, model))
    assessment = assess_dose_stacking(amount, interval_h, model, window)
    assert not assessment.will_breach_mtc


def test_dose_stacking_early_toxic():
    # Model with narrow window: MEC 10, MTC 25
    model = PKModel(t_half_h=6.0, vd_l=20.0, bioavailability=1.0, ka_per_h=2.0)
    window = TherapeuticWindow(mec=10.0, mtc=25.0)
    amount = 400.0  # Cmax around 15 ug/mL

    # Redosing immediately after peak (e.g. t = 2 h)
    assessment = assess_dose_stacking(amount, 2.0, model, window)
    assert assessment.stacked_peak > window.mtc
    assert assessment.will_breach_mtc is True
    assert assessment.risk_level in ("critical_toxic", "high")
    assert not assessment.safe_to_redose


def test_dose_stacking_late_safe():
    model = PKModel(t_half_h=2.0, vd_l=30.0, bioavailability=0.8, ka_per_h=1.5)
    window = TherapeuticWindow(mec=5.0, mtc=40.0)
    amount = 200.0

    # Wait 5 half-lives (10 hours)
    assessment = assess_dose_stacking(amount, 10.0, model, window)
    assert assessment.will_breach_mtc is False
    assert assessment.risk_level == "safe"
    assert assessment.safe_to_redose is True


def test_guidelines_content():
    guidelines = get_harm_reduction_guidelines("ritalin")
    assert len(guidelines) >= 5
    for g in guidelines:
        assert "title_cs" in g and "title_en" in g
        assert "desc_cs" in g and "desc_en" in g
        assert len(g["title_cs"]) > 5
        assert len(g["desc_cs"]) > 10


def test_route_timelines_no_double_multiplication():
    drug = get_drug("ritalin")
    amount = drug.amount_from_dose_mg(10.0)
    # Ritalin has custom route config in drugs.json: insufflated ka = 4.5, bio = 0.7
    nasal_model = drug.to_pk_model(70.0, route=Route.INSUFFLATED)
    assert nasal_model.ka_per_h == 4.5
    assert nasal_model.bioavailability == 0.7

    tl = get_route_timeline(Route.INSUFFLATED, amount, nasal_model)
    # Should use 4.5 and 0.7 directly, NOT 4.5 * 2.5 = 11.25
    expected_tmax = tmax_single_dose(nasal_model)
    assert tl.tmax_h == pytest.approx(expected_tmax, rel=1e-4)
    assert tl.bioavailability_pct == 70.0


def test_safe_redosing_interval_consistency_with_stacking():
    # Model where single peak is above MEC
    drug = get_drug("ibuprofen")
    model = drug.to_pk_model(70.0)
    amount = drug.amount_from_dose_mg(400.0)
    window = TherapeuticWindow(drug.mec, drug.mtc)

    interval_h = compute_safe_redosing_interval_h(amount, model, window)
    assessment = assess_dose_stacking(amount, interval_h, model, window)

    # At the computed safe interval, it must actually be safe to redose
    assert assessment.safe_to_redose is True
    assert not assessment.will_breach_mtc
    assert assessment.risk_level == "safe"


def test_phases_monotonicity():
    from core.harm_reduction import estimate_phases

    # Normal case
    model = PKModel(4.0, 40.0, bioavailability=0.8, ka_per_h=1.5)
    window = TherapeuticWindow(5.0, 25.0)
    p1 = estimate_phases(200.0, model, window)
    assert p1.onset_h <= p1.peak_start_h <= p1.tmax_h <= p1.peak_end_h <= p1.offset_h <= p1.comedown_end_h

    # Extreme low MEC case (previously inverted offset_h and comedown_end_h)
    p2 = estimate_phases(100.0, PKModel(4.0, 10.0, 1.0, 2.0), TherapeuticWindow(0.01, 50.0))
    assert p2.onset_h <= p2.peak_start_h <= p2.tmax_h <= p2.peak_end_h <= p2.offset_h <= p2.comedown_end_h

    # IV Bolus case
    p_iv = estimate_phases(100.0, PKModel(4.0, 20.0, route=Route.IV_BOLUS))
    assert p_iv.onset_h == 0.0
    assert p_iv.tmax_h == 0.0
    assert p_iv.peak_end_h <= p_iv.offset_h <= p_iv.comedown_end_h


def test_safe_redosing_interval_large_redose():
    model = PKModel(t_half_h=4.0, vd_l=40.0, bioavailability=0.8, ka_per_h=1.5)
    window = TherapeuticWindow(mec=5.0, mtc=30.0)
    # Base dose: 800 mg (peak is ~16 ug/mL, ceiling is 25.5 ug/mL)
    # Moderate next dose: 400 mg (peak ~8 ug/mL) -> max_residual 17.5 > MEC (5.0), so MEC governs
    interval_moderate = compute_safe_redosing_interval_h(800.0, model, window, redose_amount=400.0)
    # Very large next dose: 1500 mg (peak ~22.6 ug/mL) -> max_residual 2.9 < MEC (5.0), so stacking ceiling governs
    interval_large = compute_safe_redosing_interval_h(800.0, model, window, redose_amount=1500.0)
    assert interval_large > interval_moderate


def test_dose_stacking_time_above_mtc():
    model = PKModel(t_half_h=6.0, vd_l=20.0, bioavailability=1.0, ka_per_h=2.0)
    window = TherapeuticWindow(mec=10.0, mtc=25.0)
    assessment = assess_dose_stacking(400.0, 2.0, model, window)
    assert assessment.will_breach_mtc is True
    assert assessment.time_above_mtc_h > 0.0
    assert "Předpokládaná doba v toxickém pásmu" in assessment.risk_explanation


def test_pk_model_route_aliases():
    m_iv = PKModel(4.0, 40.0, route="iv")
    assert m_iv.route == Route.IV_BOLUS
    assert m_iv.f_eff == 1.0

    m_nasal = PKModel(4.0, 40.0, route="nasal")
    assert m_nasal.route == Route.INSUFFLATED

    m_inhaled = PKModel(4.0, 40.0, route="inhaled")
    assert m_inhaled.route == Route.SMOKED



def test_mdma_stacking_uses_saturable_elimination():
    drug = get_drug("mdma")
    model = drug.to_pk_model(70.0)
    assert model.is_saturable
    window = TherapeuticWindow(drug.mec, drug.mtc)
    a = assess_dose_stacking(100.0, 2.0, model, window, redose_amount=100.0)
    t = np.linspace(0.0, 24.0, 2401)
    superposed = single_dose_concentration(t + 2.0, 100.0, model) + single_dose_concentration(t, 100.0, model)
    assert a.stacked_peak > superposed.max()


def test_mdma_safe_interval_not_shorter_than_linear_estimate():
    drug = get_drug("mdma")
    model = drug.to_pk_model(70.0)
    linear = PKModel(model.t_half_h, model.vd_l, model.bioavailability, model.ka_per_h, model.route)
    window = TherapeuticWindow(drug.mec, drug.mtc)
    assert compute_safe_redosing_interval_h(100.0, model, window, redose_amount=60.0) >= \
        compute_safe_redosing_interval_h(100.0, linear, window, redose_amount=60.0)
