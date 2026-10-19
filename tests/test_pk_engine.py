import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from dataclasses import replace

from core.pk_engine import (
    DoseEvent,
    Metabolite,
    PKModel,
    Route,
    accumulation_index,
    apply_missed_dose,
    auc_trapezoid,
    clearance,
    cmax_single_dose,
    concentration_profile,
    elimination_rate_constant,
    regular_schedule,
    single_dose_concentration,
    steady_state_concentration,
    steady_state_summary,
    time_grid,
    time_to_steady_state,
    time_to_steady_state_range,
    tmax_single_dose,
)

ORAL = PKModel(t_half_h=4.0, vd_l=40.0, bioavailability=0.8, ka_per_h=1.2)
IV = PKModel(t_half_h=4.0, vd_l=40.0, route=Route.IV_BOLUS)
DOSE = 500.0


def test_ke_and_clearance():
    assert elimination_rate_constant(2.0) == pytest.approx(math.log(2) / 2)
    assert clearance(0.1, 30.0) == pytest.approx(3.0)
    assert ORAL.clearance_l_per_h == pytest.approx(ORAL.ke * 40.0)


def test_bateman_matches_ode_solution():
    ke, ka = ORAL.ke, ORAL.ka_per_h

    def rhs(_, y):  # y = [gut amount, central amount]
        return [-ka * y[0], ka * y[0] - ke * y[1]]

    t = np.linspace(0, 48, 200)
    sol = solve_ivp(rhs, (0, 48), [ORAL.f_eff * DOSE, 0.0], t_eval=t, rtol=1e-10, atol=1e-12)
    expected = sol.y[1] / ORAL.vd_l
    np.testing.assert_allclose(single_dose_concentration(t, DOSE, ORAL), expected, rtol=1e-6, atol=1e-9)


def test_bateman_zero_at_t0_and_before():
    assert float(single_dose_concentration(0.0, DOSE, ORAL)) == 0.0
    assert float(single_dose_concentration(-3.0, DOSE, ORAL)) == 0.0


def test_tmax_and_cmax_analytic():
    ka, ke = ORAL.ka_per_h, ORAL.ke
    assert tmax_single_dose(ORAL) == pytest.approx(math.log(ka / ke) / (ka - ke))
    t = np.linspace(0, 24, 200001)
    c = single_dose_concentration(t, DOSE, ORAL)
    assert t[np.argmax(c)] == pytest.approx(tmax_single_dose(ORAL), abs=1e-3)
    assert c.max() == pytest.approx(cmax_single_dose(DOSE, ORAL), rel=1e-8)


def test_flip_flop_limit_is_continuous():
    ke = ORAL.ke
    equal = PKModel(4.0, 40.0, 0.8, ka_per_h=ke)
    near = PKModel(4.0, 40.0, 0.8, ka_per_h=ke * (1 + 1e-5))
    t = np.linspace(0, 30, 50)
    np.testing.assert_allclose(
        single_dose_concentration(t, DOSE, equal), single_dose_concentration(t, DOSE, near), rtol=1e-4
    )
    assert tmax_single_dose(equal) == pytest.approx(1 / ke)


def test_iv_bolus_monoexponential():
    c0 = DOSE / IV.vd_l
    assert float(single_dose_concentration(0.0, DOSE, IV)) == pytest.approx(c0)
    assert float(single_dose_concentration(4.0, DOSE, IV)) == pytest.approx(c0 / 2)
    assert IV.f_eff == 1.0 and tmax_single_dose(IV) == 0.0


def test_iv_ignores_bioavailability():
    iv_low_f = PKModel(4.0, 40.0, bioavailability=0.1, route=Route.IV_BOLUS)
    assert float(single_dose_concentration(1.0, DOSE, iv_low_f)) == pytest.approx(
        float(single_dose_concentration(1.0, DOSE, IV))
    )


@pytest.mark.parametrize("model", [ORAL, IV])
def test_auc_inf_equals_fd_over_cl(model):
    t = np.linspace(0, 400, 400001)
    c = single_dose_concentration(t, DOSE, model)
    expected = model.f_eff * DOSE / model.clearance_l_per_h
    assert auc_trapezoid(t, c) == pytest.approx(expected, rel=1e-4)


def test_auc_trapezoid_bounds_and_exactness():
    t = np.array([0.0, 1.0, 2.0, 3.0])
    c = np.array([0.0, 2.0, 2.0, 0.0])
    assert auc_trapezoid(t, c) == pytest.approx(4.0)
    assert auc_trapezoid(t, c, 0.5, 2.5) == pytest.approx(0.5 * (1 + 2) / 2 + 2 + 0.5 * (2 + 1) / 2)
    with pytest.raises(ValueError):
        auc_trapezoid(t, c, 2.0, 1.0)


def test_regular_schedule():
    assert [e.time_h for e in regular_schedule(DOSE, 8, 24)] == [0, 8, 16]
    assert len(regular_schedule(DOSE, 7, 24)) == 4
    assert len(regular_schedule(DOSE, 24, 24 * 5)) == 5
    with pytest.raises(ValueError):
        regular_schedule(DOSE, 0, 24)


def test_superposition_is_linear():
    t = time_grid(48)
    two = concentration_profile(t, [DoseEvent(0, DOSE), DoseEvent(12, DOSE)], ORAL)
    expected = single_dose_concentration(t, DOSE, ORAL) + single_dose_concentration(t - 12, DOSE, ORAL)
    np.testing.assert_allclose(two, expected)


@pytest.mark.parametrize("model", [ORAL, IV])
@pytest.mark.parametrize("tau", [6.0, 12.0, 24.0])
def test_steady_state_matches_simulation(model, tau):
    duration = (math.ceil(40 * model.t_half_h / tau) + 2) * tau  # whole intervals, >> 5 half-lives
    events = regular_schedule(DOSE, tau, duration)
    t = np.linspace(0, duration, 200001)
    c = concentration_profile(t, events, model)
    last = t >= events[-1].time_h
    ss = steady_state_summary(DOSE, tau, model)
    assert c[last].max() == pytest.approx(ss.cmax_ss, rel=2e-4)
    assert c[-1] == pytest.approx(ss.cmin_ss, rel=1e-6)  # trough = end of last interval
    # Ave over a full interval equals F*D/(CL*tau)
    n_last = len(events) - 1
    start = events[n_last - 1].time_h
    auc = auc_trapezoid(t, c, start, start + tau)
    # IV bolus jumps at every dose; trapezoid over the grid step straddling it adds a tiny error.
    tol = 1e-3 if model.route is Route.IV_BOLUS else 1e-5
    assert auc / tau == pytest.approx(ss.cavg_ss, rel=tol)


def test_steady_state_curve_is_periodic():
    tau = 8.0
    c0 = steady_state_concentration(0.0, DOSE, tau, ORAL)
    ctau = steady_state_concentration(tau, DOSE, tau, ORAL)
    assert float(c0) == pytest.approx(float(ctau))


def test_accumulation_index():
    ke = ORAL.ke
    assert accumulation_index(ke, 4.0) == pytest.approx(2.0)  # tau == t1/2 -> R = 2
    ss = steady_state_summary(DOSE, 4.0, IV)
    assert ss.accumulation_index == pytest.approx(2.0)
    # for IV bolus Cmax,ss = R * C0 and Cmin,ss = Cmax,ss - C0
    c0 = DOSE / IV.vd_l
    assert ss.cmax_ss == pytest.approx(2 * c0)
    assert ss.cmin_ss == pytest.approx(c0)


def test_time_to_steady_state():
    assert time_to_steady_state(10.0) == pytest.approx(45.0)
    assert time_to_steady_state_range(10.0) == (40.0, 50.0)
    # after 4.32 half-lives ~95 % of steady state is reached
    ss = steady_state_summary(DOSE, 24.0, IV)
    assert ss.time_to_ss_h == pytest.approx(4.5 * IV.t_half_h)


def test_ordering_of_steady_state_levels():
    ss = steady_state_summary(DOSE, 8.0, ORAL)
    assert ss.cmin_ss < ss.cavg_ss < ss.cmax_ss
    assert 0 < ss.tmax_ss_h < 8.0
    assert ss.fluctuation_pct > 0


def test_missed_dose_skip_and_delay():
    events = regular_schedule(DOSE, 8, 48)
    skipped = apply_missed_dose(events, 2)
    assert len(skipped) == len(events) - 1 and 16 not in [e.time_h for e in skipped]
    delayed = apply_missed_dose(events, 2, delay_h=3)
    assert [e.time_h for e in delayed][2] == 19 and len(delayed) == len(events)
    with pytest.raises(IndexError):
        apply_missed_dose(events, 99)


def test_missed_dose_lowers_then_recovers():
    tau = 8.0
    events = regular_schedule(DOSE, tau, 200)
    t = time_grid(200)
    ref = concentration_profile(t, events, ORAL)
    skipped = concentration_profile(t, apply_missed_dose(events, 10), ORAL)
    after = t > events[10].time_h
    assert skipped[after].min() <= ref[after].min()
    assert np.all(skipped <= ref + 1e-12)
    np.testing.assert_allclose(skipped[t < events[10].time_h], ref[t < events[10].time_h])


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(t_half_h=0, vd_l=10),
        dict(t_half_h=2, vd_l=-1),
        dict(t_half_h=2, vd_l=10, bioavailability=1.5),
        dict(t_half_h=2, vd_l=10, ka_per_h=0),
        dict(t_half_h=2, vd_l=10, ka_per_h=0, route=Route.INSUFFLATED),
        dict(t_half_h=2, vd_l=10, ka_per_h=0, route=Route.SMOKED),
    ],
)
def test_model_validation(kwargs):
    with pytest.raises(ValueError):
        PKModel(**kwargs)


def test_routes_pk_kinetics_ordering():
    dose = 200.0
    t_half = 3.0
    vd = 30.0
    oral = PKModel(t_half_h=t_half, vd_l=vd, bioavailability=0.6, ka_per_h=1.2, route=Route.ORAL)
    nasal = PKModel(t_half_h=t_half, vd_l=vd, bioavailability=0.75, ka_per_h=3.5, route=Route.INSUFFLATED)
    smoked = PKModel(t_half_h=t_half, vd_l=vd, bioavailability=0.55, ka_per_h=8.0, route=Route.SMOKED)
    iv = PKModel(t_half_h=t_half, vd_l=vd, route=Route.IV_BOLUS)

    tmax_iv = tmax_single_dose(iv)
    tmax_smk = tmax_single_dose(smoked)
    tmax_nas = tmax_single_dose(nasal)
    tmax_or = tmax_single_dose(oral)

    assert tmax_iv == 0.0
    assert 0.0 < tmax_smk < tmax_nas < tmax_or

    # AUC conservation: AUC_inf = F * Dose / CL
    t = np.linspace(0, 300, 300001)
    for model in (nasal, smoked):
        c = single_dose_concentration(t, dose, model)
        expected_auc = model.f_eff * dose / model.clearance_l_per_h
        assert auc_trapezoid(t, c) == pytest.approx(expected_auc, rel=1e-3)



# --------------------------------------------------------- saturable elimination
SAT = PKModel(t_half_h=7.5, vd_l=390.0, bioavailability=0.8, ka_per_h=1.5, km=1.2, saturable_fraction=0.76)


def _unsaturated(model: PKModel) -> PKModel:
    """Same model with km so large that elimination is effectively linear."""
    return PKModel(model.t_half_h, model.vd_l, model.bioavailability, model.ka_per_h, model.route,
                   km=1e9, saturable_fraction=1.0)


@pytest.mark.parametrize("model", [ORAL, IV])
def test_saturable_path_matches_analytic_model_when_unsaturated(model):
    sat = _unsaturated(model)
    t = time_grid(72.0)
    events = regular_schedule(DOSE, 12.0, 72.0)
    np.testing.assert_allclose(concentration_profile(t, events, sat), concentration_profile(t, events, model),
                               rtol=1e-5, atol=1e-8)
    assert cmax_single_dose(DOSE, sat) == pytest.approx(cmax_single_dose(DOSE, model), rel=1e-4)
    assert tmax_single_dose(sat, DOSE) == pytest.approx(tmax_single_dose(model), abs=0.01)
    ss_sat, ss_lin = steady_state_summary(DOSE, 12.0, sat), steady_state_summary(DOSE, 12.0, model)
    assert ss_sat.cavg_ss == pytest.approx(ss_lin.cavg_ss, rel=1e-3)
    assert ss_sat.cmin_ss == pytest.approx(ss_lin.cmin_ss, rel=1e-3)


def test_saturable_exposure_more_than_dose_proportional():
    t = np.linspace(0.0, 150.0, 6001)
    auc = {d: auc_trapezoid(t, single_dose_concentration(t, d, SAT)) for d in (50.0, 100.0, 200.0)}
    assert auc[100.0] / auc[50.0] > 2.0
    assert auc[200.0] / auc[100.0] > auc[100.0] / auc[50.0]
    assert tmax_single_dose(SAT, 200.0) > tmax_single_dose(SAT, 50.0)


def test_saturable_second_dose_exceeds_superposition():
    t = np.linspace(0.0, 24.0, 2401)
    stacked = concentration_profile(t, [DoseEvent(0.0, 100.0), DoseEvent(3.0, 100.0)], SAT)
    superposed = single_dose_concentration(t, 100.0, SAT) + single_dose_concentration(t - 3.0, 100.0, SAT)
    assert stacked.max() > superposed.max()
    assert auc_trapezoid(t, stacked) > 1.02 * auc_trapezoid(t, superposed)


def test_saturable_steady_state_is_consistent_with_profile():
    ss = steady_state_summary(100.0, 24.0, SAT)
    assert ss.cmin_ss < ss.cavg_ss < ss.cmax_ss
    t = np.linspace(0.0, 24.0, 97)
    np.testing.assert_allclose(steady_state_concentration(t + 24.0, 100.0, 24.0, SAT),
                               steady_state_concentration(t, 100.0, 24.0, SAT), rtol=1e-5)


def test_saturable_iv_steady_state_dose_time_is_trough():
    iv = PKModel(7.5, 390.0, route=Route.IV_BOLUS, km=1.2, saturable_fraction=0.76)
    ss = steady_state_summary(100.0, 12.0, iv)
    c = steady_state_concentration(np.array([0.0, 12.0, 24.0]), 100.0, 12.0, iv)
    assert c[0] == pytest.approx(ss.cmax_ss) and c[1] == pytest.approx(ss.cmin_ss) == c[2]


def test_with_half_life_keeps_saturable_parameters():
    m = SAT.with_half_life(10.0)
    assert (m.km, m.saturable_fraction, m.t_half_h) == (SAT.km, SAT.saturable_fraction, 10.0)


@pytest.mark.parametrize("kwargs", [dict(km=0.0, saturable_fraction=0.5), dict(km=1.0, saturable_fraction=0.0),
                                    dict(km=1.0, saturable_fraction=1.5)])
def test_saturable_validation(kwargs):
    with pytest.raises(ValueError):
        PKModel(4.0, 40.0, **kwargs)


# ------------------------------------------------- two compartments, lag, inactivation
TWO = PKModel(t_half_h=6.0, vd_l=60.0, bioavailability=0.8, ka_per_h=1.5, v1_fraction=0.3, q_l_per_h=20.0)


def _ode_two_compartment(model: PKModel, dose: float, t: np.ndarray) -> np.ndarray:
    k10, k12, k21 = model.micro_constants
    ka = 0.0 if model.route is Route.IV_BOLUS else model.ka_per_h

    def rhs(_t, y):
        return [-ka * y[0], ka * y[0] - (k10 + k12) * y[1] + k21 * y[2], k12 * y[1] - k21 * y[2]]
    y0 = [0.0, dose, 0.0] if model.route is Route.IV_BOLUS else [model.f_eff * dose, 0.0, 0.0]
    sol = solve_ivp(rhs, (0, t[-1]), y0, t_eval=t, rtol=1e-10, atol=1e-12)
    return sol.y[1] / model.v1_l


@pytest.mark.parametrize("route", [Route.ORAL, Route.IV_BOLUS])
def test_two_compartment_closed_form_matches_ode(route):
    model = PKModel(6.0, 60.0, 0.8, 1.5, route, v1_fraction=0.3, q_l_per_h=20.0)
    t = np.linspace(0.0, 48.0, 481)
    np.testing.assert_allclose(single_dose_concentration(t, DOSE, model), _ode_two_compartment(model, DOSE, t),
                               rtol=1e-6, atol=1e-9)


def test_two_compartment_keeps_clearance_and_auc_but_raises_early_peak():
    one = PKModel(6.0, 60.0, 0.8, 1.5)
    t = np.linspace(0.0, 400.0, 400001)
    assert auc_trapezoid(t, single_dose_concentration(t, DOSE, TWO)) == pytest.approx(0.8 * DOSE / one.clearance_l_per_h,
                                                                                      rel=1e-3)
    assert cmax_single_dose(DOSE, TWO) > cmax_single_dose(DOSE, one)
    assert TWO.terminal_half_life_h > one.t_half_h
    a, b = TWO.disposition_rates
    assert a > b > 0


def test_two_compartment_steady_state_matches_long_superposition():
    tau = 12.0
    events = regular_schedule(DOSE, tau, 50 * tau)
    t_last = np.linspace(49 * tau, 50 * tau, 241)
    np.testing.assert_allclose(steady_state_concentration(t_last - 49 * tau, DOSE, tau, TWO),
                               concentration_profile(t_last, events, TWO), rtol=1e-6)
    ss = steady_state_summary(DOSE, tau, TWO)
    assert ss.cavg_ss == pytest.approx(0.8 * DOSE / (TWO.clearance_l_per_h * tau))
    assert ss.time_to_ss_h == pytest.approx(4.5 * TWO.terminal_half_life_h)


@pytest.mark.parametrize("model", [ORAL, TWO])
def test_lag_shifts_the_whole_profile(model):
    lagged = PKModel(**{**model.__dict__, "tlag_h": 0.5})
    t = np.linspace(0.0, 24.0, 241)
    np.testing.assert_allclose(single_dose_concentration(t, DOSE, lagged), single_dose_concentration(t - 0.5, DOSE, model),
                               atol=1e-12)
    assert tmax_single_dose(lagged) == pytest.approx(tmax_single_dose(model) + 0.5, abs=1e-3)
    assert cmax_single_dose(DOSE, lagged) == pytest.approx(cmax_single_dose(DOSE, model), rel=1e-6)
    ss, ss_lag = steady_state_summary(DOSE, 12.0, model), steady_state_summary(DOSE, 12.0, lagged)
    assert ss_lag.cavg_ss == pytest.approx(ss.cavg_ss) and ss_lag.cmin_ss == pytest.approx(ss.cmin_ss, rel=1e-6)


def test_lag_ignored_for_iv_bolus():
    iv = PKModel(4.0, 40.0, route=Route.IV_BOLUS, tlag_h=1.0)
    assert float(single_dose_concentration(0.0, DOSE, iv)) == pytest.approx(DOSE / 40.0)


MBI = PKModel(t_half_h=2.0, vd_l=300.0, bioavailability=0.8, ka_per_h=1.5, inactivation_fraction=0.8,
              kinact_per_h=1.0, ki=0.15, kdeg_per_h=math.log(2) / 46.6)


def test_enzyme_inactivation_and_recovery():
    from core.pk_engine import enzyme_activity
    ev = [DoseEvent(0.0, 100.0)]
    e = enzyme_activity(np.array([0.0, 6.0, 46.6 + 24.0, 240.0]), ev, MBI)
    assert e[0] == pytest.approx(1.0) and e[1] < 0.3
    assert e[1] < e[2] < e[3] <= 1.0  # recovers with the enzyme turnover half-life
    assert MBI.terminal_half_life_h == pytest.approx(2.0 / (1 - 0.8))


def test_enzyme_inactivation_makes_second_dose_larger():
    t = np.linspace(0.0, 48.0, 4801)
    c = concentration_profile(t, [DoseEvent(0.0, 100.0), DoseEvent(24.0, 100.0)], MBI)
    first, second = t < 24, t >= 24
    assert auc_trapezoid(t[second], c[second]) > 1.3 * auc_trapezoid(t[first], c[first])


def test_enzyme_inactivation_validation():
    with pytest.raises(ValueError):
        PKModel(2.0, 300.0, inactivation_fraction=0.8, kinact_per_h=0.0, ki=0.1, kdeg_per_h=0.01)
    with pytest.raises(ValueError):
        PKModel(2.0, 300.0, km=1.0, saturable_fraction=0.5, inactivation_fraction=0.8, kinact_per_h=1.0, ki=0.1,
                kdeg_per_h=0.01)


# ------------------------------------------------------------- metabolites
PARENT_ONLY = PKModel(44.0, 83.0, 0.93, 4.0)
NORDAZEPAM = Metabolite("nordazepam", t_half_h=65.0, vd_l=87.0, formation_fraction=0.5,
                        molar_ratio=1.052, reported="parent")
WITH_MET = PKModel(44.0, 83.0, 0.93, 4.0, metabolite=NORDAZEPAM)


def test_metabolite_needs_integration_without_being_nonlinear():
    # A linear parent plus metabolite is still a sum of exponentials; it just cannot be superposed.
    assert WITH_MET.needs_integration and not WITH_MET.is_nonlinear
    assert not PARENT_ONLY.needs_integration
    assert hash(WITH_MET)  # must stay hashable: the engine lru_caches on the model


def test_reporting_the_parent_is_numerically_a_no_op():
    t = np.linspace(0.0, 400.0, 2001)
    np.testing.assert_allclose(single_dose_concentration(t, 10.0, WITH_MET),
                               single_dose_concentration(t, 10.0, PARENT_ONLY), atol=1e-7)
    assert cmax_single_dose(10.0, WITH_MET) == pytest.approx(cmax_single_dose(10.0, PARENT_ONLY), rel=1e-5)
    ss_a = steady_state_summary(10.0, 24.0, WITH_MET)
    ss_b = steady_state_summary(10.0, 24.0, PARENT_ONLY)
    assert ss_a.cavg_ss == pytest.approx(ss_b.cavg_ss, rel=1e-3)


def test_metabolite_mass_balance():
    """Everything formed must be eliminated: AUC_m = (f / molar) * CL_p * AUC_p / (ke_m * V_m)."""
    from core.pk_engine import metabolite_profile
    t = np.linspace(0.0, 2000.0, 20001)
    parent = single_dose_concentration(t, 10.0, WITH_MET)
    met = metabolite_profile(t, [DoseEvent(0.0, 10.0)], WITH_MET)
    expected = (NORDAZEPAM.formation_fraction / NORDAZEPAM.molar_ratio
                * WITH_MET.clearance_l_per_h * auc_trapezoid(t, parent)
                / (NORDAZEPAM.ke * NORDAZEPAM.vd_l))
    assert auc_trapezoid(t, met) == pytest.approx(expected, rel=1e-3)


def test_metabolite_is_formed_in_its_own_mass_units():
    """A metabolite half the parent's molecular weight forms half the mass, never double.

    The engine once multiplied by MW_parent / MW_metabolite, which is neither the
    metabolite's mass nor parent-equivalents: nordazepam came out 11 % high and
    9-hydroxyrisperidone 7 % low against the ng/mL values labs report.
    """
    from core.pk_engine import metabolite_profile
    light = replace(NORDAZEPAM, formation_fraction=1.0, molar_ratio=2.0)
    model = PKModel(44.0, 83.0, 0.93, 4.0, metabolite=light)
    t = np.linspace(0.0, 2000.0, 20001)
    formed = light.ke * light.vd_l * auc_trapezoid(t, metabolite_profile(t, [DoseEvent(0.0, 10.0)], model))
    eliminated = model.clearance_l_per_h * auc_trapezoid(t, single_dose_concentration(t, 10.0, model))
    assert formed / eliminated == pytest.approx(0.5, rel=1e-3)


def test_a_slow_metabolite_peaks_long_after_the_parent():
    from core.pk_engine import metabolite_profile
    t = np.linspace(0.0, 400.0, 4001)
    met = metabolite_profile(t, [DoseEvent(0.0, 10.0)], WITH_MET)
    assert t[int(met.argmax())] > 10 * tmax_single_dose(PARENT_ONLY)


@pytest.mark.parametrize("reported, expected_t_half", [("parent", 44.0), ("metabolite", 65.0), ("sum", 65.0)])
def test_terminal_half_life_follows_the_reported_quantity(reported, expected_t_half):
    m = PKModel(44.0, 83.0, 0.93, 4.0, metabolite=replace(NORDAZEPAM, reported=reported))
    assert m.terminal_half_life_h == pytest.approx(expected_t_half)
    assert m.reported_quantity == reported


def test_sum_exceeds_parent_and_tracks_potency():
    t = np.linspace(0.0, 400.0, 2001)
    as_sum = PKModel(44.0, 83.0, 0.93, 4.0, metabolite=replace(NORDAZEPAM, reported="sum"))
    half_potency = replace(as_sum, metabolite=replace(NORDAZEPAM, reported="sum", potency_ratio=0.5))
    c_parent = single_dose_concentration(t, 10.0, PARENT_ONLY)
    c_sum = single_dose_concentration(t, 10.0, as_sum)
    c_half = single_dose_concentration(t, 10.0, half_potency)
    assert auc_trapezoid(t, c_sum) > auc_trapezoid(t, c_half) > auc_trapezoid(t, c_parent)


@pytest.mark.parametrize("kwargs", [
    dict(t_half_h=0.0), dict(vd_l=-1.0), dict(formation_fraction=0.0),
    dict(formation_fraction=1.5), dict(molar_ratio=0.0), dict(reported="nonsense")])
def test_metabolite_validation(kwargs):
    with pytest.raises(ValueError):
        replace(NORDAZEPAM, **kwargs)


# ------------------------------------------------- coupled (shared-enzyme) profiles
from core.pk_engine import (
    INDUCTION,
    REVERSIBLE,
    ClampedPerpetrator,
    CoDrug,
    CoupledProfiles,
    EnzymeEffect,
    concentration_profile as _cp,
    coupled_profiles,
)
from core.enzymes import derived_exposure_factor

_VICTIM = PKModel(t_half_h=20.0, vd_l=90.0, bioavailability=1.0, ka_per_h=1.0)  # CYP2D6 substrate


def test_coupled_leaves_uncoupled_drugs_identical():
    # Two drugs that share no inactivated enzyme get exactly their standalone curves.
    t = time_grid(48.0)
    a = CoDrug(PARENT_ONLY, tuple(regular_schedule(100.0, 12.0, 48.0)))
    b = CoDrug(_VICTIM, tuple(regular_schedule(50.0, 24.0, 48.0)), {"CYP2D6": 0.6})
    res = coupled_profiles(t, [a, b])
    assert np.allclose(res.concentrations[0], _cp(t, a.events, PARENT_ONLY))
    assert np.allclose(res.concentrations[1], _cp(t, b.events, _VICTIM))
    assert res.enzyme_activity == {}


def test_coupled_inactivator_matches_the_single_drug_engine():
    from core.pk_engine import enzyme_activity
    t = time_grid(72.0)
    ev = tuple(regular_schedule(100.0, 24.0, 72.0))
    res = coupled_profiles(t, [CoDrug(MBI, ev, inactivates="CYP2D6")])
    assert np.allclose(res.concentrations[0], _cp(t, ev, MBI), rtol=1e-5, atol=1e-9)
    assert np.allclose(res.enzyme_activity["CYP2D6"], enzyme_activity(t, ev, MBI), rtol=1e-5, atol=1e-6)


def test_coupled_victim_follows_the_inactivator():
    # A CYP2D6 substrate taken with an inactivator is cleared more slowly -> higher exposure.
    t = time_grid(72.0)
    ev_v = tuple(regular_schedule(50.0, 24.0, 72.0))
    victim = CoDrug(_VICTIM, ev_v, {"CYP2D6": 0.6})
    perp = CoDrug(MBI, tuple(regular_schedule(100.0, 24.0, 72.0)), inactivates="CYP2D6")
    alone = _cp(t, ev_v, _VICTIM)
    with_perp = coupled_profiles(t, [victim, perp]).concentrations[0]
    assert auc_trapezoid(t, with_perp) > 1.1 * auc_trapezoid(t, alone)
    # A drug with no CYP2D6 share is untouched by the same perpetrator.
    inert = CoDrug(_VICTIM, ev_v, {})
    assert np.allclose(coupled_profiles(t, [inert, perp]).concentrations[0], alone)


def test_coupled_validates_inactivates_and_shares():
    with pytest.raises(ValueError):
        CoDrug(_VICTIM, (), inactivates="CYP2D6")  # names an enzyme but the model inactivates none
    with pytest.raises(ValueError):
        CoDrug(MBI, (), inactivates="")  # inactivates one but no enzyme named
    with pytest.raises(ValueError):
        CoDrug(_VICTIM, (), {"CYP2D6": 0.7, "CYP3A4": 0.5})  # shares sum > 1


def _last_interval_exposure(c, t, tau_h):
    """AUC over the final dosing interval - a steady-state exposure proxy."""
    mask = t >= t[-1] - tau_h
    return auc_trapezoid(t[mask], c[mask])


def test_clamped_reversible_inhibitor_matches_the_static_model():
    # A constant reversible inhibitor level must reproduce the steady-state static model
    # 1 / (fm*residual + 1-fm), with residual = 1/(1 + I/ki). This ties the dynamic
    # reversible term to core.enzymes.derived_exposure_factor.
    fm, ki, level = 0.6, 5.0, 5.0            # I = ki, so residual availability = 0.5
    tau, dur = 24.0, 24.0 * 25
    t = time_grid(dur)
    ev = tuple(regular_schedule(50.0, tau, dur))
    victim = CoDrug(_VICTIM, ev, {"CYP2D6": fm})
    perp = ClampedPerpetrator(level, (EnzymeEffect("CYP2D6", REVERSIBLE, ki=ki),))
    base = coupled_profiles(t, [victim]).concentrations[0]
    hit = coupled_profiles(t, [victim], clamped=[perp]).concentrations[0]
    ratio = _last_interval_exposure(hit, t, tau) / _last_interval_exposure(base, t, tau)
    residual = 1.0 / (1.0 + level / ki)
    assert ratio == pytest.approx(derived_exposure_factor(fm, residual), rel=5e-3)


def test_clamped_inducer_matches_the_static_model():
    # A constant inducer level must reproduce the same static relation with residual > 1.
    fm, emax, ec50, level = 0.6, 4.0, 5.0, 5.0   # induction term d*emax*I/(ec50+I) = 2 -> residual 3
    tau, dur = 24.0, 24.0 * 25
    t = time_grid(dur)
    ev = tuple(regular_schedule(50.0, tau, dur))
    victim = CoDrug(_VICTIM, ev, {"CYP3A4": fm})
    perp = ClampedPerpetrator(level, (EnzymeEffect("CYP3A4", INDUCTION, emax=emax, ec50=ec50),))
    turnover = {"CYP3A4": 0.0193}
    base = coupled_profiles(t, [victim], enzyme_turnover=turnover).concentrations[0]
    hit = coupled_profiles(t, [victim], clamped=[perp], enzyme_turnover=turnover).concentrations[0]
    ratio = _last_interval_exposure(hit, t, tau) / _last_interval_exposure(base, t, tau)
    residual = 1.0 + emax * level / (ec50 + level)
    assert ratio == pytest.approx(derived_exposure_factor(fm, residual), rel=5e-3)


def test_dynamic_reversible_perpetrator_raises_a_victim_and_spares_a_non_substrate():
    # A co-dosed drug that reversibly inhibits CYP2D6 raises a CYP2D6 victim's exposure,
    # and leaves a drug with no CYP2D6 share untouched.
    t = time_grid(72.0)
    ev_v = tuple(regular_schedule(50.0, 24.0, 72.0))
    perp = CoDrug(PARENT_ONLY, tuple(regular_schedule(100.0, 24.0, 72.0)),
                  effects=(EnzymeEffect("CYP2D6", REVERSIBLE, ki=2.0),))
    victim = CoDrug(_VICTIM, ev_v, {"CYP2D6": 0.6})
    alone = _cp(t, ev_v, _VICTIM)
    with_perp = coupled_profiles(t, [victim, perp]).concentrations[0]
    assert auc_trapezoid(t, with_perp) > 1.1 * auc_trapezoid(t, alone)
    inert = CoDrug(_VICTIM, ev_v, {})
    assert np.allclose(coupled_profiles(t, [inert, perp]).concentrations[0], alone)


def test_induction_without_a_turnover_is_rejected():
    t = time_grid(48.0)
    victim = CoDrug(_VICTIM, tuple(regular_schedule(50.0, 24.0, 48.0)), {"CYP3A4": 0.6})
    perp = ClampedPerpetrator(5.0, (EnzymeEffect("CYP3A4", INDUCTION, emax=4.0, ec50=5.0),))
    with pytest.raises(ValueError):
        coupled_profiles(t, [victim], clamped=[perp])  # no enzyme_turnover supplied


def test_enzyme_effect_and_clamp_validation():
    with pytest.raises(ValueError):
        EnzymeEffect("CYP2D6", REVERSIBLE, ki=0.0)
    with pytest.raises(ValueError):
        EnzymeEffect("CYP2D6", INDUCTION, emax=4.0, ec50=0.0)
    with pytest.raises(ValueError):
        EnzymeEffect("CYP2D6", "nonsense")
    with pytest.raises(ValueError):
        ClampedPerpetrator(5.0, ())
