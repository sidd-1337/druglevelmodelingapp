import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from core.pk_engine import (
    DoseEvent,
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
