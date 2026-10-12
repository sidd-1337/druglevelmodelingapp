import math

import numpy as np
import pytest

from core.drug_database import get_drug
from core.individualize import Observation, _apply_eta, fit_individual, individual_band
from core.pk_engine import concentration_profile, regular_schedule, time_grid

LI = get_drug("lithium")
POP = LI.to_pk_model(70)
EVENTS = regular_schedule(LI.amount_from_dose_mg(400), 12, 168)


def _measure(true_model, times):
    c = concentration_profile(np.asarray(times, dtype=float), EVENTS, true_model)
    return [Observation(t, float(v)) for t, v in zip(times, c)]


def test_level_matching_population_returns_population_model():
    fit = fit_individual(POP, EVENTS, _measure(POP, [156.0]), LI.cv_cl)
    assert fit.eta == pytest.approx((0.0, 0.0), abs=1e-4)
    assert fit.model.t_half_h == pytest.approx(POP.t_half_h, rel=1e-3)


def test_estimate_moves_toward_truth_and_shrinks_to_prior():
    true = _apply_eta(POP, [math.log(0.7), 0.0])
    fit = fit_individual(POP, EVENTS, _measure(true, [156.0]), LI.cv_cl)
    assert 0.7 < fit.cl_ratio < 1.0
    assert fit.model.t_half_h > POP.t_half_h


def test_precise_rich_data_recovers_true_parameters():
    true = _apply_eta(POP, [math.log(0.6), math.log(1.2)])
    times = [2.0, 6.0, 12.0, 60.0, 146.0, 150.0, 156.0, 166.0]
    fit = fit_individual(POP, EVENTS, _measure(true, times), LI.cv_cl, residual_cv=0.01)
    assert fit.cl_ratio == pytest.approx(0.6, rel=0.03)
    assert fit.v_ratio == pytest.approx(1.2, rel=0.05)
    np.testing.assert_allclose(fit.predicted, [o.concentration for o in fit.observations], rtol=0.02)


def test_more_samples_narrow_the_posterior():
    true = _apply_eta(POP, [math.log(0.8), 0.0])
    one = fit_individual(POP, EVENTS, _measure(true, [156.0]), LI.cv_cl)
    three = fit_individual(POP, EVENTS, _measure(true, [134.0, 146.0, 156.0]), LI.cv_cl)
    assert three.cov[0][0] < one.cov[0][0] < math.log1p(LI.cv_cl ** 2)


def test_band_brackets_individual_profile():
    true = _apply_eta(POP, [math.log(0.8), 0.0])
    fit = fit_individual(POP, EVENTS, _measure(true, [156.0]), LI.cv_cl)
    t = time_grid(168)
    lo, hi = individual_band(fit, t, EVENTS)
    c = concentration_profile(t, EVENTS, fit.model)
    late = t > 24
    assert np.all(lo[late] < c[late]) and np.all(c[late] < hi[late])


def test_saturable_model_can_be_individualised():
    d = get_drug("mdma")
    pop = d.to_pk_model(70)
    events = regular_schedule(100.0, 24, 24)
    true = _apply_eta(pop, [math.log(0.7), 0.0])
    obs = [Observation(t, float(c)) for t, c in zip((4.0, 10.0), concentration_profile(np.array([4.0, 10.0]), events, true))]
    fit = fit_individual(pop, events, obs, d.cv_cl)
    assert fit.model.is_saturable and fit.cl_ratio < 1.0


@pytest.mark.parametrize("obs", [[], [Observation(10.0, 0.0)], [Observation(0.0, 0.5)]])
def test_invalid_observations_rejected(obs):
    with pytest.raises(ValueError):
        fit_individual(POP, EVENTS, obs, LI.cv_cl)
