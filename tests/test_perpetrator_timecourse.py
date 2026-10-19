"""Time course of an interaction factor: onset, offset, induction lag, calibration."""
import math

import numpy as np
import pytest

from core.perpetrator_timecourse import (
    ActingFactor,
    PerpetratorCourse,
    apply_factor_course,
    combined_factor_at,
    effect_fraction,
    enzyme_activity_curve,
)

DAY = 24.0
T = np.linspace(0.0, 21 * DAY, 21 * 96 + 1)


def _at(curve, day):
    return float(curve[T >= day * DAY][0])


def test_inhibitor_fraction_builds_from_zero_to_one():
    c = PerpetratorCourse("CYP2D6", 0.07, t_half_h=21.0, induces=False)
    g = effect_fraction(c, T)
    assert g[0] == pytest.approx(0.0, abs=1e-6)
    assert _at(g, 14) == pytest.approx(1.0, abs=0.01)  # ~16 half-lives in: fully accumulated


def test_inhibitor_starts_late_and_stops_early():
    started = effect_fraction(PerpetratorCourse("CYP2D6", 0.07, 21.0, False, start_h=7 * DAY), T)
    assert _at(started, 5) == pytest.approx(0.0, abs=1e-6)  # nothing before it is taken
    assert _at(started, 14) > 0.7
    stopped = effect_fraction(PerpetratorCourse("CYP2D6", 0.07, 21.0, False, stop_h=10 * DAY), T)
    assert _at(stopped, 10) == pytest.approx(1.0, abs=0.02)
    assert _at(stopped, 20) < 0.05  # washed out ~10 half-lives after stopping


def test_enzyme_activity_curve_reaches_the_residual():
    c = PerpetratorCourse("CYP2D6", 0.07, 21.0, False)
    a = enzyme_activity_curve(c, T)
    assert a[0] == pytest.approx(1.0)
    assert _at(a, 14) == pytest.approx(0.07, abs=0.005)


def test_induction_lags_the_inducer_and_needs_turnover():
    c = PerpetratorCourse("CYP3A4", 26.0, t_half_h=3.5, induces=True)
    kdeg = math.log(2) / 36.0  # enzyme turnover much slower than the inducer's 3.5 h
    g = effect_fraction(c, T, turnover_per_h=kdeg)
    # The inducer itself is at steady state within a day, but the effect is still building.
    assert _at(g, 1) < 0.6
    assert _at(g, 2) < _at(g, 5) < _at(g, 14)
    assert _at(g, 14) == pytest.approx(1.0, abs=0.03)
    with pytest.raises(ValueError):
        effect_fraction(c, T)  # inducer with no turnover rate


def test_induction_recovers_slowly_after_stopping():
    c = PerpetratorCourse("CYP3A4", 26.0, 3.5, True, stop_h=14 * DAY)
    t = np.linspace(0.0, 30 * DAY, 30 * 96 + 1)
    g = effect_fraction(c, t, turnover_per_h=math.log(2) / 36.0)
    at = lambda d: float(g[t >= d * DAY][0])
    assert at(14) == pytest.approx(1.0, abs=0.03)
    assert at(16) > 0.3  # still induced two days after stopping - the slow part
    assert at(29) < 0.05


def test_apply_factor_course_reproduces_the_factor_at_full_effect():
    base = np.full_like(T, 2.0)
    full = np.ones_like(T)
    assert np.allclose(apply_factor_course(base, [ActingFactor(1.6, full)]), base * 1.6)
    assert np.allclose(apply_factor_course(base, [ActingFactor(1.6, np.zeros_like(T))]), base)
    # Two perpetrators multiply.
    both = apply_factor_course(base, [ActingFactor(1.6, full), ActingFactor(0.5, full)])
    assert np.allclose(both, base * 1.6 * 0.5)


def test_combined_factor_at_is_the_product():
    full = np.ones_like(T)
    assert combined_factor_at([ActingFactor(1.6, full), ActingFactor(2.0, full)]) == pytest.approx(3.2)
    assert combined_factor_at([ActingFactor(1.6, np.zeros_like(T))]) == pytest.approx(1.0)


@pytest.mark.parametrize("kwargs", [
    {"residual": 1.5, "induces": False},   # inhibitor cannot raise activity
    {"residual": 0.5, "induces": True},    # inducer cannot lower it
    {"residual": -0.1, "induces": False},  # negative
    {"t_half_h": 0.0, "induces": False},   # non-positive half-life
])
def test_validation_rejects_inconsistent_courses(kwargs):
    base = {"enzyme": "CYP2D6", "residual": 0.5, "t_half_h": 10.0, "induces": False}
    with pytest.raises(ValueError):
        PerpetratorCourse(**{**base, **kwargs})


def test_stop_before_start_rejected():
    with pytest.raises(ValueError):
        PerpetratorCourse("CYP2D6", 0.5, 10.0, False, start_h=48.0, stop_h=24.0)


def test_already_taking_means_full_effect_from_the_start():
    inhib = PerpetratorCourse("CYP2D6", 0.07, 120.0, False, start_h=None)
    assert effect_fraction(inhib, T)[0] == pytest.approx(1.0)
    induc = PerpetratorCourse("CYP3A4", 26.0, 3.5, True, start_h=None)
    g = effect_fraction(induc, T, turnover_per_h=math.log(2) / 36.0)
    assert g[0] == pytest.approx(1.0) and _at(g, 10) == pytest.approx(1.0, abs=1e-6)


def test_offset_half_life_slows_only_the_fade():
    short = PerpetratorCourse("CYP2D6", 0.07, 21.0, False, start_h=None, stop_h=2 * DAY)
    long = PerpetratorCourse("CYP2D6", 0.07, 21.0, False, start_h=None, stop_h=2 * DAY, offset_t_half_h=220.0)
    gs, gl = effect_fraction(short, T), effect_fraction(long, T)
    assert _at(gs, 1) == _at(gl, 1) == pytest.approx(1.0)  # identical while still taken
    assert _at(gs, 9) < 0.01 and _at(gl, 9) > 0.55  # norfluoxetine-like tail: still largely active a week later
