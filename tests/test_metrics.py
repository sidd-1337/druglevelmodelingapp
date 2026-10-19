import numpy as np
import pytest

from core.metrics import (
    TherapeuticWindow,
    assess_tdm,
    compute_metrics,
    first_crossing_time,
    fluctuation_percent,
    peak_trough_in_interval,
    recovery_time,
    terminal_half_life,
)
from core.pk_engine import (
    PKModel,
    apply_missed_dose,
    concentration_profile,
    regular_schedule,
    single_dose_concentration,
    time_grid,
    tmax_single_dose,
)

MODEL = PKModel(t_half_h=3.0, vd_l=50.0, bioavailability=0.9, ka_per_h=1.5)
WINDOW = TherapeuticWindow(mec=2.0, mtc=8.0)


def test_compute_metrics_single_dose():
    t = time_grid(48)
    c = single_dose_concentration(t, 400.0, MODEL)
    m = compute_metrics(t, c, MODEL)
    assert m.tmax_h == pytest.approx(tmax_single_dose(MODEL), abs=0.05)
    assert m.clearance_l_per_h == pytest.approx(MODEL.ke * 50.0)
    assert m.t_half_h == 3.0
    assert 0 < m.auc_0_24 < MODEL.f_eff * 400.0 / MODEL.clearance_l_per_h


def test_auc_window_clipped_to_profile_length():
    t = time_grid(12)
    c = single_dose_concentration(t, 400.0, MODEL)
    assert compute_metrics(t, c, MODEL).auc_0_24 > 0


def test_terminal_half_life_recovers_input():
    t = np.linspace(0, 60, 600)
    c = single_dose_concentration(t, 400.0, MODEL)
    assert terminal_half_life(t, c, n_points=100) == pytest.approx(3.0, rel=1e-3)
    with pytest.raises(ValueError):
        terminal_half_life(t[:5], np.zeros(5))


def test_fluctuation_and_peak_trough():
    assert fluctuation_percent(10, 4, 6) == pytest.approx(100.0)
    t = np.linspace(0, 10, 1001)
    c = 5 + 3 * np.sin(2 * np.pi * t / 10)
    peak, trough = peak_trough_in_interval(t, c, 0, 10)
    assert peak == pytest.approx(8, abs=1e-3) and trough == pytest.approx(2, abs=1e-3)


def test_window_classification():
    zones = WINDOW.classify(np.array([0.5, 2.0, 5.0, 8.0, 9.0]))
    assert list(zones) == [0, 1, 1, 1, 2]
    with pytest.raises(ValueError):
        TherapeuticWindow(5, 5)


def test_first_crossing_interpolates():
    t = np.array([0.0, 1.0, 2.0])
    c = np.array([0.0, 4.0, 8.0])
    assert first_crossing_time(t, c, 6.0) == pytest.approx(1.5)
    assert first_crossing_time(t, c, 100.0) is None


def test_assess_toxic_reports_time():
    t = np.linspace(0, 10, 1001)
    c = np.where(t > 6.5, 10.0, 3.0)
    res = assess_tdm(t, c, WINDOW)
    assert res.status == "toxic"
    assert res.first_toxic_time_h == pytest.approx(6.5, abs=0.02)
    assert res.time_above_mtc_h == pytest.approx(3.5, abs=0.02)


def test_assess_never_effective_and_trough_and_ok():
    t = np.linspace(0, 10, 1001)
    assert assess_tdm(t, np.full_like(t, 1.0), WINDOW).status == "never_effective"
    saw = np.where(t < 5, 5.0, 1.0)
    assert assess_tdm(t, saw, WINDOW).status == "trough_below_mec"
    # trough only evaluated after assess_from_h
    assert assess_tdm(t, saw, WINDOW, assess_from_h=0.0).status == "trough_below_mec"
    good = np.full_like(t, 4.0)
    res = assess_tdm(t, good, WINDOW)
    assert res.status == "within_window"
    assert res.time_in_window_h == pytest.approx(10.0, abs=0.02)


def test_assess_ignores_initial_ramp_when_assess_from_given():
    t = np.linspace(0, 20, 2001)
    c = np.where(t < 2, 0.0, 4.0)
    assert assess_tdm(t, c, WINDOW).status == "trough_below_mec"
    assert assess_tdm(t, c, WINDOW, assess_from_h=5.0).status == "within_window"


def test_recovery_time_after_skipped_dose():
    tau = 8.0
    events = regular_schedule(400.0, tau, 400)
    t = time_grid(400)
    ref = concentration_profile(t, events, MODEL)
    alt = concentration_profile(t, apply_missed_dose(events, 20), MODEL)
    rec = recovery_time(t, ref, alt, after_h=events[20].time_h, rel_tol=0.05)
    assert rec is not None and rec > 0
    # identical profiles recover instantly
    assert recovery_time(t, ref, ref, after_h=0.0) == 0.0
    # never recovers within a too-short horizon
    short = time_grid(events[20].time_h + 2)
    ref_s = concentration_profile(short, events, MODEL)
    alt_s = concentration_profile(short, apply_missed_dose(events, 20), MODEL)
    assert recovery_time(short, ref_s, alt_s, after_h=events[20].time_h) is None


def test_window_position_is_zero_at_mec_and_one_at_mtc():
    import numpy as np
    from core.metrics import TherapeuticWindow
    w = TherapeuticWindow(0.6, 1.5)
    assert w.position(0.6) == pytest.approx(0.0)
    assert w.position(1.5) == pytest.approx(1.0)
    assert w.position(0.3) < 0.0  # sub-therapeutic is negative
    assert np.allclose(w.position(np.array([0.6, 1.05, 1.5])), [0.0, 0.5, 1.0])
