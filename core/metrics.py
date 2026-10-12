"""Derived clinical metrics and TDM (therapeutic drug monitoring) validation."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from core.pk_engine import PKModel, auc_trapezoid


@dataclass(frozen=True)
class ClinicalMetrics:
    cmax: float
    tmax_h: float
    t_half_h: float
    clearance_l_per_h: float
    auc_0_24: float


def compute_metrics(
    t: np.ndarray, c: np.ndarray, model: PKModel, auc_window_h: tuple[float, float] = (0.0, 24.0)
) -> ClinicalMetrics:
    """Cmax/Tmax over the whole profile, AUC over ``auc_window_h`` (default 0-24 h)."""
    t = np.asarray(t, dtype=float)
    c = np.asarray(c, dtype=float)
    i = int(np.argmax(c))
    hi = min(auc_window_h[1], float(t[-1]))
    return ClinicalMetrics(
        cmax=float(c[i]),
        tmax_h=float(t[i]),
        t_half_h=model.terminal_half_life_h,
        clearance_l_per_h=model.clearance_l_per_h,
        auc_0_24=auc_trapezoid(t, c, auc_window_h[0], hi),
    )


def peak_trough_in_interval(
    t: np.ndarray, c: np.ndarray, start_h: float, end_h: float
) -> tuple[float, float]:
    """(peak, trough) of the profile restricted to [start_h, end_h]."""
    grid = np.linspace(start_h, end_h, 2001)
    seg = np.interp(grid, t, c)
    return float(seg.max()), float(seg.min())


def fluctuation_percent(cmax: float, cmin: float, cavg: float) -> float:
    """Peak-trough fluctuation, 100 * (Cmax - Cmin) / Cavg."""
    return 100.0 * (cmax - cmin) / cavg


def terminal_half_life(t: np.ndarray, c: np.ndarray, n_points: int = 20) -> float:
    """Half-life from log-linear regression over the last ``n_points`` positive samples.

    Intended for single-dose profiles sampled in the elimination phase.
    """
    t = np.asarray(t, dtype=float)
    c = np.asarray(c, dtype=float)
    mask = c > 0
    tt, cc = t[mask][-n_points:], c[mask][-n_points:]
    if len(tt) < 3:
        raise ValueError("not enough positive samples for regression")
    slope = stats.linregress(tt, np.log(cc)).slope
    if slope >= 0:
        raise ValueError("profile is not in a declining phase")
    return float(np.log(2.0) / -slope)


# ----------------------------------------------------------------------- TDM
@dataclass(frozen=True)
class TherapeuticWindow:
    mec: float
    mtc: float

    def __post_init__(self) -> None:
        if not 0 < self.mec < self.mtc:
            raise ValueError("require 0 < MEC < MTC")

    def classify(self, c: np.ndarray | float) -> np.ndarray:
        """0 = sub-therapeutic, 1 = therapeutic, 2 = toxic."""
        c = np.asarray(c, dtype=float)
        return np.where(c > self.mtc, 2, np.where(c >= self.mec, 1, 0))


ZONE_LABELS = {0: "subtherapeutic", 1: "therapeutic", 2: "toxic"}


@dataclass(frozen=True)
class TDMAssessment:
    status: str  # "toxic" | "never_effective" | "trough_below_mec" | "within_window"
    cmax: float
    t_cmax_h: float
    first_toxic_time_h: float | None
    time_above_mtc_h: float
    time_in_window_h: float
    time_below_mec_h: float
    assessed_from_h: float
    trough: float


def _time_where(t: np.ndarray, mask: np.ndarray) -> float:
    """Total time spent where ``mask`` holds (left-point rule on the grid)."""
    return float(np.sum(np.diff(t)[mask[:-1]]))


def first_crossing_time(t: np.ndarray, c: np.ndarray, threshold: float) -> float | None:
    """First time the profile exceeds ``threshold`` (linear interpolation)."""
    idx = np.flatnonzero(c > threshold)
    if idx.size == 0:
        return None
    j = int(idx[0])
    if j == 0:
        return float(t[0])
    c0, c1 = c[j - 1], c[j]
    return float(t[j - 1] + (threshold - c0) / (c1 - c0) * (t[j] - t[j - 1]))


def assess_tdm(
    t: np.ndarray, c: np.ndarray, window: TherapeuticWindow, assess_from_h: float = 0.0
) -> TDMAssessment:
    """Evaluate a profile against the therapeutic window.

    Toxicity is checked over the entire profile; the trough check only from
    ``assess_from_h`` (typically the start of steady state).
    """
    t = np.asarray(t, dtype=float)
    c = np.asarray(c, dtype=float)
    i = int(np.argmax(c))
    sel = t >= assess_from_h
    ts, cs = t[sel], c[sel]
    trough = float(cs.min())
    first_tox = first_crossing_time(t, c, window.mtc)

    if first_tox is not None:
        status = "toxic"
    elif c[i] < window.mec:
        status = "never_effective"
    elif trough < window.mec:
        status = "trough_below_mec"
    else:
        status = "within_window"

    return TDMAssessment(
        status=status,
        cmax=float(c[i]),
        t_cmax_h=float(t[i]),
        first_toxic_time_h=first_tox,
        time_above_mtc_h=_time_where(ts, cs > window.mtc),
        time_in_window_h=_time_where(ts, (cs >= window.mec) & (cs <= window.mtc)),
        time_below_mec_h=_time_where(ts, cs < window.mec),
        assessed_from_h=float(ts[0]),
        trough=trough,
    )


def recovery_time(
    t: np.ndarray,
    c_ref: np.ndarray,
    c_alt: np.ndarray,
    after_h: float,
    rel_tol: float = 0.05,
) -> float | None:
    """Hours after ``after_h`` until the altered profile is back within ``rel_tol``
    of the reference for good. ``None`` if it never recovers within the horizon."""
    t = np.asarray(t, dtype=float)
    ref = np.asarray(c_ref, dtype=float)
    alt = np.asarray(c_alt, dtype=float)
    mask = t >= after_h
    dev = np.abs(alt - ref) > rel_tol * np.maximum(ref, 1e-12)
    dev &= mask
    idx = np.flatnonzero(dev)
    if idx.size == 0:
        return 0.0
    last = int(idx[-1])
    if last >= len(t) - 1:
        return None
    return float(t[last + 1] - after_h)
