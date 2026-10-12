"""Pure pharmacokinetic mathematics (no UI, no I/O).

Conventions
-----------
* Time in hours, volumes in litres, rate constants in 1/h.
* ``amount`` is the administered dose expressed in the *amount unit of the
  reported concentration* (e.g. mg for ug/mL, mmol for mmol/L). Converting a
  dose in mg to that unit is the caller's job (see ``Drug.amount_from_dose_mg``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from enum import Enum
from functools import lru_cache
from typing import Iterable, Sequence

import numpy as np
from scipy.integrate import solve_ivp

LN2 = math.log(2.0)
# Relative closeness of ka and ke below which the flip-flop limit is used.
_KA_KE_REL_TOL = 1e-9
TSS_HALF_LIVES = 4.5  # midpoint of the textbook 4-5 half-lives rule


class Route(str, Enum):
    ORAL = "oral"
    IV_BOLUS = "iv_bolus"
    INSUFFLATED = "insufflated"
    SMOKED = "smoked"


@dataclass(frozen=True)
class PKModel:
    """One-compartment model parameters.

    With ``km`` set, elimination is partly saturable (Michaelis-Menten):
    CL(C) = (1 - f) * CL + f * CL * km / (km + C), f = ``saturable_fraction``.
    ``t_half_h`` and ``clearance_l_per_h`` then describe the low-concentration
    limit, and profiles are integrated numerically instead of superposed.
    """

    t_half_h: float
    vd_l: float
    bioavailability: float = 1.0
    ka_per_h: float = 1.0
    route: Route = Route.ORAL
    km: float | None = None
    saturable_fraction: float = 0.0

    def __post_init__(self) -> None:
        r = self.route
        if isinstance(r, str) and not isinstance(r, Route):
            r_str = r.strip().lower()
            if r_str in ("inhaled", "smoked"):
                r = Route.SMOKED
            elif r_str in ("iv", "iv_bolus", "intravenous"):
                r = Route.IV_BOLUS
            elif r_str in ("insufflated", "nasal", "snorted"):
                r = Route.INSUFFLATED
            elif r_str in ("oral", "po", "peroral"):
                r = Route.ORAL
            else:
                r = Route(r_str)
            object.__setattr__(self, "route", r)

        if self.t_half_h <= 0:
            raise ValueError("t_half_h must be > 0")
        if self.vd_l <= 0:
            raise ValueError("vd_l must be > 0")
        if not 0 < self.bioavailability <= 1:
            raise ValueError("bioavailability must be in (0, 1]")
        if self.route != Route.IV_BOLUS and self.ka_per_h <= 0:
            raise ValueError(f"ka_per_h must be > 0 for {self.route.value} administration")
        if self.km is not None:
            if self.km <= 0:
                raise ValueError("km must be > 0")
            if not 0 < self.saturable_fraction <= 1:
                raise ValueError("saturable_fraction must be in (0, 1] when km is set")

    @property
    def ke(self) -> float:
        return elimination_rate_constant(self.t_half_h)

    @property
    def f_eff(self) -> float:
        """Fraction reaching systemic circulation (1 for IV bolus)."""
        return 1.0 if self.route is Route.IV_BOLUS else self.bioavailability

    @property
    def clearance_l_per_h(self) -> float:
        return clearance(self.ke, self.vd_l)

    @property
    def is_saturable(self) -> bool:
        return self.km is not None

    def with_half_life(self, t_half_h: float) -> "PKModel":
        return replace(self, t_half_h=t_half_h)


@dataclass(frozen=True)
class DoseEvent:
    time_h: float
    amount: float


# ---------------------------------------------------------------- basic relations
def elimination_rate_constant(t_half_h: float) -> float:
    """ke = ln2 / t1/2."""
    return LN2 / t_half_h


def half_life_from_ke(ke: float) -> float:
    return LN2 / ke


def clearance(ke: float, vd_l: float) -> float:
    """CL = ke * Vd (L/h)."""
    return ke * vd_l


def _same_rate(ka: float, ke: float) -> bool:
    return abs(ka - ke) <= _KA_KE_REL_TOL * max(ka, ke)


# ------------------------------------------------------------------ single dose
def single_dose_concentration(t: np.ndarray | float, amount: float, model: PKModel) -> np.ndarray:
    """Concentration after one dose given at t = 0 (zero for t < 0).

    Oral: Bateman function; IV bolus: C0 * exp(-ke t) with C0 = D / Vd.
    """
    if model.is_saturable:
        return _saturable_profile(t, (DoseEvent(0.0, amount),), model)
    t = np.asarray(t, dtype=float)
    tp = np.clip(t, 0.0, None)
    ke = model.ke
    if model.route is Route.IV_BOLUS:
        c = amount / model.vd_l * np.exp(-ke * tp)
    else:
        ka = model.ka_per_h
        scale = model.f_eff * amount / model.vd_l
        if _same_rate(ka, ke):  # flip-flop limit ka -> ke
            c = scale * ka * tp * np.exp(-ke * tp)
        else:
            # exp(-ke t) - exp(-ka t) written with expm1 to avoid cancellation.
            c = scale * ka * np.exp(-ke * tp) * (-np.expm1(-(ka - ke) * tp)) / (ka - ke)
    return np.where(t >= 0, c, 0.0)


def tmax_single_dose(model: PKModel, amount: float | None = None) -> float:
    """Time of the peak after a single dose.

    Saturable models peak later at higher doses, so pass ``amount`` for them;
    without it the low-concentration (linear) limit is returned.
    """
    if model.route is Route.IV_BOLUS:
        return 0.0
    if model.is_saturable and amount is not None:
        return _saturable_single_peak(amount, model)[0]
    ka, ke = model.ka_per_h, model.ke
    if _same_rate(ka, ke):
        return 1.0 / ke
    return math.log(ka / ke) / (ka - ke)


def cmax_single_dose(amount: float, model: PKModel) -> float:
    if model.is_saturable:
        return _saturable_single_peak(amount, model)[1]
    return float(single_dose_concentration(tmax_single_dose(model), amount, model))


# -------------------------------------------------------------------- schedules
def regular_schedule(amount: float, interval_h: float, duration_h: float) -> list[DoseEvent]:
    """Doses at 0, tau, 2 tau, ... strictly before ``duration_h``."""
    if interval_h <= 0 or duration_h <= 0:
        raise ValueError("interval_h and duration_h must be > 0")
    n = max(1, math.ceil(duration_h / interval_h - 1e-9))
    return [DoseEvent(k * interval_h, amount) for k in range(n)]


def apply_missed_dose(
    events: Sequence[DoseEvent], index: int, delay_h: float | None = None
) -> list[DoseEvent]:
    """Skip dose ``index`` (``delay_h is None``) or take it ``delay_h`` hours late."""
    if not 0 <= index < len(events):
        raise IndexError("dose index out of range")
    if delay_h is not None and delay_h < 0:
        raise ValueError("delay_h must be >= 0")
    out = list(events)
    if delay_h is None:
        del out[index]
    else:
        out[index] = DoseEvent(events[index].time_h + delay_h, events[index].amount)
    return sorted(out, key=lambda e: e.time_h)


# ------------------------------------------------------------- multiple dosing
def concentration_profile(
    t: np.ndarray, events: Iterable[DoseEvent], model: PKModel
) -> np.ndarray:
    """Superposition of single-dose profiles (valid: linear, time-invariant PK).

    Saturable models are integrated numerically over the whole schedule instead.
    """
    if model.is_saturable:
        return _saturable_profile(t, tuple(events), model)
    t = np.asarray(t, dtype=float)
    total = np.zeros_like(t)
    for ev in events:
        total += single_dose_concentration(t - ev.time_h, ev.amount, model)
    return total


def time_grid(duration_h: float, max_step_h: float = 0.05, max_points: int = 20000) -> np.ndarray:
    n = int(min(max_points, max(500, math.ceil(duration_h / max_step_h))))
    return np.linspace(0.0, duration_h, n + 1)


# ------------------------------------------------------------------------- AUC
def auc_trapezoid(
    t: np.ndarray, c: np.ndarray, t_start: float | None = None, t_end: float | None = None
) -> float:
    """Linear trapezoidal AUC over [t_start, t_end] (bounds are interpolated)."""
    t = np.asarray(t, dtype=float)
    c = np.asarray(c, dtype=float)
    lo = t[0] if t_start is None else t_start
    hi = t[-1] if t_end is None else t_end
    if hi < lo:
        raise ValueError("t_end must be >= t_start")
    inside = (t > lo) & (t < hi)
    tt = np.concatenate(([lo], t[inside], [hi]))
    cc = np.concatenate(([np.interp(lo, t, c)], c[inside], [np.interp(hi, t, c)]))
    return float(np.sum((cc[1:] + cc[:-1]) * 0.5 * np.diff(tt)))


# ----------------------------------------------------------------- steady state
def accumulation_index(ke: float, tau_h: float) -> float:
    """R = 1 / (1 - exp(-ke tau))."""
    return 1.0 / (1.0 - math.exp(-ke * tau_h))


def time_to_steady_state(t_half_h: float, n_half_lives: float = TSS_HALF_LIVES) -> float:
    return n_half_lives * t_half_h


def time_to_steady_state_range(t_half_h: float) -> tuple[float, float]:
    """The classical 4-5 half-lives window."""
    return 4.0 * t_half_h, 5.0 * t_half_h


def _nudged_ka(ka: float, ke: float) -> float:
    return ke * (1 + 1e-6) if abs(ka - ke) <= 1e-6 * max(ka, ke) else ka


def steady_state_concentration(
    t_since_dose_h: np.ndarray | float, amount: float, tau_h: float, model: PKModel
) -> np.ndarray:
    """Concentration at time t after a dose once steady state is reached."""
    t = np.asarray(t_since_dose_h, dtype=float)
    if model.is_saturable:
        ss = _saturable_steady_state(amount, tau_h, model)
        phase = np.mod(t, tau_h)
        phase = np.where((phase == 0) & (t > 0), tau_h, phase)  # t = k*tau is the trough
        return np.interp(phase, ss.t, ss.c)
    ke = model.ke
    if model.route is Route.IV_BOLUS:
        return amount / model.vd_l * np.exp(-ke * t) / (1.0 - math.exp(-ke * tau_h))
    ka = _nudged_ka(model.ka_per_h, ke)
    scale = model.f_eff * amount * ka / (model.vd_l * (ka - ke))
    return scale * (
        np.exp(-ke * t) / (1.0 - math.exp(-ke * tau_h))
        - np.exp(-ka * t) / (1.0 - math.exp(-ka * tau_h))
    )


@dataclass(frozen=True)
class SteadyState:
    cmax_ss: float
    tmax_ss_h: float
    cmin_ss: float
    cavg_ss: float
    accumulation_index: float
    fluctuation_pct: float
    time_to_ss_h: float
    time_to_ss_range_h: tuple[float, float]


def steady_state_summary(amount: float, tau_h: float, model: PKModel) -> SteadyState:
    """Analytic steady-state descriptors for repeated dosing every tau hours."""
    if tau_h <= 0:
        raise ValueError("tau_h must be > 0")
    if model.is_saturable:
        return _saturable_steady_state(amount, tau_h, model).summary
    ke = model.ke
    if model.route is Route.IV_BOLUS:
        tmax = 0.0
        cmax = float(steady_state_concentration(0.0, amount, tau_h, model))
    else:
        ka = _nudged_ka(model.ka_per_h, ke)
        ratio = ka * (1 - math.exp(-ke * tau_h)) / (ke * (1 - math.exp(-ka * tau_h)))
        tmax = math.log(ratio) / (ka - ke)
        cmax = float(steady_state_concentration(tmax, amount, tau_h, model))
    cmin = float(steady_state_concentration(tau_h, amount, tau_h, model))
    cavg = model.f_eff * amount / (model.clearance_l_per_h * tau_h)
    return SteadyState(
        cmax_ss=cmax,
        tmax_ss_h=tmax,
        cmin_ss=cmin,
        cavg_ss=cavg,
        accumulation_index=accumulation_index(ke, tau_h),
        fluctuation_pct=100.0 * (cmax - cmin) / cavg,
        time_to_ss_h=time_to_steady_state(model.t_half_h),
        time_to_ss_range_h=time_to_steady_state_range(model.t_half_h),
    )


# ------------------------------------------------- saturable (Michaelis-Menten)
def _saturable_rhs(model: PKModel):
    v = model.vd_l
    cl = model.clearance_l_per_h
    f = model.saturable_fraction
    cl_lin, vmax, km = (1.0 - f) * cl, f * cl * model.km, model.km
    ka = 0.0 if model.route is Route.IV_BOLUS else model.ka_per_h

    def rhs(_t, y):
        c = max(y[1], 0.0) / v
        return [-ka * y[0], ka * y[0] - cl_lin * c - vmax * c / (km + c)]

    return rhs


def _dose_into(state: np.ndarray, amount: float, model: PKModel) -> np.ndarray:
    """Add a dose to [depot, central] amounts (absorbed fraction for extravascular routes)."""
    out = state.copy()
    if model.route is Route.IV_BOLUS:
        out[1] += amount
    else:
        out[0] += model.f_eff * amount
    return out


def _integrate(rhs, y0: np.ndarray, t0: float, t1: float, t_eval: np.ndarray, atol: float):
    """Integrate [t0, t1]; return (state at t1, states at sorted ``t_eval`` within the span)."""
    if t1 <= t0:
        return y0, np.tile(y0[:, None], (1, len(t_eval)))
    ends_at_t1 = len(t_eval) > 0 and t_eval[-1] == t1
    pts = t_eval if ends_at_t1 else np.append(t_eval, t1)
    sol = solve_ivp(rhs, (t0, t1), y0, method="LSODA", t_eval=pts, rtol=1e-7, atol=atol)
    if not sol.success:
        raise RuntimeError(f"saturable PK integration failed: {sol.message}")
    return sol.y[:, -1], (sol.y if ends_at_t1 else sol.y[:, :-1])


def _saturable_profile(t, events: tuple[DoseEvent, ...], model: PKModel) -> np.ndarray:
    """Concentration over ``t`` for a dose schedule, integrating the nonlinear ODE."""
    t_arr = np.asarray(t, dtype=float)
    flat = t_arr.reshape(-1)
    out = np.zeros_like(flat)
    evs = sorted((e for e in events if e.amount > 0), key=lambda e: e.time_h)
    if not evs or flat.size == 0:
        return out.reshape(t_arr.shape)
    order = np.argsort(flat)
    ts = flat[order]
    vals = np.zeros_like(ts)
    rhs = _saturable_rhs(model)
    atol = 1e-10 * max(e.amount for e in evs)
    times = sorted({e.time_h for e in evs})
    y = np.zeros(2)
    for i, t0 in enumerate(times):
        for e in evs:
            if e.time_h == t0:
                y = _dose_into(y, e.amount, model)
        t1 = times[i + 1] if i + 1 < len(times) else max(ts[-1], t0)
        last = i + 1 == len(times)
        mask = (ts >= t0) & ((ts <= t1) if last else (ts < t1))
        y, ys = _integrate(rhs, y, t0, t1, ts[mask], atol)
        vals[mask] = ys[1] / model.vd_l
    out[order] = np.clip(vals, 0.0, None)
    return out.reshape(t_arr.shape)


@lru_cache(maxsize=256)
def _saturable_single_peak(amount: float, model: PKModel) -> tuple[float, float]:
    """(tmax, cmax) of a single dose, located numerically."""
    if model.route is Route.IV_BOLUS:
        return 0.0, amount / model.vd_l
    lin_tmax = tmax_single_dose(replace(model, km=None, saturable_fraction=0.0))
    coarse = np.linspace(0.0, max(4.0 * lin_tmax, 1.0), 801)
    c = _saturable_profile(coarse, (DoseEvent(0.0, amount),), model)
    i = int(np.argmax(c))
    lo, hi = coarse[max(i - 1, 0)], coarse[min(i + 1, len(coarse) - 1)]
    fine = np.linspace(lo, hi, 401)
    cf = _saturable_profile(fine, (DoseEvent(0.0, amount),), model)
    j = int(np.argmax(cf))
    return float(fine[j]), float(cf[j])


@dataclass(frozen=True)
class _SaturableSS:
    t: np.ndarray
    c: np.ndarray
    summary: SteadyState


@lru_cache(maxsize=128)
def _saturable_steady_state(amount: float, tau_h: float, model: PKModel) -> _SaturableSS:
    """Repeat the dosing interval until the trough converges (max 2000 doses)."""
    rhs = _saturable_rhs(model)
    atol = 1e-10 * amount
    grid = np.linspace(0.0, tau_h, 401)
    y = np.zeros(2)
    troughs: list[float] = []
    first_auc = None
    for _ in range(2000):
        y = _dose_into(y, amount, model)
        y, ys = _integrate(rhs, y, 0.0, tau_h, grid, atol)
        c = np.clip(ys[1] / model.vd_l, 0.0, None)
        if first_auc is None:
            first_auc = auc_trapezoid(grid, c)
        troughs.append(float(c[-1]))
        if len(troughs) > 3 and abs(troughs[-1] - troughs[-2]) <= 1e-6 * troughs[-1]:
            break
    i = int(np.argmax(c))
    cmin = float(c[-1])
    cavg = auc_trapezoid(grid, c) / tau_h
    tr = np.asarray(troughs)

    def reached(frac: float) -> float:
        k = int(np.argmax(tr >= frac * cmin))
        return (k + 1) * tau_h

    summary = SteadyState(
        cmax_ss=float(c[i]),
        tmax_ss_h=float(grid[i]),
        cmin_ss=cmin,
        cavg_ss=cavg,
        accumulation_index=cavg * tau_h / first_auc,
        fluctuation_pct=100.0 * (float(c[i]) - cmin) / cavg,
        time_to_ss_h=reached(1.0 - 2.0 ** -TSS_HALF_LIVES),
        time_to_ss_range_h=(reached(1.0 - 2.0 ** -4), reached(1.0 - 2.0 ** -5)),
    )
    return _SaturableSS(grid, c, summary)
