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
    """Compartmental model parameters.

    ``t_half_h`` and ``vd_l`` define clearance, CL = ln2 * vd_l / t_half_h, and the
    steady-state volume; together they fix AUC and the average steady-state level.

    Two compartments (``v1_fraction`` < 1): the central volume is v1_fraction * vd_l,
    the rest is peripheral, exchanging with inter-compartmental clearance ``q_l_per_h``.
    The terminal half-life is then longer than t_half_h (see ``terminal_half_life_h``).

    Nonlinear elimination (profiles integrated numerically instead of superposed):
    * ``km``: saturable (Michaelis-Menten) share f = ``saturable_fraction`` of CL,
      CL(C) = (1 - f) * CL + f * CL * km / (km + C).
    * ``inactivation_fraction`` fm > 0: mechanism-based enzyme inactivation, the drug
      destroys the enzyme carrying fm of its clearance. Active enzyme E (1 = baseline)
      follows dE/dt = kdeg * (1 - E) - kinact * C / (KI + C) * E and CL = CL * (1 - fm + fm * E).
    CL then describes the low-concentration, fresh-enzyme limit.

    ``tlag_h`` delays absorption of extravascular doses (e.g. gastric emptying).
    """

    t_half_h: float
    vd_l: float
    bioavailability: float = 1.0
    ka_per_h: float = 1.0
    route: Route = Route.ORAL
    km: float | None = None
    saturable_fraction: float = 0.0
    v1_fraction: float = 1.0
    q_l_per_h: float = 0.0
    inactivation_fraction: float = 0.0
    kinact_per_h: float = 0.0
    ki: float = 1.0
    kdeg_per_h: float = 0.0
    tlag_h: float = 0.0  # absorption lag (extravascular routes only)

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
        if not 0 < self.v1_fraction <= 1:
            raise ValueError("v1_fraction must be in (0, 1]")
        if self.v1_fraction < 1 and self.q_l_per_h <= 0:
            raise ValueError("q_l_per_h must be > 0 for a two-compartment model")
        if self.tlag_h < 0:
            raise ValueError("tlag_h must be >= 0")
        if self.inactivation_fraction:
            if not 0 < self.inactivation_fraction <= 1:
                raise ValueError("inactivation_fraction must be in (0, 1]")
            if min(self.kinact_per_h, self.ki, self.kdeg_per_h) <= 0:
                raise ValueError("kinact_per_h, ki and kdeg_per_h must be > 0 for enzyme inactivation")
            if self.km is not None:
                raise ValueError("use either saturable (km) or enzyme-inactivation elimination, not both")

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

    @property
    def is_nonlinear(self) -> bool:
        """Needs numerical integration (no superposition)."""
        return self.km is not None or self.inactivation_fraction > 0

    @property
    def is_two_compartment(self) -> bool:
        return self.v1_fraction < 1

    @property
    def v1_l(self) -> float:
        return self.vd_l * self.v1_fraction

    @property
    def micro_constants(self) -> tuple[float, float, float]:
        """(k10, k12, k21) in 1/h; k12 = k21 = 0 for one compartment."""
        k10 = self.clearance_l_per_h / self.v1_l
        if not self.is_two_compartment:
            return k10, 0.0, 0.0
        return k10, self.q_l_per_h / self.v1_l, self.q_l_per_h / (self.vd_l - self.v1_l)

    @property
    def disposition_rates(self) -> tuple[float, ...]:
        """Exponents of the IV-bolus curve: (ke,) or (alpha, beta), alpha > beta."""
        k10, k12, k21 = self.micro_constants
        if not self.is_two_compartment:
            return (k10,)
        s = k10 + k12 + k21
        root = math.sqrt(s * s - 4.0 * k10 * k21)
        return 0.5 * (s + root), 0.5 * (s - root)

    @property
    def lag_h(self) -> float:
        """Absorption lag actually applied (none for IV bolus)."""
        return 0.0 if self.route is Route.IV_BOLUS else self.tlag_h

    @property
    def terminal_half_life_h(self) -> float:
        """Half-life of the elimination phase. With enzyme inactivation this is the
        inactivated-enzyme value, which is what follows any usual dose."""
        rate = min(self.disposition_rates)
        if self.inactivation_fraction:
            rate *= 1.0 - self.inactivation_fraction
        return LN2 / rate

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
    if model.is_nonlinear:
        return _nonlinear_profile(t, (DoseEvent(0.0, amount),), model)
    t = np.asarray(t, dtype=float) - model.lag_h
    tp = np.clip(t, 0.0, None)
    if model.is_two_compartment:
        coefs, rates = _exp_terms(model, amount)
        c = np.sum(coefs[:, None] * np.exp(-np.outer(rates, tp.reshape(-1))), axis=0).reshape(tp.shape)
        return np.where(t >= 0, np.clip(c, 0.0, None), 0.0)
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
    if model.is_nonlinear and amount is not None:
        return _nonlinear_single_peak(amount, model)[0]
    if model.is_two_compartment:
        return _linear_single_peak(_linearized(model))[0]
    ka, ke = model.ka_per_h, model.ke
    if _same_rate(ka, ke):
        return model.lag_h + 1.0 / ke
    return model.lag_h + math.log(ka / ke) / (ka - ke)


def cmax_single_dose(amount: float, model: PKModel) -> float:
    if model.is_nonlinear:
        return _nonlinear_single_peak(amount, model)[1]
    if model.is_two_compartment:
        return amount * _linear_single_peak(model)[1]
    return float(single_dose_concentration(tmax_single_dose(model), amount, model))


def _exp_terms(model: PKModel, amount: float) -> tuple[np.ndarray, np.ndarray]:
    """Single-dose curve of a linear model as sum(coef_i * exp(-rate_i * t))."""
    k10, k12, k21 = model.micro_constants
    if not model.is_two_compartment:
        rates = (k10,)
        if model.route is Route.IV_BOLUS:
            return np.array([amount / model.vd_l]), np.array(rates)
    else:
        rates = model.disposition_rates
    if model.route is Route.IV_BOLUS:
        a, b = rates
        return amount / model.v1_l * np.array([(a - k21) / (a - b), (k21 - b) / (a - b)]), np.array(rates)
    ka = model.ka_per_h
    for r in rates:  # keep ka distinct from the disposition rates
        if abs(ka - r) <= 1e-6 * max(ka, r):
            ka = r * (1 + 1e-5)
    all_rates = (*rates, ka)
    coefs = []
    for i, r in enumerate(all_rates):
        num = k21 - r if model.is_two_compartment else 1.0
        den = np.prod([o - r for j, o in enumerate(all_rates) if j != i])
        coefs.append(num / den)
    scale = model.f_eff * amount * ka / model.v1_l
    return scale * np.array(coefs), np.array(all_rates)


def _linearized(model: PKModel) -> PKModel:
    return replace(model, km=None, saturable_fraction=0.0, inactivation_fraction=0.0)


@lru_cache(maxsize=256)
def _linear_single_peak(model: PKModel) -> tuple[float, float]:
    """(tmax, cmax per unit dose) of a linear model, located numerically."""
    if model.route is Route.IV_BOLUS:
        return 0.0, float(single_dose_concentration(0.0, 1.0, model))
    horizon = model.lag_h + 10.0 / model.ka_per_h + 2.0 * model.terminal_half_life_h
    coarse = np.linspace(0.0, horizon, 4001)
    c = single_dose_concentration(coarse, 1.0, model)
    i = int(np.argmax(c))
    fine = np.linspace(coarse[max(i - 1, 0)], coarse[min(i + 1, len(coarse) - 1)], 801)
    cf = single_dose_concentration(fine, 1.0, model)
    j = int(np.argmax(cf))
    return float(fine[j]), float(cf[j])


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
    if model.is_nonlinear:
        return _nonlinear_profile(t, tuple(events), model)
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
    if model.lag_h:
        # At steady state the lagged profile is the unlagged one shifted by the lag.
        unlagged = replace(model, tlag_h=0.0)
        return steady_state_concentration(np.mod(t - model.lag_h, tau_h), amount, tau_h, unlagged)
    if model.is_nonlinear:
        ss = _nonlinear_steady_state(amount, tau_h, model)
        phase = np.mod(t, tau_h)
        phase = np.where((phase == 0) & (t > 0), tau_h, phase)  # t = k*tau is the trough
        return np.interp(phase, ss.t, ss.c)
    if model.is_two_compartment:
        coefs, rates = _exp_terms(model, amount)
        geo = coefs / (1.0 - np.exp(-rates * tau_h))
        return np.sum(geo[:, None] * np.exp(-np.outer(rates, t.reshape(-1))), axis=0).reshape(t.shape)
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
    if model.lag_h:
        ss = steady_state_summary(amount, tau_h, replace(model, tlag_h=0.0))
        return replace(ss, tmax_ss_h=(ss.tmax_ss_h + model.lag_h) % tau_h)
    if model.is_nonlinear:
        return _nonlinear_steady_state(amount, tau_h, model).summary
    if model.is_two_compartment:
        grid = np.linspace(0.0, tau_h, 2001)
        c = steady_state_concentration(grid, amount, tau_h, model)
        i = int(np.argmax(c))
        cmin, cavg = float(c.min()), model.f_eff * amount / (model.clearance_l_per_h * tau_h)
        t_term = model.terminal_half_life_h
        return SteadyState(
            cmax_ss=float(c[i]), tmax_ss_h=float(grid[i]), cmin_ss=cmin, cavg_ss=cavg,
            accumulation_index=accumulation_index(LN2 / t_term, tau_h),
            fluctuation_pct=100.0 * (float(c[i]) - cmin) / cavg,
            time_to_ss_h=time_to_steady_state(t_term),
            time_to_ss_range_h=time_to_steady_state_range(t_term),
        )
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


# ------------------------------------------------------ nonlinear (numerical)
# State: [depot amount, central amount, peripheral amount, active enzyme fraction].
def _nonlinear_rhs(model: PKModel):
    v1 = model.v1_l
    cl = model.clearance_l_per_h
    _, k12, k21 = model.micro_constants
    ka = 0.0 if model.route is Route.IV_BOLUS else model.ka_per_h
    f_sat, km = model.saturable_fraction, model.km
    fm, kinact, ki, kdeg = model.inactivation_fraction, model.kinact_per_h, model.ki, model.kdeg_per_h

    def rhs(_t, y):
        c = max(y[1], 0.0) / v1
        if km is not None:
            elim = ((1.0 - f_sat) + f_sat * km / (km + c)) * cl * c
        else:
            elim = ((1.0 - fm) + fm * y[3]) * cl * c
        d_enzyme = kdeg * (1.0 - y[3]) - kinact * c / (ki + c) * y[3] if fm else 0.0
        return [-ka * y[0], ka * y[0] - elim - k12 * y[1] + k21 * y[2], k12 * y[1] - k21 * y[2], d_enzyme]

    return rhs


def _initial_state() -> np.ndarray:
    return np.array([0.0, 0.0, 0.0, 1.0])


def _dose_into(state: np.ndarray, amount: float, model: PKModel) -> np.ndarray:
    """Add a dose (absorbed fraction for extravascular routes) to the state vector."""
    out = state.copy()
    if model.route is Route.IV_BOLUS:
        out[1] += amount
    else:
        out[0] += model.f_eff * amount
    return out


def _atol(amount_scale: float) -> np.ndarray:
    return np.array([1e-10 * amount_scale] * 3 + [1e-10])


def _integrate(rhs, y0: np.ndarray, t0: float, t1: float, t_eval: np.ndarray, atol):
    """Integrate [t0, t1]; return (state at t1, states at sorted ``t_eval`` within the span)."""
    if t1 <= t0:
        return y0, np.tile(y0[:, None], (1, len(t_eval)))
    ends_at_t1 = len(t_eval) > 0 and t_eval[-1] == t1
    pts = t_eval if ends_at_t1 else np.append(t_eval, t1)
    sol = solve_ivp(rhs, (t0, t1), y0, method="LSODA", t_eval=pts, rtol=1e-7, atol=atol)
    if not sol.success:
        raise RuntimeError(f"nonlinear PK integration failed: {sol.message}")
    return sol.y[:, -1], (sol.y if ends_at_t1 else sol.y[:, :-1])


def _nonlinear_states(t, events: tuple[DoseEvent, ...], model: PKModel) -> tuple[np.ndarray, np.ndarray]:
    """(sorted unique-safe times, state matrix 4 x n) for a dose schedule; zeros before the first dose."""
    ts = np.sort(np.asarray(t, dtype=float).reshape(-1))
    states = np.tile(_initial_state()[:, None], (1, ts.size))
    lag = model.lag_h
    evs = sorted((DoseEvent(e.time_h + lag, e.amount) for e in events if e.amount > 0), key=lambda e: e.time_h)
    if not evs or ts.size == 0:
        return ts, states
    rhs = _nonlinear_rhs(model)
    atol = _atol(max(e.amount for e in evs))
    times = sorted({e.time_h for e in evs})
    y = _initial_state()
    for i, t0 in enumerate(times):
        for e in evs:
            if e.time_h == t0:
                y = _dose_into(y, e.amount, model)
        t1 = times[i + 1] if i + 1 < len(times) else max(ts[-1], t0)
        last = i + 1 == len(times)
        mask = (ts >= t0) & ((ts <= t1) if last else (ts < t1))
        y, ys = _integrate(rhs, y, t0, t1, ts[mask], atol)
        states[:, mask] = ys
    return ts, states


def _nonlinear_profile(t, events: tuple[DoseEvent, ...], model: PKModel) -> np.ndarray:
    """Concentration over ``t`` for a dose schedule, integrating the nonlinear ODE."""
    t_arr = np.asarray(t, dtype=float)
    flat = t_arr.reshape(-1)
    order = np.argsort(flat)
    _, states = _nonlinear_states(flat, events, model)
    out = np.empty_like(flat)
    out[order] = np.clip(states[1] / model.v1_l, 0.0, None)
    return out.reshape(t_arr.shape)


def enzyme_activity(t, events: Iterable[DoseEvent], model: PKModel) -> np.ndarray:
    """Active fraction of the inactivatable enzyme over ``t`` (1 = baseline)."""
    t_arr = np.asarray(t, dtype=float)
    if not model.inactivation_fraction:
        return np.ones_like(t_arr)
    flat = t_arr.reshape(-1)
    order = np.argsort(flat)
    _, states = _nonlinear_states(flat, tuple(events), model)
    out = np.empty_like(flat)
    out[order] = states[3]
    return out.reshape(t_arr.shape)


@lru_cache(maxsize=256)
def _nonlinear_single_peak(amount: float, model: PKModel) -> tuple[float, float]:
    """(tmax, cmax) of a single dose, located numerically."""
    if model.route is Route.IV_BOLUS:
        return 0.0, amount / model.v1_l
    lin_tmax = tmax_single_dose(_linearized(model))
    coarse = np.linspace(0.0, max(4.0 * lin_tmax, model.lag_h + 1.0), 801)
    c = _nonlinear_profile(coarse, (DoseEvent(0.0, amount),), model)
    i = int(np.argmax(c))
    lo, hi = coarse[max(i - 1, 0)], coarse[min(i + 1, len(coarse) - 1)]
    fine = np.linspace(lo, hi, 401)
    cf = _nonlinear_profile(fine, (DoseEvent(0.0, amount),), model)
    j = int(np.argmax(cf))
    return float(fine[j]), float(cf[j])


@dataclass(frozen=True)
class _NonlinearSS:
    t: np.ndarray
    c: np.ndarray
    summary: SteadyState


@lru_cache(maxsize=128)
def _nonlinear_steady_state(amount: float, tau_h: float, model: PKModel) -> _NonlinearSS:
    """Repeat the dosing interval until the trough converges (max 2000 doses)."""
    rhs = _nonlinear_rhs(model)
    atol = _atol(amount)
    grid = np.linspace(0.0, tau_h, 401)
    y = _initial_state()
    troughs: list[float] = []
    first_auc = None
    for _ in range(2000):
        y = _dose_into(y, amount, model)
        y, ys = _integrate(rhs, y, 0.0, tau_h, grid, atol)
        c = np.clip(ys[1] / model.v1_l, 0.0, None)
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
    return _NonlinearSS(grid, c, summary)
