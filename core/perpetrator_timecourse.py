"""Time course of a drug-drug interaction: when the perpetrator's effect builds and fades.

The static model in :mod:`core.interactions` multiplies a victim's exposure by a
constant factor - the full effect of the perpetrator's stated regimen, applied the
instant the perpetrator is added. This module keeps that factor's *magnitude* (the
measured one where a study exists, the mechanistic estimate otherwise) and gives it
a *time course*, so a perpetrator started or stopped part-way through a simulation
is right, and a slow one does not act before it has accumulated.

The effect is written as a factor on the victim's exposure, factor(t) = 1 + (F - 1)
* g(t), where F is the steady-state exposure factor - the measured AUC ratio where a
study exists, the mechanistic estimate otherwise - and g(t) in [0, 1] is how far the
perpetrator's effect has developed. The factor scales the victim's reported
concentration curve: an exposure factor is an empirical ratio of exposure, so it is
applied to exposure, which reproduces the measured F exactly for any victim,
including a drug reported as a parent-plus-metabolite moiety where scaling a single
clearance would not. The victim's own re-equilibration half-life is not modelled
separately; the curve follows the perpetrator's time course, which for the slow
perpetrators that matter here (SSRIs, inducers) is the dominant one.

Two mechanisms set g(t):

* **Reversible inhibition.** The effect tracks the inhibitor's own concentration,
  with no lag of its own, so g(t) is the inhibitor's exposure as a fraction of its
  steady-state level: it builds as 1 - exp(-ke (t - start)) while the drug is taken
  and decays as exp(-ke (t - stop)) after it is stopped, ke = ln2 / t_half. The
  within-interval ripple is left out deliberately - it is far faster than any
  enzyme effect. Onset and offset are the inhibitor's own half-life.

* **Induction.** The inducer raises enzyme *synthesis*, but the active enzyme pool
  builds and decays only at its own turnover rate, so induction lags the inducer by
  days. The enzyme activity E (1 = baseline) follows dE/dt = kdeg (S(t) - E) with
  S(t) = 1 + (r - 1) * h(t), r the enzyme's induced activity and h the inducer's
  exposure fraction; g(t) = (E(t) - 1) / (r - 1) is that same build-up and decay
  carried onto the victim's factor. The slow approach, and the equally slow return
  after the inducer stops, are the clinically important part a static factor misses.

The magnitude F is supplied by the caller (measured study or mechanistic estimate);
this module only shapes its timing from the perpetrator's pharmacokinetics.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

LN2 = math.log(2.0)


@dataclass(frozen=True)
class PerpetratorCourse:
    """How one perpetrator's effect on one enzyme develops over a stated window.

    ``residual`` is the enzyme activity its regimen leaves at steady state (< 1 for
    an inhibitor, > 1 for an inducer); it shapes the induction lag and the display
    curve, not the victim magnitude. ``t_half_h`` is the perpetrator's own
    elimination half-life, which sets how fast its concentration - and so its
    effect - builds. ``offset_t_half_h`` sets how fast the effect fades after it is
    stopped, where that is slower than its own clearance (an active metabolite, a
    tissue reservoir, or an enzyme it destroyed having to be resynthesised); it
    defaults to ``t_half_h``. ``start_h`` is when it is first taken, ``None`` meaning
    it was already being taken at steady state before the simulation begins;
    ``stop_h`` is when it is stopped (``None`` = throughout). ``interval_h`` records
    the regimen the residual was measured at.
    """

    enzyme: str
    residual: float
    t_half_h: float
    induces: bool
    interval_h: float = 24.0
    start_h: float | None = 0.0
    stop_h: float | None = None
    offset_t_half_h: float | None = None

    def __post_init__(self) -> None:
        if self.residual < 0:
            raise ValueError("residual must be >= 0")
        if self.induces and self.residual <= 1:
            raise ValueError("an inducer must have residual > 1")
        if not self.induces and self.residual > 1:
            raise ValueError("a reversible inhibitor must have residual <= 1")
        if self.t_half_h <= 0 or self.interval_h <= 0:
            raise ValueError("t_half_h and interval_h must be > 0")
        if self.offset_t_half_h is not None and self.offset_t_half_h <= 0:
            raise ValueError("offset_t_half_h must be > 0")
        if self.stop_h is not None and self.start_h is not None and self.stop_h < self.start_h:
            raise ValueError("stop_h must be >= start_h")


def _exposure_envelope(course: PerpetratorCourse, t: np.ndarray) -> np.ndarray:
    """The perpetrator's exposure over ``t`` as a fraction of its steady-state level (0..1)."""
    t = np.asarray(t, dtype=float)
    ke = LN2 / course.t_half_h
    if course.start_h is None:  # already at steady state
        env = np.ones_like(t)
        at_stop = 1.0
    else:
        env = np.where(t >= course.start_h, 1.0 - np.exp(-ke * np.clip(t - course.start_h, 0.0, None)), 0.0)
        at_stop = None if course.stop_h is None else 1.0 - math.exp(-ke * max(course.stop_h - course.start_h, 0.0))
    if course.stop_h is not None:
        ke_off = LN2 / (course.offset_t_half_h or course.t_half_h)
        after = t >= course.stop_h
        env = np.where(after, at_stop * np.exp(-ke_off * np.clip(t - course.stop_h, 0.0, None)), env)
    return env


def _induction_pool(t: np.ndarray, h: np.ndarray, residual: float, kdeg_per_h: float) -> np.ndarray:
    """Active enzyme fraction over ``t`` for synthesis target S(t) = 1 + (residual-1)*h(t)."""
    if kdeg_per_h <= 0:
        raise ValueError("enzyme turnover kdeg must be > 0")
    target = 1.0 + (residual - 1.0) * np.asarray(h, dtype=float)
    e = np.empty_like(t)
    # Start at the pool's steady value for the exposure at t[0]: baseline for a
    # perpetrator not yet taken, fully induced for one already taken at steady state.
    e[0] = target[0]
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        decay = math.exp(-kdeg_per_h * dt)
        # Exact over a step for a piecewise-linear target (trapezoidal source).
        e[i] = e[i - 1] * decay + 0.5 * (target[i] + target[i - 1]) * (1.0 - decay)
    return e


def effect_fraction(course: PerpetratorCourse, t: np.ndarray, turnover_per_h: float | None = None) -> np.ndarray:
    """How far the perpetrator's effect has developed over ``t``, in [0, 1].

    0 = not yet acting, 1 = full stated effect. For an inducer this follows the
    enzyme's turnover (``turnover_per_h`` required); for a reversible inhibitor it
    follows the inhibitor's own concentration envelope.
    """
    env = _exposure_envelope(course, np.asarray(t, dtype=float))
    if not course.induces:
        return env
    if turnover_per_h is None:
        raise ValueError(f"{course.enzyme}: an inducer needs the enzyme turnover rate")
    pool = _induction_pool(np.asarray(t, dtype=float), env, course.residual, turnover_per_h)
    return (pool - 1.0) / (course.residual - 1.0)


def exposure_fraction(course: PerpetratorCourse, t: np.ndarray) -> np.ndarray:
    """The perpetrator's own exposure over ``t`` as a fraction of its steady state (no enzyme lag)."""
    return _exposure_envelope(course, np.asarray(t, dtype=float))


def activity_from_fraction(course: PerpetratorCourse, g: np.ndarray) -> np.ndarray:
    """Enzyme activity (1 = baseline) for an effect developed to ``g``; ``residual`` at g = 1.

    A reversible inhibitor acts as 1 / (1 + beta g) with beta = (1 - r) / r, the
    concentration-response of competitive inhibition; an inducer's pool is linear in g.
    """
    g = np.asarray(g, dtype=float)
    if course.induces:
        return 1.0 + (course.residual - 1.0) * g
    beta = (1.0 - course.residual) / max(course.residual, 1e-9)
    return 1.0 / (1.0 + beta * g)


def enzyme_activity_curve(course: PerpetratorCourse, t: np.ndarray,
                          turnover_per_h: float | None = None) -> np.ndarray:
    """The enzyme's activity over ``t`` (1 = baseline). Reaches ``residual`` at full effect."""
    return activity_from_fraction(course, effect_fraction(course, t, turnover_per_h))


def combined_activity(courses: list[PerpetratorCourse], t: np.ndarray,
                      turnover_per_h: float | None = None) -> np.ndarray:
    """One enzyme's activity over ``t`` with several perpetrators acting on it at once.

    Reversible inhibitors compete for the same enzyme, so their occupancies add rather
    than their activities multiplying: A = E / (1 + sum_j beta_j g_j), beta_j as in
    :func:`activity_from_fraction`. Inducers raise the same enzyme's synthesis, so they
    share one pool, dE/dt = kdeg (1 + sum_j (r_j - 1) h_j(t) - E). Inhibition then acts
    on the induced pool, the product of the two terms in the net-effect model (Fahmi
    2008). With one perpetrator this is exactly :func:`enzyme_activity_curve`.
    """
    t = np.asarray(t, dtype=float)
    if len({c.enzyme for c in courses}) > 1:
        raise ValueError("combined_activity takes perpetrators of one enzyme")
    occupancy = np.zeros_like(t)
    synthesis = np.zeros_like(t)
    for c in courses:
        if c.induces:
            synthesis = synthesis + (c.residual - 1.0) * _exposure_envelope(c, t)
        else:
            occupancy = occupancy + (1.0 - c.residual) / max(c.residual, 1e-9) * _exposure_envelope(c, t)
    pool = np.ones_like(t)
    if any(c.induces for c in courses):
        if turnover_per_h is None:
            raise ValueError(f"{courses[0].enzyme}: an inducer needs the enzyme turnover rate")
        pool = _induction_pool(t, synthesis, 2.0, turnover_per_h)
    return pool / (1.0 + occupancy)


@dataclass(frozen=True)
class ActingFactor:
    """A perpetrator's time-varying exposure factor on one victim.

    ``factor`` is the steady-state exposure factor (the measured AUC ratio where a
    study exists, else the mechanistic estimate); ``fraction`` is ``g(t)`` from
    :func:`effect_fraction` on the simulation grid. A static interaction is the case
    ``fraction`` identically 1.
    """

    factor: float
    fraction: np.ndarray

    def __post_init__(self) -> None:
        if self.factor <= 0:
            raise ValueError("factor must be > 0")


def apply_factor_course(c_alone: np.ndarray, factors: list[ActingFactor]) -> np.ndarray:
    """The victim's curve with the interaction factors applied over time.

    Each factor multiplies the exposure by ``1 + (factor - 1) * g(t)``, so at full
    development (g = 1) the curve is scaled by the steady-state factor and at g = 0 it
    is unchanged. The product over several perpetrators matches their combined static
    factor once all have accumulated.
    """
    c = np.asarray(c_alone, dtype=float).copy()
    for f in factors:
        c = c * (1.0 + (f.factor - 1.0) * np.asarray(f.fraction, dtype=float))
    return c


def combined_factor_at(factors: list[ActingFactor], index: int = -1) -> float:
    """The product of the factors at one grid point (default the last), for display."""
    product = 1.0
    for f in factors:
        product *= 1.0 + (f.factor - 1.0) * float(np.asarray(f.fraction)[index])
    return product
