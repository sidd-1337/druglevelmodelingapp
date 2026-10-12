"""Bayesian (MAP) individualisation from measured drug levels.

The population model is the prior: log CL and log V of the patient are normal
around the preset values with the published between-patient SDs. Measured levels
carry a proportional (log-normal) error. The maximum a posteriori estimate is the
individual model; a Laplace approximation of the posterior gives its uncertainty.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Sequence

import numpy as np
from scipy.optimize import minimize

from core.pk_engine import DoseEvent, PKModel, Route, concentration_profile
from core.population import Z_90, log_sd

MIN_PRIOR_SD = 0.05  # keeps the prior proper when a CV of 0 is entered
DEFAULT_CV_V = 0.25
DEFAULT_RESIDUAL_CV = 0.2  # assay + model misspecification
_TINY = 1e-300


@dataclass(frozen=True)
class Observation:
    time_h: float  # since the start of the dosing schedule
    concentration: float


@dataclass(frozen=True)
class IndividualFit:
    model: PKModel
    population_model: PKModel
    eta: tuple[float, float]  # (log CL_i/CL_pop, log V_i/V_pop)
    cov: tuple[tuple[float, float], tuple[float, float]]  # posterior covariance of eta
    observations: tuple[Observation, ...]
    predicted: tuple[float, ...]
    population_predicted: tuple[float, ...]

    @property
    def cl_ratio(self) -> float:
        return math.exp(self.eta[0])

    @property
    def v_ratio(self) -> float:
        return math.exp(self.eta[1])


def _apply_eta(pop: PKModel, eta: Sequence[float]) -> PKModel:
    cl = pop.clearance_l_per_h * math.exp(eta[0])
    v = pop.vd_l * math.exp(eta[1])
    return replace(pop, vd_l=v, t_half_h=math.log(2.0) * v / cl)


def _hessian(f, x: np.ndarray, h: float = 1e-3) -> np.ndarray:
    n = len(x)
    H = np.empty((n, n))
    for i in range(n):
        for j in range(i, n):
            ei, ej = np.eye(n)[i] * h, np.eye(n)[j] * h
            H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h * h)
    return H


def fit_individual(
    population_model: PKModel,
    events: Sequence[DoseEvent],
    observations: Sequence[Observation],
    cv_cl: float,
    cv_f: float = 0.0,
    cv_v: float = DEFAULT_CV_V,
    residual_cv: float = DEFAULT_RESIDUAL_CV,
) -> IndividualFit:
    """MAP estimate of the patient's clearance and volume given measured levels.

    F is not estimated separately: its variability is folded into the prior of
    apparent clearance (CL/F), which is what oral levels inform.
    """
    obs = tuple(sorted(observations, key=lambda o: o.time_h))
    if not obs:
        raise ValueError("at least one observation is required")
    if not events:
        raise ValueError("at least one dose is required")
    first_dose = min(e.time_h for e in events)
    if any(o.concentration <= 0 for o in obs):
        raise ValueError("measured concentrations must be > 0")
    if any(o.time_h <= first_dose for o in obs):
        raise ValueError("measurements must be taken after the first dose")
    if residual_cv <= 0:
        raise ValueError("residual_cv must be > 0")

    w_f = 0.0 if population_model.route is Route.IV_BOLUS else log_sd(cv_f)
    w_cl = max(math.hypot(log_sd(cv_cl), w_f), MIN_PRIOR_SD)
    w_v = max(log_sd(cv_v), MIN_PRIOR_SD)
    sigma = log_sd(residual_cv)
    t_obs = np.array([o.time_h for o in obs])
    log_c = np.log([o.concentration for o in obs])

    def predict(eta) -> np.ndarray:
        return concentration_profile(t_obs, events, _apply_eta(population_model, eta))

    def objective(eta) -> float:
        # -2 log posterior up to a constant
        resid = (log_c - np.log(np.maximum(predict(eta), _TINY))) / sigma
        return float(resid @ resid + (eta[0] / w_cl) ** 2 + (eta[1] / w_v) ** 2)

    res = minimize(objective, np.zeros(2), method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-10, maxiter=2000))
    eta = np.asarray(res.x, dtype=float)
    H = _hessian(objective, eta)
    try:
        cov = 2.0 * np.linalg.inv(H)
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        cov = np.diag([w_cl ** 2, w_v ** 2])

    return IndividualFit(
        model=_apply_eta(population_model, eta),
        population_model=population_model,
        eta=(float(eta[0]), float(eta[1])),
        cov=((float(cov[0, 0]), float(cov[0, 1])), (float(cov[1, 0]), float(cov[1, 1]))),
        observations=obs,
        predicted=tuple(float(v) for v in predict(eta)),
        population_predicted=tuple(float(v) for v in predict(np.zeros(2))),
    )


def individual_band(
    fit: IndividualFit, t: np.ndarray, events: Sequence[DoseEvent], z: float = Z_90,
    n_samples: int | None = None, seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """5th-95th percentile (for z = 1.645) of the patient's profile given the fit uncertainty."""
    n = n_samples or (60 if fit.model.is_nonlinear else 300)
    rng = np.random.default_rng(seed)
    etas = rng.multivariate_normal(np.array(fit.eta), np.array(fit.cov), size=n)
    curves = np.array([concentration_profile(t, events, _apply_eta(fit.population_model, e)) for e in etas])
    q = 100.0 * float(0.5 * math.erfc(z / math.sqrt(2.0)))
    return np.percentile(curves, q, axis=0), np.percentile(curves, 100.0 - q, axis=0)
