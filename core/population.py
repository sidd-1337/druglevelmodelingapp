"""Patient covariates and inter-individual variability (no UI, no I/O).

Clearance model
---------------
CL = CL_ref * [(1 - fe) * (WT / 70) ** 0.75 + fe * CrCl / 100]

``CL_ref`` is the label clearance of a 70 kg adult with normal renal function
(taken as CrCl 100 mL/min), ``fe`` the fraction of it that tracks creatinine
clearance. Without a CrCl the whole clearance scales allometrically with weight.
Volume scales linearly with weight, so t1/2 = ln2 * V / CL follows from both.
"""
from __future__ import annotations

import math
from dataclasses import replace

from core.pk_engine import PKModel, Route

REFERENCE_WEIGHT_KG = 70.0
REFERENCE_CRCL_ML_MIN = 100.0
ALLOMETRIC_EXPONENT = 0.75
Z_90 = 1.6448536269514722  # 5th-95th percentile


def cockcroft_gault(age_y: float, weight_kg: float, scr_umol_l: float, female: bool) -> float:
    """Creatinine clearance [mL/min] from serum creatinine in umol/L."""
    if age_y <= 0 or weight_kg <= 0 or scr_umol_l <= 0:
        raise ValueError("age, weight and serum creatinine must be > 0")
    crcl = (140.0 - age_y) * weight_kg * 1.23 / scr_umol_l
    return max(crcl * (0.85 if female else 1.0), 1.0)


def clearance_factor(weight_kg: float, crcl_ml_min: float | None = None, fe_renal: float = 0.0) -> float:
    """Patient clearance relative to the 70 kg, CrCl 100 mL/min reference."""
    if weight_kg <= 0:
        raise ValueError("weight_kg must be > 0")
    if not 0.0 <= fe_renal <= 1.0:
        raise ValueError("fe_renal must be in [0, 1]")
    size = (weight_kg / REFERENCE_WEIGHT_KG) ** ALLOMETRIC_EXPONENT
    if crcl_ml_min is None or fe_renal == 0.0:
        return size
    if crcl_ml_min <= 0:
        raise ValueError("crcl_ml_min must be > 0")
    return (1.0 - fe_renal) * size + fe_renal * crcl_ml_min / REFERENCE_CRCL_ML_MIN


def log_sd(cv: float) -> float:
    """SD of log(X) for a log-normal X with coefficient of variation ``cv``."""
    if cv < 0:
        raise ValueError("cv must be >= 0")
    return math.sqrt(math.log1p(cv * cv))


def exposure_band_models(
    model: PKModel, cv_cl: float, cv_f: float = 0.0, z: float = Z_90
) -> tuple[PKModel, PKModel]:
    """(low, high) exposure models at the z-quantiles of log-normal CL and F variability.

    The z-shift of total log exposure is split between CL and F in proportion to
    their variances, i.e. each model is the most likely (CL, F) pair for a patient
    at that exposure percentile. CL moves at constant V, so t1/2 moves with it.
    """
    w_cl = log_sd(cv_cl)
    w_f = 0.0 if model.route is Route.IV_BOLUS else log_sd(cv_f)
    w_tot = math.hypot(w_cl, w_f)
    if w_tot == 0.0:
        return model, model
    s_cl = z * w_cl * w_cl / w_tot
    s_f = z * w_f * w_f / w_tot

    def shifted(sign: float) -> PKModel:
        f = min(1.0, model.bioavailability * math.exp(sign * s_f))
        return replace(model, t_half_h=model.t_half_h * math.exp(sign * s_cl), bioavailability=f)

    return shifted(-1.0), shifted(+1.0)
