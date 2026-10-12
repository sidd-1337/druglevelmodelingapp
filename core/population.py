"""Patient covariates and inter-individual variability (no UI, no I/O).

Clearance model
---------------
Clearance is split between elimination pathways, each carrying a fraction ``f_p``
of the reference clearance and scaling with its own patient factor ``a_p``::

    CL / CL_ref = sum_p f_p * a_p  +  (1 - sum_p f_p) * size

    size     = (WT / 70) ** 0.75
    a_renal  = CrCl / 100, or ``size`` when no CrCl is known
    a_enzyme = size * enzyme activity (genotype; 1.0 for a normal metaboliser)

The unassigned remainder is treated as size-scaled and unaffected by genotype.
``CL_ref`` is the clearance of a 70 kg adult with normal renal function (taken as
CrCl 100 mL/min) and normal enzyme activity, so that reference patient gives 1.0.
Volume scales linearly with weight, so t1/2 = ln2 * V / CL follows from both.
"""
from __future__ import annotations

import math
from dataclasses import replace
from typing import Mapping

from core.pk_engine import PKModel, Route

REFERENCE_WEIGHT_KG = 70.0
REFERENCE_CRCL_ML_MIN = 100.0
ALLOMETRIC_EXPONENT = 0.75
Z_90 = 1.6448536269514722  # 5th-95th percentile
RENAL = "renal"  # the one pathway that scales with CrCl instead of body size


def cockcroft_gault(age_y: float, weight_kg: float, scr_umol_l: float, female: bool) -> float:
    """Creatinine clearance [mL/min] from serum creatinine in umol/L."""
    if age_y <= 0 or weight_kg <= 0 or scr_umol_l <= 0:
        raise ValueError("age, weight and serum creatinine must be > 0")
    crcl = (140.0 - age_y) * weight_kg * 1.23 / scr_umol_l
    return max(crcl * (0.85 if female else 1.0), 1.0)


def pathway_clearance_factor(
    weight_kg: float,
    pathways: Mapping[str, float] | None = None,
    crcl_ml_min: float | None = None,
    enzyme_activity: Mapping[str, float] | None = None,
) -> float:
    """Patient clearance relative to the 70 kg, CrCl 100 mL/min, normal-metaboliser reference.

    ``pathways`` maps a pathway name to its fraction of reference clearance; the
    unassigned remainder scales with body size only. Renal clearance tracks CrCl
    and gets no allometric term (CrCl already carries body size), so with an
    unknown CrCl it falls back to ``size`` and the whole expression reduces to it.
    """
    if weight_kg <= 0:
        raise ValueError("weight_kg must be > 0")
    if crcl_ml_min is not None and crcl_ml_min <= 0:
        raise ValueError("crcl_ml_min must be > 0")
    size = (weight_kg / REFERENCE_WEIGHT_KG) ** ALLOMETRIC_EXPONENT
    renal_rel = size if crcl_ml_min is None else crcl_ml_min / REFERENCE_CRCL_ML_MIN

    assigned = total = 0.0
    for name, fraction in (pathways or {}).items():
        if not 0.0 <= fraction <= 1.0:
            raise ValueError(f"pathway {name!r}: fraction must be in [0, 1]")
        activity = renal_rel if name == RENAL else size * (enzyme_activity or {}).get(name, 1.0)
        assigned += fraction
        total += fraction * activity
    if assigned > 1.0 + 1e-9:
        raise ValueError("pathway fractions must sum to <= 1")
    return total + (1.0 - assigned) * size


def clearance_factor(weight_kg: float, crcl_ml_min: float | None = None, fe_renal: float = 0.0) -> float:
    """Renal-only special case of :func:`pathway_clearance_factor`."""
    if not 0.0 <= fe_renal <= 1.0:
        raise ValueError("fe_renal must be in [0, 1]")
    return pathway_clearance_factor(weight_kg, {RENAL: fe_renal} if fe_renal else None, crcl_ml_min)


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
