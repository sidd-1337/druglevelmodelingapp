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


def pathway_shares(
    weight_kg: float,
    pathways: Mapping[str, float] | None = None,
    crcl_ml_min: float | None = None,
    enzyme_activity: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Each pathway's fraction of THIS patient's clearance.

    The reference fractions in ``pathways`` shift with the patient: a poor
    metaboliser clears less through the enzyme, so it carries less of what is left.
    """
    total = pathway_clearance_factor(weight_kg, pathways, crcl_ml_min, enzyme_activity)
    size = (weight_kg / REFERENCE_WEIGHT_KG) ** ALLOMETRIC_EXPONENT
    renal_rel = size if crcl_ml_min is None else crcl_ml_min / REFERENCE_CRCL_ML_MIN
    return {name: fraction * (renal_rel if name == RENAL else size * (enzyme_activity or {}).get(name, 1.0)) / total
            for name, fraction in (pathways or {}).items()}


def clearance_factor(weight_kg: float, crcl_ml_min: float | None = None, fe_renal: float = 0.0) -> float:
    """Renal-only special case of :func:`pathway_clearance_factor`."""
    if not 0.0 <= fe_renal <= 1.0:
        raise ValueError("fe_renal must be in [0, 1]")
    return pathway_clearance_factor(weight_kg, {RENAL: fe_renal} if fe_renal else None, crcl_ml_min)


def phenotype_log_sd(pathway_fraction: float, activities_with_frequencies: Mapping[float, float]) -> float:
    """Log-scale SD of clearance from an unknown phenotype on one pathway.

    ``activities_with_frequencies`` maps an activity multiplier to how common it is.
    A pathway carrying ``pathway_fraction`` of clearance turns activity ``a`` into a
    relative clearance of ``1 - f + f * a``; the spread of that mixture is what the
    band must widen by.
    """
    if not 0 <= pathway_fraction <= 1:
        raise ValueError("pathway_fraction must be in [0, 1]")
    pairs = [(1.0 - pathway_fraction + pathway_fraction * a, w)
             for a, w in activities_with_frequencies.items() if w > 0]
    total = sum(w for _, w in pairs)
    if not pairs or total <= 0:
        return 0.0
    mean = sum(cl * w for cl, w in pairs) / total
    if mean <= 0:
        return 0.0
    variance = sum(w * (cl - mean) ** 2 for cl, w in pairs) / total
    return log_sd(math.sqrt(variance) / mean)


def log_sd(cv: float) -> float:
    """SD of log(X) for a log-normal X with coefficient of variation ``cv``."""
    if cv < 0:
        raise ValueError("cv must be >= 0")
    return math.sqrt(math.log1p(cv * cv))


def exposure_band_models(
    model: PKModel, cv_cl: float, cv_f: float = 0.0, z: float = Z_90,
    extra_cl_log_sd: float = 0.0,
) -> tuple[PKModel, PKModel]:
    """(low, high) exposure models at the z-quantiles of log-normal CL and F variability.

    The z-shift of total log exposure is split between CL and F in proportion to
    their variances, i.e. each model is the most likely (CL, F) pair for a patient
    at that exposure percentile. CL moves at constant V, so t1/2 moves with it.

    ``extra_cl_log_sd`` adds a further source of clearance spread in quadrature -
    used for an unknown metaboliser phenotype. It widens the band without moving
    its centre, so the prediction the band is drawn around is unchanged.
    """
    w_cl = math.hypot(log_sd(cv_cl), extra_cl_log_sd)
    w_f = 0.0 if model.route is Route.IV_BOLUS else log_sd(cv_f)
    w_tot = math.hypot(w_cl, w_f)
    if w_tot == 0.0:
        return model, model
    s_cl = z * w_cl * w_cl / w_tot
    s_f = z * w_f * w_f / w_tot

    def shifted(sign: float) -> PKModel:
        f = min(1.0, model.bioavailability * math.exp(sign * s_f))
        scale = math.exp(sign * s_cl)
        met = model.metabolite
        if met is not None and met.reported != "parent":
            # When the reported quantity includes the metabolite, the metabolite's own
            # clearance is most of it. Shifting only the parent would leave the band
            # almost flat, which is exactly the buffering that makes the sum insensitive
            # to the parent's clearance - true of genotype, but not of a patient who
            # simply eliminates this drug faster or slower overall.
            met = replace(met, t_half_h=met.t_half_h * scale)
        return replace(model, t_half_h=model.t_half_h * scale, bioavailability=f, metabolite=met)

    return shifted(-1.0), shifted(+1.0)
