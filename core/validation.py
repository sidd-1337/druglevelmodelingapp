"""Compare preset predictions with published PK data.

Three independent back-tests, all printed by ``python -m core.validation``:

* ``validate()`` - summary metrics (Cmax, Tmax, AUC, steady state) against
  ``data/reference_pk.json``. Answers "is the exposure right?".
* ``validate_interactions()`` - measured interaction factors against the ones
  fm x residual enzyme activity predicts. Answers "is the mechanism right?".
* ``validate_curves()`` - whole concentration-time curves against
  ``data/reference_curves/``. Answers "is the *shape* right?", which the
  summary metrics cannot: a model with the right Cmax and AUC can still have
  the wrong absorption rate or a missing distribution phase.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from core.drug_database import DEFAULT_DB_PATH, get_drug
from core.pk_engine import (
    DoseEvent,
    Route,
    auc_trapezoid,
    cmax_single_dose,
    concentration_profile,
    steady_state_summary,
    tmax_single_dose,
)

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
REFERENCE_PATH = DATA_DIR / "reference_pk.json"
REFERENCE_CURVE_DIR = DATA_DIR / "reference_curves"
SINGLE_DOSE_METRICS = ("cmax", "tmax_h", "auc_inf", "auc_0_24")
STEADY_STATE_METRICS = ("cavg_ss", "cmin_ss")
# Two equal doses ``interval_h`` apart: the combined profile versus superposition of two
# single doses (0-24 h), or the second dosing interval versus the first.
REDOSING_METRICS = ("auc_vs_superposition", "cmax_vs_superposition", "auc_second_vs_first", "cmax_second_vs_first")


@dataclass(frozen=True)
class ValidationRow:
    drug_id: str
    metric: str
    dose_mg: float
    interval_h: float | None
    observed: float
    predicted: float
    max_fold: float
    source: str
    verified: bool
    note: str

    @property
    def ratio(self) -> float:
        """Predicted / observed."""
        return self.predicted / self.observed

    @property
    def fold_error(self) -> float:
        """Symmetric fold error, >= 1."""
        r = self.ratio
        return max(r, 1.0 / r)

    @property
    def within(self) -> bool:
        return self.fold_error <= self.max_fold


def predict(drug_id: str, metric: str, dose_mg: float, interval_h: float | None = None,
            body_weight_kg: float = 70.0, db_path: str | Path = DEFAULT_DB_PATH) -> float:
    """Model prediction of one reference metric for a 70 kg adult (oral route)."""
    drug = get_drug(drug_id, db_path)
    model = drug.to_pk_model(body_weight_kg)
    amount = drug.amount_from_dose_mg(dose_mg)
    if metric in STEADY_STATE_METRICS:
        if interval_h is None:
            raise ValueError(f"{metric} needs interval_h")
        ss = steady_state_summary(amount, interval_h, model)
        return ss.cavg_ss if metric == "cavg_ss" else ss.cmin_ss
    if metric == "cmax":
        return cmax_single_dose(amount, model)
    if metric == "tmax_h":
        return tmax_single_dose(model, amount)
    if metric in REDOSING_METRICS:
        if interval_h is None:
            raise ValueError(f"{metric} needs interval_h (time between the doses)")
        return _redosing_ratio(metric, amount, interval_h, model)
    if metric == "auc_inf":
        if not model.needs_integration:
            return model.f_eff * amount / model.clearance_l_per_h
        horizon = 20.0 * model.terminal_half_life_h
    elif metric == "auc_0_24":
        horizon = 24.0
    else:
        raise ValueError(f"unknown metric {metric!r}")
    t = np.linspace(0.0, horizon, 4001)
    c = concentration_profile(t, (DoseEvent(0.0, amount),), model)
    return auc_trapezoid(t, c)


def _redosing_ratio(metric: str, amount: float, gap_h: float, model) -> float:
    t = np.linspace(0.0, 24.0 + gap_h, 4001)
    both = concentration_profile(t, (DoseEvent(0.0, amount), DoseEvent(gap_h, amount)), model)
    if metric.endswith("_vs_superposition"):
        single = concentration_profile(t, (DoseEvent(0.0, amount),), model)
        superposed = single + np.interp(t - gap_h, t, single, left=0.0)
        w = t <= 24.0
        if metric.startswith("auc"):
            return auc_trapezoid(t[w], both[w]) / auc_trapezoid(t[w], superposed[w])
        return float(both[w].max() / superposed[w].max())
    first, second = t < gap_h, t >= gap_h
    if metric.startswith("auc"):
        return auc_trapezoid(t[second], both[second]) / auc_trapezoid(t[first], both[first])
    return float(both[second].max() / both[first].max())


@dataclass(frozen=True)
class InteractionValidationRow:
    """A measured interaction factor next to the one the mechanism predicts."""

    victim_id: str
    perpetrator: str
    enzyme: str
    fm: float
    residual: float
    observed: float
    predicted: float
    source: str
    max_fold: float = 2.0

    @property
    def ratio(self) -> float:
        return self.predicted / self.observed

    @property
    def fold_error(self) -> float:
        r = self.ratio
        return max(r, 1.0 / r)

    @property
    def within(self) -> bool:
        return self.fold_error <= self.max_fold


def validate_interactions(
    interactions_path: str | Path | None = None, db_path: str | Path = DEFAULT_DB_PATH,
) -> list[InteractionValidationRow]:
    """Back-test: predict each MEASURED interaction factor from fm x residual activity.

    Only pairs the mechanism can even address are covered - a quantified clearance
    effect, a perpetrator in the enzyme registry, and an enzyme the victim is
    actually cleared by. Renal handling, transporter and first-pass interactions are
    out of scope by construction, not by omission.
    """
    from core.enzymes import load_perpetrators
    from core.interactions import CLEARANCE, load_interactions, mechanistic_exposure_factor

    by_name = {p.name: p for p in load_perpetrators()}
    kwargs = {"path": interactions_path} if interactions_path else {}
    rows = []
    for measured in load_interactions(**kwargs):
        if measured.qualitative or measured.site != CLEARANCE or measured.estimated:
            continue
        perpetrator = by_name.get(measured.perpetrator)
        if perpetrator is None:
            continue
        drug = get_drug(measured.victim_id, db_path)
        fm = (drug.pathways or {}).get(perpetrator.enzyme, 0.0)
        if fm <= 0:
            continue
        residual = perpetrator.residual_activity
        predicted, _ = mechanistic_exposure_factor(drug, perpetrator.enzyme, residual)
        rows.append(InteractionValidationRow(
            victim_id=measured.victim_id, perpetrator=measured.perpetrator,
            enzyme=perpetrator.enzyme, fm=fm, residual=residual,
            observed=measured.exposure_factor, predicted=predicted, source=measured.source))
    return rows


# --------------------------------------------------------------- curve validation
# Published curves are quoted in whatever unit the paper used. Factors convert that
# unit into the drug's own ``concentration_unit``; an unlisted pair is an error rather
# than a guess, because a silent 1000-fold slip would look like a model failure.
UNIT_FACTORS: dict[tuple[str, str], float] = {
    ("ug/mL", "ug/mL"): 1.0,
    ("mg/L", "ug/mL"): 1.0,  # same thing, both spellings appear in print
    ("ng/mL", "ug/mL"): 1e-3,
    ("ug/L", "ug/mL"): 1e-3,
    ("mmol/L", "mmol/L"): 1.0,
    ("mEq/L", "mmol/L"): 1.0,  # lithium is monovalent, so mEq and mmol coincide
    ("umol/L", "mmol/L"): 1e-3,
}
# "observed": concentrations measured in subjects - the real test.
# "published_model": a curve generated from someone else's fitted compartmental
# parameters. Not a measurement, but still an external check of shape, because
# their model was fitted to data this one has never seen. Scored separately.
CURVE_KINDS = ("observed", "published_model")


# Mass to molar needs the molar mass of the species the paper weighed, which the
# unit alone never says: lithium is printed as ug/mL of Li+ (6.94 g/mol), and reading
# it as lithium carbonate (73.89 g/mol) would be a silent tenfold error. So the curve
# file must state it, and without it the pair still raises.
MASS_UNITS_IN_MG_PER_L = {"mg/L": 1.0, "ug/mL": 1.0, "ng/mL": 1e-3, "ug/L": 1e-3}
MOLAR_UNITS_IN_MMOL_PER_L = {"mmol/L": 1.0, "umol/L": 1e-3}


def convert_concentration(values, from_unit: str, to_unit: str, molar_mass_g_per_mol: float | None = None):
    """Concentrations in ``from_unit`` expressed in ``to_unit``."""
    values = np.asarray(values, dtype=float)
    if (from_unit, to_unit) in UNIT_FACTORS:
        return values * UNIT_FACTORS[(from_unit, to_unit)]
    if molar_mass_g_per_mol and from_unit in MASS_UNITS_IN_MG_PER_L and to_unit in MOLAR_UNITS_IN_MMOL_PER_L:
        # mg/L divided by g/mol is mmol/L
        return values * MASS_UNITS_IN_MG_PER_L[from_unit] / molar_mass_g_per_mol / MOLAR_UNITS_IN_MMOL_PER_L[to_unit]
    hint = " (a mass unit needs molar_mass_g_per_mol)" if to_unit in MOLAR_UNITS_IN_MMOL_PER_L else ""
    raise ValueError(f"no conversion from {from_unit!r} to {to_unit!r}{hint}")


@dataclass(frozen=True)
class CurveFit:
    """How well a predicted curve matches an observed one, point by point."""

    times_h: np.ndarray
    observed: np.ndarray
    predicted: np.ndarray

    def __post_init__(self) -> None:
        if not (self.times_h.shape == self.observed.shape == self.predicted.shape):
            raise ValueError("times, observed and predicted must have the same shape")
        if self.observed.size == 0:
            raise ValueError("no points to score")
        if np.any(self.observed <= 0) or np.any(self.predicted <= 0):
            raise ValueError("fold error is undefined at a zero or negative concentration")

    @property
    def log_errors(self) -> np.ndarray:
        return np.log10(self.predicted / self.observed)

    @property
    def afe(self) -> float:
        """Average fold error: >1 means the model runs high, <1 means low. Bias."""
        return float(10.0 ** np.mean(self.log_errors))

    @property
    def aafe(self) -> float:
        """Absolute average fold error, >= 1. Precision, regardless of direction.

        Known in PBPK papers as GMFE - the geometric mean of the per-point fold
        errors - and the same quantity :class:`ValidationRow` reports for a single
        metric, so curve and summary scores are directly comparable.
        """
        return float(10.0 ** np.mean(np.abs(self.log_errors)))

    @property
    def fraction_within_2fold(self) -> float:
        return float(np.mean(np.abs(self.log_errors) <= np.log10(2.0)))


def score_curve(times_h, observed, predicted) -> CurveFit:
    """Score two concentration vectors in the same unit at the same times."""
    return CurveFit(*(np.asarray(a, dtype=float) for a in (times_h, observed, predicted)))


def predict_curve(
    drug_id: str,
    times_h,
    doses: list[tuple[float, float]],
    body_weight_kg: float = 70.0,
    route: Route = Route.ORAL,
    db_path: str | Path = DEFAULT_DB_PATH,
):
    """Predicted concentrations at ``times_h`` for a (time_h, dose_mg) schedule."""
    drug = get_drug(drug_id, db_path)
    model = drug.to_pk_model(body_weight_kg, route)
    events = tuple(DoseEvent(t, drug.amount_from_dose_mg(mg)) for t, mg in doses)
    return concentration_profile(np.asarray(times_h, dtype=float), events, model)


@dataclass(frozen=True)
class CurveValidationRow:
    """One published curve against the model, with the fit already scored."""

    curve_id: str
    drug_id: str
    label: str
    kind: str
    fit: CurveFit
    max_aafe: float
    source: str
    verified: bool
    note: str

    @property
    def n_points(self) -> int:
        return int(self.fit.observed.size)

    @property
    def within(self) -> bool:
        return self.fit.aafe <= self.max_aafe


def load_reference_curves(directory: str | Path = REFERENCE_CURVE_DIR) -> list[dict]:
    """Every curve file in ``directory``, in filename order. Missing directory: none."""
    directory = Path(directory)
    if not directory.is_dir():
        return []
    curves = []
    for path in sorted(directory.glob("*.json")):
        with open(path, encoding="utf-8") as fh:
            curve = json.load(fh)
        curve.setdefault("curve_id", path.stem)
        curves.append(curve)
    return curves


def validate_curve(curve: dict, db_path: str | Path = DEFAULT_DB_PATH) -> CurveValidationRow:
    """Score one curve definition (see ``data/reference_curves/README.md``)."""
    kind = curve.get("kind", "observed")
    if kind not in CURVE_KINDS:
        raise ValueError(f"{curve.get('curve_id')}: kind must be one of {CURVE_KINDS}")
    points = sorted((float(t), float(c)) for t, c in curve["points"])
    drug = get_drug(curve["drug_id"], db_path)
    molar_mass = curve.get("molar_mass_g_per_mol")
    obs = convert_concentration([c for _, c in points], curve["unit"], drug.concentration_unit, molar_mass)
    times = np.array([t for t, _ in points])
    # Below the assay's limit of quantification a "concentration" is a rounding
    # artefact; scoring it would let one unquantifiable tail point dominate the fit.
    keep = obs > convert_concentration(curve.get("lloq", 0.0), curve["unit"], drug.concentration_unit, molar_mass)
    if not keep.any():
        raise ValueError(f"{curve.get('curve_id')}: every point is at or below the LLOQ")
    times, obs = times[keep], obs[keep]
    doses = [(float(d.get("time_h", 0.0)), float(d["dose_mg"])) for d in curve["doses"]]
    predicted = predict_curve(
        curve["drug_id"], times, doses,
        body_weight_kg=float(curve.get("body_weight_kg", 70.0)),
        route=Route(curve.get("route", "oral")), db_path=db_path,
    )
    if np.any(predicted <= 0):
        # Almost always a sample taken before the preset's absorption lag, where the
        # model says there is no drug at all and the fold error is infinite. That is a
        # real structural disagreement, not a rounding issue: either the early points
        # belong below ``lloq``, or this preset's ``tlag_h`` is too long.
        early = times[predicted <= 0]
        raise ValueError(
            f"{curve.get('curve_id')}: the model predicts zero at t = "
            f"{', '.join(f'{t:g}' for t in early)} h, so those points cannot be scored. "
            f"Exclude them with \"lloq\" or correct the preset's absorption lag."
        )
    fit = CurveFit(times, obs, predicted)
    return CurveValidationRow(
        curve_id=curve.get("curve_id", curve["drug_id"]), drug_id=curve["drug_id"],
        label=curve.get("label", curve["drug_id"]), kind=kind, fit=fit,
        max_aafe=float(curve.get("max_aafe", 1.5)), source=curve["source"],
        verified=bool(curve.get("verified", False)), note=curve.get("note", ""),
    )


def validate_curves(
    directory: str | Path = REFERENCE_CURVE_DIR, db_path: str | Path = DEFAULT_DB_PATH,
) -> list[CurveValidationRow]:
    return [validate_curve(c, db_path) for c in load_reference_curves(directory)]


def load_references(path: str | Path = REFERENCE_PATH) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)["references"]


def validate(path: str | Path = REFERENCE_PATH, db_path: str | Path = DEFAULT_DB_PATH) -> list[ValidationRow]:
    rows = []
    for ref in load_references(path):
        interval = ref.get("interval_h")
        rows.append(ValidationRow(
            drug_id=ref["drug_id"],
            metric=ref["metric"],
            dose_mg=float(ref["dose_mg"]),
            interval_h=float(interval) if interval is not None else None,
            observed=float(ref["value"]),
            predicted=predict(ref["drug_id"], ref["metric"], float(ref["dose_mg"]), interval, db_path=db_path),
            max_fold=float(ref.get("max_fold", 1.5)),
            source=ref["source"],
            verified=bool(ref.get("verified", False)),
            note=ref.get("note", ""),
        ))
    return rows


METRIC_GROUPS = {
    "exposure": ("auc_inf", "auc_0_24", *STEADY_STATE_METRICS),
    "cmax": ("cmax",),
    "tmax": ("tmax_h",),
    "redosing": REDOSING_METRICS,
}


@dataclass(frozen=True)
class GroupSummary:
    n: int
    within: int
    geometric_mean_fold_error: float


def summarize(rows: list[ValidationRow]) -> dict[str, GroupSummary]:
    """Per metric group: count, number within tolerance, geometric mean fold error."""
    out = {}
    for name, metrics in METRIC_GROUPS.items():
        sel = [r for r in rows if r.metric in metrics]
        if sel:
            gm = float(np.exp(np.mean([np.log(r.fold_error) for r in sel])))
            out[name] = GroupSummary(len(sel), sum(r.within for r in sel), gm)
    return out


def regimen_label(row: ValidationRow) -> str:
    if row.metric in REDOSING_METRICS:
        return f"2 x {row.dose_mg:g} mg, {row.interval_h:g} h apart"
    return f"{row.dose_mg:g} mg" + (f" q{row.interval_h:g}h" if row.interval_h else "")


def main() -> None:
    rows = validate()
    print(f"{'drug':14} {'metric':20} {'regimen':>26} {'observed':>10} {'predicted':>10} {'pred/obs':>8}  ok")
    for r in rows:
        flag = "yes" if r.within else f"NO (>{r.max_fold:g}x)"
        print(f"{r.drug_id:14} {r.metric:20} {regimen_label(r):>26} {r.observed:>10.4g} {r.predicted:>10.4g} "
              f"{r.ratio:>8.2f}  {flag}")
    print()
    for name, g in summarize(rows).items():
        print(f"{name:9} {g.within}/{g.n} within tolerance, geometric mean fold error {g.geometric_mean_fold_error:.2f}")

    _print_interaction_gate()
    _print_curve_table()


def _print_interaction_gate() -> None:
    ddi = validate_interactions()
    print("\n--- mechanistic back-test of measured interaction factors ---")
    if not ddi:
        print("no measured pair is addressable by the mechanism (needs fm and a known perpetrator)")
        return
    print(f"{'victim':14} {'perpetrator':30} {'enzyme':9} {'fm':>5} {'resid':>6} "
          f"{'observed':>9} {'predicted':>10} {'pred/obs':>9}  ok")
    for r in ddi:
        print(f"{r.victim_id:14} {r.perpetrator:30} {r.enzyme:9} {r.fm:>5.2f} {r.residual:>6.3f} "
              f"{r.observed:>9.2f} {r.predicted:>10.2f} {r.ratio:>9.2f}  {'yes' if r.within else 'NO'}")
    gm = float(np.exp(np.mean([np.log(r.fold_error) for r in ddi])))
    print(f"\n{sum(r.within for r in ddi)}/{len(ddi)} within 2x, geometric mean fold error {gm:.2f}")


def _print_curve_table() -> None:
    curves = validate_curves()
    print("\n--- whole-curve fit against published concentration-time data ---")
    if not curves:
        print(f"no reference curves in {REFERENCE_CURVE_DIR}/ "
              f"(see its README for the file format)")
        return
    print(f"{'curve':34} {'kind':16} {'n':>3} {'AFE':>6} {'AAFE':>6} {'<2x':>5}  ok")
    for r in curves:
        flag = "yes" if r.within else f"NO (>{r.max_aafe:g}x)"
        print(f"{r.curve_id:34} {r.kind:16} {r.n_points:>3} {r.fit.afe:>6.2f} "
              f"{r.fit.aafe:>6.2f} {r.fit.fraction_within_2fold:>5.0%}  {flag}")
    for kind in CURVE_KINDS:
        sel = [r for r in curves if r.kind == kind]
        if sel:
            gm = float(10.0 ** np.mean([np.log10(r.fit.aafe) for r in sel]))
            print(f"{kind:16} {sum(r.within for r in sel)}/{len(sel)} within tolerance, "
                  f"geometric mean AAFE {gm:.2f}")


if __name__ == "__main__":
    main()
