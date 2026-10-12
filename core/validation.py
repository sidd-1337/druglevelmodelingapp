"""Compare preset predictions with published PK data (data/reference_pk.json).

Run ``python -m core.validation`` for a table of predicted / observed ratios.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from core.drug_database import DEFAULT_DB_PATH, get_drug
from core.pk_engine import (
    DoseEvent,
    auc_trapezoid,
    cmax_single_dose,
    concentration_profile,
    steady_state_summary,
    tmax_single_dose,
)

REFERENCE_PATH = Path(__file__).resolve().parent.parent / "data" / "reference_pk.json"
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
    from core.enzymes import derived_exposure_factor, load_perpetrators
    from core.interactions import CLEARANCE, load_interactions

    by_name = {p.name: p for p in load_perpetrators()}
    kwargs = {"path": interactions_path} if interactions_path else {}
    rows = []
    for measured in load_interactions(**kwargs):
        if measured.qualitative or measured.site != CLEARANCE or measured.estimated:
            continue
        perpetrator = by_name.get(measured.perpetrator)
        if perpetrator is None:
            continue
        fm = (get_drug(measured.victim_id, db_path).pathways or {}).get(perpetrator.enzyme, 0.0)
        if fm <= 0:
            continue
        residual = perpetrator.residual_activity
        rows.append(InteractionValidationRow(
            victim_id=measured.victim_id, perpetrator=measured.perpetrator,
            enzyme=perpetrator.enzyme, fm=fm, residual=residual,
            observed=measured.exposure_factor,
            predicted=derived_exposure_factor(fm, residual), source=measured.source))
    return rows


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

    ddi = validate_interactions()
    print(f"\n--- mechanistic back-test of measured interaction factors ---")
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


if __name__ == "__main__":
    main()
