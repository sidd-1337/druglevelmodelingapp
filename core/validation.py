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
    if metric == "auc_inf":
        if not model.is_saturable:
            return model.f_eff * amount / model.clearance_l_per_h
        horizon = 20.0 * model.t_half_h
    elif metric == "auc_0_24":
        horizon = 24.0
    else:
        raise ValueError(f"unknown metric {metric!r}")
    t = np.linspace(0.0, horizon, 4001)
    c = concentration_profile(t, (DoseEvent(0.0, amount),), model)
    return auc_trapezoid(t, c)


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


def main() -> None:
    rows = validate()
    print(f"{'drug':14} {'metric':9} {'regimen':>12} {'observed':>10} {'predicted':>10} {'pred/obs':>8}  ok")
    for r in rows:
        regimen = f"{r.dose_mg:g} mg" + (f" q{r.interval_h:g}h" if r.interval_h else "")
        flag = "yes" if r.within else f"NO (>{r.max_fold:g}x)"
        print(f"{r.drug_id:14} {r.metric:9} {regimen:>12} {r.observed:>10.4g} {r.predicted:>10.4g} {r.ratio:>8.2f}  {flag}")
    exposure = [r for r in rows if r.metric in ("auc_inf", "auc_0_24", *STEADY_STATE_METRICS)]
    gm = float(np.exp(np.mean([np.log(r.fold_error) for r in exposure])))
    print(f"\nexposure metrics: {sum(r.within for r in exposure)}/{len(exposure)} within tolerance, "
          f"geometric mean fold error {gm:.2f}")


if __name__ == "__main__":
    main()
