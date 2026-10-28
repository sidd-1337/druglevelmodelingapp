"""Print what a preset predicts for a 70 kg adult next to its published references.

Usage (from the repository root):
    python .claude/skills/drug-analysis/scripts/preset_report.py <drug_id>
        [--formulation pr] [--dose 300] [--interval 12]
        [--set vd_l_per_kg=8.5 --set metabolite.t_half_h=12 ...]
        [--ref cmax:1:0.0091 --ref cavg_ss:3:0.0207:24 ...]

--set edits a temporary copy of data/drugs.json, so you can try fits without touching
the file. --ref scores a candidate reference (metric:dose_mg:value[:interval_h]) without
adding it to data/reference_pk.json. References and curves are always scored against
the default formulation, which is what core.validation does.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))

from core.drug_database import DEFAULT_DB_PATH, get_drug  # noqa: E402
import numpy as np  # noqa: E402

from core.enzymes import load_enzymes  # noqa: E402
from core.pk_engine import (  # noqa: E402
    DoseEvent, auc_trapezoid, cmax_single_dose, concentration_profile, metabolite_profile, regular_schedule,
    steady_state_summary, tmax_single_dose,
)
from core.validation import predict, validate, validate_curves  # noqa: E402


def _apply_sets(drug_id: str, sets: list[str]) -> Path:
    raw = json.loads(DEFAULT_DB_PATH.read_text(encoding="utf-8"))
    entry = next((d for d in raw["drugs"] if d["id"] == drug_id), None)
    if entry is None:
        sys.exit(f"unknown drug id {drug_id!r}")
    for item in sets:
        key, _, value = item.partition("=")
        *parents, leaf = key.split(".")
        node = entry
        for p in parents:
            node = node.setdefault(p, {})
        node[leaf] = json.loads(value)
    tmp = Path(tempfile.mkdtemp()) / "drugs.json"
    tmp.write_text(json.dumps(raw), encoding="utf-8")
    return tmp


def _metabolite_ratios(model, amount: float, tau: float) -> None:
    """Steady-state metabolite/parent ratios over the last of ten dosing intervals (mass units)."""
    n = 10
    t = np.linspace(0.0, n * tau, n * 400 + 1)
    doses = regular_schedule(amount, tau, n * tau)
    last = t >= (n - 1) * tau
    # concentration_profile gives the REPORTED quantity; switch it to the parent alone.
    parent_only = replace(model, metabolite=replace(model.metabolite, reported="parent"))
    parent = concentration_profile(t, doses, parent_only)[last]
    met = metabolite_profile(t, doses, model)[last]
    tl = t[last]
    print(f"  metabolite/parent at steady state: AUC {auc_trapezoid(tl, met) / auc_trapezoid(tl, parent):.3f}  "
          f"Cmax {met.max() / parent.max():.3f} (mass; x{model.metabolite.molar_ratio:.4g} for molar)")


def _phenotype_ratios(drug, amount: float, tau: float) -> None:
    """Steady-state Cavg of the reported quantity per metaboliser phenotype, relative to normal."""
    enzymes = load_enzymes()
    base = steady_state_summary(amount, tau, drug.to_pk_model(70)).cavg_ss
    for name in drug.pathways or {}:
        labels = [p for p in getattr(enzymes.get(name), "genotype_activity", {}) if p != "normal"]
        if not labels:
            continue
        ratios = "  ".join(
            f"{p} {steady_state_summary(amount, tau, drug.to_pk_model(70, genotypes={name: p})).cavg_ss / base:.2f}"
            for p in labels)
        print(f"  {name} phenotype / normal (Cavg): {ratios}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("drug_id")
    ap.add_argument("--formulation")
    ap.add_argument("--dose", type=float)
    ap.add_argument("--interval", type=float)
    ap.add_argument("--set", action="append", default=[], metavar="FIELD=JSON")
    ap.add_argument("--ref", action="append", default=[], metavar="METRIC:DOSE:VALUE[:INTERVAL]")
    args = ap.parse_args()

    db = _apply_sets(args.drug_id, args.set) if args.set else DEFAULT_DB_PATH
    drug = get_drug(args.drug_id, db).with_formulation(args.formulation)
    model = drug.to_pk_model(70)
    unit = drug.concentration_unit
    dose = args.dose or drug.default_dose_mg
    tau = args.interval or drug.default_interval_h
    amount = drug.amount_from_dose_mg(dose)

    print(f"{drug.name} ({drug.id}{', ' + args.formulation if args.formulation else ''}), 70 kg")
    print(f"  t1/2 {model.t_half_h:.3g} h  V {model.vd_l:.4g} L  CL {model.clearance_l_per_h:.4g} L/h  "
          f"CL/F {model.clearance_l_per_h / model.f_eff:.4g} L/h  F {model.f_eff:.3g}  ka {model.ka_per_h:.3g}/h")
    if model.metabolite:
        m = model.metabolite
        print(f"  metabolite {m.name}: t1/2 {m.t_half_h:g} h  V {m.vd_l:.4g} L  fm {m.formation_fraction:g}  "
              f"reported={m.reported}")
    print(f"  terminal t1/2 of the reported quantity {model.terminal_half_life_h:.3g} h")
    t = np.linspace(0.0, 20.0 * model.terminal_half_life_h, 8001)
    auc = auc_trapezoid(t, concentration_profile(t, (DoseEvent(0.0, amount),), model))
    print(f"\nSingle {dose:g} mg: Cmax {cmax_single_dose(amount, model):.4g} {unit}  "
          f"Tmax {tmax_single_dose(model, amount):.3g} h  AUCinf {auc:.4g} {unit}*h")
    ss = steady_state_summary(amount, tau, model)
    print(f"Steady state {dose:g} mg q{tau:g}h: Cmin {ss.cmin_ss:.4g}  Cavg {ss.cavg_ss:.4g}  Cmax {ss.cmax_ss:.4g} "
          f"{unit}  fluctuation {ss.fluctuation_pct:.0f} %  window {drug.mec:g}-{drug.mtc:g}")

    if model.metabolite:
        _metabolite_ratios(model, amount, tau)
    _phenotype_ratios(drug, amount, tau)

    for item in args.ref:
        metric, d_mg, value, *iv = item.split(":")
        interval = float(iv[0]) if iv else None
        pred = predict(drug.id, metric, float(d_mg), interval, db_path=db)
        print(f"  candidate {metric} {float(d_mg):g} mg{f' q{interval:g}h' if interval else ''}: obs {float(value):.4g} "
              f"pred {pred:.4g} ratio {pred / float(value):.2f}")

    rows = [r for r in validate(db_path=db) if r.drug_id == drug.id]
    print(f"\nReferences ({len(rows)}):")
    for r in rows:
        mark = "ok  " if r.within else "MISS"
        iv = f" q{r.interval_h:g}h" if r.interval_h else ""
        print(f"  {mark} {r.metric:<10} {r.dose_mg:>6g} mg{iv:<6} obs {r.observed:<9.4g} pred {r.predicted:<9.4g} "
              f"ratio {r.ratio:.2f} (limit {r.max_fold:g})  {'verified' if r.verified else 'unverified'}")
    try:
        curves = [c for c in validate_curves(db_path=db) if c.drug_id == drug.id]
    except ValueError as exc:
        print(f"\nCurves: scoring stopped: {exc}")
        return
    print(f"\nCurves ({len(curves)}):")
    for c in curves:
        print(f"  {'ok  ' if c.within else 'MISS'} {c.curve_id}: AAFE {c.fit.aafe:.2f}  AFE {c.fit.afe:.2f}  "
              f"n={c.n_points}")


if __name__ == "__main__":
    main()
