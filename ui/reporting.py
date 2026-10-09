"""Exports: tabular CSV data and a human-readable dosing protocol (Markdown)."""
from __future__ import annotations

from datetime import date
from typing import Callable, Sequence

import numpy as np
import pandas as pd

from core.drug_database import Drug
from core.metrics import ZONE_LABELS, ClinicalMetrics, TherapeuticWindow
from core.pk_engine import DoseEvent, PKModel, SteadyState


def build_dataframe(
    t: np.ndarray,
    c: np.ndarray,
    window: TherapeuticWindow,
    unit: str,
    band: tuple[np.ndarray, np.ndarray] | None = None,
    reference: np.ndarray | None = None,
) -> pd.DataFrame:
    df = pd.DataFrame({"time_h": t, "time_day": t / 24.0, f"concentration_{unit}": c})
    if band is not None:
        df[f"conc_t_half_slow_{unit}"] = band[1]
        df[f"conc_t_half_fast_{unit}"] = band[0]
    if reference is not None:
        df[f"concentration_no_missed_dose_{unit}"] = reference
    df["zone"] = [ZONE_LABELS[int(z)] for z in window.classify(c)]
    return df.round(6)


def build_protocol(
    T: Callable[..., str],
    drug: Drug,
    model: PKModel,
    events: Sequence[DoseEvent],
    dose_mg: float,
    tau_h: float,
    days: int,
    weight_kg: float,
    window: TherapeuticWindow,
    metrics: ClinicalMetrics,
    ss: SteadyState,
    assessment_text: str,
    missed_text: str | None = None,
) -> str:
    unit = drug.concentration_unit
    shown = events[:12]
    times = ", ".join(f"{e.time_h:g}" for e in shown) + (" ..." if len(events) > len(shown) else "")
    lines = [
        f"# {T('protocol_title')}",
        f"_{T('protocol_generated', d=date.today().isoformat())}_",
        "",
        f"## {T('protocol_patient_regimen')}",
        f"- **{T('drug')}:** {drug.name} ({drug.category})" if drug.category else f"- **{T('drug')}:** {drug.name}",
        f"- **{T('body_weight')}:** {weight_kg:g} kg",
        f"- **{T('route')}:** {T('route_' + model.route.value)}",
        f"- **{T('regimen')}:** {dose_mg:g} mg ({drug.dose_basis}) {T('every')} {tau_h:g} h, {days} {T('days_unit')} "
        f"({len(events)} {T('doses_unit')})",
        f"- **{T('dose_times')} [h]:** {times}",
    ]
    if missed_text:
        lines.append(f"- **{T('missed_scenario')}:** {missed_text}")
    lines += [
        "",
        f"## {T('protocol_params')}",
        f"| {T('parameter')} | {T('value')} |",
        "|---|---|",
        f"| t½ | {model.t_half_h:g} h |",
        f"| k_e | {model.ke:.4f} 1/h |",
        f"| V_d | {model.vd_l:.2f} L |",
        f"| F | {model.f_eff:g} |",
        f"| k_a | {model.ka_per_h:g} 1/h |" if model.route.value == "oral" else f"| k_a | – |",
        f"| CL | {model.clearance_l_per_h:.3f} L/h |",
        f"| MEC / MTC | {window.mec:g} / {window.mtc:g} {unit} |",
        "",
        f"## {T('protocol_results')}",
        f"| {T('parameter')} | {T('value')} |",
        "|---|---|",
        f"| C_max (simulation) | {metrics.cmax:.4g} {unit} at {metrics.tmax_h:.2f} h |",
        f"| AUC 0–24 h | {metrics.auc_0_24:.4g} {unit}·h |",
        f"| C_ss,max / C_ss,min / C_ss,avg | {ss.cmax_ss:.4g} / {ss.cmin_ss:.4g} / {ss.cavg_ss:.4g} {unit} |",
        f"| {T('accumulation')} R | {ss.accumulation_index:.3f} |",
        f"| {T('fluctuation')} | {ss.fluctuation_pct:.1f} % |",
        f"| {T('time_to_ss')} | {ss.time_to_ss_h:.1f} h ({ss.time_to_ss_range_h[0]:.0f}–{ss.time_to_ss_range_h[1]:.0f} h) |",
        "",
        f"## {T('protocol_assessment')}",
        assessment_text,
        "",
        "---",
        f"_{T('disclaimer')}_",
        "",
    ]
    return "\n".join(lines)
