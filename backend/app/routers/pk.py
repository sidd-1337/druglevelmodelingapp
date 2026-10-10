from __future__ import annotations

import dataclasses
import math
import numpy as np
from fastapi import APIRouter, HTTPException

from backend.app.schemas.pk import (
    ClinicalMetricsResponse,
    DoseStackingResponse,
    DrugPreset,
    HarmReductionRequest,
    HarmReductionResponse,
    PhaseTimingsResponse,
    RouteTimelineResponse,
    SimulateRequest,
    SimulateResponse,
    SteadyStateResponse,
    TDMAssessmentResponse,
)
from core.drug_database import get_drug, load_drugs
from core.harm_reduction import (
    assess_dose_stacking,
    compute_safe_redosing_interval_h,
    estimate_phases,
    get_harm_reduction_guidelines,
    get_route_timeline,
)
from core.metrics import (
    TherapeuticWindow,
    assess_tdm,
    compute_metrics,
)
from core.pk_engine import (
    PKModel,
    Route,
    apply_missed_dose,
    concentration_profile,
    regular_schedule,
    single_dose_concentration,
    steady_state_summary,
    time_grid,
    tmax_single_dose,
)

router = APIRouter()


def _resolve_route(route_str: str) -> Route:
    clean = route_str.strip().lower()
    if clean in ("inhaled", "smoked"):
        return Route.SMOKED
    if clean in ("iv", "iv_bolus", "intravenous"):
        return Route.IV_BOLUS
    if clean in ("insufflated", "nasal", "snorted"):
        return Route.INSUFFLATED
    if clean in ("oral", "po", "peroral"):
        return Route.ORAL
    try:
        return Route(clean)
    except ValueError:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid route: {route_str!r}. Supported routes: {', '.join(r.value for r in Route)} or 'inhaled'.",
        )


@router.get("/presets", response_model=list[DrugPreset])
def list_presets() -> list[DrugPreset]:
    """Retrieve all built-in clinical drug models with pharmacokinetic parameters and route configs."""
    presets = load_drugs()
    return [
        DrugPreset(
            id=d.id,
            name=d.name,
            category=d.category,
            t_half_h=d.t_half_h,
            t_half_range_h=d.t_half_range_h,
            vd_l_per_kg=d.vd_l_per_kg,
            bioavailability=d.bioavailability,
            ka_per_h=d.ka_per_h,
            mec=d.mec,
            mtc=d.mtc,
            concentration_unit=d.concentration_unit,
            dose_conversion=d.dose_conversion,
            default_dose_mg=d.default_dose_mg,
            default_interval_h=d.default_interval_h,
            dose_basis=d.dose_basis,
            notes=d.notes,
            routes=d.routes,
        )
        for d in presets
    ]


@router.post("/simulate", response_model=SimulateResponse)
def simulate_regimen(req: SimulateRequest) -> SimulateResponse:
    """Run one-compartment pharmacokinetic simulation supporting all routes of administration."""
    try:
        drug = get_drug(req.drug_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Drug with id {req.drug_id!r} not found in presets.")

    route = _resolve_route(req.route)

    # Base model from drug
    model = drug.to_pk_model(req.body_weight_kg, route=route)

    # Apply custom overrides if provided
    if req.custom_pk is not None:
        cp = req.custom_pk
        t_half = cp.t_half_h if cp.t_half_h is not None else model.t_half_h
        vd = cp.vd_l_per_kg * req.body_weight_kg if cp.vd_l_per_kg is not None else model.vd_l
        bio = cp.bioavailability if cp.bioavailability is not None else model.bioavailability
        ka = cp.ka_per_h if cp.ka_per_h is not None else model.ka_per_h
        model = PKModel(
            t_half_h=t_half,
            vd_l=vd,
            bioavailability=bio,
            ka_per_h=ka,
            route=route,
        )

    mec = req.custom_pk.mec if req.custom_pk and req.custom_pk.mec is not None else drug.mec
    mtc = req.custom_pk.mtc if req.custom_pk and req.custom_pk.mtc is not None else drug.mtc
    unit = req.custom_pk.concentration_unit if req.custom_pk and req.custom_pk.concentration_unit else drug.concentration_unit
    window = TherapeuticWindow(mec=mec, mtc=mtc)

    amount = drug.amount_from_dose_mg(req.dose_mg)
    events = regular_schedule(amount, req.interval_h, req.duration_h)
    t = time_grid(req.duration_h)
    c = concentration_profile(t, events, model)

    # Classify zones
    zone_codes = window.classify(c)
    zone_names = [("toxic" if z == 2 else ("therapeutic" if z == 1 else "subtherapeutic")) for z in zone_codes]

    # Clinical metrics & steady state
    metrics = compute_metrics(t, c, model)
    ss = steady_state_summary(amount, req.interval_h, model)
    tdm = assess_tdm(t, c, window)

    missed_c_list = None
    if req.missed_dose is not None:
        md = req.missed_dose
        if 0 <= md.dose_index < len(events):
            alt_events = apply_missed_dose(events, md.dose_index, delay_h=md.delay_h)
            c_alt = concentration_profile(t, alt_events, model)
            missed_c_list = [round(float(v), 5) for v in c_alt]

    return SimulateResponse(
        drug_id=drug.id,
        drug_name=drug.name,
        route=route.value,
        dose_mg=req.dose_mg,
        interval_h=req.interval_h,
        duration_h=req.duration_h,
        unit=unit,
        time_grid=[round(float(x), 4) for x in t],
        concentrations=[round(float(y), 5) for y in c],
        zones=zone_names,
        dose_times=[round(float(e.time_h), 3) for e in events],
        metrics=ClinicalMetricsResponse(
            cmax=round(metrics.cmax, 4),
            tmax_h=round(metrics.tmax_h, 3),
            t_half_h=round(metrics.t_half_h, 3),
            clearance_l_per_h=round(metrics.clearance_l_per_h, 4),
            auc_0_24=round(metrics.auc_0_24, 4),
        ),
        steady_state=SteadyStateResponse(
            cmax_ss=round(ss.cmax_ss, 4),
            tmax_ss_h=round(ss.tmax_ss_h, 3),
            cmin_ss=round(ss.cmin_ss, 4),
            cavg_ss=round(ss.cavg_ss, 4),
            accumulation_index=round(ss.accumulation_index, 3),
            fluctuation_pct=round(ss.fluctuation_pct, 2),
            time_to_ss_h=round(ss.time_to_ss_h, 2),
            time_to_ss_range_h=(round(ss.time_to_ss_range_h[0], 2), round(ss.time_to_ss_range_h[1], 2)),
        ),
        tdm=TDMAssessmentResponse(
            status=tdm.status,
            cmax=round(tdm.cmax, 4),
            t_cmax_h=round(tdm.t_cmax_h, 3),
            first_toxic_time_h=round(tdm.first_toxic_time_h, 3) if tdm.first_toxic_time_h is not None else None,
            time_above_mtc_h=round(tdm.time_above_mtc_h, 3),
            time_in_window_h=round(tdm.time_in_window_h, 3),
            time_below_mec_h=round(tdm.time_below_mec_h, 3),
            trough=round(tdm.trough, 4),
        ),
        missed_concentrations=missed_c_list,
    )


@router.post("/harm-reduction", response_model=HarmReductionResponse)
def calculate_harm_reduction(req: HarmReductionRequest) -> HarmReductionResponse:
    """Calculate harm reduction indicators: safe redosing interval, stacking toxicity, onset/peak/offset timelines, and educational guidelines."""
    try:
        drug = get_drug(req.drug_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Drug with id {req.drug_id!r} not found in presets.")

    route = _resolve_route(req.route)
    model = drug.to_pk_model(req.body_weight_kg, route=route)

    # Custom PK overrides
    if req.custom_pk is not None:
        cp = req.custom_pk
        t_half = cp.t_half_h if cp.t_half_h is not None else model.t_half_h
        vd = cp.vd_l_per_kg * req.body_weight_kg if cp.vd_l_per_kg is not None else model.vd_l
        bio = cp.bioavailability if cp.bioavailability is not None else model.bioavailability
        ka = cp.ka_per_h if cp.ka_per_h is not None else model.ka_per_h
        model = PKModel(
            t_half_h=t_half,
            vd_l=vd,
            bioavailability=bio,
            ka_per_h=ka,
            route=route,
        )

    mec = req.custom_pk.mec if req.custom_pk and req.custom_pk.mec is not None else drug.mec
    mtc = req.custom_pk.mtc if req.custom_pk and req.custom_pk.mtc is not None else drug.mtc
    unit = req.custom_pk.concentration_unit if req.custom_pk and req.custom_pk.concentration_unit else drug.concentration_unit
    window = TherapeuticWindow(mec=mec, mtc=mtc)

    amount = drug.amount_from_dose_mg(req.dose_mg)
    next_dose_mg = req.next_dose_mg if req.next_dose_mg is not None else req.dose_mg
    redose_amount = drug.amount_from_dose_mg(next_dose_mg)

    # 1. Safe redosing interval
    safe_interval_h = compute_safe_redosing_interval_h(
        amount, model, window, redose_amount=redose_amount
    )

    # 2. Dose stacking assessment
    assessment = assess_dose_stacking(
        amount, req.elapsed_h, model, window, unit=unit, redose_amount=redose_amount
    )

    # 3. Subjective and PK phase estimations
    phases = estimate_phases(amount, model, window)

    # 4. Route timeline comparisons for all major routes
    route_timelines: list[RouteTimelineResponse] = []
    for r in [Route.ORAL, Route.IV_BOLUS, Route.INSUFFLATED, Route.SMOKED]:
        tl = get_route_timeline(r, amount, model, drug=drug)
        route_timelines.append(
            RouteTimelineResponse(
                route=tl.route.value,
                route_name_en=tl.route_name_en,
                route_name_cs=tl.route_name_cs,
                onset_minutes_range=tl.onset_minutes_range,
                tmax_h=round(tl.tmax_h, 3),
                duration_hours_range=tl.duration_hours_range,
                cmax=round(tl.cmax, 4),
                bioavailability_pct=tl.bioavailability_pct,
                description_en=tl.description_en,
                description_cs=tl.description_cs,
            )
        )

    # 5. Harm reduction guidelines
    guidelines = get_harm_reduction_guidelines(drug.id)

    # 6. Curves for visualization with fine resolution around absorption peaks
    horizon = max(24.0, req.elapsed_h + 4.0 * model.t_half_h)
    t_peak = tmax_single_dose(model)
    fine_near_0 = np.linspace(0.0, max(2.0, 4.0 * t_peak), 150)
    fine_near_redose = np.linspace(req.elapsed_h, req.elapsed_h + max(2.0, 4.0 * t_peak), 150)
    grid_full = np.linspace(0.0, horizon, 300)
    t_curve = np.unique(np.concatenate((fine_near_0, fine_near_redose, grid_full)))
    t_curve = t_curve[(t_curve >= 0.0) & (t_curve <= horizon)]
    t_curve.sort()

    c_single = single_dose_concentration(t_curve, amount, model)
    single_curve_data = [{"time_h": round(float(th), 3), "concentration": round(float(ch), 5)} for th, ch in zip(t_curve, c_single)]

    c_post_redose = single_dose_concentration(t_curve - req.elapsed_h, redose_amount, model)
    c_stacked = c_single + c_post_redose
    stacked_curve_data = [{"time_h": round(float(th), 3), "concentration": round(float(ch), 5)} for th, ch in zip(t_curve, c_stacked)]

    return HarmReductionResponse(
        drug_id=drug.id,
        drug_name=drug.name,
        route=route.value,
        dose_mg=req.dose_mg,
        next_dose_mg=next_dose_mg,
        unit=unit,
        safe_redosing_interval_h=round(safe_interval_h, 2),
        phases=PhaseTimingsResponse(
            onset_h=phases.onset_h,
            peak_start_h=phases.peak_start_h,
            tmax_h=phases.tmax_h,
            peak_end_h=phases.peak_end_h,
            offset_h=phases.offset_h,
            comedown_end_h=phases.comedown_end_h,
        ),
        stacking_assessment=DoseStackingResponse(
            elapsed_h=round(assessment.elapsed_h, 2),
            current_level=round(assessment.current_level, 4),
            single_cmax=round(assessment.single_cmax, 4),
            stacked_peak=round(assessment.stacked_peak, 4),
            mtc=assessment.mtc,
            mec=assessment.mec,
            unit=assessment.unit,
            risk_level=assessment.risk_level,
            risk_headline=assessment.risk_headline,
            risk_explanation=assessment.risk_explanation,
            will_breach_mtc=assessment.will_breach_mtc,
            safe_to_redose=assessment.safe_to_redose,
            time_above_mtc_h=assessment.time_above_mtc_h,
        ),
        route_timelines=route_timelines,
        guidelines=guidelines,
        single_curve=single_curve_data,
        stacked_curve=stacked_curve_data,
    )
