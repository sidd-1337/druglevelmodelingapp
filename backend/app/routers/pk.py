from __future__ import annotations

import dataclasses
import math
import numpy as np
from fastapi import APIRouter, HTTPException

from backend.app.schemas.pk import (
    ClinicalMetricsResponse,
    CustomPKParams,
    DoseStackingResponse,
    DrugPreset,
    HarmReductionRequest,
    HarmReductionResponse,
    IndividualizeRequest,
    IndividualizeResponse,
    ObservationFitResponse,
    PatientIn,
    PhaseTimingsResponse,
    RouteTimelineResponse,
    SimulateRequest,
    SimulateResponse,
    SteadyStateResponse,
    AppliedInteractionResponse,
    PooledEnzymeResponse,
    ClashResponse,
    InteractionResponse,
    MultiDrugRequest,
    MultiDrugResponse,
    MultiDrugSeriesResponse,
    PerpetratorResponse,
    TDMAssessmentResponse,
    ValidationGroupResponse,
    ValidationResponse,
    ValidationRowResponse,
)
from core.drug_database import Drug, get_drug, load_drugs
from core.population import cockcroft_gault
from core.validation import summarize, validate
from core.harm_reduction import (
    assess_dose_stacking,
    compute_safe_redosing_interval_h,
    estimate_phases,
    get_harm_reduction_guidelines,
    get_route_timeline,
)
from core.individualize import Observation, fit_individual, individual_band
from core.enzymes import load_perpetrators, perpetrator_pk
from core.interactions import interactions_for, load_interactions
from core.regimen import AddedSubstance, RegimenDrug, simulate_regimen as run_regimen
from core.metrics import (
    TDMAssessment,
    TherapeuticWindow,
    assess_tdm,
    compute_metrics,
)
from core.pk_engine import (
    DoseEvent,
    PKModel,
    Route,
    SteadyState,
    apply_missed_dose,
    concentration_profile,
    regular_schedule,
    single_dose_concentration,
    steady_state_summary,
    time_grid,
    tmax_single_dose,
)

router = APIRouter()


def _with_formulation(drug: Drug, formulation: str | None) -> Drug:
    try:
        return drug.with_formulation(formulation)
    except KeyError:
        raise HTTPException(status_code=400,
                            detail=f"Unknown formulation {formulation!r} for {drug.id!r}; "
                                   f"available: {drug.formulation_ids or 'none'}")


def _resolve_crcl(patient: PatientIn | None, body_weight_kg: float) -> float | None:
    if patient is None:
        return None
    if patient.crcl_ml_min is not None:
        return patient.crcl_ml_min
    if patient.age_y is not None and patient.serum_creatinine_umol_l is not None:
        return cockcroft_gault(patient.age_y, body_weight_kg, patient.serum_creatinine_umol_l, patient.female)
    return None


def _patient_model(
    drug: Drug, route: Route, body_weight_kg: float, crcl_ml_min: float | None,
    custom: CustomPKParams | None, genotypes: dict[str, str] | None = None,
) -> PKModel:
    """Preset model for this patient, with explicit custom overrides applied last."""
    try:
        model = drug.to_pk_model(body_weight_kg, route=route, crcl_ml_min=crcl_ml_min, genotypes=genotypes)
    except KeyError as exc:
        raise HTTPException(status_code=400, detail=f"Unknown enzyme or phenotype: {exc}")
    if custom is None:
        return model
    changes = {}
    if custom.t_half_h is not None:
        changes["t_half_h"] = custom.t_half_h
    if custom.vd_l_per_kg is not None:
        changes["vd_l"] = custom.vd_l_per_kg * body_weight_kg
    if custom.bioavailability is not None:
        changes["bioavailability"] = custom.bioavailability
    if custom.ka_per_h is not None:
        changes["ka_per_h"] = custom.ka_per_h
    return dataclasses.replace(model, **changes)


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


def _steady_state_response(ss: SteadyState) -> SteadyStateResponse:
    return SteadyStateResponse(
        cmax_ss=round(ss.cmax_ss, 4),
        tmax_ss_h=round(ss.tmax_ss_h, 3),
        cmin_ss=round(ss.cmin_ss, 4),
        cavg_ss=round(ss.cavg_ss, 4),
        accumulation_index=round(ss.accumulation_index, 3),
        fluctuation_pct=round(ss.fluctuation_pct, 2),
        time_to_ss_h=round(ss.time_to_ss_h, 2),
        time_to_ss_range_h=(round(ss.time_to_ss_range_h[0], 2), round(ss.time_to_ss_range_h[1], 2)),
    )


def _tdm_response(tdm: TDMAssessment) -> TDMAssessmentResponse:
    return TDMAssessmentResponse(
        status=tdm.status,
        cmax=round(tdm.cmax, 4),
        t_cmax_h=round(tdm.t_cmax_h, 3),
        first_toxic_time_h=round(tdm.first_toxic_time_h, 3) if tdm.first_toxic_time_h is not None else None,
        time_above_mtc_h=round(tdm.time_above_mtc_h, 3),
        time_in_window_h=round(tdm.time_in_window_h, 3),
        time_below_mec_h=round(tdm.time_below_mec_h, 3),
        trough=round(tdm.trough, 4),
    )


def _round_list(arr, ndigits: int = 5) -> list[float]:
    return [round(float(v), ndigits) for v in arr]


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
            formulations=d.formulations,
            cv_cl=d.cv_cl,
            cv_f=d.cv_f,
            fe_renal=d.fe_renal,
            pathways=d.pathways,
            fu=d.fu,
            km=d.km,
            saturable_fraction=d.saturable_fraction,
            v1_fraction=d.v1_fraction,
            q_l_per_h=d.q_l_per_h,
            tlag_h=d.tlag_h,
            inactivation_fraction=d.inactivation_fraction,
            kinact_per_h=d.kinact_per_h,
            ki=d.ki,
            kdeg_per_h=d.kdeg_per_h,
            inactivated_enzyme=d.inactivated_enzyme,
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

    drug = _with_formulation(drug, req.formulation)
    route = _resolve_route(req.route)
    crcl = _resolve_crcl(req.patient, req.body_weight_kg)
    genotypes = req.patient.phenotypes if req.patient else None
    model = _patient_model(drug, route, req.body_weight_kg, crcl, req.custom_pk, genotypes)

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

    band_models = drug.exposure_band(model, genotypes)
    band = [concentration_profile(t, events, m) for m in band_models]

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
        steady_state=_steady_state_response(ss),
        tdm=_tdm_response(tdm),
        missed_concentrations=missed_c_list,
        crcl_ml_min=round(crcl, 1) if crcl is not None else None,
        band_low=[round(float(v), 5) for v in np.minimum(*band)],
        band_high=[round(float(v), 5) for v in np.maximum(*band)],
    )


@router.post("/harm-reduction", response_model=HarmReductionResponse)
def calculate_harm_reduction(req: HarmReductionRequest) -> HarmReductionResponse:
    """Calculate harm reduction indicators: safe redosing interval, stacking toxicity, onset/peak/offset timelines, and educational guidelines."""
    try:
        drug = get_drug(req.drug_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Drug with id {req.drug_id!r} not found in presets.")

    drug = _with_formulation(drug, req.formulation)
    route = _resolve_route(req.route)
    crcl = _resolve_crcl(req.patient, req.body_weight_kg)
    model = _patient_model(drug, route, req.body_weight_kg, crcl, req.custom_pk)

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
    horizon = max(24.0, req.elapsed_h + 4.0 * model.terminal_half_life_h)
    t_peak = tmax_single_dose(model, amount)
    fine_near_0 = np.linspace(0.0, max(2.0, 4.0 * t_peak), 150)
    fine_near_redose = np.linspace(req.elapsed_h, req.elapsed_h + max(2.0, 4.0 * t_peak), 150)
    grid_full = np.linspace(0.0, horizon, 300)
    t_curve = np.unique(np.concatenate((fine_near_0, fine_near_redose, grid_full)))
    t_curve = t_curve[(t_curve >= 0.0) & (t_curve <= horizon)]
    t_curve.sort()

    c_single = single_dose_concentration(t_curve, amount, model)
    single_curve_data = [{"time_h": round(float(th), 3), "concentration": round(float(ch), 5)} for th, ch in zip(t_curve, c_single)]

    stack_events = [DoseEvent(0.0, amount), DoseEvent(req.elapsed_h, redose_amount)]
    c_stacked = concentration_profile(t_curve, stack_events, model)
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


@router.post("/individualize", response_model=IndividualizeResponse)
def individualize_regimen(req: IndividualizeRequest) -> IndividualizeResponse:
    """Fit the patient's clearance and volume to measured levels (Bayesian MAP, population prior)."""
    try:
        drug = get_drug(req.drug_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Drug with id {req.drug_id!r} not found in presets.")

    drug = _with_formulation(drug, req.formulation)
    route = _resolve_route(req.route)
    crcl = _resolve_crcl(req.patient, req.body_weight_kg)
    pop = _patient_model(drug, route, req.body_weight_kg, crcl, req.custom_pk)
    mec = req.custom_pk.mec if req.custom_pk and req.custom_pk.mec is not None else drug.mec
    mtc = req.custom_pk.mtc if req.custom_pk and req.custom_pk.mtc is not None else drug.mtc
    unit = req.custom_pk.concentration_unit if req.custom_pk and req.custom_pk.concentration_unit else drug.concentration_unit
    window = TherapeuticWindow(mec=mec, mtc=mtc)

    amount = drug.amount_from_dose_mg(req.dose_mg)
    if req.dose_history:
        events = sorted((DoseEvent(d.time_h, drug.amount_from_dose_mg(d.dose_mg)) for d in req.dose_history),
                        key=lambda e: e.time_h)
    else:
        events = regular_schedule(amount, req.interval_h, req.duration_h)
        if req.missed_dose is not None and 0 <= req.missed_dose.dose_index < len(events):
            events = apply_missed_dose(events, req.missed_dose.dose_index, delay_h=req.missed_dose.delay_h)
    observations = [Observation(o.time_h, o.concentration) for o in req.observations]
    try:
        fit = fit_individual(pop, events, observations, drug.cv_cl if drug.cv_cl is not None else 0.3,
                             drug.cv_f, residual_cv=req.residual_cv)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))

    t = time_grid(max(req.duration_h, max(o.time_h for o in observations), events[-1].time_h))
    c_pop = concentration_profile(t, events, pop)
    c_ind = concentration_profile(t, events, fit.model)
    lo, hi = individual_band(fit, t, events)
    ss = steady_state_summary(amount, req.interval_h, fit.model)
    return IndividualizeResponse(
        drug_id=drug.id,
        unit=unit,
        time_grid=_round_list(t, 4),
        population_concentrations=_round_list(c_pop),
        individual_concentrations=_round_list(c_ind),
        band_low=_round_list(lo),
        band_high=_round_list(hi),
        dose_times=[round(float(e.time_h), 3) for e in events],
        observations=[
            ObservationFitResponse(time_h=o.time_h, observed=o.concentration,
                                   population_predicted=round(p, 5), individual_predicted=round(i, 5))
            for o, p, i in zip(fit.observations, fit.population_predicted, fit.predicted)
        ],
        clearance_population_l_per_h=round(pop.clearance_l_per_h, 4),
        clearance_individual_l_per_h=round(fit.model.clearance_l_per_h, 4),
        vd_population_l=round(pop.vd_l, 3),
        vd_individual_l=round(fit.model.vd_l, 3),
        t_half_population_h=round(pop.terminal_half_life_h, 3),
        t_half_individual_h=round(fit.model.terminal_half_life_h, 3),
        steady_state=_steady_state_response(ss),
        tdm=_tdm_response(assess_tdm(t, c_ind, window)),
    )


@router.get("/validation", response_model=ValidationResponse)
def validation_report() -> ValidationResponse:
    """Preset predictions versus published PK data (data/reference_pk.json)."""
    rows = validate()
    return ValidationResponse(
        summary={k: ValidationGroupResponse(n=g.n, within=g.within,
                                            geometric_mean_fold_error=round(g.geometric_mean_fold_error, 3))
                 for k, g in summarize(rows).items()},
        rows=[ValidationRowResponse(
            drug_id=r.drug_id, metric=r.metric, dose_mg=r.dose_mg, interval_h=r.interval_h,
            published=r.observed, predicted=round(r.predicted, 6), ratio=round(r.ratio, 3), max_fold=r.max_fold,
            within_tolerance=r.within, source=r.source, verified=r.verified, note=r.note) for r in rows],
    )


def _interaction_response(i) -> InteractionResponse:
    return InteractionResponse(
        victim_id=i.victim_id, perpetrator=i.perpetrator, perpetrator_cs=i.perpetrator_cs, site=i.site,
        exposure_factor=i.exposure_factor, direction=i.direction, onset_h=i.onset_h, qualitative=i.qualitative,
        source=i.source, verified=i.verified, note=i.note, estimated=i.estimated, basis=i.basis)


@router.get("/interactions/{drug_id}", response_model=list[InteractionResponse])
def drug_interactions(drug_id: str) -> list[InteractionResponse]:
    """Known quantified interactions where this drug is the victim."""
    try:
        get_drug(drug_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Drug with id {drug_id!r} not found in presets.")
    return [_interaction_response(i) for i in interactions_for(drug_id)]


@router.get("/perpetrators", response_model=list[PerpetratorResponse])
def list_perpetrators() -> list[PerpetratorResponse]:
    """Substances that can be added to a /multi-drug regimen without being simulated themselves.

    Every perpetrator named in data/interactions.json plus the mechanistic ones in
    data/enzymes.json, minus the preset drugs (those are added by drug_id).
    """
    presets = {d.name for d in load_drugs()}
    names: dict[str, tuple[str, set[str]]] = {}
    for i in load_interactions():
        cs, victims = names.get(i.perpetrator, (i.perpetrator_cs, set()))
        names[i.perpetrator] = (cs or i.perpetrator_cs, victims | {i.victim_id})
    for p in load_perpetrators():
        names.setdefault(p.name, (p.name_cs, set()))
    out = []
    for name, (cs, victims) in sorted(names.items()):
        if name in presets:
            continue
        p = perpetrator_pk(name)
        out.append(PerpetratorResponse(
            name=name, name_cs=cs or name, timed=bool(p and p.t_half_h), enzyme=p.enzyme if p else None,
            induces=bool(p and p.induces), victims=sorted(victims)))
    return out


@router.post("/multi-drug", response_model=MultiDrugResponse)
def simulate_multi_drug(req: MultiDrugRequest) -> MultiDrugResponse:
    """Simulate several drugs on one time grid, with their interactions over time and PD clashes.

    Same model as the Streamlit co-medication tab (core.regimen): a shared enzyme a
    drug destroys is simulated as a pool, an added perpetrator with known PK builds
    and fades from start_h / stop_h (omit start_h when it is already being taken),
    and drugs or perpetrators sharing a PD effect class are flagged while both are
    active.
    """
    picked, added = [], []
    for item in req.drugs:
        if item.drug_id is None:
            if not item.perpetrator:
                raise HTTPException(status_code=422, detail="each entry needs drug_id or perpetrator")
            if item.stop_h is not None and item.start_h is not None and item.stop_h <= item.start_h:
                raise HTTPException(status_code=422, detail=f"{item.perpetrator}: stop_h must be after start_h")
            added.append(AddedSubstance(item.perpetrator, item.start_h, item.stop_h))
            continue
        try:
            picked.append((item, _with_formulation(get_drug(item.drug_id), item.formulation)))
        except KeyError:
            raise HTTPException(status_code=404, detail=f"Drug with id {item.drug_id!r} not found in presets.")
    if not picked:
        raise HTTPException(status_code=422, detail="at least one entry must name a preset drug_id")

    crcl = _resolve_crcl(req.patient, req.body_weight_kg)
    t = time_grid(req.duration_h)
    entries = []
    for item, drug in picked:
        model = drug.to_pk_model(req.body_weight_kg, route=_resolve_route(item.route), crcl_ml_min=crcl)
        dose_mg = item.dose_mg if item.dose_mg is not None else drug.default_dose_mg
        interval_h = item.interval_h if item.interval_h is not None else drug.default_interval_h
        events = tuple(regular_schedule(drug.amount_from_dose_mg(dose_mg), interval_h, req.duration_h))
        entries.append(RegimenDrug(drug, model, events, interval_h, drug.enzyme_shares(req.body_weight_kg, crcl)))
    res = run_regimen(t, entries, added)
    out_t = t if req.output_points is None else np.linspace(0.0, float(t[-1]), req.output_points)

    def resample(a) -> list[float]:
        return _round_list(a if req.output_points is None else np.interp(out_t, t, a))

    series = []
    for (item, drug), r in zip(picked, res.drugs):
        dose_mg = item.dose_mg if item.dose_mg is not None else drug.default_dose_mg
        amount, f = drug.amount_from_dose_mg(dose_mg), r.exposure_factor
        alone_ss = steady_state_summary(amount, r.entry.interval_h, r.entry.model)
        ss = dataclasses.replace(alone_ss, cmax_ss=alone_ss.cmax_ss * f, cmin_ss=alone_ss.cmin_ss * f,
                                 cavg_ss=alone_ss.cavg_ss * f)
        window = TherapeuticWindow(mec=drug.mec, mtc=drug.mtc)
        assess_from = min(ss.time_to_ss_h, max(0.0, req.duration_h - r.entry.interval_h))
        series.append(MultiDrugSeriesResponse(
            drug_id=drug.id, drug_name=drug.name, unit=drug.concentration_unit, dose_mg=dose_mg,
            interval_h=r.entry.interval_h, concentrations=resample(r.concentrations),
            concentrations_alone=resample(r.alone), mec=drug.mec, mtc=drug.mtc,
            exposure_factor=round(f, 4), t_half_h=round(r.entry.model.terminal_half_life_h, 3),
            steady_state=_steady_state_response(ss),
            tdm=_tdm_response(assess_tdm(t, r.concentrations, window, assess_from)),
            interactions=[_interaction_response(a.interaction) for a in r.applied],
            applied=[AppliedInteractionResponse(
                perpetrator=a.interaction.perpetrator, site=a.interaction.site,
                factor=round(a.factor, 4) if a.factor is not None else None, timed=a.timed,
                via_enzyme=a.via_enzyme, skipped=a.skipped, pooled=a.pooled) for a in r.applied],
            drivers=r.drivers,
            pooled=[PooledEnzymeResponse(enzyme=e, multiplied=round(m, 4), combined=round(c, 4))
                    for e, (m, c) in r.pooled.items()],
        ))
    return MultiDrugResponse(
        time_grid=_round_list(out_t, 4), crcl_ml_min=round(crcl, 1) if crcl is not None else None, series=series,
        unmatched_perpetrators=sorted(res.unmatched),
        enzyme_activity={e: resample(a) for e, a in res.enzyme_activity.items()},
        enzyme_drivers=res.enzyme_drivers,
        clashes=[ClashResponse(effect=c.effect, participants=list(c.drug_ids),
                               names=[res.labels[k][0] for k in c.drug_ids], overlap_h=round(c.overlap_h, 2),
                               first_h=round(c.first_h, 2) if c.first_h is not None else None)
                 for c in res.clashes],
        pd_unassessed=res.unassessed,
    )
