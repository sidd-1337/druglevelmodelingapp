from __future__ import annotations

from typing import Literal
from pydantic import BaseModel, ConfigDict, Field


class RouteConfig(BaseModel):
    bioavailability: float
    ka_per_h: float


class DrugPreset(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: str
    name: str
    category: str = ""
    t_half_h: float
    t_half_range_h: tuple[float, float] | None = None
    vd_l_per_kg: float
    bioavailability: float
    ka_per_h: float
    mec: float
    mtc: float
    concentration_unit: str = "ug/mL"
    dose_conversion: float = 1.0
    default_dose_mg: float = 100.0
    default_interval_h: float = 12.0
    dose_basis: str = "mg"
    notes: str = ""
    routes: dict[str, dict[str, float]] | None = None
    formulations: dict[str, dict] | None = None
    cv_cl: float | None = None
    cv_f: float = 0.0
    fe_renal: float = 0.0
    pathways: dict[str, float] | None = None
    fu: float = 1.0
    km: float | None = None
    saturable_fraction: float = 0.0
    v1_fraction: float = 1.0
    q_l_per_h: float = 0.0
    tlag_h: float = 0.0
    inactivation_fraction: float = 0.0
    kinact_per_h: float = 0.0
    ki: float = 1.0
    kdeg_per_h: float = 0.0
    inactivated_enzyme: str = ""


class PatientIn(BaseModel):
    """Renal function: give crcl_ml_min directly, or age + serum creatinine (Cockcroft-Gault)."""

    crcl_ml_min: float | None = Field(default=None, gt=0)
    age_y: float | None = Field(default=None, gt=0, lt=120)
    female: bool = False
    serum_creatinine_umol_l: float | None = Field(default=None, gt=0)
    phenotypes: dict[str, str] | None = Field(
        default=None,
        description="Metaboliser phenotype per enzyme, e.g. {'CYP2D6': 'poor'}. An enzyme left out is "
                    "treated as unknown: clearance stays at the normal-metaboliser value but the "
                    "between-patient band widens to cover the phenotype mix.")


class CustomPKParams(BaseModel):
    t_half_h: float | None = None
    vd_l_per_kg: float | None = None
    bioavailability: float | None = None
    ka_per_h: float | None = None
    mec: float | None = None
    mtc: float | None = None
    concentration_unit: str | None = None


class MissedDoseRequest(BaseModel):
    dose_index: int = Field(ge=0, description="0-based index of the affected dose")
    delay_h: float | None = Field(default=None, ge=0, description="Delay in hours (None = skip dose completely)")


class SimulateRequest(BaseModel):
    drug_id: str = "ibuprofen"
    dose_mg: float = Field(gt=0, default=400.0)
    interval_h: float = Field(gt=0, default=8.0)
    duration_h: float = Field(gt=0, default=72.0)
    body_weight_kg: float = Field(gt=0, default=70.0)
    route: str = Field(default="oral", description="Route of administration: oral, iv_bolus, insufflated, smoked, inhaled")
    formulation: str | None = Field(default=None, description="Formulation id from the preset (e.g. 'pr')")
    patient: PatientIn | None = None
    custom_pk: CustomPKParams | None = None
    missed_dose: MissedDoseRequest | None = None


class ClinicalMetricsResponse(BaseModel):
    cmax: float
    tmax_h: float
    t_half_h: float
    clearance_l_per_h: float
    auc_0_24: float


class SteadyStateResponse(BaseModel):
    cmax_ss: float
    tmax_ss_h: float
    cmin_ss: float
    cavg_ss: float
    accumulation_index: float
    fluctuation_pct: float
    time_to_ss_h: float
    time_to_ss_range_h: tuple[float, float]


class TDMAssessmentResponse(BaseModel):
    status: str
    cmax: float
    t_cmax_h: float
    first_toxic_time_h: float | None = None
    time_above_mtc_h: float
    time_in_window_h: float
    time_below_mec_h: float
    trough: float


class SimulateResponse(BaseModel):
    drug_id: str
    drug_name: str
    route: str
    dose_mg: float
    interval_h: float
    duration_h: float
    unit: str
    time_grid: list[float]
    concentrations: list[float]
    zones: list[str]
    dose_times: list[float]
    metrics: ClinicalMetricsResponse
    steady_state: SteadyStateResponse
    tdm: TDMAssessmentResponse
    missed_concentrations: list[float] | None = None
    crcl_ml_min: float | None = None
    band_low: list[float] = Field(description="5th percentile of inter-individual exposure")
    band_high: list[float] = Field(description="95th percentile of inter-individual exposure")


class HarmReductionRequest(BaseModel):
    drug_id: str = "ritalin"
    dose_mg: float = Field(gt=0, default=10.0)
    body_weight_kg: float = Field(gt=0, default=70.0)
    route: str = Field(default="oral", description="Route: oral, iv_bolus, insufflated, smoked, inhaled")
    elapsed_h: float = Field(ge=0, default=2.0, description="Hours elapsed since first dose")
    formulation: str | None = None
    next_dose_mg: float | None = Field(default=None, gt=0, description="Planned next dose in mg")
    patient: PatientIn | None = None
    custom_pk: CustomPKParams | None = None


class PhaseTimingsResponse(BaseModel):
    onset_h: float
    peak_start_h: float
    tmax_h: float
    peak_end_h: float
    offset_h: float
    comedown_end_h: float


class RouteTimelineResponse(BaseModel):
    route: str
    route_name_en: str
    route_name_cs: str
    onset_minutes_range: tuple[int, int]
    tmax_h: float
    duration_hours_range: tuple[float, float]
    cmax: float
    bioavailability_pct: float
    description_en: str
    description_cs: str


class DoseStackingResponse(BaseModel):
    elapsed_h: float
    current_level: float
    single_cmax: float
    stacked_peak: float
    mtc: float
    mec: float
    unit: str
    risk_level: str
    risk_headline: str
    risk_explanation: str
    will_breach_mtc: bool
    safe_to_redose: bool
    time_above_mtc_h: float = 0.0


class HarmReductionResponse(BaseModel):
    drug_id: str
    drug_name: str
    route: str
    dose_mg: float
    next_dose_mg: float
    unit: str
    safe_redosing_interval_h: float
    phases: PhaseTimingsResponse
    stacking_assessment: DoseStackingResponse
    route_timelines: list[RouteTimelineResponse]
    guidelines: list[dict[str, str]]
    single_curve: list[dict[str, float]]
    stacked_curve: list[dict[str, float]]


class ObservationIn(BaseModel):
    time_h: float = Field(gt=0, description="Hours since the first dose of the schedule")
    concentration: float = Field(gt=0, description="Measured level in the drug's concentration unit")


class DoseIn(BaseModel):
    time_h: float = Field(ge=0, description="Hours since the start of the record")
    dose_mg: float = Field(gt=0)


class IndividualizeRequest(SimulateRequest):
    observations: list[ObservationIn] = Field(min_length=1, max_length=20)
    dose_history: list[DoseIn] | None = Field(
        default=None, max_length=1000,
        description="Doses actually taken; replaces dose_mg/interval_h/missed_dose for the fit when given")
    residual_cv: float = Field(default=0.2, gt=0, le=1.0, description="Assay + model error (proportional CV)")


class ObservationFitResponse(BaseModel):
    time_h: float
    observed: float
    population_predicted: float
    individual_predicted: float


class IndividualizeResponse(BaseModel):
    drug_id: str
    unit: str
    time_grid: list[float]
    population_concentrations: list[float]
    individual_concentrations: list[float]
    band_low: list[float]
    band_high: list[float]
    dose_times: list[float]
    observations: list[ObservationFitResponse]
    clearance_population_l_per_h: float
    clearance_individual_l_per_h: float
    vd_population_l: float
    vd_individual_l: float
    t_half_population_h: float
    t_half_individual_h: float
    steady_state: SteadyStateResponse
    tdm: TDMAssessmentResponse


class ValidationRowResponse(BaseModel):
    drug_id: str
    metric: str
    dose_mg: float
    interval_h: float | None = None
    published: float
    predicted: float
    ratio: float = Field(description="predicted / published")
    max_fold: float
    within_tolerance: bool
    source: str
    verified: bool
    note: str = ""


class ValidationGroupResponse(BaseModel):
    n: int
    within: int
    geometric_mean_fold_error: float


class ValidationResponse(BaseModel):
    summary: dict[str, ValidationGroupResponse]
    rows: list[ValidationRowResponse]


class CoMedicationIn(BaseModel):
    """A drug in the regimen; preset drugs are simulated, perpetrator names only act on others."""

    drug_id: str | None = Field(default=None, description="Preset id; omit for a perpetrator-only entry")
    perpetrator: str | None = Field(default=None, description="Exact perpetrator name from data/interactions.json")
    dose_mg: float | None = Field(default=None, gt=0)
    interval_h: float | None = Field(default=None, gt=0)
    route: str = "oral"
    formulation: str | None = None
    start_h: float | None = Field(default=None, ge=0, description=(
        "Perpetrator entries: hour it is first taken. Omit when it is already being taken at steady state."))
    stop_h: float | None = Field(default=None, gt=0, description="Perpetrator entries: hour it is stopped")


class PerpetratorResponse(BaseModel):
    name: str = Field(description="Exact perpetrator name to send as CoMedicationIn.perpetrator")
    name_cs: str
    timed: bool = Field(description="Has its own PK, so start_h / stop_h shape its effect")
    enzyme: str | None = None
    induces: bool = False
    victims: list[str] = Field(description="Preset ids it has a documented interaction with")


class AppliedInteractionResponse(BaseModel):
    perpetrator: str
    site: str
    factor: float | None = Field(description="Steady-state exposure factor applied; null if not applied as a factor")
    timed: bool = Field(description="Builds and fades with the perpetrator's own pharmacokinetics")
    via_enzyme: bool = Field(description="Applied through a shared enzyme pool instead of as a factor")
    skipped: str = Field(default="", description="'qualitative' or 'iv' when it changes nothing")
    pooled: str = Field(default="", description="Enzyme it acts on together with other added substances")


class PooledEnzymeResponse(BaseModel):
    enzyme: str
    multiplied: float = Field(description="The separate factors multiplied, at the end of the simulation")
    combined: float = Field(description="Combined in the enzyme instead (what the curve uses)")


class ClashResponse(BaseModel):
    effect: str
    participants: list[str] = Field(description="Drug ids, or 'substance:<name>' for an added perpetrator")
    names: list[str]
    overlap_h: float = Field(description="Time at least two of them are active together")
    first_h: float | None = None


class InteractionResponse(BaseModel):
    victim_id: str
    perpetrator: str
    perpetrator_cs: str
    site: str
    exposure_factor: float | None = None
    direction: str
    onset_h: float
    qualitative: bool
    source: str
    verified: bool
    note: str = ""
    estimated: bool = Field(default=False, description="Derived from the mechanism, not measured for this pair")
    basis: str = ""


class MultiDrugSeriesResponse(BaseModel):
    drug_id: str
    drug_name: str
    unit: str
    dose_mg: float
    interval_h: float
    concentrations: list[float]
    concentrations_alone: list[float]
    mec: float
    mtc: float
    exposure_factor: float = Field(description="Exposure over the last dosing interval, with co-medication / alone")
    t_half_h: float = Field(description="The drug's own terminal half-life; interactions scale exposure, not it")
    steady_state: SteadyStateResponse = Field(description="Steady state alone, levels scaled by exposure_factor")
    tdm: TDMAssessmentResponse
    interactions: list[InteractionResponse]
    applied: list[AppliedInteractionResponse] = Field(default_factory=list)
    drivers: list[str] = Field(default_factory=list, description="What makes its interaction effect time-varying")
    pooled: list[PooledEnzymeResponse] = Field(
        default_factory=list, description="Enzymes several added substances act on together for this drug")


class MultiDrugRequest(BaseModel):
    drugs: list[CoMedicationIn] = Field(min_length=1, max_length=8)
    duration_h: float = Field(gt=0, default=168.0)
    body_weight_kg: float = Field(gt=0, default=70.0)
    patient: PatientIn | None = None
    output_points: int | None = Field(default=None, ge=50, le=5000, description=(
        "Return curves resampled to this many points (simulation still runs on the fine grid); "
        "omit for the full grid"))


class MultiDrugResponse(BaseModel):
    time_grid: list[float]
    crcl_ml_min: float | None = None
    series: list[MultiDrugSeriesResponse]
    unmatched_perpetrators: list[str] = Field(
        default_factory=list, description="Named perpetrators with no interaction data for any simulated drug")
    enzyme_activity: dict[str, list[float]] = Field(
        default_factory=dict, description="Active enzyme as a fraction of baseline over time_grid")
    enzyme_drivers: dict[str, list[str]] = Field(default_factory=dict)
    clashes: list[ClashResponse] = Field(default_factory=list)
    pd_unassessed: list[str] = Field(default_factory=list, description="Participants not assessed for PD effects")
