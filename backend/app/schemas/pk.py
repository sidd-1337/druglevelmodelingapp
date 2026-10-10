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


class HarmReductionRequest(BaseModel):
    drug_id: str = "ritalin"
    dose_mg: float = Field(gt=0, default=10.0)
    body_weight_kg: float = Field(gt=0, default=70.0)
    route: str = Field(default="oral", description="Route: oral, iv_bolus, insufflated, smoked, inhaled")
    elapsed_h: float = Field(ge=0, default=2.0, description="Hours elapsed since first dose")
    next_dose_mg: float | None = Field(default=None, gt=0, description="Planned next dose in mg")
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
