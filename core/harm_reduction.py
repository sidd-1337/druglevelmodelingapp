"""Harm reduction pharmacokinetics module.

Provides:
- Safe redosing interval calculations.
- Dose stacking warnings (predicting cumulative peak and toxic window breaches).
- Route of administration onset/peak/offset timelines (Oral, IV Bolus, Insufflated, Smoked).
- Practical harm reduction advice and risk mitigation protocols.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from core.metrics import TherapeuticWindow
from core.pk_engine import (
    PKModel,
    Route,
    cmax_single_dose,
    single_dose_concentration,
    tmax_single_dose,
)

if TYPE_CHECKING:
    from core.drug_database import Drug


@dataclass(frozen=True)
class DoseStackingAssessment:
    """Assessment of risk when taking an additional dose at elapsed time t."""

    elapsed_h: float
    current_level: float
    single_cmax: float
    stacked_peak: float
    mtc: float
    mec: float
    unit: str
    risk_level: Literal["safe", "moderate", "high", "critical_toxic"]
    risk_headline: str
    risk_explanation: str
    will_breach_mtc: bool
    safe_to_redose: bool
    time_above_mtc_h: float = 0.0


@dataclass(frozen=True)
class RouteTimeline:
    """Clinical and harm-reduction timeline for a specific route of administration."""

    route: Route
    route_name_en: str
    route_name_cs: str
    onset_minutes_range: tuple[int, int]
    tmax_h: float
    duration_hours_range: tuple[float, float]
    cmax: float
    bioavailability_pct: float
    description_en: str
    description_cs: str


def compute_safe_redosing_interval_h(
    amount: float,
    model: PKModel,
    window: TherapeuticWindow,
    safety_margin_factor: float = 0.85,
    redose_amount: float | None = None,
) -> float:
    """Compute the minimum time in hours to wait before a second dose is safe.

    A redose is deemed pharmacologically safe when the residual blood level plus
    the peak of the new dose does not exceed the safety ceiling (safety_margin_factor * MTC),
    is below the minimum effective concentration (MEC) so the primary dose is cleared,
    and is past the peak of the primary dose.
    """
    new_dose = amount if redose_amount is None else redose_amount
    peak = cmax_single_dose(new_dose, model)
    ceiling = safety_margin_factor * window.mtc
    max_residual = max(0.0, ceiling - peak)

    # If even the new single dose exceeds the ceiling, recommend at least 5 half-lives
    if max_residual <= 0:
        return max(5.0 * model.t_half_h, 6.0)

    # Target residual: must satisfy residual + peak <= ceiling AND residual <= MEC
    # (if primary peak >= MEC; if primary peak < MEC, dose was sub-therapeutic so ceiling governs)
    primary_peak = cmax_single_dose(amount, model)
    target_residual = min(max_residual, window.mec) if primary_peak >= window.mec else max_residual

    # Search on a time grid for residual <= target_residual after tmax
    t_peak = tmax_single_dose(model)
    horizon = max(48.0, 10.0 * model.t_half_h)
    t_test = np.linspace(t_peak, horizon, 3000)
    c_test = single_dose_concentration(t_test, amount, model)

    below = np.where((t_test >= t_peak) & (c_test <= target_residual))[0]
    if len(below) > 0:
        idx = below[0]
        if idx > 0:
            t0, t1 = t_test[idx - 1], t_test[idx]
            c0, c1 = c_test[idx - 1], c_test[idx]
            if c0 != c1:
                safe_time = float(t0 + (target_residual - c0) / (c1 - c0) * (t1 - t0))
            else:
                safe_time = float(t_test[idx])
        else:
            safe_time = float(t_test[idx])
    else:
        safe_time = 5.0 * model.t_half_h

    # Ensure it's at least past the active peak window
    min_bound = 0.5 if model.route == Route.IV_BOLUS else t_peak + 0.5
    return max(safe_time, min_bound)


def assess_dose_stacking(
    amount: float,
    elapsed_h: float,
    model: PKModel,
    window: TherapeuticWindow,
    unit: str = "ug/mL",
    redose_amount: float | None = None,
) -> DoseStackingAssessment:
    """Evaluate what happens if a dose is redosed at `elapsed_h` hours after the first dose."""
    new_dose = amount if redose_amount is None else redose_amount
    current_level = float(single_dose_concentration(elapsed_h, amount, model))
    single_cmax = cmax_single_dose(new_dose, model)

    # Multi-resolution sampling: dense near absorption peaks to resolve fast routes (smoked, insufflated)
    t_peak = tmax_single_dose(model)
    t_near = np.linspace(0.0, max(2.0, 4.0 * t_peak), 400)
    t_far = np.linspace(0.0, max(24.0, 6.0 * model.t_half_h), 800)
    t_post = np.unique(np.concatenate(([0.0, t_peak], t_near, t_far)))
    t_post.sort()

    c_residual = single_dose_concentration(elapsed_h + t_post, amount, model)
    c_new = single_dose_concentration(t_post, new_dose, model)
    stacked_curve = c_residual + c_new
    stacked_peak = float(np.max(stacked_curve))

    mtc = window.mtc
    mec = window.mec
    will_breach_mtc = stacked_peak >= mtc

    time_above_mtc_h = 0.0
    if will_breach_mtc:
        above_mask = stacked_curve >= mtc
        dt = np.diff(t_post)
        time_above_mtc_h = float(np.sum(dt[above_mask[:-1]]))

    if will_breach_mtc:
        risk_level = "critical_toxic"
        headline = "Kritické riziko předávkování / toxicity (překročení MTC)"
        explanation = (
            f"Při podání další dávky v čase {elapsed_h:.1f} h stoupne kumulativní hladina na "
            f"{stacked_peak:.3g} {unit}, což překračuje toxickou hranici MTC ({mtc:.3g} {unit})! "
            f"Předpokládaná doba v toxickém pásmu: {time_above_mtc_h:.1f} h. "
            f"Předchozí dávka je stále aktivní (aktuální zbytková hladina: {current_level:.3g} {unit})."
        )
        safe_to_redose = False
    elif stacked_peak > 0.85 * mtc:
        risk_level = "high"
        headline = "Vysoké riziko – hladina se nebezpečně blíží toxickému pásmu"
        explanation = (
            f"Kumulativní vrchol {stacked_peak:.3g} {unit} dosahuje přes 85 % MTC ({mtc:.3g} {unit}). "
            f"Vlivem biologické variability hrozí nežádoucí účinky a předávkování."
        )
        safe_to_redose = False
    elif current_level > mec:
        risk_level = "moderate"
        headline = "Střední riziko vrstvení – látka stále plně účinkuje"
        explanation = (
            f"Zbytková hladina ({current_level:.3g} {unit}) je stále nad minimální účinnou mezí MEC ({mec:.3g} {unit}). "
            f"Časné podání zesílí účinky a prodlouží celkovou dobu působení (akumulace)."
        )
        safe_to_redose = False
    else:
        risk_level = "safe"
        headline = "Nízké riziko kumulace – předchozí dávka je z větší části odbourána"
        explanation = (
            f"Zbytková hladina klesla na {current_level:.3g} {unit} (pod MEC). "
            f"Předpokládaný nový vrchol {stacked_peak:.3g} {unit} zůstane bezpečně v terapeutickém okně."
        )
        safe_to_redose = True

    return DoseStackingAssessment(
        elapsed_h=elapsed_h,
        current_level=current_level,
        single_cmax=single_cmax,
        stacked_peak=stacked_peak,
        mtc=mtc,
        mec=mec,
        unit=unit,
        risk_level=risk_level,
        risk_headline=headline,
        risk_explanation=explanation,
        will_breach_mtc=will_breach_mtc,
        safe_to_redose=safe_to_redose,
        time_above_mtc_h=round(time_above_mtc_h, 2),
    )


def get_route_timeline(
    route: Route,
    amount: float,
    model: PKModel,
    drug: Drug | None = None,
) -> RouteTimeline:
    """Generate kinetic timeline data for a route of administration."""
    if drug is not None:
        bw = model.vd_l / drug.vd_l_per_kg if drug.vd_l_per_kg > 0 else 70.0
        route_base = drug.to_pk_model(body_weight_kg=bw, route=route)
        route_model = PKModel(
            t_half_h=model.t_half_h,
            vd_l=model.vd_l,
            bioavailability=route_base.bioavailability,
            ka_per_h=route_base.ka_per_h,
            route=route,
        )
    elif model.route == route:
        route_model = model
    else:
        ka = model.ka_per_h
        bio = model.bioavailability
        if route == Route.IV_BOLUS:
            ka = 1.0
            bio = 1.0
        elif route == Route.INSUFFLATED:
            ka = max(model.ka_per_h * 2.5, 3.5) if model.route == Route.ORAL else model.ka_per_h
            bio = min(1.0, max(model.bioavailability, 0.70)) if model.route == Route.ORAL else model.bioavailability
        elif route in (Route.SMOKED, getattr(Route, "INHALED", Route.SMOKED)):
            ka = max(model.ka_per_h * 5.0, 8.0) if model.route == Route.ORAL else model.ka_per_h
            bio = min(1.0, max(model.bioavailability * 0.9, 0.60)) if model.route == Route.ORAL else model.bioavailability
        elif route == Route.ORAL:
            if model.route in (Route.INSUFFLATED, Route.SMOKED):
                ka = max(1.0, model.ka_per_h / 2.5)
                bio = min(1.0, model.bioavailability * 0.8)

        route_model = PKModel(
            t_half_h=model.t_half_h,
            vd_l=model.vd_l,
            bioavailability=bio,
            ka_per_h=ka,
            route=route,
        )

    tmax = tmax_single_dose(route_model)
    cmax = cmax_single_dose(amount, route_model)
    t_half = model.t_half_h

    if route == Route.ORAL:
        onset = (20, 45)
        duration = (round(2.5 * t_half, 1), round(4.5 * t_half, 1))
        en_desc = "First-pass metabolism; gradual absorption and sustained systemic exposure."
        cs_desc = "První průchod játry (first-pass); pozvolný nástup a déle trvající expozice."
    elif route == Route.IV_BOLUS:
        onset = (0, 1)
        duration = (round(2.0 * t_half, 1), round(3.5 * t_half, 1))
        en_desc = "100% bioavailability, instant systemic concentration, no absorption delay."
        cs_desc = "100% biologická dostupnost, okamžitý nástup do krevního oběhu bez zpoždění."
    elif route == Route.INSUFFLATED:
        onset = (5, 15)
        duration = (round(2.0 * t_half, 1), round(3.8 * t_half, 1))
        en_desc = "Rapid absorption through nasal mucous membranes, bypasses primary gut transit."
        cs_desc = "Rychlá absorpce nosní sliznicí, obchází žaludeční průchod, rychlejší nástup maxima."
    elif route == Route.SMOKED:
        onset = (1, 5)
        duration = (round(1.5 * t_half, 1), round(3.2 * t_half, 1))
        en_desc = "Ultra-rapid alveolar uptake into pulmonary capillary bed, peak within minutes."
        cs_desc = "Extrémně rychlý přestup plicními sklípky do krve a mozku, maximum během minut."
    else:
        onset = (15, 30)
        duration = (round(2.0 * t_half, 1), round(4.0 * t_half, 1))
        en_desc = "Extravascular systemic absorption."
        cs_desc = "Extravaskulární systémová absorpce."

    names = {
        Route.ORAL: ("Oral", "Perorální (spolknutí)"),
        Route.IV_BOLUS: ("IV Bolus", "Intravenózní (injekce)"),
        Route.INSUFFLATED: ("Insufflated", "Nazální (šňupání)"),
        Route.SMOKED: ("Smoked / Inhaled", "Inhalační (kouření/vaping)"),
    }
    name_en, name_cs = names.get(route, (route.value, route.value))

    return RouteTimeline(
        route=route,
        route_name_en=name_en,
        route_name_cs=name_cs,
        onset_minutes_range=onset,
        tmax_h=tmax,
        duration_hours_range=duration,
        cmax=cmax,
        bioavailability_pct=round(route_model.f_eff * 100, 1),
        description_en=en_desc,
        description_cs=cs_desc,
    )


@dataclass(frozen=True)
class PhaseTimings:
    """Pharmacokinetic and subjective phase timing estimations (in hours)."""

    onset_h: float
    peak_start_h: float
    tmax_h: float
    peak_end_h: float
    offset_h: float
    comedown_end_h: float


def estimate_phases(
    amount: float,
    model: PKModel,
    window: TherapeuticWindow | None = None,
) -> PhaseTimings:
    """Estimate onset, peak plateau, offset, and comedown phases with monotonic progression."""
    t_half = model.t_half_h
    t_peak = tmax_single_dose(model)
    horizon = max(24.0, 8.0 * t_half)

    # Dense time grid ensuring peak and phase boundaries are accurately captured
    t_dense = np.unique(np.concatenate((
        [0.0, t_peak],
        np.linspace(0.0, max(1.0, 3.0 * t_peak), 500),
        np.linspace(0.0, horizon, 3000),
    )))
    t_dense.sort()
    c = single_dose_concentration(t_dense, amount, model)
    cmax = float(np.max(c))

    if cmax <= 0:
        return PhaseTimings(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    if model.route == Route.IV_BOLUS:
        onset_h = 0.0
        peak_start_h = 0.0
        tmax_h = 0.0
        after_peak = t_dense >= 0.0
    else:
        tmax_h = t_peak
        # Onset: ~20% of Cmax
        onset_idx = np.where(c >= 0.20 * cmax)[0]
        onset_h = float(t_dense[onset_idx[0]]) if len(onset_idx) > 0 else 0.1
        # Peak start: ~80% of Cmax rising
        rising_80 = np.where((t_dense <= t_peak) & (c >= 0.80 * cmax))[0]
        peak_start_h = float(t_dense[rising_80[0]]) if len(rising_80) > 0 else max(onset_h, 0.7 * t_peak)
        after_peak = t_dense >= t_peak

    # Peak end: ~80% of Cmax falling
    falling_80 = np.where(after_peak & (c <= 0.80 * cmax))[0]
    peak_end_h = float(t_dense[falling_80[0]]) if len(falling_80) > 0 else tmax_h + 0.5 * t_half
    peak_end_h = max(peak_end_h, tmax_h)

    # Offset: ~50% of Cmax falling (primary subjective plateau wearing off)
    falling_50 = np.where(after_peak & (c <= 0.50 * cmax))[0]
    offset_h = float(t_dense[falling_50[0]]) if len(falling_50) > 0 else tmax_h + 1.5 * t_half
    offset_h = max(offset_h, peak_end_h)

    # Comedown end: ~15% of Cmax (or below MEC if MEC is in reasonable range [0.05*Cmax, 0.40*Cmax])
    comedown_thresh = 0.15 * cmax
    if window is not None and 0.05 * cmax <= window.mec <= 0.40 * cmax:
        comedown_thresh = window.mec
    falling_comedown = np.where(after_peak & (c <= comedown_thresh))[0]
    comedown_end_h = float(t_dense[falling_comedown[0]]) if len(falling_comedown) > 0 else tmax_h + 3.5 * t_half
    comedown_end_h = max(comedown_end_h, offset_h + 0.1)

    return PhaseTimings(
        onset_h=round(onset_h, 2),
        peak_start_h=round(peak_start_h, 2),
        tmax_h=round(tmax_h, 2),
        peak_end_h=round(peak_end_h, 2),
        offset_h=round(offset_h, 2),
        comedown_end_h=round(comedown_end_h, 2),
    )


def get_harm_reduction_guidelines(drug_id: str = "") -> list[dict[str, str]]:
    """Return core harm reduction guidelines in English and Czech."""
    return [
        {
            "id": "test_before_ingest",
            "title_cs": "1. Otestujte látku před užitím (Drug Checking)",
            "title_en": "1. Test before you ingest (Drug Checking)",
            "desc_cs": "Vždy ověřte přítomnost nebezpečných příměsí (fentanyl, nitazeny, syntetické kanabinoidy) testovacími proužky nebo činidly (Marquis, Mecke, Mandelin).",
            "desc_en": "Always test for dangerous adulterants (fentanyl, nitazenes, synthetic adulterants) using rapid test strips or chemical reagents.",
        },
        {
            "id": "start_low_go_slow",
            "title_cs": "2. Začněte s minimální dávkou (Start low, go slow)",
            "title_en": "2. Start low and go slow",
            "desc_cs": "Nejprve užijte malou testovací dávku (např. 1/4 obvyklého množství) k vyloučení alergické reakce a neočekávané síly. Vyčkejte na plný vrchol (Tmax) před zvažováním další dávky.",
            "desc_en": "Take an initial allergy/threshold test dose. Wait past the peak concentration (Tmax) before considering any further intake.",
        },
        {
            "id": "prevent_dose_stacking",
            "title_cs": "3. Zabraňte vrstvení dávek (Dose Stacking)",
            "title_en": "3. Prevent dose stacking",
            "desc_cs": "Pokud látka 'ještě nezabrala' nebo účinek slábne, neznamená to, že zmizela z těla. Předčasné přídavky vedou k nekontrolované akumulaci a akutní toxicitě.",
            "desc_en": "If effects feel delayed or fading, active substance remains in circulation. Redosing too early causes dangerous stacking and acute toxicity.",
        },
        {
            "id": "polysubstance_warning",
            "title_cs": "4. Vyvarujte se kombinací tlumivých látek (Polydrug risk)",
            "title_en": "4. Avoid combining CNS depressants",
            "desc_cs": "Kombinace benzodiazepinů (např. Rivotril, Diazepam), alkoholu a opioidů vede k synergickému útlumu dechového centra a je nejčastější příčinou úmrtí z předávkování.",
            "desc_en": "Combining benzodiazepines, alcohol, and opioids causes synergistic respiratory depression and is the primary driver of fatal overdoses.",
        },
        {
            "id": "hydration_and_care",
            "title_cs": "5. Hydratace, chlazení a nikdy nebuďte sami",
            "title_en": "5. Hydration, cooling, and never use alone",
            "desc_cs": "Při užití stimulantů doplňujte tekutiny elektrolyty (cca 250–500 ml/h, nepřeplňujte se vodou), dělejte přestávky v chladu a zajistěte přítomnost střízlivého dohledu.",
            "desc_en": "Sip electrolyte fluids (250–500 ml/h), take cool rest breaks with stimulants, and never use alone without someone knowing.",
        },
    ]
