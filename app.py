"""Streamlit app: interactive PK simulation & Therapeutic Drug Monitoring (TDM)."""
from __future__ import annotations

import math

import numpy as np
import streamlit as st

from core.drug_database import CUSTOM_ID, custom_drug, formulation_from_dosage_form, load_drugs
from core.enzymes import load_enzymes, perpetrator_pk, polymorphic_enzymes
from core.individualize import Observation, fit_individual, individual_band
from core.interactions import load_interactions
from core.regimen import SUBSTANCE, AddedSubstance, RegimenDrug, simulate_regimen
from core.validation import regimen_label, summarize, validate
from core.population import clearance_factor, cockcroft_gault
from core.metrics import (
    TherapeuticWindow,
    assess_tdm,
    compute_metrics,
    recovery_time,
)
from core.pk_engine import (
    DoseEvent,
    Route,
    apply_missed_dose,
    metabolite_profile,
    concentration_profile,
    regular_schedule,
    steady_state_concentration,
    steady_state_summary,
    time_grid,
)
from core.harm_reduction import (
    assess_dose_stacking,
    compute_safe_redosing_interval_h,
    get_harm_reduction_guidelines,
    get_route_timeline,
)
from core.sukl_client import (
    PRESET_DEFAULT_SUKL_CODES,
    fetch_sukl_detail,
    query_sukl_api,
)
from ui.charts import SERIES as OVERLAY_SERIES, build_enzyme_figure, build_figure, build_overlay_figure
from ui.i18n import LANGUAGES, make_translator
from ui.reporting import build_dataframe, build_protocol

st.set_page_config(page_title="PK/PD Simulator & TDM", page_icon="💊", layout="wide")


@st.cache_data(ttl=3600, show_spinner=False)
def search_sukl(query: str):
    return query_sukl_api(query)


@st.cache_data(ttl=3600, show_spinner=False)
def get_sukl_detail(kod_sukl: str):
    return fetch_sukl_detail(kod_sukl)


@st.cache_data(show_spinner=False)
def validation_rows():
    return validate()


def assessment_message(T, res, window: TherapeuticWindow, unit: str) -> tuple[str, str]:
    """Map a TDM assessment to an (alert level, text) pair."""
    if res.status == "toxic":
        return "error", T("alert_toxic", t=res.first_toxic_time_h, mtc=window.mtc, unit=unit,
                          hours=res.time_above_mtc_h, cmax=res.cmax)
    if res.status == "never_effective":
        return "warning", T("alert_never", mec=window.mec, unit=unit, cmax=res.cmax)
    if res.status == "trough_below_mec":
        return "warning", T("alert_trough", trough=res.trough, mec=window.mec, unit=unit,
                            t0=res.assessed_from_h)
    return "success", T("alert_ok")


def show_alert(level: str, text: str) -> None:
    {"error": st.error, "warning": st.warning, "success": st.success}[level](text)


# Parse incoming URL query parameters (e.g. from Next.js portal)
qp = st.query_params
qp_drug = str(qp.get("drug") or "").lower()
qp_sukl = str(qp.get("sukl") or "")
try:
    qp_dose = float(qp.get("dose")) if qp.get("dose") else None
except (ValueError, TypeError):
    qp_dose = None
try:
    qp_interval = float(qp.get("interval")) if qp.get("interval") else None
except (ValueError, TypeError):
    qp_interval = None
try:
    qp_weight = int(float(qp.get("weight"))) if qp.get("weight") else 90
except (ValueError, TypeError):
    qp_weight = 90
qp_route = str(qp.get("route") or "").lower()

# ------------------------------------------------------------------ sidebar
with st.sidebar:
    lang = st.radio("Language / Jazyk", list(LANGUAGES), horizontal=True)
    T = make_translator(LANGUAGES[lang])

    st.header(T("sb_drug"))

    init_sukl = qp_sukl if "sukl_q" not in st.session_state else st.session_state.get("sukl_q", "")
    sukl_query = st.text_input(
        T("sukl_search_label"),
        value=init_sukl,
        placeholder=T("sukl_search_placeholder"),
        key="sukl_q",
    )
    sukl_med = None
    if sukl_query.strip():
        sukl_results = search_sukl(sukl_query.strip())
        if sukl_results:
            st.caption(T("sukl_found_count", n=len(sukl_results)))
            chosen_sukl_idx = st.selectbox(
                T("sukl_select_pkg"),
                range(len(sukl_results)),
                format_func=lambda i: sukl_results[i].display_name,
                key="sukl_sel",
            )
            sukl_med = sukl_results[chosen_sukl_idx]
            st.markdown(f"[{T('sukl_link')}]({sukl_med.sukl_url})")
            if sukl_med.matched_preset_id:
                st.success(T("sukl_matched_preset", name=sukl_med.nazev))
            else:
                st.info(T("sukl_no_preset"))

    presets = load_drugs() + [custom_drug()]
    by_id = {d.id: d for d in presets}
    preset_ids = list(by_id)

    default_idx = 0
    if qp_drug and qp_drug in by_id:
        default_idx = preset_ids.index(qp_drug)
    elif sukl_med and sukl_med.matched_preset_id and sukl_med.matched_preset_id in by_id:
        default_idx = preset_ids.index(sukl_med.matched_preset_id)
    elif sukl_med and not sukl_med.matched_preset_id:
        default_idx = preset_ids.index(CUSTOM_ID)

    drug_id = st.selectbox(
        T("drug"),
        preset_ids,
        index=default_idx,
        format_func=lambda i: by_id[i].name,
        key=f"drug_sel_{sukl_med.kod_sukl if sukl_med else 'none'}",
    )
    base = by_id[drug_id]
    k = drug_id  # per-drug widget keys: defaults reset when the drug changes
    if base.formulation_ids:
        sukl_form = " ".join(filter(None, (getattr(sukl_med, "lekova_forma", ""), getattr(sukl_med, "doplnek", "")))) \
            if sukl_med else ""
        auto_form = formulation_from_dosage_form(base, sukl_form)
        form_ids = base.formulation_ids
        form_id = st.radio(
            T("formulation"), form_ids, index=form_ids.index(auto_form) if auto_form else 0,
            format_func=lambda f: base.formulation_name(f, LANGUAGES[lang]), horizontal=True, key=f"form_{k}")
        if auto_form and sukl_med:
            st.caption(T("formulation_auto", form=base.formulation_name(auto_form, LANGUAGES[lang])))
        base = base.with_formulation(form_id)
    if base.notes:
        st.caption(base.notes)

    init_weight = qp_weight if 30 <= qp_weight <= 150 else 70
    weight = st.slider(T("body_weight") + " [kg]", 30, 150, init_weight, key="weight")

    crcl = None
    if st.checkbox(T("renal_enable"), False, key="renal_on"):
        rc = st.columns(2)
        age = rc[0].number_input(T("age") + f" [{T('years')}]", 18, 110, 50, key="age")
        female = rc[1].checkbox(T("sex_female"), False, key="female")
        scr = st.number_input(T("creatinine") + " [µmol/L]", 20.0, 1500.0, 80.0, 1.0, key="scr")
        crcl = cockcroft_gault(age, weight, scr, female)
        if base.fe_renal > 0:
            st.caption(T("crcl_result", crcl=crcl, factor=clearance_factor(weight, crcl, base.fe_renal)))
        else:
            st.caption(T("renal_not_used", name=base.name))

    # Only enzymes this drug is actually cleared by, and only those with a phenotype system.
    genotype_enzymes = sorted((set(base.pathways or {}) & polymorphic_enzymes()))
    genotypes: dict[str, str] = {}
    if genotype_enzymes:
        with st.expander(T("sb_genotype")):
            st.caption(T("genotype_help"))
            all_enzymes = load_enzymes()
            for enzyme in genotype_enzymes:
                labels = [T("phenotype_unknown")] + list(all_enzymes[enzyme].genotype_activity)
                choice = st.selectbox(
                    T("genotype_of", enzyme=enzyme), labels,
                    format_func=lambda p: p if p == T("phenotype_unknown") else T("phenotype_" + p),
                    key=f"gt_{enzyme}_{k}")
                if choice != T("phenotype_unknown"):
                    genotypes[enzyme] = choice
                    st.caption(T("genotype_activity",
                                 pct=100 * all_enzymes[enzyme].genotype_activity[choice]))

    with st.expander(T("sb_params")):
        t_half = st.number_input("t½ [h]", 0.1, 500.0, base.t_half_h, 0.1, key=f"th_{k}")
        cv_cl = st.number_input(T("cv_cl") + " [%]", 0.0, 200.0, 100.0 * (base.cv_cl or 0.3), 1.0,
                                key=f"cvcl_{k}") / 100.0
        cv_f = st.number_input(T("cv_f") + " [%]", 0.0, 200.0, 100.0 * base.cv_f, 1.0, key=f"cvf_{k}") / 100.0
        vd_kg = st.number_input("V_d [L/kg]", 0.01, 20.0, base.vd_l_per_kg, 0.01, key=f"vd_{k}")
        f_bio = st.slider("F", 0.05, 1.0, base.bioavailability, 0.01, key=f"f_{k}")
        ka = st.number_input("k_a [1/h]", 0.05, 20.0, base.ka_per_h, 0.05, key=f"ka_{k}")
        unit = st.text_input(T("conc_unit"), base.concentration_unit, key=f"u_{k}") if drug_id == CUSTOM_ID \
            else base.concentration_unit
        mec = st.number_input(f"MEC [{unit}]", 0.0001, 10000.0, base.mec, format="%g", key=f"mec_{k}")
        mtc = st.number_input(f"MTC [{unit}]", 0.0001, 10000.0, base.mtc, format="%g", key=f"mtc_{k}")

    st.header(T("sb_regimen"))
    route_choices = list(Route)
    route_vals = [r.value for r in route_choices]
    if qp_route in route_vals:
        default_route_idx = route_vals.index(qp_route)
    elif sukl_med and any(t in (sukl_med.cesta_podani or "").lower() for t in ["intraven", "injek", "infuz", "i.v."]):
        default_route_idx = 1
    else:
        default_route_idx = 0
    route = st.radio(
        T("route"),
        route_choices,
        index=default_route_idx,
        format_func=lambda r: T("route_" + r.value),
        horizontal=True,
        key=f"route_{k}_{sukl_med.kod_sukl if sukl_med else 'none'}",
    )
    if qp_dose is not None:
        suggested_dose = float(qp_dose)
    elif sukl_med and sukl_med.parsed_dose_mg is not None:
        suggested_dose = float(sukl_med.parsed_dose_mg)
    else:
        suggested_dose = float(base.default_dose_mg)
    dose_mg = st.number_input(
        T("dose") + " [mg]",
        0.1,
        100000.0,
        suggested_dose,
        format="%g",
        key=f"dose_{k}_{sukl_med.kod_sukl if sukl_med else 'none'}",
    )
    if sukl_med and sukl_med.parsed_dose_mg is not None:
        st.caption(T("sukl_dose_set", dose=suggested_dose))
    else:
        st.caption(base.dose_basis)
    presets_tau = [6, 8, 12, 24]
    default_tau = qp_interval if qp_interval is not None else base.default_interval_h
    tau_choice = st.selectbox(
        T("interval"), presets_tau + ["custom"],
        index=presets_tau.index(default_tau) if default_tau in presets_tau else len(presets_tau),
        format_func=lambda v: T("custom") if v == "custom" else f"{v} h", key=f"tau_{k}")
    tau = float(st.number_input(T("interval") + " [h]", 0.5, 168.0, float(default_tau), 0.5, key=f"tauc_{k}")) \
        if tau_choice == "custom" else float(tau_choice)
    default_days = int(min(30, max(3, math.ceil(5 * base.t_half_h / 24) + 1)))
    days = st.slider(T("duration") + f" [{T('days_unit')}]", 1, 30, default_days, key=f"days_{k}")

    st.header(T("sb_display"))
    log_scale = st.radio(T("axis_scale"), [False, True], format_func=lambda v: T("semilog") if v else T("linear"),
                         horizontal=True)
    show_band = st.checkbox(T("show_band"), True)
    show_zones = st.checkbox(T("show_zones"), True)

    duration_h = days * 24.0
    events = regular_schedule(1.0, tau, duration_h)  # only for counting doses here
    st.header(T("sb_missed"))
    missed_on = st.checkbox(T("missed_enable"), False)
    missed_index, delay_h = 0, None
    if missed_on:
        if len(events) > 1:
            missed_index = st.slider(T("missed_which"), 1, len(events), min(3, len(events)), key=f"mi_{k}") - 1
        st.caption(T("missed_at", t=events[missed_index].time_h))
        mode = st.radio(T("missed_mode"), ["skip", "delay"], format_func=lambda m: T("mode_" + m), horizontal=True)
        if mode == "delay":
            delay_h = st.slider(T("missed_delay") + " [h]", 0.5, float(max(1.0, tau)), min(2.0, float(tau)), 0.5)

# ---------------------------------------------------------------- simulation
try:
    drug = base.with_overrides(
        t_half_h=t_half, t_half_range_h=None, cv_cl=cv_cl, cv_f=cv_f, vd_l_per_kg=vd_kg,
        bioavailability=f_bio, ka_per_h=ka, mec=mec, mtc=mtc, concentration_unit=unit)
    window = TherapeuticWindow(mec, mtc)
except ValueError:
    st.error(T("err_window"))
    st.stop()

amount = drug.amount_from_dose_mg(dose_mg)
model = drug.to_pk_model(weight, route, crcl_ml_min=crcl, genotypes=genotypes or None)
ref_events = regular_schedule(amount, tau, duration_h)
active_events = apply_missed_dose(ref_events, missed_index, delay_h) if missed_on else ref_events
t = time_grid(duration_h)
c_ref = concentration_profile(t, ref_events, model)
c = concentration_profile(t, active_events, model) if missed_on else c_ref
curves = [concentration_profile(t, active_events, m)
          for m in drug.exposure_band(model, genotypes or None)]
band = (np.minimum(*curves), np.maximum(*curves))

metrics = compute_metrics(t, c, model)
ss = steady_state_summary(amount, tau, model)
assess_from = max(0.0, min(ss.time_to_ss_h, duration_h - tau))
assessment = assess_tdm(t, c, window, assess_from)
level, message = assessment_message(T, assessment, window, unit)

# ----------------------------------------------------------------- main page
st.title("💊 " + T("title"))
st.caption(T("subtitle"))

if qp_sukl or qp_drug:
    ref_code = qp_sukl or (active_sukl_code if "active_sukl_code" in locals() else "")
    st.info(
        f"🔗 **Propojeno z Lékového Portálu ČR** (SÚKL: `{ref_code}`, Lék: `{base.name}`). "
        f"[← Zpět na moderní portál](http://localhost:3000/leky/{ref_code})"
    )

active_sukl_code = sukl_med.kod_sukl if sukl_med else PRESET_DEFAULT_SUKL_CODES.get(drug_id)
drug_detail = get_sukl_detail(active_sukl_code) if active_sukl_code else None

if drug_detail:
    c_disp, c_reimb, c_supp = st.columns(3)
    disp_icon = "🔵" if drug_detail.zpusob_vydeje_kod == "blue_stripe" else ("💊" if drug_detail.zpusob_vydeje_kod == "rx" else "🟢")
    disp_text = T("sukl_status_" + drug_detail.zpusob_vydeje_kod)
    c_disp.info(f"{disp_icon} **{T('card_dispensing')}:** {disp_text}")

    reimb_icon = "🟢" if drug_detail.uhrada_rezim == "full" else ("🟡" if drug_detail.uhrada_rezim == "partial" else "⚪")
    reimb_text = T("sukl_reimb_" + drug_detail.uhrada_rezim)
    if drug_detail.uhrada_rezim == "partial" and drug_detail.max_doplatek is not None:
        reimb_detail = f" (doplatek {drug_detail.max_doplatek:.0f} Kč)"
    elif drug_detail.uhrada_rezim == "none" and drug_detail.max_cena_lek is not None:
        reimb_detail = f" (cena cca {drug_detail.max_cena_lek:.0f} Kč)"
    else:
        reimb_detail = ""
    c_reimb.info(f"{reimb_icon} **{T('card_reimbursement')}:** {reimb_text}{reimb_detail}")

    supp_icon = "✅" if drug_detail.dostupnost_status == "active" else ("⚠️" if drug_detail.dostupnost_status == "interrupted" else "❌")
    supp_text = T("sukl_shortage_" + drug_detail.dostupnost_status)
    c_supp.info(f"{supp_icon} **{T('card_supply_status')}:** {supp_text}")

cols = st.columns(6)
cols[0].metric("C_max", f"{metrics.cmax:.3g} {unit}", help=T("help_cmax"))
cols[1].metric("T_max", f"{metrics.tmax_h:.2f} h", help=T("help_tmax"))
cols[2].metric("AUC₀₋₂₄", f"{metrics.auc_0_24:.3g} {unit}·h", help=T("help_auc"))
cols[3].metric(T("clearance"), f"{metrics.clearance_l_per_h:.2f} L/h", help=T("help_clearance"))
cols[4].metric("t½", f"{metrics.t_half_h:.3g} h", help=T("help_thalf"))
cols[5].metric(T("time_to_ss"), f"{ss.time_to_ss_h:.0f} h", help=T("tss_help"))

show_alert(level, message)
if duration_h < ss.time_to_ss_h:
    st.info(T("info_short_sim", tss=ss.time_to_ss_h, days=math.ceil(ss.time_to_ss_h / 24)))

missed_time = ref_events[missed_index].time_h if missed_on else None
fig = build_figure(
    t, c, window, unit, T, log_scale=log_scale, band=band if show_band else None,
    reference=c_ref if missed_on else None, dose_times=[e.time_h for e in active_events],
    missed_time=missed_time, show_zones=show_zones,
    extra=((metabolite_profile(t, active_events, model), T("metabolite_curve", name=model.metabolite.name))
           if model.metabolite is not None and not missed_on else None))
st.plotly_chart(fig, width="stretch")
if model.metabolite is not None:
    st.caption(T("metabolite_caption", name=model.metabolite.name,
                 t_half=model.metabolite.t_half_h, reported=T("metabolite_reported_" + model.reported_quantity)))

tab_ss, tab_fit, tab_multi, tab_missed, tab_hr, tab_sukl, tab_glossary, tab_data, tab_model = st.tabs(
    [T("tab_ss"), T("tab_fit"), T("tab_multi"), T("tab_missed"), T("tab_hr"), T("tab_sukl"), T("tab_glossary"),
     T("tab_data"), T("tab_model")])

with tab_ss:
    r = st.columns(4)
    r[0].metric("C_ss,max", f"{ss.cmax_ss:.3g} {unit}", help=f"{T('help_css_max')} (t = {ss.tmax_ss_h:.2f} h)")
    r[1].metric("C_ss,min (trough)", f"{ss.cmin_ss:.3g} {unit}", help=T("help_css_min"))
    r[2].metric("C_ss,avg", f"{ss.cavg_ss:.3g} {unit}", help=T("help_css_avg"))
    r[3].metric(T("accumulation") + " R", f"{ss.accumulation_index:.2f}", help=T("help_r"))
    r = st.columns(4)
    r[0].metric(T("fluctuation"), f"{ss.fluctuation_pct:.0f} %", help=T("help_fluct"))
    r[1].metric("AUC_τ,ss", f"{ss.cavg_ss * tau:.3g} {unit}·h", help=T("help_auc"))
    r[2].metric(T("time_to_ss"), f"{ss.time_to_ss_h:.0f} h",
                help=f"{ss.time_to_ss_range_h[0]:.0f}–{ss.time_to_ss_range_h[1]:.0f} h (4–5 × t½)")
    r[3].metric("τ", f"{tau:g} h", help=T("help_tau"))
    ss_t = np.linspace(0.0, tau, 400)
    ss_c = steady_state_concentration(ss_t, amount, tau, model)
    ss_assess = assess_tdm(ss_t, ss_c, window)
    show_alert(*assessment_message(T, ss_assess, window, unit))
    st.plotly_chart(build_figure(ss_t, ss_c, window, unit, T, log_scale=log_scale, dose_times=[0.0],
                                 show_zones=show_zones, height=380), width="stretch")
    st.caption(T("ss_caption"))

with tab_fit:
    st.markdown(T("fit_intro"))
    if st.checkbox(T("fit_enable"), False, key=f"fit_on_{k}"):
        with st.expander(T("fit_dose_log")):
            st.caption(T("fit_dose_log_help"))
            log = st.data_editor(
                [{"t [h]": e.time_h, T("dose") + " [mg]": e.amount / drug.dose_conversion} for e in active_events],
                num_rows="dynamic", hide_index=True, key=f"fit_log_{k}")
        fit_events = sorted(
            (DoseEvent(float(row["t [h]"]), drug.amount_from_dose_mg(float(row[T("dose") + " [mg]"])))
             for row in log if row.get("t [h]") is not None and row.get(T("dose") + " [mg]")),
            key=lambda e: e.time_h)
        if not fit_events:
            st.error(T("fit_no_doses"))
        else:
            c_log = concentration_profile(t, fit_events, model)
            n_obs = int(st.number_input(T("fit_n"), 1, 5, 1, key=f"fit_n_{k}"))
            # Default: troughs just before the last doses, prefilled with the population prediction.
            last_doses = [e.time_h for e in fit_events if e.time_h > 0][-n_obs:] or [duration_h]
            default_times = ([last_doses[0]] * n_obs + last_doses)[-n_obs:]
            observations = []
            for i in range(n_obs):
                oc = st.columns(2)
                t_i = oc[0].number_input(T("fit_time", i=i + 1), 0.1, float(duration_h),
                                         float(default_times[i]), 0.5, key=f"fit_t{i}_{k}")
                c_default = float(np.interp(t_i, t, c_log))
                c_i = oc[1].number_input(T("fit_level", i=i + 1, unit=unit), 0.0, 1e6,
                                         float(f"{c_default:.3g}"), format="%g", key=f"fit_c{i}_{k}")
                observations.append(Observation(t_i, c_i))
            try:
                fit = fit_individual(model, fit_events, observations, drug.cv_cl, drug.cv_f)
            except ValueError as exc:
                st.error(str(exc))
            else:
                fit_ss = steady_state_summary(amount, tau, fit.model)
                fc = st.columns(4)
                fc[0].metric(T("fit_cl"), f"{fit.model.clearance_l_per_h:.3g} L/h",
                             f"{100 * (fit.cl_ratio - 1):+.0f} %", delta_color="off",
                             help=T("fit_cl_help", pop=model.clearance_l_per_h))
                fc[1].metric(T("fit_thalf"), f"{fit.model.terminal_half_life_h:.3g} h",
                             help=T("fit_thalf_help", pop=model.terminal_half_life_h))
                fc[2].metric(T("fit_css_avg"), f"{fit_ss.cavg_ss:.3g} {unit}")
                fc[3].metric(T("fit_css_min"), f"{fit_ss.cmin_ss:.3g} {unit}")
                c_fit = concentration_profile(t, fit_events, fit.model)
                show_alert(*assessment_message(T, assess_tdm(t, c_fit, window, assess_from), window, unit))
                st.plotly_chart(build_figure(
                    t, c_fit, window, unit, T, log_scale=log_scale, band=individual_band(fit, t, fit_events),
                    reference=c_log, dose_times=[e.time_h for e in fit_events], show_zones=show_zones,
                    observations=[(o.time_h, o.concentration) for o in fit.observations],
                    conc_label=T("fit_individual"), reference_label=T("fit_population"), band_label=T("fit_band"),
                    height=440), width="stretch")
                st.dataframe(
                    [{"t [h]": o.time_h, T("fit_table_obs"): o.concentration, T("fit_table_pop"): round(p_, 4),
                      T("fit_table_ind"): round(i_, 4)}
                     for o, p_, i_ in zip(fit.observations, fit.population_predicted, fit.predicted)],
                    hide_index=True)
                st.caption(T("fit_caption"))

with tab_multi:
    st.markdown(T("multi_intro"))
    lang_code = LANGUAGES[lang]
    other_ids = [d.id for d in presets if d.id not in (drug_id, CUSTOM_ID)]
    co_ids = st.multiselect(T("multi_select"), other_ids, format_func=lambda i: by_id[i].name,
                            max_selections=len(OVERLAY_SERIES) - 1, key=f"multi_{k}")
    perp_names = sorted({i.perpetrator for i in load_interactions()} - {by_id[i].name for i in co_ids})
    co_perps = st.multiselect(T("multi_perp"), perp_names,
                              format_func=lambda n: next(i.name(lang_code) for i in load_interactions()
                                                         if i.perpetrator == n), key=f"multiperp_{k}")
    # A substance with known PK can be started and stopped; by default it is already being taken.
    perp_label = {i.perpetrator: i.name(lang_code) for i in load_interactions()}
    windows: dict[str, tuple[float | None, float | None]] = {}
    timed = [n for n in co_perps if (p := perpetrator_pk(n)) is not None and p.t_half_h]
    if timed:
        with st.expander(T("multi_timing")):
            st.caption(T("multi_timing_help"))
            for n in timed:
                wc = st.columns(3)
                already = wc[0].checkbox(T("multi_already", name=perp_label[n]), True, key=f"pa_{n}_{k}")
                start_h = None if already else 24.0 * wc[1].number_input(
                    T("multi_start", name=perp_label[n]), 0, int(days), 0, key=f"ps_{n}_{k}")
                stop_d = wc[2].number_input(T("multi_stop", name=perp_label[n]), 0, int(days), 0, key=f"pe_{n}_{k}")
                stop_h = 24.0 * stop_d if stop_d else None
                if stop_h is not None and start_h is not None and stop_h <= start_h:
                    st.warning(T("multi_stop_before_start", name=perp_label[n]))
                    stop_h = None
                windows[n] = (start_h, stop_h)
    if not co_ids and not co_perps:
        st.info(T("multi_none"))
    else:
        entries = [RegimenDrug(drug, model, tuple(regular_schedule(amount, tau, duration_h)), tau,
                               drug.enzyme_shares(weight, crcl, genotypes or None))]
        for cid in co_ids:
            co = by_id[cid]
            mc = st.columns(2)
            co_dose = mc[0].number_input(T("multi_dose", name=co.name), 0.1, 100000.0, float(co.default_dose_mg),
                                         format="%g", key=f"md_{cid}_{k}")
            co_tau = mc[1].number_input(T("multi_interval", name=co.name), 0.5, 168.0, float(co.default_interval_h),
                                        0.5, key=f"mi_{cid}_{k}")
            entries.append(RegimenDrug(
                co, co.to_pk_model(weight, Route.ORAL, crcl_ml_min=crcl, genotypes=genotypes or None),
                tuple(regular_schedule(co.amount_from_dose_mg(co_dose), co_tau, duration_h)), co_tau,
                co.enzyme_shares(weight, crcl, genotypes or None)))
        res = simulate_regimen(t, entries, [AddedSubstance(n, *windows.get(n, (None, None))) for n in co_perps])

        def label(key: str) -> str:
            return res.labels[key][1 if lang_code == "cs" else 0]

        def shown(name: str) -> str:  # a driver is a regimen drug's name or an added substance's
            return label(SUBSTANCE + name) if SUBSTANCE + name in res.labels else name

        st.markdown(f"#### {T('multi_overlay_title')}")
        st.plotly_chart(build_overlay_figure(
            t, [(r.entry.drug.name, r.concentrations, TherapeuticWindow(r.entry.drug.mec, r.entry.drug.mtc),
                 r.entry.drug.concentration_unit) for r in res.drugs],
            T, show_zones=show_zones), width="stretch")
        st.caption(T("multi_overlay_caption"))
        if res.enzyme_activity:
            st.markdown(f"#### {T('multi_enzyme_title')}")
            st.plotly_chart(build_enzyme_figure(t, res.enzyme_activity, T), width="stretch")
            st.caption(T("multi_enzyme_caption", pairs="; ".join(
                T("multi_enzyme_driver", enzyme=e, drugs=", ".join(shown(n) for n in res.enzyme_drivers[e]))
                for e in res.enzyme_activity)))

        st.markdown(f"#### {T('pd_title')}")
        for clash in res.clashes:
            who = ", ".join(label(i) for i in clash.drug_ids)
            if clash.overlap_h > 0:
                st.warning(T("pd_clash", effect=T("pd_effect_" + clash.effect), drugs=who,
                             hours=clash.overlap_h, first=clash.first_h))
            else:
                st.info(T("pd_clash_apart", effect=T("pd_effect_" + clash.effect), drugs=who))
        if not res.clashes and len(res.tags) >= 2:
            st.caption(T("pd_none"))
        if res.unassessed:
            st.caption(T("pd_unassessed", drugs=", ".join(label(k) for k in res.unassessed)))
        pd_rows = [{T("drug"): label(key), T("pd_effect"): T("pd_effect_" + tag.effect),
                    T("validation_source"): tag.source, T("validation_note"): tag.note}
                   for key, tags in res.tags.items() for tag in tags]
        if pd_rows:
            with st.expander(T("pd_sources")):
                st.dataframe(pd_rows, hide_index=True, width="stretch",
                             column_config={T("validation_source"): st.column_config.LinkColumn(
                                 display_text="source")})

        rows = []
        for r in res.drugs:
            d = r.entry.drug
            win_i = TherapeuticWindow(d.mec, d.mtc)
            acts = [a for a in r.applied if a.factor is not None or a.via_enzyme]
            st.markdown(f"**{d.name}**")
            if acts or r.drivers:
                st.metric(T("multi_factor_sim") if r.drivers else T("multi_factor"), f"{r.exposure_factor:.2f}x",
                          f"{100 * (r.exposure_factor - 1):+.0f} %", delta_color="off")
                if r.drivers:
                    st.caption(T("multi_follows", drugs=", ".join(shown(n) for n in r.drivers)))
                for enz, (multiplied, combined) in r.pooled.items():
                    st.caption(T("multi_pooled", enzyme=enz, multiplied=multiplied, combined=combined,
                                 drugs=", ".join(shown(a.interaction.perpetrator) for a in r.applied
                                                 if a.pooled == enz)))
            else:
                st.caption(T("multi_no_interaction"))
            show_alert(*assessment_message(T, assess_tdm(t, r.concentrations, win_i, assess_from), win_i,
                                           d.concentration_unit))
            with st.expander(T("multi_curve", name=d.name)):
                st.plotly_chart(build_figure(
                    t, r.concentrations, win_i, d.concentration_unit, T, log_scale=log_scale,
                    reference=r.alone if acts or r.drivers else None, dose_times=[e.time_h for e in r.entry.events],
                    show_zones=show_zones, height=320, conc_label=T("multi_with", name=d.name),
                    reference_label=T("multi_alone", name=d.name)), width="stretch")
            for a in r.applied:
                i = a.interaction
                if a.skipped == "iv":
                    effect = T("multi_skip_iv")
                elif a.via_enzyme:
                    effect = T("multi_via_enzyme")
                elif a.factor is None:
                    effect = T("multi_qualitative")
                else:
                    effect = f"x{a.factor:.2f} ({100 * (a.factor - 1):+.0f} %)"
                    if a.pooled:
                        effect += f", {T('multi_pooled_row', enzyme=a.pooled)}"
                onset = f"{i.onset_h / 24:.0f} d" if i.onset_h >= 24 else f"{i.onset_h:.0f} h"
                rows.append({T("drug"): d.name, T("multi_perpetrator"): i.name(lang_code),
                             T("multi_effect"): effect,
                             T("multi_evidence"): T("multi_estimated") if i.estimated else T("multi_measured"),
                             T("multi_site"): T("multi_site_" + i.site),
                             T("multi_onset"): T("multi_timed") if a.timed else onset,
                             T("validation_source"): i.source,
                             T("validation_note"): (i.basis or i.note) if i.estimated else i.note})
        if any(i for i in rows if i[T("multi_evidence")] == T("multi_estimated")):
            st.info(T("multi_estimated_help"))
        if rows:
            st.markdown(f"#### {T('multi_table_title')}")
            st.dataframe(rows, hide_index=True, width="stretch",
                         column_config={T("validation_source"): st.column_config.LinkColumn(display_text="source")})
        st.caption(T("multi_caption"))

with tab_missed:
    if not missed_on:
        st.info(T("missed_hint"))
    else:
        t_m = ref_events[missed_index].time_h
        win = (max(0.0, t_m - 1.5 * tau), min(duration_h, t_m + 4 * tau))
        ref_a = assess_tdm(t, c_ref, window, t_m)
        alt_a = assess_tdm(t, c, window, t_m)
        rec = recovery_time(t, c_ref, c, t_m)
        show_alert(*assessment_message(T, alt_a, window, unit))
        a = st.columns(4)
        a[0].metric(T("min_after"), f"{alt_a.trough:.3g} {unit}", f"{alt_a.trough - ref_a.trough:+.3g}")
        a[1].metric(T("time_below_mec"), f"{alt_a.time_below_mec_h:.1f} h",
                    f"{alt_a.time_below_mec_h - ref_a.time_below_mec_h:+.1f} h", delta_color="inverse")
        a[2].metric(T("time_above_mtc"), f"{alt_a.time_above_mtc_h:.1f} h",
                    f"{alt_a.time_above_mtc_h - ref_a.time_above_mtc_h:+.1f} h", delta_color="inverse")
        a[3].metric(T("recovery"), f"{rec:.1f} h" if rec is not None else T("not_recovered"))
        st.plotly_chart(build_figure(
            t, c, window, unit, T, log_scale=log_scale, reference=c_ref,
            dose_times=[e.time_h for e in active_events], missed_time=t_m, x_range=win,
            show_zones=show_zones, height=420), width="stretch")

with tab_hr:
    st.markdown(f"### 🛡️ {T('hr_heading')}")
    st.caption("Interaktivní kalkulačka bezpečného rozestupu mezi dávkami a prevence kumulativního předávkování (dose stacking).")

    col_hr1, col_hr2 = st.columns(2)
    with col_hr1:
        elapsed_val = st.slider(
            T("hr_elapsed_slider"),
            0.0,
            float(max(24.0, math.ceil(model.terminal_half_life_h * 5))),
            float(round(tau, 1)),
            0.5,
            key=f"hr_elapsed_{k}",
        )
    with col_hr2:
        redose_mg = st.number_input(
            T("hr_redose_amount"),
            0.1,
            100000.0,
            float(dose_mg),
            format="%g",
            key=f"hr_redose_{k}",
        )
        redose_amount = drug.amount_from_dose_mg(redose_mg)

    safe_interval_h = compute_safe_redosing_interval_h(
        amount, model, window, redose_amount=redose_amount
    )

    assessment = assess_dose_stacking(
        amount, elapsed_val, model, window, unit=unit, redose_amount=redose_amount
    )

    h_cols = st.columns(4)
    h_cols[0].metric(T("hr_safe_interval_metric"), f"{safe_interval_h:.1f} h", help="Minimální doporučený odstup pro zamezení kumulace")
    h_cols[1].metric(T("hr_current_level_metric"), f"{assessment.current_level:.3g} {unit}", help="Zbytková koncentrace z předchozí dávky")
    h_cols[2].metric(T("hr_stacking_metric"), f"{assessment.stacked_peak:.3g} {unit}", help="Předpokládaný součet při podání další dávky")
    remaining_to_safe = max(0.0, safe_interval_h - elapsed_val)
    h_cols[3].metric("Do bezpečného užití", f"{remaining_to_safe:.1f} h" if remaining_to_safe > 0 else "0.0 h (Bezpečné)")

    if assessment.risk_level == "critical_toxic":
        st.error(f"🚨 **{assessment.risk_headline}**\n\n{assessment.risk_explanation}")
    elif assessment.risk_level == "high":
        st.warning(f"⚠️ **{assessment.risk_headline}**\n\n{assessment.risk_explanation}")
    elif assessment.risk_level == "moderate":
        st.info(f"ℹ️ **{assessment.risk_headline}**\n\n{assessment.risk_explanation}")
    else:
        st.success(f"✅ **{assessment.risk_headline}**\n\n{assessment.risk_explanation}")

    # Simulated curve of stacking
    t_hr = np.linspace(0.0, max(24.0, elapsed_val + 3 * model.terminal_half_life_h), 500)
    c_stacked = concentration_profile(
        t_hr, [DoseEvent(0.0, amount), DoseEvent(elapsed_val, redose_amount)], model)
    fig_hr = build_figure(
        t_hr, c_stacked, window, unit, T,
        log_scale=log_scale,
        dose_times=[0.0, elapsed_val],
        show_zones=show_zones,
        height=380,
    )
    st.plotly_chart(fig_hr, width="stretch")

    st.markdown(f"#### ⏱️ {T('hr_routes_title')}")
    r_cols = st.columns(4)
    for idx, r_val in enumerate([Route.ORAL, Route.IV_BOLUS, Route.INSUFFLATED, Route.SMOKED]):
        tl = get_route_timeline(r_val, amount, model, drug=drug)
        with r_cols[idx]:
            with st.container(border=True):
                st.markdown(f"**{tl.route_name_cs if lang == 'Čeština' else tl.route_name_en}**")
                st.write(f"⚡ **Nástup:** {tl.onset_minutes_range[0]}–{tl.onset_minutes_range[1]} min")
                st.write(f"📈 **Vrchol (Tmax):** {tl.tmax_h:.2f} h")
                st.write(f"⏳ **Trvání:** {tl.duration_hours_range[0]:.1f}–{tl.duration_hours_range[1]:.1f} h")
                st.write(f"💧 **Biolog. dostupnost (F):** {tl.bioavailability_pct:.0f} %")
                st.caption(tl.description_cs if lang == 'Čeština' else tl.description_en)

    st.markdown(f"#### 📋 {T('hr_principles_title')}")
    guidelines = get_harm_reduction_guidelines(drug.id)
    g_cols = st.columns(len(guidelines))
    for idx, g in enumerate(guidelines):
        with g_cols[idx]:
            with st.container(border=True):
                title = g["title_cs"] if lang == "Čeština" else g["title_en"]
                desc = g["desc_cs"] if lang == "Čeština" else g["desc_en"]
                st.markdown(f"**{title}**")
                st.caption(desc)

with tab_sukl:
    if not drug_detail:
        st.info(T("no_sukl_data"))
    else:
        st.markdown(f"### 🏛️ {T('sukl_detail_heading')}")
        st.caption(f"**{drug_detail.nazev} {drug_detail.doplnek}** [SÚKL: `{drug_detail.kod_sukl}`]")
        if drug_detail.drzitel:
            st.write(f"🏢 **{T('sukl_holder')}:** {drug_detail.drzitel} | 📝 **{T('sukl_reg_no')}:** {drug_detail.registracni_cislo}")

        col_p1, col_p2 = st.columns(2)
        with col_p1:
            with st.container(border=True):
                st.markdown(f"#### 💰 {T('card_reimbursement')}")
                c1, c2 = st.columns(2)
                uhrada_str = f"{drug_detail.uhrada_zp:.2f} Kč" if drug_detail.uhrada_zp is not None else "0.00 Kč (0 %)"
                c1.metric(T("reimb_zp"), uhrada_str)
                doplatek_str = f"{drug_detail.max_doplatek:.2f} Kč" if drug_detail.max_doplatek is not None else "—"
                c2.metric(T("copay_patient"), doplatek_str)

                c3, c4 = st.columns(2)
                max_cena_str = f"{drug_detail.max_cena_lek:.2f} Kč" if drug_detail.max_cena_lek is not None else "—"
                c3.metric(T("max_price_pharmacy"), max_cena_str)
                limit_str = f"{drug_detail.zapocitatelny_doplatek:.2f} Kč" if drug_detail.zapocitatelny_doplatek is not None else "0.00 Kč"
                c4.metric(T("reimb_limit_copay"), limit_str)

        with col_p2:
            with st.container(border=True):
                st.markdown(f"#### 🚨 {T('card_supply_status')} & {T('card_dispensing')}")
                st.write(f"**{T('card_supply_status')}:** {T('sukl_shortage_' + drug_detail.dostupnost_status)} – *{drug_detail.dostupnost_text}*")
                st.write(f"**{T('card_dispensing')}:** {T('sukl_status_' + drug_detail.zpusob_vydeje_kod)} – *{drug_detail.zpusob_vydeje_text}*")
                st.divider()
                doc_cols = st.columns(2)
                if drug_detail.spc_url:
                    doc_cols[0].markdown(f"📄 [{T('download_spc')}]({drug_detail.spc_url})")
                if drug_detail.pil_url:
                    doc_cols[1].markdown(f"📦 [{T('download_pil')}]({drug_detail.pil_url})")
                st.markdown(f"🔗 [{T('sukl_link')}]({drug_detail.sukl_url})")

        if drug_detail.preskripcni_omezeni or drug_detail.indikacni_omezeni:
            with st.expander(f"🩺 {T('reimb_indications')} & {T('reimb_doctor')}", expanded=True):
                if drug_detail.preskripcni_omezeni:
                    st.markdown(f"**{T('reimb_doctor')}:** {drug_detail.preskripcni_omezeni}")
                if drug_detail.indikacni_omezeni:
                    st.markdown(f"**{T('reimb_indications')}:** {drug_detail.indikacni_omezeni}")

with tab_glossary:
    st.markdown(f"### 📖 {T('tab_glossary')}")
    st.caption(T("glossary_intro"))

    col_g1, col_g2 = st.columns(2)
    with col_g1:
        with st.container(border=True):
            st.markdown(f"**{T('glossary_cmax_title')}**")
            st.write(T("glossary_cmax_desc"))
        with st.container(border=True):
            st.markdown(f"**{T('glossary_tmax_title')}**")
            st.write(T("glossary_tmax_desc"))
        with st.container(border=True):
            st.markdown(f"**{T('glossary_thalf_title')}**")
            st.write(T("glossary_thalf_desc"))
        with st.container(border=True):
            st.markdown(f"**{T('glossary_cl_title')}**")
            st.write(T("glossary_cl_desc"))
        with st.container(border=True):
            st.markdown(f"**{T('glossary_vd_title')}**")
            st.write(T("glossary_vd_desc"))

    with col_g2:
        with st.container(border=True):
            st.markdown(f"**{T('glossary_css_title')}**")
            st.write(T("glossary_css_desc"))
        with st.container(border=True):
            st.markdown(f"**{T('glossary_window_title')}**")
            st.write(T("glossary_window_desc"))
        with st.container(border=True):
            st.markdown(f"**{T('glossary_auc_title')}**")
            st.write(T("glossary_auc_desc"))
        with st.container(border=True):
            st.markdown(f"**{T('glossary_r_title')}**")
            st.write(T("glossary_r_desc"))

with tab_data:
    df = build_dataframe(t, c, window, unit, band=band, reference=c_ref if missed_on else None)
    missed_text = (T("mode_skip") if delay_h is None else f"{T('mode_delay')} +{delay_h:g} h")
    missed_text = f"#{missed_index + 1} (t = {ref_events[missed_index].time_h:g} h): {missed_text}" if missed_on else None
    protocol = build_protocol(T, drug, model, active_events, dose_mg, tau, days, weight, window, metrics, ss,
                              message, missed_text)
    d = st.columns(2)
    d[0].download_button("⬇ " + T("download_csv"), df.to_csv(index=False).encode("utf-8"),
                         file_name=f"pk_{drug.id}.csv", mime="text/csv", width="stretch")
    d[1].download_button("⬇ " + T("download_protocol"), protocol.encode("utf-8"),
                         file_name=f"protocol_{drug.id}.md", mime="text/markdown", width="stretch")
    with st.expander(T("protocol_preview")):
        st.markdown(protocol)
    st.dataframe(df.iloc[:: max(1, len(df) // 500)], width="stretch", hide_index=True)

with tab_model:
    st.markdown(T("model_intro"))
    st.latex(r"C(t)=\frac{F\,D\,k_a}{V_d\,(k_a-k_e)}\left(e^{-k_e t}-e^{-k_a t}\right)")
    st.latex(r"C_{IV}(t)=\frac{D}{V_d}\,e^{-k_e t},\qquad k_e=\frac{\ln 2}{t_{1/2}},\quad CL=k_e V_d")
    st.latex(r"R=\frac{1}{1-e^{-k_e\tau}},\qquad C_{ss,avg}=\frac{F\,D}{CL\,\tau}")
    st.markdown(T("model_notes"))

    st.markdown(f"#### {T('validation_title')}")
    st.caption(T("validation_intro"))
    val_rows = validation_rows()
    val_cols = st.columns(3)
    for col, (group, g) in zip(val_cols, summarize(val_rows).items()):
        col.metric(T("validation_" + group), f"{g.within}/{g.n}",
                   help=T("validation_gm", gm=g.geometric_mean_fold_error))
    st.dataframe(
        [{T("drug"): r.drug_id, T("validation_metric"): r.metric, T("regimen"): regimen_label(r),
          T("validation_published"): r.observed, T("validation_predicted"): round(r.predicted, 5),
          T("validation_ratio"): round(r.ratio, 2), T("validation_ok"): "✅" if r.within else f"❌ >{r.max_fold:g}×",
          T("validation_source"): r.source, T("validation_note"): r.note}
         for r in val_rows],
        hide_index=True, width="stretch",
        column_config={T("validation_source"): st.column_config.LinkColumn(display_text="source")},
    )

st.divider()
st.caption("⚕️ " + T("disclaimer"))
