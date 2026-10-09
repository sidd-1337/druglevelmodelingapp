"""Streamlit app: interactive PK simulation & Therapeutic Drug Monitoring (TDM)."""
from __future__ import annotations

import math

import numpy as np
import streamlit as st

from core.drug_database import CUSTOM_ID, custom_drug, load_drugs
from core.metrics import (
    TherapeuticWindow,
    assess_tdm,
    compute_metrics,
    recovery_time,
)
from core.pk_engine import (
    Route,
    apply_missed_dose,
    concentration_profile,
    regular_schedule,
    steady_state_concentration,
    steady_state_summary,
    time_grid,
)
from core.sukl_client import (
    PRESET_DEFAULT_SUKL_CODES,
    fetch_sukl_detail,
    query_sukl_api,
)
from ui.charts import build_figure
from ui.i18n import LANGUAGES, make_translator
from ui.reporting import build_dataframe, build_protocol

st.set_page_config(page_title="PK/PD Simulator & TDM", page_icon="💊", layout="wide")


@st.cache_data(ttl=3600, show_spinner=False)
def search_sukl(query: str):
    return query_sukl_api(query)


@st.cache_data(ttl=3600, show_spinner=False)
def get_sukl_detail(kod_sukl: str):
    return fetch_sukl_detail(kod_sukl)


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


# ------------------------------------------------------------------ sidebar
with st.sidebar:
    lang = st.radio("Language / Jazyk", list(LANGUAGES), horizontal=True)
    T = make_translator(LANGUAGES[lang])

    st.header(T("sb_drug"))

    sukl_query = st.text_input(
        T("sukl_search_label"),
        value="",
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
    if sukl_med and sukl_med.matched_preset_id and sukl_med.matched_preset_id in by_id:
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
    if base.notes:
        st.caption(base.notes)

    weight = st.slider(T("body_weight") + " [kg]", 30, 150, 70, key="weight")

    with st.expander(T("sb_params")):
        t_half = st.number_input("t½ [h]", 0.1, 500.0, base.t_half_h, 0.1, key=f"th_{k}")
        band_lo = st.number_input(T("t_half_min") + " [h]", 0.1, 500.0, base.half_life_band_h[0], 0.1, key=f"tl_{k}")
        band_hi = st.number_input(T("t_half_max") + " [h]", 0.1, 500.0, base.half_life_band_h[1], 0.1, key=f"tu_{k}")
        vd_kg = st.number_input("V_d [L/kg]", 0.01, 20.0, base.vd_l_per_kg, 0.01, key=f"vd_{k}")
        f_bio = st.slider("F", 0.05, 1.0, base.bioavailability, 0.01, key=f"f_{k}")
        ka = st.number_input("k_a [1/h]", 0.05, 20.0, base.ka_per_h, 0.05, key=f"ka_{k}")
        unit = st.text_input(T("conc_unit"), base.concentration_unit, key=f"u_{k}") if drug_id == CUSTOM_ID \
            else base.concentration_unit
        mec = st.number_input(f"MEC [{unit}]", 0.0001, 10000.0, base.mec, format="%g", key=f"mec_{k}")
        mtc = st.number_input(f"MTC [{unit}]", 0.0001, 10000.0, base.mtc, format="%g", key=f"mtc_{k}")

    st.header(T("sb_regimen"))
    default_route_idx = 0
    if sukl_med and any(t in (sukl_med.cesta_podani or "").lower() for t in ["intraven", "injek", "infuz", "i.v."]):
        default_route_idx = 1
    route = st.radio(
        T("route"),
        list(Route),
        index=default_route_idx,
        format_func=lambda r: T("route_" + r.value),
        horizontal=True,
        key=f"route_{k}_{sukl_med.kod_sukl if sukl_med else 'none'}",
    )
    suggested_dose = float(sukl_med.parsed_dose_mg) if sukl_med and sukl_med.parsed_dose_mg is not None else float(base.default_dose_mg)
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
    default_tau = base.default_interval_h
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
        t_half_h=t_half, t_half_range_h=(min(band_lo, t_half), max(band_hi, t_half)), vd_l_per_kg=vd_kg,
        bioavailability=f_bio, ka_per_h=ka, mec=mec, mtc=mtc, concentration_unit=unit)
    window = TherapeuticWindow(mec, mtc)
except ValueError:
    st.error(T("err_window"))
    st.stop()

amount = drug.amount_from_dose_mg(dose_mg)
model = drug.to_pk_model(weight, route)
ref_events = regular_schedule(amount, tau, duration_h)
active_events = apply_missed_dose(ref_events, missed_index, delay_h) if missed_on else ref_events
t = time_grid(duration_h)
c_ref = concentration_profile(t, ref_events, model)
c = concentration_profile(t, active_events, model) if missed_on else c_ref
lo_t, hi_t = drug.half_life_band_h
curves = [concentration_profile(t, active_events, model.with_half_life(h)) for h in (lo_t, hi_t)]
band = (np.minimum(*curves), np.maximum(*curves))

metrics = compute_metrics(t, c, model)
ss = steady_state_summary(amount, tau, model)
assess_from = max(0.0, min(ss.time_to_ss_h, duration_h - tau))
assessment = assess_tdm(t, c, window, assess_from)
level, message = assessment_message(T, assessment, window, unit)

# ----------------------------------------------------------------- main page
st.title("💊 " + T("title"))
st.caption(T("subtitle"))

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
cols[4].metric("t½", f"{metrics.t_half_h:g} h", help=T("help_thalf"))
cols[5].metric(T("time_to_ss"), f"{ss.time_to_ss_h:.0f} h", help=T("tss_help"))

show_alert(level, message)
if duration_h < ss.time_to_ss_h:
    st.info(T("info_short_sim", tss=ss.time_to_ss_h, days=math.ceil(ss.time_to_ss_h / 24)))

missed_time = ref_events[missed_index].time_h if missed_on else None
fig = build_figure(
    t, c, window, unit, T, log_scale=log_scale, band=band if show_band else None,
    reference=c_ref if missed_on else None, dose_times=[e.time_h for e in active_events],
    missed_time=missed_time, show_zones=show_zones)
st.plotly_chart(fig, width="stretch")

tab_ss, tab_missed, tab_sukl, tab_glossary, tab_data, tab_model = st.tabs(
    [T("tab_ss"), T("tab_missed"), T("tab_sukl"), T("tab_glossary"), T("tab_data"), T("tab_model")])

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

st.divider()
st.caption("⚕️ " + T("disclaimer"))
