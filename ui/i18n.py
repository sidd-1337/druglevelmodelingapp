"""Minimal EN/CS translation layer for the UI."""
from __future__ import annotations

from typing import Callable

LANGUAGES = {"English": "en", "Čeština": "cs"}

_TEXT: dict[str, dict[str, str]] = {
    "en": {
        "title": "Drug Level Modeling – PK/PD Simulation & TDM",
        "subtitle": "One-compartment pharmacokinetics with oral absorption (Bateman), multiple dosing, "
                    "steady state and therapeutic drug monitoring.",
        "disclaimer": "Medical disclaimer: educational and modelling tool only – not a certified medical device. "
                      "Do not use it to make dosing decisions for real patients.",
        # sidebar
        "sb_drug": "Drug & patient", "sb_params": "PK parameters (editable)", "sb_regimen": "Dosing regimen",
        "sb_display": "Display", "sb_missed": "Missed dose analysis",
        "drug": "Drug", "body_weight": "Body weight", "t_half_min": "t½ min (variability)",
        "t_half_max": "t½ max (variability)", "conc_unit": "Concentration unit",
        "route": "Route", "route_oral": "Oral", "route_iv_bolus": "IV bolus",
        "dose": "Dose", "interval": "Dosing interval τ", "custom": "Custom", "duration": "Duration",
        "days_unit": "days", "doses_unit": "doses", "every": "every",
        "axis_scale": "Concentration axis", "linear": "Linear", "semilog": "Semi-log (log C)",
        "show_band": "Show t½ variability band", "show_zones": "Show therapeutic zones",
        "missed_enable": "Simulate a missed dose", "missed_which": "Affected dose #",
        "missed_at": "Scheduled at t = {t:g} h", "missed_mode": "Patient behaviour",
        "mode_skip": "Skips the dose", "mode_delay": "Takes it late", "missed_delay": "Delay",
        # chart
        "zone_sub": "Sub-therapeutic (< MEC)", "zone_ther": "Therapeutic window", "zone_tox": "Toxic (> MTC)",
        "band": "t½ variability", "reference": "Planned regimen", "conc": "Concentration",
        "doses": "Doses", "dose_marker": "Dose", "missed_marker": "Missed / delayed dose",
        "axis_time": "Time [h]", "axis_conc": "Concentration",
        # kpis / tabs
        "clearance": "Clearance", "time_to_ss": "Time to steady state",
        "tss_help": "≈ 4–5 × t½", "accumulation": "Accumulation index", "fluctuation": "Peak–trough fluctuation",
        "tab_ss": "Steady state", "tab_missed": "Missed dose", "tab_data": "Data & export", "tab_model": "Model",
        "ss_caption": "Analytical steady-state interval (superposition limit for t → ∞).",
        "missed_hint": "Enable “Simulate a missed dose” in the sidebar to compare against the planned regimen.",
        "min_after": "Lowest level after event", "time_below_mec": "Time below MEC",
        "time_above_mtc": "Time above MTC", "recovery": "Recovery to plan (±5 %)", "not_recovered": "not within horizon",
        # alerts
        "alert_ok": "✅ The regimen keeps the concentration stable inside the therapeutic window.",
        "alert_toxic": "⚠️ Warning: the estimated level exceeds the safe limit MTC ({mtc:g} {unit}) at T = {t:.1f} h! "
                       "Time above MTC: {hours:.1f} h, peak {cmax:.3g} {unit}.",
        "alert_never": "⚠️ The concentration never reaches MEC ({mec:g} {unit}); peak is only {cmax:.3g} {unit} "
                       "– the regimen is likely ineffective.",
        "alert_trough": "⚠️ Trough levels ({trough:.3g} {unit}) fall below MEC ({mec:g} {unit}) "
                        "(assessed from t = {t0:.0f} h) – risk of sub-therapeutic exposure.",
        "info_short_sim": "The simulation is shorter than the time to steady state (≈ {tss:.0f} h). "
                          "Extend the duration to at least {days} days to evaluate steady state.",
        "err_window": "Invalid parameters: MEC must be lower than MTC and all values must be positive.",
        # exports
        "download_csv": "Download CSV", "download_protocol": "Download dosing protocol",
        "protocol_preview": "Preview dosing protocol", "protocol_title": "Dosing protocol (simulation summary)",
        "protocol_generated": "Generated {d} by the PK/PD simulator", "protocol_patient_regimen": "Patient & regimen",
        "protocol_params": "PK parameters", "protocol_results": "Results", "protocol_assessment": "TDM assessment",
        "regimen": "Regimen", "dose_times": "Dose times", "missed_scenario": "Missed-dose scenario",
        "parameter": "Parameter", "value": "Value",
        # model
        "model_intro": "**One-compartment model, first-order absorption and elimination.**",
        "model_notes": "Multiple doses use the superposition principle (linear, time-invariant kinetics). "
                       "IV bolus has F = 1 and no absorption phase. Concentrations in mg/L equal µg/mL. "
                       "If ka ≈ ke the flip-flop limit C = F·D·ka·t·exp(−ke·t)/Vd is used.",
    },
    "cs": {
        "title": "Modelování hladin léčiv – PK/PD simulace a TDM",
        "subtitle": "Jednokompartmentová farmakokinetika s perorální absorpcí (Batemanova funkce), vícenásobným "
                    "dávkováním, ustáleným stavem a terapeutickým monitorováním léčiv.",
        "disclaimer": "Lékařské upozornění: pouze výukový a modelovací nástroj – nejde o certifikovaný zdravotnický "
                      "prostředek. Nepoužívejte jej pro rozhodování o dávkování skutečných pacientů.",
        "sb_drug": "Léčivo a pacient", "sb_params": "PK parametry (upravitelné)", "sb_regimen": "Dávkovací schéma",
        "sb_display": "Zobrazení", "sb_missed": "Analýza vynechané dávky",
        "drug": "Léčivo", "body_weight": "Tělesná hmotnost", "t_half_min": "t½ min (variabilita)",
        "t_half_max": "t½ max (variabilita)", "conc_unit": "Jednotka koncentrace",
        "route": "Cesta podání", "route_oral": "Perorální", "route_iv_bolus": "IV bolus",
        "dose": "Dávka", "interval": "Dávkovací interval τ", "custom": "Vlastní", "duration": "Délka",
        "days_unit": "dnů", "doses_unit": "dávek", "every": "každých",
        "axis_scale": "Osa koncentrace", "linear": "Lineární", "semilog": "Semilogaritmická (log C)",
        "show_band": "Zobrazit pásmo variability t½", "show_zones": "Zobrazit terapeutická pásma",
        "missed_enable": "Simulovat vynechanou dávku", "missed_which": "Dotčená dávka č.",
        "missed_at": "Plánováno v čase t = {t:g} h", "missed_mode": "Chování pacienta",
        "mode_skip": "Dávku vynechá", "mode_delay": "Vezme ji později", "missed_delay": "Zpoždění",
        "zone_sub": "Subterapeutické (< MEC)", "zone_ther": "Terapeutické okno", "zone_tox": "Toxické (> MTC)",
        "band": "Variabilita t½", "reference": "Plánované schéma", "conc": "Koncentrace",
        "doses": "Dávky", "dose_marker": "Dávka", "missed_marker": "Vynechaná / opožděná dávka",
        "axis_time": "Čas [h]", "axis_conc": "Koncentrace",
        "clearance": "Clearance", "time_to_ss": "Čas do ustáleného stavu",
        "tss_help": "≈ 4–5 × t½", "accumulation": "Index akumulace", "fluctuation": "Fluktuace peak–trough",
        "tab_ss": "Ustálený stav", "tab_missed": "Vynechaná dávka", "tab_data": "Data a export", "tab_model": "Model",
        "ss_caption": "Analytický interval v ustáleném stavu (limita superpozice pro t → ∞).",
        "missed_hint": "Zapněte „Simulovat vynechanou dávku“ v postranním panelu a porovnejte s plánovaným schématem.",
        "min_after": "Nejnižší hladina po události", "time_below_mec": "Čas pod MEC",
        "time_above_mtc": "Čas nad MTC", "recovery": "Návrat k plánu (±5 %)", "not_recovered": "v horizontu ne",
        "alert_ok": "✅ Dávkování udržuje stabilní koncentraci v terapeutickém okně.",
        "alert_toxic": "⚠️ Pozor: Odhadovaná hladina překračuje bezpečný limit MTC ({mtc:g} {unit}) v čase "
                       "T={t:.1f} h! Čas nad MTC: {hours:.1f} h, maximum {cmax:.3g} {unit}.",
        "alert_never": "⚠️ Koncentrace nikdy nedosáhne MEC ({mec:g} {unit}); maximum je jen {cmax:.3g} {unit} "
                       "– schéma je pravděpodobně neúčinné.",
        "alert_trough": "⚠️ Údolní hladiny ({trough:.3g} {unit}) klesají pod MEC ({mec:g} {unit}) "
                        "(hodnoceno od t = {t0:.0f} h) – riziko subterapeutické expozice.",
        "info_short_sim": "Simulace je kratší než čas do ustáleného stavu (≈ {tss:.0f} h). "
                          "Pro vyhodnocení ustáleného stavu prodlužte délku alespoň na {days} dnů.",
        "err_window": "Neplatné parametry: MEC musí být menší než MTC a všechny hodnoty musí být kladné.",
        "download_csv": "Stáhnout CSV", "download_protocol": "Stáhnout dávkovací protokol",
        "protocol_preview": "Náhled dávkovacího protokolu", "protocol_title": "Dávkovací protokol (souhrn simulace)",
        "protocol_generated": "Vygenerováno {d} PK/PD simulátorem", "protocol_patient_regimen": "Pacient a schéma",
        "protocol_params": "PK parametry", "protocol_results": "Výsledky", "protocol_assessment": "Hodnocení TDM",
        "regimen": "Schéma", "dose_times": "Časy dávek", "missed_scenario": "Scénář vynechané dávky",
        "parameter": "Parametr", "value": "Hodnota",
        "model_intro": "**Jednokompartmentový model, absorpce a eliminace 1. řádu.**",
        "model_notes": "Vícenásobné dávkování využívá princip superpozice (lineární, časově invariantní kinetika). "
                       "IV bolus má F = 1 a žádnou absorpční fázi. Koncentrace v mg/L odpovídá µg/mL. "
                       "Pro ka ≈ ke se používá limita C = F·D·ka·t·exp(−ke·t)/Vd.",
    },
}


def make_translator(lang: str) -> Callable[..., str]:
    """Return ``T(key, **fmt)``; falls back to English for missing keys."""
    table = _TEXT.get(lang, _TEXT["en"])

    def T(key: str, **fmt) -> str:
        text = table.get(key, _TEXT["en"].get(key))
        if text is None:
            raise KeyError(f"missing translation key: {key}")
        return text.format(**fmt) if fmt else text

    return T
