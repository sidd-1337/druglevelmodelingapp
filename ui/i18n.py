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
        "drug": "Drug", "body_weight": "Body weight", "conc_unit": "Concentration unit",
        "cv_cl": "Between-patient CV of clearance", "cv_f": "Between-patient CV of bioavailability",
        "renal_enable": "Enter renal function", "age": "Age", "years": "years", "sex_female": "Female",
        "creatinine": "Serum creatinine",
        "crcl_result": "CrCl (Cockcroft-Gault) {crcl:.0f} mL/min → clearance ×{factor:.2f} vs. reference adult.",
        "renal_not_used": "{name} is not cleared renally in this model; renal function does not change its levels.",
        "route": "Route", "route_oral": "Oral", "route_iv_bolus": "IV bolus",
        "route_insufflated": "Insufflated (nasal)", "route_smoked": "Smoked / Inhaled", "route_inhaled": "Inhaled",
        "dose": "Dose", "interval": "Dosing interval τ", "custom": "Custom", "duration": "Duration",
        "days_unit": "days", "doses_unit": "doses", "every": "every",
        "axis_scale": "Concentration axis", "linear": "Linear", "semilog": "Semi-log (log C)",
        "show_band": "Show 90 % between-patient band", "show_zones": "Show therapeutic zones",
        "missed_enable": "Simulate a missed dose", "missed_which": "Affected dose #",
        "missed_at": "Scheduled at t = {t:g} h", "missed_mode": "Patient behaviour",
        "mode_skip": "Skips the dose", "mode_delay": "Takes it late", "missed_delay": "Delay",
        # sukl
        "sukl_search_label": "Search Czech SÚKL registry",
        "sukl_search_placeholder": "e.g. Livizux, Rivotril, Paralen, Diazepam...",
        "sukl_found_count": "Found {n} registered packages in SÚKL",
        "sukl_select_pkg": "Select strength / package",
        "sukl_link": "Open official SÚKL registry detail ↗",
        "sukl_matched_preset": "Matched to clinical preset: {name}",
        "sukl_no_preset": "Package found in SÚKL; using custom model (verify t½ and Vd in SPC section 5.2).",
        "sukl_dose_set": "Dose set from package: {dose:g} mg",
        "tab_sukl": "SÚKL & Pricing",
        "card_dispensing": "Dispensing mode",
        "card_reimbursement": "Reimbursement",
        "card_supply_status": "Supply status",
        "sukl_status_rx": "Prescription (Rx)",
        "sukl_status_blue_stripe": "Prescription (Controlled)",
        "sukl_status_otc": "Over-the-counter (OTC)",
        "sukl_reimb_full": "Fully covered",
        "sukl_reimb_partial": "Partially covered",
        "sukl_reimb_none": "Not covered (Self-pay)",
        "sukl_shortage_active": "Supplied in CZ",
        "sukl_shortage_interrupted": "Supply interrupted",
        "sukl_shortage_ended": "Discontinued",
        "sukl_shortage_unknown": "No supply report",
        "sukl_detail_heading": "Official SÚKL Registry, Pricing & Reimbursement",
        "sukl_holder": "Authorization holder",
        "sukl_reg_no": "Registration no.",
        "reimb_zp": "Insurance reimbursement",
        "copay_patient": "Estimated copay",
        "max_price_pharmacy": "Max price in pharmacy",
        "reimb_limit_copay": "Copay to annual limit",
        "reimb_indications": "Reimbursement conditions",
        "reimb_doctor": "Authorized medical specialty",
        "download_spc": "Download official SPC (for physicians) ↗",
        "download_pil": "Download Patient Leaflet (PIL) ↗",
        "no_sukl_data": "No official SÚKL DLP detail available for this product.",
        # chart
        "zone_sub": "Sub-therapeutic (< MEC)", "zone_ther": "Therapeutic window", "zone_tox": "Toxic (> MTC)",
        "band": "90 % of patients", "reference": "Planned regimen", "conc": "Concentration",
        "doses": "Doses", "dose_marker": "Dose", "missed_marker": "Missed / delayed dose",
        "axis_time": "Time [h]", "axis_conc": "Concentration",
        # kpis / tabs
        "clearance": "Clearance", "time_to_ss": "Time to steady state",
        "tss_help": "≈ 4–5 × t½", "accumulation": "Accumulation index", "fluctuation": "Peak–trough fluctuation",
        "tab_ss": "Steady state", "tab_missed": "Missed dose", "tab_hr": "Harm reduction",
        "tab_glossary": "Glossary & Variables", "tab_data": "Data & export", "tab_model": "Model",
        "tab_fit": "Measured levels (TDM fit)", "measured": "Measured level",
        "sb_genotype": "Metaboliser status", "genotype_of": "{enzyme} phenotype",
        "phenotype_unknown": "Unknown", "phenotype_poor": "Poor metaboliser",
        "phenotype_intermediate": "Intermediate metaboliser", "phenotype_normal": "Normal metaboliser",
        "phenotype_rapid": "Rapid metaboliser", "phenotype_ultrarapid": "Ultrarapid metaboliser",
        "genotype_activity": "Enzyme activity {pct:.0f} % of normal.",
        "genotype_help": "Only enzymes this drug is actually cleared by are listed. Leaving one unknown "
                         "keeps the normal-metaboliser prediction but widens the between-patient band, "
                         "because such a patient could be any of the phenotypes.",
        "metabolite_curve": "{name} (metabolite)",
        "metabolite_caption": "{name} is an active metabolite with a half-life of about {t_half:g} h. "
                              "The concentration and therapeutic window above are {reported}.",
        "metabolite_reported_parent": "the parent drug only",
        "metabolite_reported_metabolite": "the metabolite",
        "metabolite_reported_sum": "the two combined (active moiety)",
        "formulation": "Formulation",
        "formulation_auto": "Selected from the SÚKL dosage form: {form}.",
        "tab_multi": "Co-medication",
        "multi_intro": "Add other drugs taken alongside. All of them are drawn on one chart, each relative to its "
                       "own therapeutic window. Where a quantified interaction is documented, the affected drug's "
                       "clearance or bioavailability is scaled; an enzyme that a drug of the regimen destroys is "
                       "simulated over time, so its effect follows that drug's levels.",
        "multi_overlay_title": "All drugs on one scale",
        "multi_overlay_caption": "Each curve is the drug's position in its own therapeutic window: 0 = MEC, 1 = MTC. "
                                 "Units and windows differ between drugs, so this is the only scale they share; "
                                 "hover for the concentration in the drug's own unit.",
        "multi_enzyme_title": "Shared enzyme activity",
        "multi_destroyed_by": "{enzyme} destroyed by {drug}",
        "multi_enzyme_caption": "Active enzyme as % of this patient's baseline ({pairs}). A destroyed enzyme comes "
                                "back only as new enzyme is made, so it recovers over days, not hours, and every "
                                "drug cleared by it is slowed accordingly at each moment.",
        "axis_window": "Position in therapeutic window (0 = MEC, 1 = MTC)",
        "axis_enzyme": "Active enzyme [% of baseline]",
        "multi_factor_sim": "Exposure vs. alone (AUC over the simulated period)",
        "multi_follows": "Time-varying: {enzymes} activity follows the levels of {drugs}.",
        "multi_curve": "{name}: curve alone and with co-medication",
        "pd_title": "Effect clashes (pharmacodynamic)",
        "pd_effect": "Effect", "pd_effect_serotonergic": "serotonergic", "pd_effect_cns_depressant": "CNS depressant",
        "pd_effect_qt_prolonging": "QT-prolonging", "pd_effect_sympathomimetic": "sympathomimetic",
        "pd_clash": "{effect}: {drugs} are above their MEC together for {hours:.1f} h (first at {first:.1f} h). "
                    "Their effects add up even though neither level changes.",
        "pd_clash_apart": "{effect}: {drugs} share this effect, but in this simulation they are not above their MEC "
                          "at the same time. Another schedule could bring them together.",
        "pd_none": "No effect class is shared by two assessed drugs of this regimen.",
        "pd_unassessed": "Not assessed for effect clashes: {drugs}.",
        "pd_sources": "Where each effect is stated",
        "multi_select": "Other preset drugs taken", "multi_perp": "Other substances (interaction sources)",
        "multi_none": "Add at least one more drug or substance to see interactions.",
        "multi_dose": "{name}: dose [mg]", "multi_interval": "{name}: interval [h]",
        "multi_factor": "Exposure vs. alone", "multi_no_interaction": "No quantified interaction documented "
                        "between these drugs. That is not proof there is none.",
        "multi_alone": "{name} alone", "multi_with": "{name} with co-medication",
        "multi_table_title": "Documented interactions",
        "multi_evidence": "Evidence", "multi_measured": "measured", "multi_estimated": "estimated",
        "multi_estimated_help": "Rows marked *estimated* have no published study for that pair. They are "
                                "calculated from the fraction of the drug cleared by the affected enzyme and "
                                "how much of that enzyme's activity the other drug leaves. Where a study does "
                                "exist, the measured value is used instead. Tested against the measured pairs "
                                "it can address, this calculation is accurate to about 1.25-fold.",
        "multi_perpetrator": "Caused by", "multi_effect": "Effect on exposure", "multi_site": "Acts on",
        "multi_onset": "Full effect after", "multi_site_clearance": "clearance", "multi_site_bioavailability":
        "bioavailability", "multi_qualitative": "direction only, no number",
        "multi_caption": "Interaction factors are the full effect of the stated perpetrator regimen and are applied "
                         "immediately; in reality induction builds over days to weeks. Only enzyme inactivation by a "
                         "drug of the regimen is simulated over time. Effect clashes are flagged from the labels "
                         "while both drugs are above their MEC, a proxy for both being active.",
        "validation_title": "Validation against published data",
        "validation_intro": "Each preset (70 kg adult, oral, normal renal function) is compared with published Cmax, "
                            "Tmax, AUC and steady-state levels from labels and clinical studies. Ratio = predicted / "
                            "published. Values were collected from search excerpts and still need checking against "
                            "the full documents.",
        "validation_exposure": "Exposure (AUC, C_ss) within tolerance", "validation_cmax": "C_max within tolerance",
        "validation_tmax": "T_max within tolerance", "validation_gm": "Geometric mean fold error {gm:.2f}×",
        "validation_metric": "Metric", "validation_published": "Published", "validation_predicted": "Model",
        "validation_ratio": "Ratio", "validation_ok": "Within tolerance", "validation_source": "Source",
        "validation_note": "Note",
        "fit_intro": "Enter this patient's measured blood levels. The model estimates the patient's own clearance "
                     "and volume (Bayesian MAP, with the population model as prior) and re-predicts the levels. "
                     "Times are hours since the first dose of the schedule above.",
        "fit_enable": "Fit the model to measured levels", "fit_n": "Number of samples",
        "fit_dose_log": "Doses actually taken", "fit_no_doses": "Enter at least one dose in the dose log.",
        "fit_dose_log_help": "Prefilled from the regimen above. Edit times and doses to match what the patient "
                             "actually took (late, skipped or extra doses); the fit uses this log.",
        "fit_time": "Sample {i}: time [h]", "fit_level": "Sample {i}: level [{unit}]",
        "fit_population": "Population model", "fit_individual": "This patient (fitted)",
        "fit_band": "90 % credible band (fitted)",
        "fit_cl": "Clearance (fitted)", "fit_cl_help": "Population: {pop:.3g} L/h",
        "fit_thalf": "t½ (fitted)", "fit_thalf_help": "Population: {pop:.3g} h",
        "fit_css_avg": "C_ss,avg (fitted)", "fit_css_min": "C_ss,min (fitted)",
        "fit_table_obs": "Measured", "fit_table_pop": "Population model", "fit_table_ind": "Fitted model",
        "fit_caption": "Assumes 20 % proportional error for assay and model. One sample mainly informs clearance; "
                       "the estimate is pulled toward the population value when the level is unusual. "
                       "Check that the dosing history above matches what the patient actually took.",
        "hr_heading": "Harm Reduction & Safe Redosing Calculator",
        "hr_safe_redosing_title": "Safe Redosing Interval (Accumulation Prevention)",
        "hr_elapsed_slider": "Time elapsed since previous dose [hours]",
        "hr_redose_amount": "Next planned dose [mg]",
        "hr_stacking_metric": "Predicted Stacked Peak",
        "hr_safe_interval_metric": "Minimum Safe Interval",
        "hr_current_level_metric": "Current Residual Level",
        "hr_routes_title": "Kinetics & Onset / Peak / Duration by Route of Administration",
        "hr_principles_title": "Harm Reduction Principles & Overdose Prevention",
        "ss_caption": "Analytical steady-state interval (superposition limit for t → ∞).",
        "missed_hint": "Enable “Simulate a missed dose” in the sidebar to compare against the planned regimen.",
        "min_after": "Lowest level after event", "time_below_mec": "Time below MEC",
        "time_above_mtc": "Time above MTC", "recovery": "Recovery to plan (±5 %)", "not_recovered": "not within horizon",
        # help tooltips
        "help_cmax": "Peak concentration reached in blood after dosing.",
        "help_tmax": "Time after dosing when peak concentration occurs.",
        "help_auc": "Total systemic drug exposure over 24 h (area under the curve).",
        "help_clearance": "Volume of blood completely cleared of drug per hour by kidneys and liver.",
        "help_thalf": "Elimination half-life: time needed for the body to eliminate 50 % of the drug.",
        "help_css_max": "Peak concentration at steady state during repeated dosing.",
        "help_css_min": "Trough concentration right before the next dose (must stay above MEC).",
        "help_css_avg": "Average concentration over the dosing interval at steady state.",
        "help_r": "Accumulation index: ratio of steady-state peak to single-dose peak.",
        "help_fluct": "Peak-to-trough swing percentage relative to the average level.",
        "help_tau": "Dosing interval (time between consecutive doses).",
        # glossary
        "glossary_intro": "Plain-language guide to pharmacokinetic variables and therapeutic drug monitoring.",
        "glossary_cmax_title": "C_max – Peak Concentration",
        "glossary_cmax_desc": "The highest drug concentration in the bloodstream after taking a dose. If C_max exceeds MTC, adverse effects or toxicity may occur.",
        "glossary_tmax_title": "T_max – Time to Peak",
        "glossary_tmax_desc": "How long it takes after swallowing the medicine for the blood concentration to reach its highest point. Rapidly absorbed drugs reach T_max in 0.5–2 hours, while extended-release formulations take several hours.",
        "glossary_thalf_title": "t½ – Elimination Half-Life",
        "glossary_thalf_desc": "The time it takes for your body (liver and kidneys) to remove 50% of the active drug. It directly determines how frequently a medication needs to be dosed.",
        "glossary_css_title": "C_ss – Steady-State Concentrations (Peak, Trough, Average)",
        "glossary_css_desc": "After taking a drug regularly for 4–5 half-lives, the amount entering the body equals the amount leaving it. At this 'steady state', levels fluctuate between C_ss,max (highest point) and C_ss,min (lowest point right before the next dose). C_ss,min must stay above MEC to maintain round-the-clock therapy.",
        "glossary_auc_title": "AUC – Area Under the Curve",
        "glossary_auc_desc": "The total area under the concentration-time curve. It represents the total systemic exposure of your body to the drug over time.",
        "glossary_cl_title": "Clearance (CL)",
        "glossary_cl_desc": "The volume of blood completely cleared of the drug per unit time (in liters per hour). Higher clearance means faster drug elimination.",
        "glossary_vd_title": "V_d – Volume of Distribution",
        "glossary_vd_desc": "A theoretical volume relating the total amount of drug in the body to its concentration in blood. A large V_d means the drug spreads deeply into body fat and tissues (e.g. brain, muscles), while a low V_d means it remains predominantly in the bloodstream.",
        "glossary_r_title": "R – Accumulation Index",
        "glossary_r_desc": "Shows how much higher the drug concentration becomes during chronic dosing compared to a single initial dose. If dosing interval τ is shorter than t½, the drug accumulates (R > 1).",
        "glossary_window_title": "Therapeutic Window: MEC & MTC",
        "glossary_window_desc": "MEC (Minimum Effective Concentration) is the threshold below which the drug loses efficacy. MTC (Maximum Tolerated Concentration) is the ceiling above which toxicity and adverse events occur. The optimal therapeutic window is the safe zone in between.",
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
                       "If ka ≈ ke the flip-flop limit C = F·D·ka·t·exp(−ke·t)/Vd is used.\n\n"
                       "**Patient covariates.** Volume scales with weight; clearance scales with weight^0.75, "
                       "and its renally cleared fraction (fe) with creatinine clearance relative to 100 mL/min: "
                       "CL = CL_ref·[(1−fe)·(WT/70)^0.75 + fe·CrCl/100]. t½ follows from V and CL.\n\n"
                       "**Variability band.** Clearance and bioavailability vary log-normally between patients "
                       "(published CVs); the band shows the 5th–95th percentile of exposure.\n\n"
                       "**Saturable elimination** (MDMA): CL(C) = (1−f)·CL + f·CL·Km/(Km + C), integrated "
                       "numerically, so repeated doses are not simple superposition.\n\n"
                       "**Validation.** Presets are checked against published Cmax, Tmax, AUC and steady-state "
                       "levels (data/reference_pk.json, `python -m core.validation`). A one-compartment model "
                       "matches exposure well but underestimates early peaks of drugs with a distribution phase.",
    },
    "cs": {
        "title": "Modelování hladin léčiv – PK/PD simulace a TDM",
        "subtitle": "Jednokompartmentová farmakokinetika s perorální absorpcí (Batemanova funkce), vícenásobným "
                    "dávkováním, ustáleným stavem a terapeutickým monitorováním léčiv.",
        "disclaimer": "Lékařské upozornění: pouze výukový a modelovací nástroj – nejde o certifikovaný zdravotnický "
                      "prostředek. Nepoužívejte jej pro rozhodování o dávkování skutečných pacientů.",
        "sb_drug": "Léčivo a pacient", "sb_params": "PK parametry (upravitelné)", "sb_regimen": "Dávkovací schéma",
        "sb_display": "Zobrazení", "sb_missed": "Analýza vynechané dávky",
        "drug": "Léčivo", "body_weight": "Tělesná hmotnost",
        "cv_cl": "Variabilita clearance mezi pacienty (CV)", "cv_f": "Variabilita biologické dostupnosti (CV)",
        "renal_enable": "Zadat funkci ledvin", "age": "Věk", "years": "let", "sex_female": "Žena",
        "creatinine": "Sérový kreatinin",
        "crcl_result": "CrCl (Cockcroft-Gault) {crcl:.0f} mL/min → clearance ×{factor:.2f} oproti referenčnímu dospělému.",
        "renal_not_used": "{name} se v tomto modelu nevylučuje ledvinami; funkce ledvin jeho hladiny nemění.",
        "conc_unit": "Jednotka koncentrace",
        "route": "Cesta podání", "route_oral": "Perorální", "route_iv_bolus": "IV bolus",
        "route_insufflated": "Nazální (šňupání)", "route_smoked": "Inhalační (kouření)", "route_inhaled": "Inhalační (vdechnutí)",
        "dose": "Dávka", "interval": "Dávkovací interval τ", "custom": "Vlastní", "duration": "Délka",
        "days_unit": "dnů", "doses_unit": "dávek", "every": "každých",
        "axis_scale": "Osa koncentrace", "linear": "Lineární", "semilog": "Semilogaritmická (log C)",
        "show_band": "Zobrazit 90% pásmo variability mezi pacienty", "show_zones": "Zobrazit terapeutická pásma",
        "missed_enable": "Simulovat vynechanou dávku", "missed_which": "Dotčená dávka č.",
        "missed_at": "Plánováno v čase t = {t:g} h", "missed_mode": "Chování pacienta",
        "mode_skip": "Dávku vynechá", "mode_delay": "Vezme ji později", "missed_delay": "Zpoždění",
        # sukl
        "sukl_search_label": "Hledat v databázi SÚKL ČR",
        "sukl_search_placeholder": "např. Livizux, Rivotril, Paralen, Diazepam...",
        "sukl_found_count": "Nalezeno {n} registrovaných balení v SÚKL",
        "sukl_select_pkg": "Vyberte sílu / balení",
        "sukl_link": "Otevřít detail přípravku v databázi SÚKL ↗",
        "sukl_matched_preset": "Přiřazen klinický model: {name}",
        "sukl_no_preset": "Balení nalezeno v SÚKL; použit vlastní model (ověřte t½ a Vd v SPC bod 5.2).",
        "sukl_dose_set": "Dávka nastavena z balení: {dose:g} mg",
        "tab_sukl": "SÚKL a úhrady",
        "card_dispensing": "Režim výdeje",
        "card_reimbursement": "Úhrada ZP",
        "card_supply_status": "Stav dodávek",
        "sukl_status_rx": "Na recept (Rx)",
        "sukl_status_blue_stripe": "Recept s modrým pruhem",
        "sukl_status_otc": "Volně prodejný (OTC)",
        "sukl_reimb_full": "Plně hrazeno pojišťovnou",
        "sukl_reimb_partial": "Částečně hrazeno",
        "sukl_reimb_none": "Nehrazeno (samoplátce)",
        "sukl_shortage_active": "Běžně dodáváno v ČR",
        "sukl_shortage_interrupted": "Přerušení dodávek (výpadek)",
        "sukl_shortage_ended": "Uvádění na trh ukončeno",
        "sukl_shortage_unknown": "Bez hlášení o dodávkách",
        "sukl_detail_heading": "Oficiální registr SÚKL, ceny a úhrady zdravotní pojišťovny",
        "sukl_holder": "Držitel registrace",
        "sukl_reg_no": "Registrační číslo",
        "reimb_zp": "Úhrada pojišťovny",
        "copay_patient": "Odhadovaný doplatek",
        "max_price_pharmacy": "Maximální cena v lékárně",
        "reimb_limit_copay": "Započitatelný doplatek do limitu",
        "reimb_indications": "Indikační omezení úhrady",
        "reimb_doctor": "Preskripční omezení (lékař)",
        "download_spc": "Stáhnout oficiální SPC (pro lékaře) ↗",
        "download_pil": "Stáhnout příbalový leták PIL ↗",
        "no_sukl_data": "Pro tento přípravek nejsou v SÚKL DLP dostupná detailní data.",
        "zone_sub": "Subterapeutické (< MEC)", "zone_ther": "Terapeutické okno", "zone_tox": "Toxické (> MTC)",
        "band": "90 % pacientů", "reference": "Plánované schéma", "conc": "Koncentrace",
        "doses": "Dávky", "dose_marker": "Dávka", "missed_marker": "Vynechaná / opožděná dávka",
        "axis_time": "Čas [h]", "axis_conc": "Koncentrace",
        "clearance": "Clearance", "time_to_ss": "Čas do ustáleného stavu",
        "tss_help": "≈ 4–5 × t½", "accumulation": "Index akumulace", "fluctuation": "Fluktuace peak–trough",
        "tab_ss": "Ustálený stav", "tab_missed": "Vynechaná dávka", "tab_hr": "Harm reduction kalkulačka",
        "tab_glossary": "Slovník veličin a pojmů", "tab_data": "Data a export", "tab_model": "Model",
        "tab_fit": "Změřené hladiny (TDM fit)", "measured": "Změřená hladina",
        "sb_genotype": "Metabolizační status", "genotype_of": "Fenotyp {enzyme}",
        "phenotype_unknown": "Neznámý", "phenotype_poor": "Pomalý metabolizátor",
        "phenotype_intermediate": "Střední metabolizátor", "phenotype_normal": "Normální metabolizátor",
        "phenotype_rapid": "Rychlý metabolizátor", "phenotype_ultrarapid": "Ultrarychlý metabolizátor",
        "genotype_activity": "Aktivita enzymu {pct:.0f} % normálu.",
        "genotype_help": "Uvedeny jsou jen enzymy, kterými se dané léčivo skutečně odbourává. Ponecháte-li "
                         "fenotyp neznámý, predikce odpovídá normálnímu metabolizátorovi, ale pásmo "
                         "variability se rozšíří, protože takový pacient může být kterýkoli z fenotypů.",
        "metabolite_curve": "{name} (metabolit)",
        "metabolite_caption": "{name} je aktivní metabolit s poločasem přibližně {t_half:g} h. "
                              "Koncentrace a terapeutické okno výše se týkají {reported}.",
        "metabolite_reported_parent": "pouze mateřské látky",
        "metabolite_reported_metabolite": "metabolitu",
        "metabolite_reported_sum": "jejich součtu (aktivní frakce)",
        "formulation": "Léková forma",
        "formulation_auto": "Vybráno podle lékové formy v SÚKL: {form}.",
        "tab_multi": "Souběžná medikace",
        "multi_intro": "Přidejte další souběžně užívaná léčiva. Všechna se vykreslí do jednoho grafu, každé vůči "
                       "vlastnímu terapeutickému oknu. Je-li doložena kvantifikovaná interakce, upraví se clearance "
                       "nebo biologická dostupnost ovlivněného léčiva; enzym, který některé léčivo režimu ničí, se "
                       "simuluje v čase, takže jeho účinek sleduje hladiny tohoto léčiva.",
        "multi_overlay_title": "Všechna léčiva na jedné škále",
        "multi_overlay_caption": "Každá křivka ukazuje polohu léčiva v jeho vlastním terapeutickém okně: 0 = MEC, "
                                 "1 = MTC. Jednotky i okna se mezi léčivy liší, proto je to jediná společná škála; "
                                 "po najetí myší se zobrazí koncentrace ve vlastní jednotce léčiva.",
        "multi_enzyme_title": "Aktivita sdíleného enzymu",
        "multi_destroyed_by": "{enzyme} ničí {drug}",
        "multi_enzyme_caption": "Aktivní enzym v % výchozí hodnoty tohoto pacienta ({pairs}). Zničený enzym se "
                                "obnovuje jen tvorbou nového, tedy během dnů, ne hodin, a každé léčivo, které ho "
                                "odbourává, se v každém okamžiku odpovídajícím způsobem zpomalí.",
        "axis_window": "Poloha v terapeutickém okně (0 = MEC, 1 = MTC)",
        "axis_enzyme": "Aktivní enzym [% výchozí hodnoty]",
        "multi_factor_sim": "Expozice oproti samostatnému užití (AUC za simulované období)",
        "multi_follows": "Časově proměnné: aktivita {enzymes} sleduje hladiny {drugs}.",
        "multi_curve": "{name}: křivka samostatně a se souběžnou medikací",
        "pd_title": "Střety účinků (farmakodynamické)",
        "pd_effect": "Účinek", "pd_effect_serotonergic": "serotonergní", "pd_effect_cns_depressant": "tlumivý na CNS",
        "pd_effect_qt_prolonging": "prodlužující QT", "pd_effect_sympathomimetic": "sympatomimetický",
        "pd_clash": "{effect}: {drugs} jsou společně nad MEC po dobu {hours:.1f} h (poprvé v {first:.1f} h). "
                    "Jejich účinky se sčítají, i když se žádná hladina nezmění.",
        "pd_clash_apart": "{effect}: {drugs} sdílejí tento účinek, ale v této simulaci nejsou nad MEC současně. "
                          "Jiné rozvržení dávek je může svést dohromady.",
        "pd_none": "Žádná dvě posouzená léčiva tohoto režimu nesdílejí stejnou třídu účinku.",
        "pd_unassessed": "Neposouzeno na střety účinků: {drugs}.",
        "pd_sources": "Kde je který účinek uveden",
        "multi_select": "Další užívaná léčiva z presetů", "multi_perp": "Další látky (zdroje interakcí)",
        "multi_none": "Přidejte alespoň jedno další léčivo nebo látku, aby se zobrazily interakce.",
        "multi_dose": "{name}: dávka [mg]", "multi_interval": "{name}: interval [h]",
        "multi_factor": "Expozice oproti samostatnému užití",
        "multi_no_interaction": "Mezi těmito léčivy není doložena kvantifikovaná interakce. To neznamená, že žádná neexistuje.",
        "multi_alone": "{name} samostatně", "multi_with": "{name} se souběžnou medikací",
        "multi_table_title": "Doložené interakce",
        "multi_evidence": "Doklad", "multi_measured": "změřeno", "multi_estimated": "odhad",
        "multi_estimated_help": "Řádky označené *odhad* nemají pro danou dvojici publikovanou studii. Počítají "
                                "se z podílu léčiva odbouraného dotčeným enzymem a z toho, kolik aktivity "
                                "tohoto enzymu druhé léčivo ponechá. Tam, kde studie existuje, se použije "
                                "změřená hodnota. Proti dvojicím, které umí popsat, je tento výpočet přesný "
                                "zhruba na 1,25násobek.",
        "multi_perpetrator": "Způsobeno", "multi_effect": "Vliv na expozici", "multi_site": "Působí na",
        "multi_onset": "Plný účinek za", "multi_site_clearance": "clearance", "multi_site_bioavailability":
        "biologickou dostupnost", "multi_qualitative": "pouze směr, bez čísla",
        "multi_caption": "Faktory interakcí odpovídají plnému účinku uvedeného režimu vyvolávajícího léčiva a aplikují "
                         "se okamžitě; ve skutečnosti se indukce rozvíjí dny až týdny. V čase se simuluje jen "
                         "inaktivace enzymu léčivem z režimu. Střety účinků se hlásí podle SPC, dokud jsou obě "
                         "léčiva nad MEC, což zastupuje stav, kdy obě působí.",
        "validation_title": "Validace proti publikovaným datům",
        "validation_intro": "Každý preset (dospělý 70 kg, perorálně, normální funkce ledvin) je porovnán s publikovanými "
                            "hodnotami Cmax, Tmax, AUC a ustálených hladin z SPC a klinických studií. Poměr = model / "
                            "publikováno. Hodnoty byly převzaty z úryvků vyhledávání a je třeba je ještě ověřit "
                            "v úplných dokumentech.",
        "validation_exposure": "Expozice (AUC, C_ss) v toleranci", "validation_cmax": "C_max v toleranci",
        "validation_tmax": "T_max v toleranci", "validation_gm": "Geometrický průměr chyby {gm:.2f}×",
        "validation_metric": "Veličina", "validation_published": "Publikováno", "validation_predicted": "Model",
        "validation_ratio": "Poměr", "validation_ok": "V toleranci", "validation_source": "Zdroj",
        "validation_note": "Poznámka",
        "fit_intro": "Zadejte změřené hladiny pacienta. Model odhadne jeho vlastní clearance a distribuční objem "
                     "(bayesovský MAP odhad, populační model jako apriorní informace) a znovu předpoví hladiny. "
                     "Časy jsou v hodinách od první dávky výše zadaného schématu.",
        "fit_enable": "Napasovat model na změřené hladiny", "fit_n": "Počet vzorků",
        "fit_dose_log": "Skutečně užité dávky", "fit_no_doses": "Zadejte do záznamu alespoň jednu dávku.",
        "fit_dose_log_help": "Předvyplněno podle schématu výše. Upravte časy a dávky podle toho, co pacient skutečně "
                             "užil (pozdě, vynecháno, navíc); fit vychází z tohoto záznamu.",
        "fit_time": "Vzorek {i}: čas [h]", "fit_level": "Vzorek {i}: hladina [{unit}]",
        "fit_population": "Populační model", "fit_individual": "Tento pacient (napasováno)",
        "fit_band": "90% interval věrohodnosti (napasováno)",
        "fit_cl": "Clearance (napasováno)", "fit_cl_help": "Populační: {pop:.3g} L/h",
        "fit_thalf": "t½ (napasováno)", "fit_thalf_help": "Populační: {pop:.3g} h",
        "fit_css_avg": "C_ss,avg (napasováno)", "fit_css_min": "C_ss,min (napasováno)",
        "fit_table_obs": "Změřeno", "fit_table_pop": "Populační model", "fit_table_ind": "Napasovaný model",
        "fit_caption": "Předpokládá 20% proporcionální chybu stanovení a modelu. Jeden vzorek informuje hlavně "
                       "o clearance; u neobvyklé hladiny je odhad tažen k populační hodnotě. Ověřte, že zadané "
                       "dávkování odpovídá tomu, co pacient skutečně užil.",
        "hr_heading": "Kalkulačka pro bezpečné dávkování a Harm Reduction",
        "hr_safe_redosing_title": "Bezpečný interval pro další dávku (Prevence akumulace)",
        "hr_elapsed_slider": "Uplynulý čas od poslední dávky [hodiny]",
        "hr_redose_amount": "Plánovaná další dávka [mg]",
        "hr_stacking_metric": "Předpokládaný kumulativní vrchol",
        "hr_safe_interval_metric": "Doporučený bezpečný odstup",
        "hr_current_level_metric": "Aktuální zbytková hladina",
        "hr_routes_title": "Kinetika a časové osy účinku podle cesty podání",
        "hr_principles_title": "Zásady Harm Reduction a prevence předávkování",
        "ss_caption": "Analytický interval v ustáleném stavu (limita superpozice pro t → ∞).",
        "missed_hint": "Zapněte „Simulovat vynechanou dávku“ v postranním panelu a porovnejte s plánovaným schématem.",
        "min_after": "Nejnižší hladina po události", "time_below_mec": "Čas pod MEC",
        "time_above_mtc": "Čas nad MTC", "recovery": "Návrat k plánu (±5 %)", "not_recovered": "v horizontu ne",
        # help tooltips
        "help_cmax": "Maximální dosažená hladina v krvi po podání dávky (vrchol křivky).",
        "help_tmax": "Doba od podání dávky do dosažení nejvyšší hladiny v krvi.",
        "help_auc": "Celková celodenní expozice organismu léku (plocha pod křivkou za 24 h).",
        "help_clearance": "Objem krve, který tělo (játra a ledviny) kompletně očistí od léku za hodinu.",
        "help_thalf": "Biologický poločas: doba, za kterou tělo vyloučí 50 % léčiva.",
        "help_css_max": "Maximální hladina (vrchol) při pravidelném užívání v ustáleném stavu.",
        "help_css_min": "Minimální (údolní) hladina těsně před další dávkou (musí zůstat nad MEC).",
        "help_css_avg": "Průměrná hladina v ustáleném stavu během dávkovacího intervalu.",
        "help_r": "Index akumulace: kolikrát je hladina v ustáleném stavu vyšší než po 1. dávce.",
        "help_fluct": "Procento kolísání mezi vrcholem a údolím vzhledem k průměru.",
        "help_tau": "Dávkovací interval (čas mezi jednotlivými dávkami).",
        # glossary
        "glossary_intro": "Přehledné vysvětlení farmakokinetických veličin a monitorování léčiv pro běžného uživatele.",
        "glossary_cmax_title": "C_max – Maximální koncentrace (vrchol)",
        "glossary_cmax_desc": "Nejvyšší dosažená hladina léku v krvi po podání dávky (vrchol křivky). Pokud C_max přesáhne horní mez MTC, hrozí nežádoucí účinky či toxicita.",
        "glossary_tmax_title": "T_max – Čas do dosažení maxima",
        "glossary_tmax_desc": "Doba od spolknutí léku, než koncentrace v krvi dosáhne svého vrcholu. Rychle vstřebatelné léky dosahují maxima za 0,5–2 hodiny, léky s řízeným uvolňováním za déle.",
        "glossary_thalf_title": "t½ – Biologický poločas eliminace",
        "glossary_thalf_desc": "Doba, za kterou tělo (játra a ledviny) vyloučí přesně 50 % účinné látky z těla. Přímo určuje, jak často se lék musí užívat.",
        "glossary_css_title": "C_ss – Koncentrace v ustáleném stavu (maximum, údolí, průměr)",
        "glossary_css_desc": "Po pravidelném užívání po dobu přibližně 4 až 5 poločasů nastane rovnováha: rychlost přísunu léku se rovná rychlosti vylučování. Hladina pravidelně osciluje mezi C_ss,max (vrchol) a C_ss,min (údolí těsně před další dávkou). C_ss,min musí zůstat nad MEC, aby lék nepřetržitě účinkoval.",
        "glossary_auc_title": "AUC – Plocha pod křivkou",
        "glossary_auc_desc": "Celková plocha pod křivkou závislosti koncentrace na čase. Vyjadřuje celkovou míru expozice organismu účinné látce.",
        "glossary_cl_title": "Clearance (CL) – Očišťovací schopnost těla",
        "glossary_cl_desc": "Objem krve, který tělo dokáže kompletně zbavit léku za jednotku času (v litrech za hodinu). Vyšší clearance znamená rychlejší odbourávání léku.",
        "glossary_vd_title": "V_d – Zdánlivý distribuční objem",
        "glossary_vd_desc": "Teoretický objem tělesné tekutiny, do kterého by se lék musel rovnoměrně rozptýlit. Velký V_d znamená, že lék prostupuje hluboko do tkání a tuku (např. do mozku), malý V_d znamená, že zůstává hlavně v krevním řečišti.",
        "glossary_r_title": "R – Index akumulace",
        "glossary_r_desc": "Ukazuje, kolikrát je hladina při dlouhodobém užívání vyšší než po jedné jediné první dávce. Pokud berete lék častěji, než je jeho poločas t½, lék se v těle hromadí (R > 1).",
        "glossary_window_title": "Terapeutické okno: MEC a MTC",
        "glossary_window_desc": "MEC (minimální účinná koncentrace) je spodní mez, pod kterou lék přestává zabírat. MTC (maximální tolerovaná koncentrace) je horní mez, nad kterou hrozí toxicita a nežádoucí účinky. Cílem je udržet hladinu bezpečně mezi nimi.",
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
                       "Pro ka ≈ ke se používá limita C = F·D·ka·t·exp(−ke·t)/Vd.\n\n"
                       "**Kovariáty pacienta.** Objem roste s hmotností; clearance s hmotností^0,75 a její "
                       "ledvinami vylučovaný podíl (fe) s clearance kreatininu vůči 100 mL/min: "
                       "CL = CL_ref·[(1−fe)·(WT/70)^0,75 + fe·CrCl/100]. t½ vyplývá z V a CL.\n\n"
                       "**Pásmo variability.** Clearance a biologická dostupnost se mezi pacienty liší log-normálně "
                       "(publikované CV); pásmo ukazuje 5.–95. percentil expozice.\n\n"
                       "**Saturovatelná eliminace** (MDMA): CL(C) = (1−f)·CL + f·CL·Km/(Km + C), řešeno numericky, "
                       "takže opakované dávky nejsou prostým součtem.\n\n"
                       "**Validace.** Presety jsou porovnány s publikovanými Cmax, Tmax, AUC a ustálenými hladinami "
                       "(data/reference_pk.json, `python -m core.validation`). Jednokompartmentový model dobře "
                       "odpovídá expozici, ale podhodnocuje časné vrcholy léčiv s distribuční fází.",
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
