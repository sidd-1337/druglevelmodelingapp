# Quetiapine (kvetiapin)

> Educational summary for the drug-level modelling app. Not clinical guidance.

## 1. Identity

| | |
|---|---|
| INN / Czech | quetiapine / kvetiapin |
| ATC | N05AH04 (diazepines, oxazepines, thiazepines) |
| Class | atypical antipsychotic, dibenzothiazepine |
| Salt | fumarate (2:1); every strength is stated as base [1, 3] |
| Molar mass | base C21H25N3O2S 383.51 g/mol; norquetiapine C17H17N3S 295.40 g/mol |
| Czech market | IR film-coated tablets 25–400 mg; prolonged-release (PR) tablets 50–400 mg (Quetiapin Teva Retard, Quetiapine Accord, Quetiapine Viatris) |
| Preset | `quetiapine`, default product QUETIAPIN NEURAXPHARM 200MG TBL FLM 60 (`0249090`) |
| Usual dose | schizophrenia 300–450 mg/day in two doses, range 150–750 [3]; PR once daily |

## 2. Mechanism of action

Quetiapine is an antagonist at many receptors at once. Its antipsychotic effect is
attributed to 5-HT2A plus D2 antagonism with higher selectivity for 5-HT2A; the label
calls the mechanism "unclear" [1, 3]. Its active metabolite norquetiapine adds
noradrenaline-transporter (NET) inhibition and 5-HT1A partial agonism, which the Czech
SPC names as the likely basis of the antidepressant effect [3, 7].

**Binding** (Ki, nM; lower = tighter). Parent and metabolite from the US label [1];
NET from Jensen 2008 [7].

| target | quetiapine | norquetiapine | action | reached at therapeutic levels? |
|---|---|---|---|---|
| H1 | 4.4 | 1.1 | antagonist | yes, even at 25–50 mg |
| α1B adrenergic | 14.6 | 46.4 | antagonist | yes |
| 5-HT2A | 38 | 2.9 | antagonist | yes |
| NET | about 1200 (100× weaker) [7] | 12 [7] | inhibitor | metabolite only |
| M1 muscarinic | 1086 | 38.3 | antagonist | metabolite only |
| D1 | 428 | 99.8 | antagonist | partly |
| D2 | 626 | 489 | antagonist, fast dissociation | transiently, around the peak |
| α2 adrenergic | 617 | 1290 | antagonist | partly |
| 5-HT1A | 1040 | 191 | partial agonist (metabolite ~10× more potent) [7] | metabolite, partly |
| benzodiazepine site | none | none | — | no |

**Rough engagement check** (inferred, not measured). At 200 mg twice daily the model's
average quetiapine level is 0.26 µg/mL ≈ 670 nM total, or about 110 nM unbound at 83 %
binding [1]. That is far above the Ki at H1, α1 and 5-HT2A, near the Ki at α2 and D2,
and below it at M1. This ignores brain/plasma partitioning, so it shows rank order only.

**Measured occupancy (PET).** After 3 weeks at 150–600 mg/day, D2 occupancy 12–14 h
after the last dose was only 0–27 %. Two to three hours after a dose it was 58–64 % [8]
(abstract only). D2 blockade is therefore transient and tracks the plasma peak. That is
the usual explanation for quetiapine's low rates of extrapyramidal symptoms and
prolactin rise. The Czech SPC reports PET evidence of 5-HT2 and D2 blockade lasting up
to 12 h, which it uses to justify twice-daily dosing [3].

**Dose dependence.** A single 25 mg dose peaks at 45 ng/mL ≈ 117 nM total, about 20 nM
unbound [4]. That is about 4× the H1 Ki and well below the D2 Ki. At low doses
quetiapine therefore acts mainly as an antihistamine and α1 blocker (inferred from Ki
and the measured peak). The antipsychotic D2 component only appears at hundreds of
mg/day.

## 3. Actions in the body

| system | effect | mechanism | when | source |
|---|---|---|---|---|
| CNS | antipsychotic, antimanic | 5-HT2A + transient D2 antagonism | 300–800 mg/day | [1, 3, 8] |
| CNS | antidepressant (bipolar depression, adjunct in MDD) | norquetiapine NET inhibition, 5-HT1A partial agonism | 150–300 mg/day | [3, 7] |
| CNS | somnolence, sedation | H1 (and α1) antagonism | from the lowest doses; onset in the first 3 days | [3] |
| CNS | cognitive and motor impairment, potentiates alcohol | CNS depression | dose-related | [1] |
| CNS | low EPS, but tardive dyskinesia and NMS possible | weak, transient D2 occupancy | | [1, 8] |
| Cardiovascular | orthostatic hypotension, dizziness, syncope, falls; tachycardia | α1 antagonism | titration phase, elderly | [1, 3] |
| Cardiovascular | QT prolongation (post-marketing, overdose); cardiomyopathy and myocarditis reported | uncertain | | [2, 3] |
| Metabolic | weight gain, hyperglycaemia, dyslipidaemia | not established; H1 antagonism is the usual explanation (inferred) | label warning (5.5) | [1] |
| Endocrine | hyperprolactinaemia possible (label warning), less than with typical antipsychotics; lower thyroid hormone levels | low trough D2 occupancy (prolactin) | | [1, 8] |
| Autonomic | dry mouth, constipation, urinary retention; severe GI hypomotility with other anticholinergics | norquetiapine M1/M3/M5 antagonism | | [1, 3, 7] |
| Blood | leukopenia, neutropenia, agranulocytosis | idiosyncratic | | [1] |
| Other | lens changes (cataracts in chronic dog studies), seizures, dysphagia, impaired temperature regulation, discontinuation symptoms | | | [1] |
| Boxed warning | increased mortality in elderly people with dementia-related psychosis; suicidal thoughts in young people | | | [1] |
| Overdose | drowsiness, tachycardia, hypotension, anticholinergic delirium, coma; QT prolongation | exaggerated pharmacology | survival reported after up to 30 g | [1] |

## 4. Pharmacokinetics

| parameter | value | definition | source | verified |
|---|---|---|---|---|
| Absorption | rapid; tablet 100 % bioavailable *relative to oral solution* | absolute F not stated in any primary document | [1] | yes |
| Tmax, IR | 1.5 h (label); 1.25 h median (single 25 mg); 1.5 h (steady state) | | [1, 4] | yes |
| Tmax, PR | about 6 h (label, SPC); 5 h at steady state (Figueroa) | | [2, 5, 6] | yes / abstract |
| Food | IR: Cmax +25 %, AUC +15 %. PR: high-fat meal Cmax +44–52 %, AUC +20–22 %; light meal none | | [1, 2] | yes |
| Cmax, AUC (healthy) | 25 mg single: Cmax 45 ng/mL, AUC 181 ng·h/mL, CL/F 138 L/h | geometric means, 12 fasted men | [4] | yes |
| Cmax, AUC (patients) | 300 mg b.i.d. steady state: Cmax 1042 ng/mL, AUCτ 4650 ng·h/mL, CL/F 65 L/h | 14 patients, mean age 44 | [4] | yes |
| IR vs PR at steady state | IR 150 b.i.d.: Cmax 568, Cmin 96.5 ng/mL; PR 300 q.d.: Cmax 495, Cmin 95.3; AUC ratio PR/IR 1.04 | patients | [6] | abstract only |
| Linearity | dose-proportional up to 800 mg/day (IR divided; PR once daily) | | [2, 5] | yes |
| t½ | quetiapine ≈ 6 h [1] / ≈ 7 h [2, 3]; norquetiapine ≈ 12 h [2, 3]; steady state within 2 days | terminal | [1, 2, 3] | yes |
| Vd | 10 ± 4 L/kg (apparent) | | [1] | yes |
| Protein binding | 83 % | | [1, 3] | yes |
| Population CL/F | 87.7 L/h, inter-individual CV 32.6 % | Japanese patients, PR, n=322 | [9] | abstract only |
| Norquetiapine/quetiapine | steady-state Cmax 21–27 %, AUC 46–56 % (label); molar Cmax 35 % (SPC) | | [2, 3] | yes |

**Metabolism.** Almost entirely hepatic: < 1 % of a dose is excreted unchanged, and 73 %
of radioactivity is recovered in urine and 20–21 % in faeces [1, 3]. The main routes are
sulfoxidation and oxidation to the parent acid, giving inactive products, and
N-dealkylation to norquetiapine. CYP3A4 drives the sulfoxide and N- and O-desalkyl
pathways; CYP2D6 contributes only the minor 7-hydroxylation [4]. Norquetiapine is both
formed and cleared by CYP3A4 [3]. Quetiapine is not a P-glycoprotein substrate [4], and
it doesn't inhibit CYPs at clinical concentrations [1, 3].

**Special populations.**
- Elderly: clearance 30–50 % lower [1, 3].
- Severe renal impairment: oral clearance 25 % lower, no dose change needed [1].
- Hepatic impairment: oral clearance 25–30 % lower, and 2 of 8 patients had threefold
  exposure [1, 3].
- No effect of sex, race or smoking [1].
- Adolescents: norquetiapine exposure higher than in adults [1, 3].

## 5. Interactions

**PK, measured.**

| perpetrator (dose) | quetiapine exposure | source | in app |
|---|---|---|---|
| ketoconazole 200 mg/day × 4 days | AUC × 6.22, Cmax × 3.35, CL/F −84 % | [4], [1] Table 17, [3] "5–8-fold" | `Ketoconazole 200 mg/day` |
| carbamazepine up to 600 mg/day | AUCτ × 0.13, Cmax × 0.20, CL/F × 7.5 | [4], [3] "down to 13 %" | `Carbamazepine` |
| phenytoin 100 mg t.i.d. | oral clearance × 5 (exposure ≈ × 0.2); SPC: + 450 % | [1], [3] | `Phenytoin 300 mg/day` |
| thioridazine 200 mg b.i.d. | oral clearance + 65 % (SPC: ~70 %) | [1, 3] | not added |
| cimetidine 400 mg t.i.d. | oral clearance − 20 % (US); "unchanged" (Czech SPC) | [1, 3] conflict | not added |
| divalproex | Cmax + 17 %, clearance unchanged | [1] | not added |
| fluoxetine, imipramine, haloperidol, risperidone | no change | [1, 3] | — |

Labels: reduce the dose to one sixth with a strong CYP3A4 inhibitor (US), or
contraindicate the combination (Czech SPC); increase up to fivefold with a strong
inducer. Avoid grapefruit juice [1, 3].

**What the app estimates for unmeasured CYP3A4 perpetrators**, using the preset's
CYP3A4 share of 0.86: clarithromycin ×4.7, diltiazem ×2.9, rifampicin ×0.04,
carbamazepine 400 mg b.i.d. ×0.30. The last is 2.3× weaker than the measured 0.13. See
section 7.

**PD (the app's clash classes).**
- *CNS depressant*: the label warns about other centrally acting drugs and alcohol [1].
- *QT-prolonging*: caution with QT-prolonging drugs [2, 3].
- *Serotonergic*: the Czech SPC warns of serotonin syndrome with SSRIs, SNRIs, MAO
  inhibitors and TCAs [3]. The US labels don't.
- Outside those classes:
  - additive hypotension with antihypertensives;
  - antagonism of levodopa and dopamine agonists;
  - additive anticholinergic GI hypomotility [1].

## 6. Pharmacodynamics and therapeutic drug monitoring

- **Range** *(unverified)*: AGNP 2017 therapeutic reference range 100–500 ng/mL for
  quetiapine, as a trough, with a laboratory alert level of 1000 ng/mL. The consensus
  full text couldn't be read in this session; the publisher blocks automated access.
  One secondary source gives 150–500 ng/mL. An AGNP update was published in July 2026
  (doi 10.1055/a-2860-7861) and couldn't be checked either. Read both before relying on
  the window.
- **Sampling time**: the AGNP range is a trough. With IR's 7 h half-life, levels
  fluctuate about 2-fold or more within a 12 h interval, so sampling time dominates the
  result.
- **Exposure–response**: occupancy follows plasma levels with little delay. D2
  occupancy falls from about 60 % at the peak to under 27 % at the trough [8].
  Antipsychotic efficacy despite low trough occupancy suggests that transient D2
  blockade is enough [8].
- **Variability**: inter-individual CV of clearance is about 33 % in a population-PK
  trial [9]. CYP3A4 inducers and inhibitors move exposure 5–8-fold in either direction
  [1, 4], more than any patient covariate.

## 7. How the app models it

**Structure.**
- One compartment, oral first-order absorption.
- Norquetiapine is simulated as a second curve. The reported quantity and the window are
  the parent's (`reported: parent`), as the label and the AGNP range define them.
- IR (default) and PR formulations.
- The CYP3A4 elimination pathway.

**The decision that could have gone the other way.** It is the clearance the model is
fitted to. Healthy young men after one 25 mg dose clear quetiapine at 138 L/h; patients
at steady state on 300 mg twice daily clear it at 65 L/h, in the same paper [4]. The
preset follows the patients:
- the app simulates therapeutic dosing, which is where the AGNP range applies;
- CL/F 65 L/h with t½ 7 h gives V/F 9.4 L/kg, which independently matches the label's
  10 ± 4 L/kg;
- it reproduces both steady-state references (0.99 and 1.08).

The cost is a 2.1× over-prediction of the healthy-volunteer AUC, recorded as a strict
known miss. The healthy-volunteer half-life of 2.6 h, against 7 h on the label, also
suggests that study missed part of the terminal phase.

| field | value | how obtained | source |
|---|---|---|---|
| `t_half_h` | 7.0 | sourced | [2, 3] |
| `bioavailability` | 1.0 | absolute F unknown, so V and CL are apparent (V/F, CL/F) | [1] |
| `vd_l_per_kg` | 9.4 | derived: CL/F 65 L/h × t½ / ln 2 / 70 kg (label 10 ± 4) | [1, 4] |
| `ka_per_h` | 2.2 | fitted to Tmax 1.25–1.5 h | [1, 4] |
| `pr.ka_per_h` | 0.25 | fitted to the PR Tmax of about 6 h (single dose 6.1 h, steady state 5.5 h) | [2, 3, 6] |
| `elimination.pathways.CYP3A4` | 0.86 | derived: (1 − 1/6.22) / (1 − 0.024), from ketoconazole and the registry's residual activity | [4] |
| `metabolite` | t½ 12 h; formation fraction 0.39, V 9.4 L/kg | t½ sourced; V assumed equal to the parent's, fraction solved to an AUC ratio of 0.51 | [2, 3] |
| `protein_binding.fu` | 0.17 | sourced (83 % bound) | [1, 3] |
| `cv_cl` | 0.33 | population-PK inter-individual CV of CL/F | [9] |
| `mec` / `mtc` | 0.1 / 1.0 µg/mL | AGNP lower limit / laboratory alert level, unverified | §6 |
| default regimen | 200 mg every 12 h (PR: 400 mg every 24 h) | middle of the SPC's usual 300–450 mg/day | [3] |

**Fit** (`preset_report.py quetiapine`):

| reference | observed | predicted | ratio |
|---|---|---|---|
| Cavg,ss 300 mg every 12 h (patients) [4] | 0.3875 | 0.3837 | 0.99 |
| Cmin,ss 150 mg every 12 h (patients) [6] | 0.0965 | 0.1046 | 1.08 |
| Cmax, 25 mg single dose (healthy) [4] | 0.045 | 0.0328 | 0.73 |
| Tmax, 25 mg single dose [4] | 1.25 h | 1.48 h | 1.18 |
| AUC, 25 mg single dose (healthy) [4] | 0.181 | 0.384 | 2.12 (known miss) |
| Cmax,ss 300 mg every 12 h [4] (test only, not a reference metric) | 1.042 | 0.576 | 0.55 (known miss) |
| norquetiapine AUC ratio at steady state [2] | 0.46–0.56 | 0.515 | calibrated |

**What the model can't represent for quetiapine.**
- **Steady-state peaks** come out at about half the measured ones (0.55× in Grimm, 0.51×
  in Figueroa), while averages and troughs match. A two-compartment model fitted to
  these numbers raises peaks to about 0.8× but stretches the terminal half-life to about
  10 h, against 7 h on the label. Without a published mean concentration table, that
  would be a fit to summary numbers rather than to a curve, so it wasn't adopted. The
  peak, where D2 occupancy and orthostatic hypotension are highest, is the part the
  model draws too low.
- **Ketoconazole** acts largely before the drug reaches the circulation (the half-life
  rose only 2.6× while AUC rose 6.2×). The co-medication tool scales the whole curve by
  the measured 6.22 and keeps the 7 h half-life, so it reproduces the exposure but not
  the slower decline. The row's `site` is `clearance` because, with F = 1, a
  bioavailability factor would be capped at 1.
- **Norquetiapine** forms from total parent clearance in the engine. Under a CYP3A4
  inhibitor or inducer its curve doesn't change, whereas in reality both its formation
  and its elimination change.
- Food, hepatic impairment, age, and CYP3A5 genotype aren't modelled.

## 8. Open items

- **Read the AGNP range** (2017, and the 2026 update) from the paper itself, then set
  `mec` (100 or 150 ng/mL?) and record the source.
- **Carbamazepine is under-predicted**: measured ×0.13 against the mechanism's ×0.30.
  The registry's carbamazepine index ratio (midazolam ×0.28, Kanefendt 2024) may
  understate induction for a high-first-pass CYP3A4 substrate. The quetiapine study also
  gave 600 mg/day rather than 800 mg/day, and was done in patients.
- **The CYP3A4 share of 0.86** is derived from the one ketoconazole study. The
  ketoconazole row is named for its real dose (`Ketoconazole 200 mg/day`) so that the
  interaction back-test doesn't count it as an independent check.
- **No published mean concentration–time table** was found for the curve harness.
  Grimm 2006 prints its curves as figures only.
- **Absolute bioavailability**: often quoted as about 9 %, but not in any primary
  document reached.

## Sources

1. SEROQUEL (quetiapine) tablets, US prescribing information, AstraZeneca, DailyMed, March 2025: sections 7, 8, 10, 11, 12.2–12.3, Tables 17–18. https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=0584dda8-bc3c-48fe-1a90-79608f78e8a0
2. SEROQUEL XR (quetiapine) extended-release tablets, US prescribing information, AstraZeneca, DailyMed, March 2025: 5.12, 12.2–12.3. https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=473a3ac4-67f4-4782-baa9-7f9bdd8761f4
3. Quetiapin Neuraxpharm, SPC (SÚKL), sections 4.2, 4.4, 4.5, 5.1, 5.2. https://prehledy.sukl.cz/prehledy/v1/dokumenty/94657 (PR: Quetiapin Teva Retard SPC 5.2, https://prehledy.sukl.cz/prehledy/v1/dokumenty/72041)
4. Grimm SW et al. Effects of cytochrome P450 3A modulators ketoconazole and carbamazepine on quetiapine pharmacokinetics. Br J Clin Pharmacol 2006;61:58–69, Tables 4–5 (full text). https://pmc.ncbi.nlm.nih.gov/articles/PMC1884989
5. Estevez-Carrizo FE et al. Single-dose relative bioavailability of a quetiapine XR formulation. Clin Ther 2011;33:738–45 (abstract only; fed, healthy volunteers: reference AUC0-inf 3546 ng·h/mL, Cmax 366 ng/mL, Tmax 5 h after 300 mg). https://pubmed.ncbi.nlm.nih.gov/21704238/
6. Figueroa C et al. Pharmacokinetic profiles of extended release quetiapine fumarate compared with quetiapine immediate release. Prog Neuropsychopharmacol Biol Psychiatry 2009;33:199–204 (abstract only). https://pubmed.ncbi.nlm.nih.gov/18948162/
7. Jensen NH et al. N-desalkylquetiapine, a potent norepinephrine reuptake inhibitor and partial 5-HT1A agonist. Neuropsychopharmacology 2008;33:2303–12 (abstract only). https://pubmed.ncbi.nlm.nih.gov/18059438/
8. Kapur S et al. A PET study of quetiapine in schizophrenia: transiently high D2 occupancy. Arch Gen Psychiatry 2000;57:553–9 (abstract only). https://pubmed.ncbi.nlm.nih.gov/10839333/
9. Fukushi R et al. Population pharmacokinetics of quetiapine extended-release in Japanese patients with bipolar depression. Clin Ther 2020;42:1067–76 (abstract only). https://pubmed.ncbi.nlm.nih.gov/32518042/
