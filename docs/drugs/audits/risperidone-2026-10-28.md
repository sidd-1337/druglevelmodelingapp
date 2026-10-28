<!--
Blind audit of the risperidone preset, written by an agent running the drug-analysis
skill in Audit mode on 2026-10-28 (develop at 599b045 plus the skill commits). It is
kept as the first measurement of the skill's accuracy. Its findings about stored data
were NOT reviewed or acted on: verify each against the primary document before changing
data/drugs.json. Paths under "scratchpad/" referred to the session's temporary files and
no longer exist. Its feedback on the skill was folded into .claude/skills/drug-analysis;
one point (#7/#9, "derived means unverified") was narrowed to the project's convention
that a plain recomputation from primary numbers can be verified, while a derivation
that adds an assumption cannot.
-->

# Risperidone: blind audit (drug-analysis skill, Audit mode)

Written 2026-10-28 BEFORE opening any stored risperidone entry (drugs.json, reference_pk.json,
interactions.json, pd_effects.json, tests, docs/model-development-log.md, README, core/sukl_client.py).
Schema learned from the `_comment` headers, the diazepam/mirtazapine/quetiapine entries, and core/*.py
(engine code only; no risperidone data). Downloads: `scratchpad/dl_risp/`.

## Step 1. Identify

| item | value | source |
|---|---|---|
| INN / Czech | risperidone / RISPERIDON | SÚKL dlprc API (`atc.nazev.cs`) |
| ATC | N05AX08 | SÚKL dlp/0201072 |
| class | atypical antipsychotic, benzisoxazole derivative ("Jiná antipsychotika") | SPC 5.1; US label 11 |
| MW | risperidone C23H27FN4O2 410.49; 9-OH-risperidone (paliperidone) 426.49 (formula +O; not read in a document) | US label 11 |
| Czech market | Risperdal (Janssen) 1/2/3/4 mg film-coated tablets + 1 mg/ml oral solution (one SPC), Risperdal Consta LAI; generics: Rispen, Risperidon Actavis/Farmax/Orion/Vipharm/Grindeks..., ODT (Mylan, Jacobsen...), oral solutions; no oral PR form | SÚKL dlprc search "risperidon"/"risperdal"/"rispen" (63+8+4 hits) |
| default SÚKL product | **0201072** RISPERDAL 1MG TBL FLM 500, ATC N05AX08, SPC doc id 90474 (SPC231514.pdf, sp. zn. sukls184255/2025) | https://prehledy.sukl.cz/prehledy/v1/dlp/0201072 |

Formulations: only IR (tablets, ODT and solution are bioequivalent, SPC 5.2 / US 12.3). No `pr` entry needed.

## Step 2. Source list

| # | document | URL | gives |
|---|---|---|---|
| S1 | SÚKL SPC Risperdal tbl + sol (Czech, 2025) | https://prehledy.sukl.cz/prehledy/v1/dokumenty/90474 | 4.2 dosing, 4.4 QT, 4.5 interactions, 5.2 F, Tmax, Vd, PB, t½, mass balance |
| S2 | US Risperdal label, DailyMed setid 7e117c7e-02fc-4343-92a1-230061dfc5e0 (v46, Oct 2026) | https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=7e117c7e-02fc-4343-92a1-230061dfc5e0 | 7.1 Table 18 DDI ratios, 12.2 Ki, 12.3 PK incl. EM/PM t½ |
| S3 | Health Canada Risperdal monograph (Janssen) | https://pdf.hres.ca/dpd_pm/00047330.PDF | Table 1.9 single 1 mg active-moiety PK (young/elderly/renal/hepatic); DDI text |
| S4 | Huang 1993 Clin Pharmacol Ther 54:257 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/7690693/ | F 66 %, active-moiety F 100 %, AUC ratio 9-OH/risp 6 oral EM, t½ |
| S5 | Snoeck 1995 Psychopharmacology 122:223 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/8748391/ | origin of S3 Table 1.9; active-moiety t½ 19 h young |
| S6 | Mannens 1993 Drug Metab Dispos 21:1134 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/7507814/ | mass balance by phenotype |
| S7 | Saito 2005 J Clin Psychopharmacol 25:527 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/16282832/ | paroxetine 10/20/40 mg -> active moiety 1.3/1.6/1.8 |
| S8 | Spina 2002 J Clin Psychopharmacol 22:419 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/12172343/ | fluoxetine 20 mg -> active moiety +75 % |
| S9 | Ono 2002 Psychopharmacology 162:50 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/12107617/ | carbamazepine 400 mg/d -> sum 39.6 -> 21.9 ng/mL |
| S10 | TV-46000 Phase 1 (CNS Drugs 2026), full text | https://pmc.ncbi.nlm.nih.gov/articles/PMC13354666/ | Table 4: steady-state ORAL risperidone 3 and 5 mg/day total active moiety Cmax/Cmin/Cavg |
| S11 | EJCP 2024 CYP2D6 haplotypes popPK (Diakonhjemmet), full text | https://pmc.ncbi.nlm.nih.gov/articles/PMC11393095/ | Table 2: CL/F by allele, BSV of 9-OH CL 49 % |
| S12 | Feng 2008 BJCP CATIE popPK (abstract only; PMC/EuropePMC full text behind bot check, not bypassed) | https://pubmed.ncbi.nlm.nih.gov/18771484/ | CL EM/IM/PM 65.4/36/12.9 L/h |
| S13 | Jukic 2019 Lancet Psychiatry (abstract only) | https://pubmed.ncbi.nlm.nih.gov/31000417/ | PM/IM dose-normalised active moiety 1.57x / 1.37x NM |
| S14 | AGNP 2017 (Hiemke 2018, Pharmacopsychiatry 51:9) - **not reached** (Thieme bot check, not attempted) | https://pubmed.ncbi.nlm.nih.gov/28910830/ | - |
| S15 | Moschny 2021 Pharmaceuticals 14:514 Table 1 (second-hand AGNP) | https://pmc.ncbi.nlm.nih.gov/articles/PMC8230242/ | "Risperidone plus 9-hydroxy-risperidone 20-60 ng/mL, alert 120 ng/mL, t½ 2-4 h / 17-23 h, level 2" + ">40 ng/mL only if no response (neurological ADRs)" |
| S16 | Taurines 2022 J Neural Transm (second-hand AGNP) | https://pmc.ncbi.nlm.nih.gov/articles/PMC9188514/ | adult RISam range 20-60 ng/mL |
| S17 | Kapur 1999 Am J Psychiatry 156:286 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/9989565/ | D2 occupancy 63-89 % at 2-12 mg/day; 5-HT2 > D2 |
| S18 | Leysen 1994 J Clin Psychiatry 55 Suppl:5 (abstract only) | https://pubmed.ncbi.nlm.nih.gov/7520908/ | 9-OH binding profile similar to parent |

Note: an AGNP **Update 2026** exists (PMID 42392224, Pharmacopsychiatry 2026); not reachable either.

## Step 3. Mechanism and actions (monograph sections 2-3, condensed)

Core: D2 + 5-HT2A antagonist; clinical effect = risperidone + 9-OH-risperidone (S2 12.1; S1 5.1).

| target | Ki (nM) | engaged at therapeutic unbound levels? | source |
|---|---|---|---|
| 5-HT2A, D2, α1, α2, H1 | 0.12-7.3 (range for the group, label does not split) | yes (unbound active moiety at 4 mg/day ss ≈ 3-9 ng/mL ≈ 7-20 nM; fu 0.10 parent / 0.23 9-OH) | S2 12.2 |
| 5-HT1A/1C/1D | 47-253 | marginal/no | S2 12.2 |
| D1, sigma | 620-800 | no | S2 12.2 |
| muscarinic, β1/β2 | none at 10 µM | no -> no anticholinergic effects | S2 12.2; S3 |
| 9-OH-risperidone | "similar binding profile"; equi-effective | same targets | S18 (abstract), S2 12.3 |

PET: D2 occupancy 63-89 % at 2-12 mg/day, 5-HT2 occupancy higher than D2 at every dose (S17, abstract).
Effects: antipsychotic/EPS/hyperprolactinaemia (D2; S2 5.6; prolactin tracks the active moiety with little hysteresis, S4); orthostatic hypotension, tachycardia, syncope (α1; S2 5.7 "probably reflecting its alpha-adrenergic antagonistic properties"); somnolence, dose-related (41 % at 16 mg/day vs 16 % placebo, S2 5.10; mechanism H1/α1 inferred); metabolic (S2 5.5; mechanism inferred); QT prolongation (S1 4.4; mechanism not stated, inferred hERG); no anticholinergic action.

## Step 4. PK/PD facts with definitions

| parameter | value | definition | source | verified |
|---|---|---|---|---|
| F (absolute) risperidone | 70 % (CV 25 %) | oral, parent | S1 5.2 "Absolutní biologická dostupnost ... 70 % (CV = 25 %)"; S2 12.3 | yes (Huang S4: 66 %) |
| F active moiety | 100 % | oral & IM vs IV, 12 healthy men | S4 abstract | yes (abstract) |
| food | none | | S1 5.2, S2 12.3 | yes |
| Tmax parent | ~1 h (S2) / 1-2 h (S1) | | S1 5.2; S2 12.3 | yes |
| Tmax 9-OH | ~3 h EM, 17 h PM | | S2 12.3 | yes |
| Tmax active moiety | 2 h (median, single 1 mg young); 1.5-2 h (ss, 2-6 mg) | | S3 Table 1.9; S10 | yes |
| t½ risperidone | ~3 h (psychotic pts, S1); 3 h CV 30 % EM, 20 h CV 40 % PM (S2) | | S1 5.2; S2 12.3 | yes |
| t½ 9-OH | 21 h CV 20 % EM, 30 h PM (S2); 24 h incl. active fraction (S1) | | S2 12.3; S1 5.2 | yes (two values) |
| t½ active moiety | 20 h (S2); 17 h median young (S3); 20 ± 2.5 h (S4); 23-27 h ss patients (S10) | | | yes |
| Vd risperidone | 1-2 L/kg | | S1 5.2; S2 12.3 | yes |
| protein binding | risperidone 90 % (S1, S2), 88 % (S3 text); 16 % unbound young (S3 Table 1.9); 9-OH 77 % | | | yes |
| linearity | dose-proportional 1-16 mg/day (S2); 2-5 mg/day linear, 6 mg less than proportional (S10) | | | yes |
| metabolism | CYP2D6 9-hydroxylation main; N-dealkylation minor; CYP3A4 lesser extent; both P-gp substrates | | S1 4.5/5.2; S2 12.3 | yes (qualitative) |
| CYP2D6 fm (parent) | 0.85 (t½ 3 vs 20 h EM/PM); 0.92 (*1/*1 CL/F 50.6 vs 4.2 L/h); 0.80 (CATIE 65.4 vs 12.9 L/h) | in vivo phenotype/genotype | S2; S11 Table 2; S12 | derived |
| mass balance | 70 % urine + 14 % faeces at 1 week; urine risp + 9-OH 35-45 % of dose; EM: unchanged 4 %, 9-OH 32 %, N-dealkyl 10-13 % | 1 mg 14C solution, 3 men | S1, S2, S6 | yes |
| active-moiety AUC 1 mg | 132 ng·h/mL (median, n=8 young healthy, age 25-35, RIA); Cmax 9.1 ng/mL; Tmax 2 h; t½ 17 h; CL/F 127 mL/min | single dose, fasting (S5) | S3 Table 1.9 | yes |
| steady-state oral TAM | 3 mg/d: Cmax 39.8, Cmin(C24) 11.86, Cavg 20.7 ng/mL; 5 mg/d: 62.2 / 18.8 / 31.9 (geometric means, CV 27-56 %) | 7 days once daily, schizophrenia pts, 81 % men, 96 % Black, mean 86.9 kg; TAM = risp + 9-OH × 410/426 (molar-corrected) | S10 Table 4 | yes |
| elderly / renal | active-moiety CL -30 % elderly, -60 % moderate-severe renal | | S2 8.6; S1 5.2 | yes |
| genotype on active moiety | "similar" EM vs PM (S1, S2); PM 1.57x, IM 1.37x NM (S13, n=725) | | | yes |
| therapeutic range | 20-60 ng/mL risperidone + 9-OH (trough, total), alert 120 ng/mL, level 2; ">40 ng/mL only when response insufficient" | AGNP 2017 | S15, S16 (second-hand) | **no** |

Measured DDIs (active moiety):

| perpetrator (dose) | risperidone dose | AUC ratio | Cmax ratio | source |
|---|---|---|---|---|
| fluoxetine 20 mg/d | 2-3 mg b.i.d. | 1.4 | 1.5 | S2 Table 18 |
| fluoxetine 20 mg/d | 4-6 mg/d | +75 % trough (1.75), n=9 | - | S8 abstract |
| fluoxetine (n=10, 3 wk) | - | "about 50 %" peak and AUC | | S3 |
| paroxetine 10 / 20 / 40 mg/d | 4 mg/d | 1.3 / 1.6 / 1.8 | - | S2 Table 18 = S7 (troughs, n=12) |
| paroxetine (n=10, 4 wk) | - | +45 % (1.45) | | S3 |
| carbamazepine 573 ± 168 mg/d | 3 mg b.i.d. | 0.51 | 0.55 | S2 Table 18 |
| carbamazepine 400 mg/d, 1 wk | 6 mg/d | 0.55 (sum 39.6 -> 21.9 ng/mL), n=11 | - | S9 abstract |
| itraconazole 200 mg/d | 2-8 mg/d | ~1.7 (concentrations) | - | S1 4.5 |
| ranitidine / cimetidine / erythromycin / amitriptyline | | 1.2 / 1.1 / 1.1 / 1.2 | 1.4 / 1.3 / 0.94 / 1.1 | S2 Table 18 |

## Step 5. Model mapping (my decisions)

- Structure: one compartment for parent + one for 9-OH (no published curve in hand shows a distribution phase that matters at q12-24 h).
- Reported quantity: **sum** (AGNP range is risperidone + 9-OH; label: clinical effect from the combined concentrations). Metabolite gets its own compartment: active, t½ 21 h vs 3 h, so not proportional to the parent within a dosing interval, and the range refers to the sum. potency_ratio 1.0 (equipotent; AGNP sum is unweighted).
- Formulations: IR only (tablet = solution = ODT).
- Decision that could have gone the other way: encoding CYP2D6. The engine keeps formation_fraction constant, so a PM (activity 0) gives a simulated sum **2.5x** EM at 4 mg/day, while the label says the sum is "similar" and Jukic 2019 finds 1.57x. I would still encode it (parent curve and inhibitor direction are right) but flag it.
- Fit (my own script with core.pk_engine.PKModel, scratchpad/scripts/fit_risp.py, grid2_risp.py): F and both t½ fixed from sources; Vparent, ka, Vmet fitted to S3 (single 1 mg) and S10 (ss 3 and 5 mg).

### Values I would store

| field | mine | how | source | where | verified |
|---|---|---|---|---|---|
| t_half_h | 3.0 | sourced | S1 | 5.2 "poločasem přibližně 3 hodiny" (psychotic pts); S2 12.3 3 h CV 30 % EM | true |
| t_half_range_h | [2.0, 4.0] | sourced second-hand | S15 | AGNP t½ 2-4 h (≈ 3 h ± CV 30 %, S2) | false (second-hand) |
| vd_l_per_kg | 1.0 | fitted (at the label's lower bound 1-2 L/kg) | S1 5.2 | | false (fitted) |
| bioavailability | 0.70 | sourced | S1 | 5.2 "Absolutní biologická dostupnost ... 70 % (CV = 25 %)" | true |
| ka_per_h | 1.5 | fitted (parent Tmax 1.5 h; sum Tmax 1.8 h single, 1.7 h ss) | S1/S3/S10 | | false (fitted) |
| elimination.pathways | {CYP2D6: 0.85} | derived 1 - 3/20 from EM/PM t½ | S2 12.3 | concordant S11 Table 2 (0.92), S12 (0.80); CYP3A4 and renal (4 % unchanged, n=1 EM, S6) left out | false (derived) |
| metabolite.formation_fraction | 0.85 | derived (= CYP2D6 share; S6 EM: 9-OH 32 % of dose in urine alone) | S2/S6 | | false |
| metabolite.t_half_h | 21 | sourced | S2 | 12.3 "about 21 hours (CV=20%) in extensive metabolizers"; S1 says 24 h | true (note S1 24 h) |
| metabolite.vd_l_per_kg | 2.3 | fitted to sum AUC/Cavg (compensates no first-pass 9-OH formation) | S3, S10 | | false (fitted) |
| metabolite.molar_ratio | 0.9625 (410.49/426.49) | derived | S2 11 | | n/a |
| metabolite.potency_ratio | 1.0 | equipotent / unweighted AGNP sum | S2 12.3 | | false (placeholder-ish) |
| metabolite.reported | sum | | S2 12.1, S15 | | - |
| protein_binding.fu | 0.10 | sourced (9-OH fu 0.23) | S1 | 5.2 "90 %, u 9-hydroxy-risperidonu 77 %"; S3 says 88 %, Table 1.9 16 % unbound | true |
| cv_cl | 0.49 | BSV of 9-OH clearance (dominates the sum), n=512 TDM popPK | S11 | Table 2 "BSV CLmet/fmet (% CV) 49%"; cross-check S10 Cavg GM CV 30-42 % | (n/a; derived from a published spread) |
| cv_f | none | active-moiety F ≈ 100 % (S4); the 25 % CV is the parent's F only | | | |
| mec | 0.020 µg/mL | AGNP lower limit (sum, trough) | S15/S16 | Table 1 | false (second-hand) |
| mtc | 0.120 µg/mL | AGNP laboratory alert level; notes must give range 20-60 and the >40 EPS caution | S15 | Table 1 | false (second-hand) |
| default regimen | 4 mg q24h | SPC: start 2 mg, day 2 4 mg, optimum 4-6 mg/day, once or twice daily | S1 | 4.2 | - |
| SÚKL | name map RISPERIDON / RISPERIDONE -> risperidone; code 0201072 (RISPERDAL 1MG TBL FLM 500); ATC N05AX08 | | SÚKL API | | |

### Exposure references (default formulation, oral)

| metric | dose | value (µg/mL, µg·h/mL, h) | definition | source / where | verified | my pred/obs |
|---|---|---|---|---|---|---|
| auc_inf | 1 mg | 0.132 | active moiety, median, n=8 young healthy, fasting, RIA | S3 Table 1.9 (orig. S5) | true | 1.21 |
| cmax | 1 mg | 0.0091 | same | S3 Table 1.9 | true | 0.87 |
| tmax_h | 1 mg | 2.0 | median | S3 Table 1.9 | true | 0.89 |
| cavg_ss | 3 mg q24h | 0.0207 | TAM GM, day 7, patients | S10 Table 4 | true | 0.97 |
| cmin_ss | 3 mg q24h | 0.01186 | C24h GM | S10 Table 4 | true | 0.97 |
| cavg_ss | 5 mg q24h | 0.0319 | TAM GM | S10 Table 4 | true | 1.04 |

(Model Cmax_ss/obs 0.86 at 3 mg, 0.92 at 5 mg. The AUC_9OH/AUC_risp ratio comes out 2.7 instead of Huang's 6, because the engine cannot form 9-OH presystemically; the parent carries the early peak instead.)

### Interactions I would store (`data/interactions.json`)

| perpetrator (registry name) | factor | site | onset_h | source / where | verified | note |
|---|---|---|---|---|---|---|
| Paroxetine 20 mg/day | 1.6 | clearance | ~168 (estimate) | S2 7.1 Table 18 "Paroxetine 20 mg/day, 4 mg/day, AUC 1.6" (= S7 troughs, not significant) | true | S3 cites another study at +45 % (1.45) |
| Fluoxetine 20 mg/day | 1.4 | clearance | ~700 (estimate; norfluoxetine) | S2 7.1 Table 18 "Fluoxetine 20 mg/day ... AUC 1.4, Cmax 1.5" | true | conflicts with S8 +75 % (troughs) and S3 "about 50 %" |
| Carbamazepine (alias of "Carbamazepine 400 mg b.i.d.") | 0.51 | clearance | 336 (S1 4.5: at least 2 weeks to full induction) | S2 7.1 Table 18 "573 ± 168 mg/day, 3 mg twice daily, AUC 0.51, Cmax 0.55" | true | measured dose 573 mg/day differs from the registry's 800 mg/day; S9 gives 0.55 at 400 mg/day |

Engine caveat (my simulation): a `clearance` factor only rescales the parent's t½, so the simulated **sum** moves 1.16 for paroxetine 1.6, 1.11 for fluoxetine 1.4 and 0.87 for carbamazepine 0.51. Only `bioavailability` reproduces a measured sum ratio (0.51 -> 0.51; 1.6 is capped by F≤1 at 1.43).

### PD clash classes (`data/pd_effects.json`)

| class | tag? | evidence | source | verified |
|---|---|---|---|---|
| serotonergic | no | 5-HT2A antagonist; no serotonin-syndrome warning or interaction in S1/S2 | - | - |
| cns_depressant | yes | S1 4.5 "opatrně v kombinaci s dalšími centrálně účinkujícími léčivými přípravky včetně alkoholu, opiátů, antihistaminik a benzodiazepinů vzhledem ke zvýšenému riziku sedace"; S2 7.2 "Centrally Acting Drugs and Alcohol" | S1 | true |
| qt_prolonging | yes | S1 4.4 "Prodloužení QT" warning + 4.5 caution with QT-prolonging drugs; S2 has no W&P QT section (overdose only, 10.1) | S1 | true |
| sympathomimetic | no | α1 antagonist (hypotension); nothing sympathomimetic | - | - |

---

# Step 4. Diff against the stored preset (written after the blind section above)

Stored entries read: data/drugs.json, reference_pk.json, interactions.json, pd_effects.json (risperidone),
core/sukl_client.py (lines 38-39, 75, 99), README.md:274 / README.cs.md:301, tests (test_drug_database,
test_population, test_regimen), docs/model-development-log.md. preset_report.py run after the blind write-up:

```
t1/2 3 h  V 108.5 L  CL 25.07 L/h  CL/F 35.81 L/h  F 0.7  ka 2/h
metabolite 9-OH: t1/2 21 h  V 113.4 L  fm 0.8  reported=sum
ok   auc_inf  1 mg  obs 0.132  pred 0.1834  ratio 1.39 (limit 1.5)  verified
```

Both parameter sets scored against every primary reference I found (scratchpad/scripts/compare_risp.py):

| reference (source) | obs | stored pred/obs | mine pred/obs |
|---|---|---|---|
| auc_inf 1 mg (S3 T1.9) | 0.132 | 1.39 | 1.21 |
| cmax 1 mg (S3 T1.9) | 0.0091 | **0.65** (miss) | 0.87 |
| tmax 1 mg (S3 T1.9) | 2.0 h | 0.97 | 0.89 |
| cavg_ss 3 mg q24h (S10 T4) | 0.0207 | 1.11 | 0.96 |
| cmin_ss 3 mg q24h (S10 T4) | 0.01186 | 1.27 | 0.97 |
| cmax_ss 3 mg q24h (S10 T4) | 0.0398 | 0.80 | 0.86 |
| cavg_ss 5 mg q24h (S10 T4) | 0.0319 | 1.20 | 1.04 |
| cmin_ss 5 mg q24h (S10 T4) | 0.0188 | 1.34 | 1.02 |
| 9-OH : parent AUC, oral EM (S4: 6) | 6 | 5.6 | **2.7** |
| CYP2D6 PM / NM sum, registry PM activity 0.12 (S13: 1.57; label "similar") | 1.57 | 1.36 | 1.80 |

## Diff table

| # | field | mine | stored | my source | stored source | verdict |
|---|---|---|---|---|---|---|
| 1 | t_half_h | 3.0 | 3.0 | SPC 5.2; US 12.3 | US label (metabolite note) | agree |
| 2 | t_half_range_h | [2, 4] | [2, 4] | AGNP second-hand (S15) | none given | agree |
| 3 | vd_l_per_kg | 1.0 (fitted) | 1.55 (solved) | fit to S3+S10 | joint solve, AUC ratio ~5 | not comparable (both fitted) |
| 4 | bioavailability | 0.70 | 0.70 | SPC 5.2 | (no block) | agree |
| 5 | ka_per_h | 1.5 (fitted) | 2.0 | fit | - | not comparable (both fitted) |
| 6 | pathways CYP2D6 | 0.85 | 0.80 | US 12.3 t½ EM/PM 3 vs 20 h | Feng 2008 CL 65.4 vs 12.9 L/h | sources conflict (S11 genotype popPK gives 0.92) |
| 7 | elimination.verified | false (derived) | **true** | SKILL ground rule | - | stored is wrong |
| 8 | metabolite.formation_fraction | 0.85 | 0.80 | = CYP2D6 share | = CYP2D6 share | sources conflict (same method, follows #6) |
| 9 | metabolite.t_half_h | 21 | 21 | US 12.3 | US 12.3 | agree (SPC says 24; both noted it) |
| 10 | metabolite.vd_l_per_kg | 2.3 (fitted) | 1.62 (solved) | fit | solve | not comparable |
| 11 | molar_ratio | 0.9625 | 0.9625 | MW | MW | agree |
| 12 | potency_ratio | 1.0 | 1.0 | equipotent | equipotent | agree |
| 13 | reported | sum | sum | AGNP / US 12.1 | same | agree |
| 14 | protein_binding.fu | 0.10 | 0.10 | SPC 5.2 | US label | agree (both note HC 88 % and 16 % unbound) |
| 15 | cv_cl | 0.49 | **1.08** | S11 Table 2: 9-OH CL BSV 49 %, n=512 | "published ~20-fold spread 1.8-36.8 nmol/L per mg/day", no source named anywhere | stored is wrong (see below) |
| 16 | mec | 0.020 | 0.020 | AGNP lower limit (second-hand) | AGNP (notes) | agree |
| 17 | mtc | 0.120 (AGNP alert) | **0.080** ("riziko EPS") | S15 Table 1 (second-hand) | none | stored is wrong (unsourced) |
| 18 | default regimen | 4 mg q24h | 2 mg q12h | SPC 4.2 | - | not comparable (same 4 mg/day; SPC allows both) |
| 19 | default SÚKL code | 0201072 RISPERDAL 1 mg | 0102979 MEDORISPER 2 mg | SÚKL API | SÚKL API (live: registered, N05AX08, SPC 94819) | not comparable (both valid) |
| 20 | name map / ATC | RISPERIDON, RISPERIDONE; N05AX08 | same | SÚKL | | agree |
| 21 | ref auc_inf 1 mg | 0.132 | 0.132 | HC T1.9 | HC T1.9 | agree |
| 22 | ref cmax 1 mg 0.0091, tmax 2 h | stored as rows | **absent** (Cmax quoted in the auc note only) | HC T1.9 | - | stored is wrong (omits rows its own source gives; the Cmax row would be a 0.65 miss -> strict xfail) |
| 23 | ref ss 3/5 mg (S10) | 3 rows | absent | S10 Table 4 | - | not comparable (source newer than the preset; stored model passes them, 0.80-1.34) |
| 24 | Paroxetine 20 mg/day | 1.6, clearance, 168 h | 1.6, clearance, 168 h | US T18 = Saito 2005 | US T18 | agree |
| 25 | Fluoxetine 20 mg/day | 1.4, clearance, ~700 h | 1.4, clearance, 168 h | US T18 | US T18 | agree (onset_h is unverified documentation; stored note itself says "accumulates over weeks") |
| 26 | Carbamazepine (alias) | 0.51, clearance, 336 h | 0.51, clearance, 336 h | US T18; SPC 4.5 ≥2 weeks | US T18 | agree |
| 27 | interaction notes | - | "after 3 weeks", "(Bondolfi 2002)", "(EM 663.0/470.0)" attributed to the label | - | US label | stored is wrong (minor: none of these is in the current DailyMed v46 Table 18 or section 7) |
| 28 | my "engine caveat" (clearance factor moves the sum only 1.16x) | claimed | regimen path scales the reported curve; test_regimen pins 1.6 | core/interactions.apply_interactions | core/regimen.py + perpetrator_timecourse.py | **mine is wrong** (apply_interactions is used only by tests) |
| 29 | my PM claim (sum 2.5x) | claimed with activity 0 | registry PM activity is 0.12 -> 1.80 (stored 1.36) | - | data/enzymes.json | **mine is wrong** (wrong activity; corrected value 1.80) |
| 30 | PD serotonergic | none | none | - | - | agree |
| 31 | PD cns_depressant | yes | yes | SPC 4.5 + US 7.2 | US 7.2 | agree |
| 32 | PD qt_prolonging | yes | yes | Risperdal SPC 4.4/4.5 | Medorisper SPC 4.4/4.5 | agree |
| 33 | PD sympathomimetic | none | none | - | - | agree |
| 34 | README rows (EN/CS) | - | t½ **20.0**, Vd **1.70**, ka 2.00 | - | - | stored is wrong (stale vs drugs.json 3.0 / 1.55) |
| 35 | tests/test_drug_database.py strict xfail | - | pins Cavg 0.028 at 2 mg q12h | - | the reference_pk note says 0.028 came from a long-acting-injectable review and was replaced | stored is wrong (stale test keeps a discredited target) |

**Counts (35 rows):** agree 18 · stored is wrong 7 (#7, 15, 17, 22, 27, 34, 35) · mine is wrong 2 (#28, 29) · sources conflict 2 (#6, 8) · not comparable 6 (#3, 5, 10, 18, 19, 23).
Restricted to the values the task listed (t½, range, Vd, F, ka, pathways, 6 metabolite fields, fu, cv_cl, MEC, MTC, regimen, exposure ref, 3 interaction factors, 4 PD classes = 25 values): agree 17, stored wrong 2 (cv_cl, mtc), sources conflict 2 (CYP2D6 share, formation fraction), not comparable 4 (Vd, ka, metabolite Vd, regimen).

## Disagreements resolved against the primaries

- **#6/#8 CYP2D6 share (0.85 vs 0.80).** All three are in-vivo and all are derived: US label 12.3 t½ 3 h (EM) vs 20 h (PM) -> 0.85 if V is unchanged; Feng 2008 abstract CL 65.4 vs 12.9 L/h -> 0.80, but its PM class is a mixture-model class (the stored note says 41 % of subjects, far above the ~7 % PM frequency the label gives), which dilutes it; EJCP 2024 Table 2 (genotype-defined, *3-*6 fixed at 0) -> 4.2 vs 50.6 L/h -> 0.92. No primary states a fraction. Either value is defensible; both must be `verified: false`. The stored PM/NM sum (1.36) sits nearer the label's "similar", mine (1.80) nearer Jukic's 1.57 from above.
- **#7 verified flag.** The stored block computes 1 - 12.9/65.4 and assumes no residual PM activity. SKILL.md: "A value that is derived ... is verified: false". Its source URL (PMC2661978) shows a reCAPTCHA page to automated access; the CL numbers are in the PubMed abstract (PMID 18771484).
- **#15 cv_cl.** Stored: CV 1.08 from a "~20-fold spread (1.8-36.8 nmol/L per mg/day)" read as a 90 % interval (log-SD 0.91). No document is named in drugs.json, the notes, the dev log or commits b980179/99ef6b0, so the range cannot be checked, and a 20-fold min-max range from a TDM set of hundreds is not a 5-95 % interval. Primaries I could open: EJCP 2024 (PMC11393095) Table 2, between-subject variability of 9-OH-risperidone clearance 49 % (53 % before covariates), parent CL/F 86 % after genotype, n=512 TDM patients; TV-46000 Phase 1 (PMC13354666) Table 4, steady-state oral active-moiety Cavg geometric CV 30.2 % (3 mg) and 41.6 % (5 mg); Snoeck 1995: the active moiety varies "much less" than risperidone. These imply a CV of ~0.4-0.5 for the sum, so the stored 1.08 (90 % band 19.9-fold at 2 mg q12h, pinned by test_population 12 < span < 30) is about twice too wide in log-SD. I believe the stored value is wrong; it is at least unsourced.
- **#17 MTC.** Stored 0.08 µg/mL "riziko EPS" has no source. The second-hand AGNP 2017 table (S15) names 120 ng/mL as the laboratory alert level and 40 ng/mL as the point above which neurological adverse effects are a concern; no reachable document names 80. Under SKILL.md ("MTC is the level a primary source names") the stored value has no basis, and the notes omit the alert level. Use 0.12 (alert) or 0.04 (EPS caution) and say which.
- **#22 missing Cmax/Tmax rows.** HC Table 1.9 gives Cmax 9.1 ng/mL and Tmax 2 h next to the AUC the preset uses; app-data.md says to add them. The stored model predicts 0.0059 (0.65), a miss that would need a KNOWN_MISSES entry, which is what the ground rule asks for. My parameters hit it (0.87) only by letting the parent carry the early peak (9-OH : parent AUC 2.7 against the primary 6): the structural limit (no first-pass 9-OH formation) forces this trade-off either way.
- **#27 note details.** The current label (setid 7e117c7e..., v46) Table 18 gives fluoxetine 1.4/1.5 and carbamazepine 0.51/0.55 with doses only; it does not say "after 3 weeks", name Bondolfi 2002, or give "EM 663.0/470.0". The values are right; the notes say more than the cited document does.

## Skill feedback (concrete edits)

1. **SKILL.md, Audit mode:** widen the blind list to "every file that names the drug": tests in any file (test_regimen, test_population pin risperidone numbers), docs/model-development-log.md, README tables, core/sukl_client.py, and code comments. Also, app-data.md tells the reader to copy `risperidone` as the metabolite template, so an audit of risperidone (or diazepam) must use the other one; say so.
2. **SKILL.md step 2:** give the SÚKL route: `GET https://prehledy.sukl.cz/prehledy/v1/dlp/<kodSUKL>` returns `dokumenty[]` with the SPC id; then `GET .../dokumenty/<id>` (PDF). I had to probe for it. Also list the open full-text route: `eutils efetch db=pmc&id=<PMCID>&rettype=xml` works for open-access articles where the PMC/Europe PMC web pages show a bot check; a stub saying "publisher does not allow downloading" means record "abstract only".
3. **SKILL.md, therapeutic ranges:** AGNP 2017 is Thieme-only and bot-checked. Say what to do: cite an open-access secondary that reproduces the AGNP table, mark `verified: false`, name it. Note that an AGNP Update 2026 exists (PMID 42392224) and say which edition governs.
4. **SKILL.md / app-data.md, MEC/MTC provenance:** drugs.json has no `source`/`verified` for `mec`/`mtc`, so an unsourced MTC (risperidone 0.08) passes every test. Either add a `window` block with source/verified/note, or require the notes to cite the document and the definition (alert level / upper limit / named adverse threshold). Pick one default: "MTC = AGNP alert level unless a primary names a different adverse threshold".
5. **SKILL.md step 5, cv_cl:** "a percentile range of dose-normalised levels" needs a warning: a min-max range is not a percentile range; state the percentile assumption in the note, and prefer a popPK BSV. For a sum-reported drug say the CV must describe the sum (dominated by the metabolite's clearance), and that the band adds unknown-phenotype spread on top, so a genotype-inclusive spread double counts.
6. **SKILL.md step 5, metabolite formed presystemically:** the engine forms the metabolite only from systemic parent. Say which target wins when sum Cmax and the metabolite:parent ratio cannot both be met (risperidone: stored kept the ratio and missed Cmax; I kept Cmax and missed the ratio), and require recording the loser as a known miss or a note.
7. **app-data.md, reference_pk:** add to the definition checklist "sum: plain mass sum (AGNP, engine) or molar-corrected total active moiety (TAM = risp + 9-OH × 410/426)". It is a 4 % difference but a definition all the same. Add to step 7: "every Cmax/Tmax the reference document prints is a row".
8. **app-data.md, interactions:** say that the live co-medication path (core/regimen.py) scales the reported concentration curve, so a measured active-moiety ratio is reproduced as stored, and `site` only decides IV skipping and the F ≤ 1 cap. core/interactions.py's module docstring and `apply_interactions` (used only in tests) describe the old parent-clearance behaviour and misled me (#28).
9. **SKILL.md ground rules, verified:** add the explicit case "a pathway fraction computed from phenotype or genotype clearances is derived: `verified: false` even when every input was read" (the stored risperidone block breaks this, and nothing catches it).
10. **preset_report.py:** for a metabolite preset, also print the metabolite:parent AUC ratio and the PM/NM ratio of the reported quantity, and accept ad-hoc candidate references on the command line, so an audit can score new primaries without editing reference_pk.json.
11. **app-data.md:** says `notes` are Czech; diazepam's are English. Pick one.
12. **SKILL.md step 1:** say whether the default product should be the originator (whose SPC usually has the fullest 5.2) or any registered product; risperidone's default is a generic, so the preset cites two different SPCs.

## Addendum: the skill changed during the run

Another session committed the skill (781bbd9) and quetiapine (59cc866) while this audit ran. The committed
SKILL.md has 209 lines, against the 174 I followed. The stored risperidone values are unchanged at HEAD.
Effect on the diff:

- **#19 default SÚKL code -> mine is wrong** under the committed step 1, which adds "currently supplied
  (`jeDodavka`)". Live API: 0201072 RISPERDAL 1 mg `jeDodavka=false`; 0046969 RISPERDAL 3 mg and the stored
  0102979 MEDORISPER 2 mg are `true`. The stored code is valid; a Risperdal choice would be 0046969.
- The new rule "fit the population the app simulates (patients at steady state)" favours S10 (patients,
  steady state) over S3 (healthy, single dose) as the fitting target. The stored model meets S10 at 0.80-1.34
  without having used it.
- The new rule "a misfit that only a more complex structure would fix: pin it as a strict xfail, don't fit
  around it" supports the stored choice (keep the 9-OH:parent ratio and accept a Cmax miss) over mine (move the
  parent to hit Cmax). The stored preset still omits that Cmax row and its xfail (#22 stands).
- Feedback items 2 (dlprc route) and 3 (AGNP 2026 update and the Thieme rule) are now partly in SKILL.md.
  What remains: the `dlp/<kod>` -> `dokumenty/<id>` route (the skill points to
  `core.sukl_client.fetch_sukl_detail`, which an Audit must not open for its own drug), and the eutils
  PMC-XML route.

**Revised counts (35 rows):** agree 18 · stored is wrong 7 · mine is wrong 3 (#19, 28, 29) · sources conflict 2 · not comparable 5 (#3, 5, 10, 18, 23).
