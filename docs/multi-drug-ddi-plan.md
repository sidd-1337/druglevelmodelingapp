# Dynamic multi-drug interaction model — design and validation plan

**Status:** plan, not yet built. **Base:** `206adf8` (tip of `develop`).

The question behind this plan: *can we simulate the levels of several drugs taken
together, rigorously, accounting for the interactions between them and the metabolic
pathways they share?* Commit `813dd19` took the first step — one shared enzyme pool,
integrated as a single ODE across co-dosed drugs. This plan generalises that step into a
mechanistic drug-drug-interaction (DDI) engine and says exactly how far the evidence lets
it go, and where it does not.

The guiding rule of this repository still holds: **a number is only as good as the primary
document behind it, and the model is only as good as what it is measured against.** Every
parameter table below carries a provenance status, and every mechanism carries a
validation target. Where the evidence does not support a mechanism, this plan says so and
keeps the existing fixed factor rather than inventing one.

---

## 1. What exists today, and what is missing

The engine already does more than its age suggests. `core/pk_engine.py` implements
**mechanism-based inactivation (MBI) dynamically** for a single drug: the active fraction
`E` of the enzyme carrying fraction `fm` of a drug's clearance follows

```
dE/dt = kdeg·(1 − E) − kinact·C/(KI + C)·E        (pk_engine.py:567)
CL(t) = CL·(1 − fm + fm·E)                          (pk_engine.py:553)
```

and `coupled_profiles` / `_coupled_states` (`pk_engine.py:786`, `:831`) already share
**one** such enzyme pool across every co-dosed drug: each drug carrying a share `s` of its
clearance on that enzyme is slowed moment-to-moment by `factor = 1 + Σ s·(E − 1)`
(`pk_engine.py:851`). The static limit of all this — `1 / (fm·residual + (1 − fm))` — lives
in `core/enzymes.py:159` (`derived_exposure_factor`) and is used today both to derive
unmeasured interaction factors and as a sanity check.

So the dynamic scaffolding is in place. What is missing, to call the model a genuine
mechanistic DDI engine:

1. **Only inactivation is dynamic.** Reversible inhibition and induction are not modelled
   at all as time-varying effects — they exist only as fixed exposure factors in
   `data/interactions.json`, or as back-derived residual activities in `data/enzymes.json`.
2. **Enzymes are not first-class.** Each drug inactivates at most one enzyme
   (`CoDrug.inactivates`, a single string). There is no registry of enzyme *pools* that
   several perpetrators and victims attach to with their own mechanisms.
3. **Perpetrators have no PK.** Paroxetine, fluvoxamine, carbamazepine, rifampicin and the
   rest exist only as static factors. A time-varying effect needs the perpetrator's own
   concentration curve `C(t)`, which they do not have.
4. **Bioavailability is a single constant.** There is no separation of gut-wall from
   hepatic first pass, so an oral CYP3A4 interaction (the common case) cannot be modelled
   mechanistically, and the risperidone first-pass gap noted in `data/drugs.json` cannot
   be closed.

---

## 2. The one decision that shapes everything: middle-out, not pure bottom-up

The brief asked for a **bottom-up** model (in-vitro constants → predicted clinical effect)
validated by a 2-fold gate. The research says a *pure* bottom-up model will fail that gate
for several of the drugs we care about, and it is worth saying plainly why, because it
changes the design.

- **Fluvoxamine** is roughly **10× more potent in vivo** on CYP1A2 and CYP2C19 than its in
  vitro Ki predicts (Yao 2003; the EMA scaled fluvoxamine's in vitro Ki down 10× to match
  observed interactions).
- **Ciprofloxacin**'s in vitro CYP1A2 IC50 (~135 µM) **cannot explain** the observed ~10×
  tizanidine AUC rise (Granfors 2004).
- **Bupropion** modelled as parent-only reversible inhibition **underpredicts** its ~5×
  desipramine effect; the validated model needs three inhibitory metabolites *and* a
  CYP2D6-downregulation term (Sager 2017).
- Microsomal and plasma binding (`fu,inc`, `fu,p`) move an apparent Ki by **one to two
  orders of magnitude** (bupropion ~100×, paroxetine ~15×). A Ki quoted without its binding
  correction is close to meaningless.

This is not a defect peculiar to us — it is why every production PBPK network (Open Systems
Pharmacology, Simcyp) is in practice **middle-out**: in-vitro mechanism and structure, but
a few sensitive constants (a Ki, an induction Emax) *anchored* to one clinical observation,
then validated against the *rest* of the clinical database. And it is, pleasingly, exactly
what `data/enzymes.json` already does — it back-derives residual activity from a clinical
index-substrate AUC ratio.

**Decision.** Build middle-out:

- Take mechanism and structure (MBI vs reversible vs induction; which metabolites matter)
  from in-vitro and PBPK primaries.
- Carry every rate constant with its binding basis, as **unbound** `Ki,u` / `KI,u` / `EC50,u`.
- Where a drug's in-vitro constant is known to mispredict (fluvoxamine, ciprofloxacin,
  bupropion), use the **in-vivo-anchored** value the OSP/Simcyp networks use, and mark it
  `anchored` in the data, distinct from `verified` (measured in vitro) and `estimated`.
- **Validate on held-out pairs.** For each perpetrator, anchor on *one* clinical pair and
  require every *other* studied pair to land within 2-fold (geometric). A perpetrator with
  only one studied pair is anchored, not validated, and is flagged as such.

This keeps the repository's honesty contract: an anchored number is not dressed up as a
measured one.

---

## 3. Scope of mechanisms — and what we are deliberately *not* building

You asked to cover reversible inhibition, induction with onset, gut-wall first pass, and
transporters/renal. The first three are in. The fourth is where the evidence stops us, and
this plan pushes back on it.

### In scope (mechanistic, time-varying)

| Mechanism | How | Drives |
|---|---|---|
| **Reversible inhibition** | activity `A = 1/(1 + Σ Iu/Ki,u)`, instantaneous | paroxetine(no), fluoxetine/norfluoxetine, terbinafine, fluvoxamine, fluconazole, ciprofloxacin, ketoconazole |
| **Mechanism-based inactivation** | `dE/dt` loss term (already built) | MDMA, paroxetine, clarithromycin, diltiazem + N-desmethyldiltiazem, (S)-norfluoxetine |
| **Induction with onset** | Emax term on enzyme *synthesis*; builds/decays over days at `kdeg` | carbamazepine (+autoinduction), rifampicin |
| **Gut-wall first pass** | separate gut enzyme pool → `FG`; oral `F = fa·FG·FH` | all oral CYP3A4 victims; fixes risperidone first-pass gap |

### Explicitly out of scope (kept as fixed factors), with the reason

A mechanistic **transporter/renal** model is *physiologically* correct but **not
defensibly parameterisable** from the available primaries, and building one would trade the
repository's honesty for false precision. Specifically:

- **Lithium–NSAID and lithium–thiazide** stay fixed factors. Lithium reabsorption is
  sodium-coupled in the proximal tubule, so a Na-coupled reabsorption fraction is the right
  *shape* — but the clinical magnitudes are sparse, perpetrator-specific (ibuprofen +15% vs
  indomethacin +30–59%), and sodium-status-dependent, and the hydrochlorothiazide magnitude
  **has no primary number at all** (the existing `×1.25` is already a class-bracketed
  estimate). Any coupling coefficient would be fitted to nothing. **One cheap, honest win:**
  honour the documented thiazide onset-over-days as a ramp on the existing factor — `onset_h`
  already exists on `Interaction` and is currently documentation only.
- **No P-gp sub-model.** For every drug in the library the clinically relevant interaction
  is already carried by CYP (risperidone/paliperidone via CYP2D6; ketamine via CYP3A4 first
  pass; carbamazepine and rifampicin as inducers). The P-gp data are in-vitro IC50 only,
  assay-dependent (ketoconazole spans 1–6 µM across assays), largely unreachable as
  primaries, and adding them risks **double-counting** the CYP effect. Low yield, high risk.
- **Amoxicillin–probenecid (×2)** is the one secretion interaction worth keeping — it has a
  primary number (Paulsen 1989) and a clean single mechanism (OAT block). A fixed clearance
  factor already captures it correctly; no sub-model needed.
- **Levetiracetam–probenecid** and **pregabalin** are out of scope and should be documented
  as deliberately excluded: the levetiracetam effect is on a metabolite the repo does not
  track, and pregabalin's label states no meaningful interaction. This matches how
  `interactions.json` already documents unquantified victims.

If you want a transporter/renal mechanism anyway, that is a research project to generate the
missing parameters, not a modelling task — I would not add it to this engine on today's data.

---

## 4. Architecture

### 4.1 Enzyme pools as first-class objects

Replace the single `CoDrug.inactivates: str` with a mechanism-carrying attachment. Each
**enzyme** becomes a pool with a turnover rate; each **drug** attaches to pools in two
roles, which are independent:

- **victim role:** a share `s_e` of the drug's clearance runs through enzyme `e` (today's
  `enzyme_shares` plus the inactivated enzyme). The drug's clearance scales with that pool's
  net activity.
- **perpetrator role:** the drug perturbs enzyme `e` by one or more of {reversible,
  inactivation, induction}, each with its own constants, driven by the drug's own `Iu(t)`.

A drug can be both (fluoxetine is a CYP2D6 perpetrator and a CYP2D6 victim), and the model
must handle that self-consistently — its inhibition of CYP2D6 slows its own clearance. The
single-enzyme ODE already couples a drug to the pool it destroys; this is the same coupling
generalised to N pools and three mechanisms.

### 4.2 The enzyme ODE, generalised

For each hepatic enzyme pool `e`, amount `E_e` (1 = this patient's baseline):

```
dE_e/dt = kdeg_e·(1 + Σ_j d·Emax_je · Iu_j/(EC50_je + Iu_j))      ← induction (synthesis up)
        − [ kdeg_e + Σ_j kinact_je · Iu_j/(KI_je + Iu_j) ]·E_e    ← baseline + MBI turnover
```

summed over perpetrators `j` present. Reversible inhibition does **not** consume enzyme, so
it is not in the ODE; it is an instantaneous activity multiplier on the pool:

```
A_e(t) = 1 / (1 + Σ_j Iu_j/Ki,u_je)
```

A victim drug `i` with clearance shares `{s_ie}` is then eliminated at

```
factor_i(t) = 1 − Σ_e s_ie·(1 − E_e(t)·A_e(t))
```

This is a drop-in generalisation of `pk_engine.py:851`, which today reads
`factor = 1 + Σ s·(E − 1)` — i.e. the special case `A_e ≡ 1`, one enzyme. **Sanity checks
that must hold:**

- No perpetrator ⇒ `E_e → 1`, `A_e = 1`, `factor = 1` ⇒ the engine reproduces today's
  single-drug curves (the existing `to ~1e-7` equivalence test must still pass).
- Induction only ⇒ steady-state `E_e = 1 + Σ d·Emax·Iu/(EC50+Iu) > 1` ⇒ clearance up,
  exposure down. Onset and offset governed by `kdeg_e`, i.e. **days**, not the inducer's
  half-life — this is the whole point of "induction with onset."
- MBI only ⇒ `E_e = kdeg/(kdeg + Σ kinact·Iu/(KI+Iu)) < 1`, reproducing the current MBI.

### 4.3 Gut-wall first pass

This is the largest structural change and should be its own phase. Add a **gut enzyme
pool** per enzyme (CYP3A4 is the one that matters), driven by the perpetrator's **gut-inlet**
concentration `Ig = fa·ka·Dose/Qent`, which is far higher than unbound plasma — the
mechanistic reason oral CYP3A4 interactions are large and why the `enzymes.json` caveat
about orally-dosed index substrates exists. Oral bioavailability becomes

```
F = fa · FG · FH,  FG = Qgut/(Qgut + fu,gut·CLint,G·E_g),  FH = QH/(QH + fu,b·CLint,H·E_h)
```

(Yang 2007 Qgut model; well-stirred liver). `Qent ≈ 18 L/h`, `QH ≈ 90 L/h`,
`fu,b = fu/(B:P)`, `fu,gut ≈ 1`.

**Pragmatic path.** The engine is parameterised by `t_half` and `Vd`, not `CLint`. Two
options, cheapest first:
1. **FG factor only (recommended first).** Keep systemic CL as the hepatic effect (the
   share mechanism already does this), and add an `FG(t)` multiplier on absorbed dose for
   oral 3A4 substrates, with `FG` driven by the gut pool. This closes most of the oral-DDI
   gap and the risperidone first-pass gap with minimal reparameterisation.
2. **Full `F = fa·FG·FH` from CLint.** Reparameterise 3A4 victims by CLint so both FG and FH
   move mechanistically. More faithful, more parameters, do it only if option 1's validation
   misses demand it.

### 4.4 Perpetrators need PK

A time-varying effect needs `C(t)` for the perpetrator, which means promoting the key
perpetrators to full presets with their own PK. Priority order (by how much dynamics
matter):

1. **fluoxetine + norfluoxetine** — norfluoxetine's ~1-week half-life makes it the
   persistent perpetrator; the pair *must* be two compartments/species.
2. **paroxetine** (MBI, accumulates), **carbamazepine** (+ 10,11-epoxide, autoinduction),
   **rifampicin** (the induction calibrator).
3. **bupropion + 3 metabolites**, **diltiazem + N-desmethyldiltiazem**, clarithromycin,
   fluvoxamine, fluconazole, ciprofloxacin, terbinafine, ketoconazole.

Where a perpetrator has no PK yet, it keeps its **static** factor (today's behaviour) — the
model degrades gracefully to the current engine for un-promoted perpetrators.

---

## 5. Parameters (from the research pass; pull exact values from the cited supplements)

Status key: `verified-ft` full text read · `verified-abs` abstract read · `anchored`
in-vivo-fitted value used by an OSP/Simcyp network · `unconfirmed` from memory, must be
confirmed before use · `PubChem` chemical fact.

### CYP2D6 perpetrators

| Drug | Mechanism | Constants (unbound where noted) | fu / B:P | MW | Source · status |
|---|---|---|---|---|---|
| Paroxetine | MBI | KI 4.85 µM (apparent, total), kinact 0.17 min⁻¹ (10.2 h⁻¹); unbound KI ~0.315 µM | 0.05 / ~1 unconf | 329.4 | Bertelsen 2003 DMD 31:289 `verified-abs`; Rüdesheim 2022 OSP `anchored` |
| Fluoxetine / norfluoxetine | reversible (no MBI on 2D6) | Ki ~0.2–1.5 µM (bufuralol, HLM); norFLX persistent | 0.05 / ~1 unconf | 309.3 / 295.3 PubChem | Stevens & Wrighton 1993 JPET 266:964 `verified-abs` |
| Bupropion + OH/threo/erythro | reversible + 2D6 downregulation | Ki (unbound, fu,mic=0.01) 21 / 13 / 5.4 / 1.7 µM | bup 0.16 / 0.82 | 239.7 / 255.7 / 241.8 / 241.8 | Pharmaceutics 2018 10(1):1 (PMC5874814) `verified-ft`; Sager 2017 |
| Terbinafine | reversible, very high affinity | Ki 22–44 nM (HLM / expressed); quasi-irreversible clinically (4-wk washout) | ~0.001 / unconf | 291.4 PubChem | Abdel-Rahman 1999 DMD 27:770 `verified-abs` |
| MDMA | MBI (already in repo) | kinact/KI fitted to own PK; kdeg = ln2/46.6 h | — | 193.24 | repo; O'Mathúna 2008 `verified` |

### CYP2C19 / CYP1A2 perpetrators

| Drug | Enzyme · mechanism | Constants | fu / B:P | MW | Source · status |
|---|---|---|---|---|---|
| Fluvoxamine | 1A2 reversible; 2C19 reversible + TDI | **1A2 Ki,u ~0.87 nM (in-vivo-anchored)**; 2C19 ~10× more potent in vivo than in vitro; (S)-norFLX 2C19 KI 7 µM, kinact 0.059 min⁻¹ | 0.23 / unconf | 318.3 | Kanacher 2020 OSP `anchored`; Yao 2003 `verified-abs` |
| Fluconazole | 2C19 / 2C9 / 3A4 reversible | 2C19 IC50 12.3 µM; 3A4 Ki 9.2 µM; 2C9 Ki 7.9 µM in vitro (refit 20.4 µM) | 0.89 / ~1 unconf | 306.3 | Niwa 2005 `verified-abs`; Simcyp refit `anchored` |
| Ciprofloxacin | 1A2 reversible, weak in vitro | IC50 ~135 µM in vitro — **underpredicts**; use in-vivo anchor to tizanidine ×10 | 0.6–0.7 / ~1 unconf | 331.3 | Granfors 2004 (clinical) `anchored` |

### CYP3A4 perpetrators, turnover, gut physiology

Full tables, provenance, and the PBPK supplements to pull exact EC50/Emax/kinact/KI from
(Hanke 2018 for clarithromycin & rifampicin; Fuhr 2021 for carbamazepine; Rowland-Yeo 2010
/ Zhang 2009 for N-desmethyldiltiazem) are in the research appendix below.

**Enzyme turnover (kdeg), the single biggest source of MBI/induction uncertainty — mandate
a sensitivity analysis (Yang 2008):**

| Enzyme | t½ (h) | kdeg (h⁻¹) | note |
|---|---|---|---|
| CYP3A4 hepatic | ~36 (10–140) | ~0.0193 | most contested |
| CYP3A4 gut | ~23 | ~0.03 | most influential for oral DDIs, least certain |
| CYP1A2 | ~39 | ~0.018 | |
| CYP2C9 | ~104 | ~0.0067 | long-lived |
| CYP2C19 | ~26 | ~0.027 | |
| CYP2D6 | ~51 (ex-vivo 70–104) | ~0.0136 | Simcyp value; ex-vivo disagrees — flag |

Gut (Qgut model): `Qent ≈ 18 L/h`, gut CYP3A4 ≈ 1% of hepatic abundance yet FG ≈ 0.4–0.6
for midazolam. Liver (well-stirred): `QH ≈ 90 L/h`.

---

## 6. Validation plan

The model is worth exactly what it is measured against, so measurement is built alongside,
in `core/validation.py` (which already runs the reference-curve and interaction back-tests).

1. **Equivalence (regression guard).** A single-drug regimen, and a co-regimen sharing no
   enzyme, must reproduce `concentration_profile` to ~1e-7. The existing coupled-equivalence
   test is extended to the N-enzyme, three-mechanism engine.
2. **Static cross-check.** At steady state with a constant inhibitor level, the dynamic
   engine must converge to `derived_exposure_factor` (`enzymes.json`'s static model) for a
   purely reversible or purely inducing perpetrator. This ties the new code to the already-
   validated static layer.
3. **Clinical DDI gate (the real test).** Predicted AUC ratio within **2-fold (geometric)**
   of observed, for a curated set of studied pairs. Anchoring pairs are excluded from the
   gate; held-out pairs are the score. Target pairs (all with primary clinical numbers):
   paroxetine→desipramine, fluoxetine→desipramine, fluvoxamine→tizanidine &
   →caffeine, ciprofloxacin→tizanidine, ketoconazole/clarithromycin/diltiazem→midazolam,
   carbamazepine→midazolam (& →risperidone, already in repo), rifampicin→midazolam,
   MDMA→risperidone (already in repo). Benchmark to beat: the OSP networks report ~93% of
   pairs within 2-fold.
4. **Time-course, qualitative.** Induction onset/offset over days (carbamazepine,
   rifampicin); MBI persistence after washout (paroxetine, MDMA, terbinafine). No clinical
   time-series to score against for most, so these are asserted directionally, not gated.
5. **Provenance test.** Every new parameter passes `tests/test_data_provenance.py`:
   `verified` ⇒ primary document cited; `anchored` ⇒ the clinical pair and PBPK source
   named; `unconfirmed` ⇒ not shippable, blocks the build.

---

## 7. Phasing

Each phase is independently shippable and leaves the app working.

- **Phase 0 — turnover + provenance.** Add the enzyme `kdeg` table with sources and the
  `anchored` provenance status. No behaviour change. Unblocks everything.
- **Phase 1 — generalise the enzyme ODE.** N enzyme pools, reversible + induction terms on
  top of the existing MBI. Victim-side only (perpetrators still static factors). Gated by
  equivalence + static cross-check.
- **Phase 2 — promote perpetrators to PK.** fluoxetine/norfluoxetine first, then paroxetine,
  carbamazepine, rifampicin. Now the clinical DDI gate can run; anchor + validate.
- **Phase 3 — gut-wall first pass (FG factor).** Closes oral-3A4 and the risperidone gap.
- **Phase 4 (optional) — full `F = fa·FG·FH` from CLint**, only if Phase 3 validation misses
  demand it; and the thiazide onset ramp on the lithium fixed factor.

Recommended stopping point: **end of Phase 3**. That gives dynamic reversible + MBI +
induction, gut first pass, a validated clinical gate, and the risperidone fix — a genuine
mechanistic multi-drug model — without the CLint reparameterisation of Phase 4.

---

## 8. Limitations to state up front (so the numbers don't flatter the model)

- It predicts a studied pair's **average** AUC ratio within ~2-fold, not an individual's.
- Turnover `kdeg` uncertainty (CYP3A4 10–140 h) propagates directly into every MBI and
  induction prediction — the sensitivity analysis is not optional.
- Anchored constants (fluvoxamine, ciprofloxacin, bupropion) are calibrated, not measured;
  they predict *other* pairs but inherit the anchor's assumptions.
- No transporter or mechanistic renal model — those interactions remain fixed factors by
  design (§3), and lithium's magnitude in particular is weakly sourced.
- Gut first pass uses the Qgut approximation, not a segmented gut model.

---

## Appendix — research provenance

The three research passes behind §5 (CYP2D6/2C19/1A2 inhibitors; CYP3A4/turnover/gut;
renal/transporters) and their primary sources are recorded in the session that produced
this plan. Key PBPK supplements to extract exact constants from, none of which are reliably
retrievable without pulling the supplementary tables:

- **Hanke 2018** CPT:PSP 7(10):647 (PMC6202474) — clarithromycin kinact/KI; rifampicin
  per-pathway EC50/Emax and OATP Ki.
- **Fuhr 2021** (OSP, Saarland repo) — carbamazepine CYP3A4 EC50/Emax + autoinduction;
  Kanefendt 2024 doi:10.1002/cpt.3151 for the midazolam AUCR.
- **Rowland-Yeo 2010 / Zhang 2009** (DMD 37:1587) — N-desmethyldiltiazem kinact/KI.
- **Kanacher 2020** Pharmaceutics 12:1191 — fluvoxamine CYP1A2/2C19 in-vivo-anchored Ki.
- **Rüdesheim 2022** Pharmaceutics 14:1734 — paroxetine CYP2D6 MBI network.
- **Pharmaceutics 2018** 10(1):1 (PMC5874814) — bupropion + 3 metabolites, with fu,mic.
- **Yang 2008** Curr Drug Metab 9:384 — enzyme turnover review; **Yang 2007** Curr Drug
  Metab 8:676 — Qgut gut-wall model.
