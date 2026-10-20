# PK model development log

What was built on `claude/pk-model-validation`, why, and what it is measured to be worth.

**Base:** `cfb2bb8` (tip of `master`) · **Head:** `d8d2e9b` · 12 commits · 43 files,
+5 864 / −256 lines · tests 97 → 235 functions (441 passed, 4 skipped)

The starting question was *"how good can predictions from public documentation be, and
can we improve the model?"* It could not be answered, because nothing measured the
model against anything. So the first thing built was the measurement, and every change
after it was justified or rejected by that measurement.

**Verification pass** (`claude/new-session-xnrjfn`, 9 commits on `22b5915`): every
search-excerpt value was checked against its primary document, the curve harness was
seeded with 23 published curves, and the live SÚKL tests now run and pass. It found an
engine bug and two absorption lags the data contradict. Its numbers supersede §1, §6 and
§8; the details are in §10.

---

## 1. Headline result

50 published reference values (Cmax, Tmax, AUC, steady-state averages), scored as
geometric mean fold error — 1.00 is perfect, 1.50 is the tolerance:

| group | n | at session start | end of build | after verification (§10) |
|---|---|---|---|---|
| exposure (AUC, Cavg,ss, Cmin,ss) | 16 | 1.56 | 1.08 | **1.16** (15/16) |
| Cmax | 16 | 1.45 | 1.14 | **1.18** (15/16) |
| Tmax | 12 | 1.41 | 1.10 | **1.28** (10/12) |
| redosing (accumulation vs superposition) | 6 → 4 | — | 1.06 | **1.07** (4/4) |

At the end of the build all 50 were within tolerance. Verification then corrected
values against their primary sources and removed two that only a secondary source
reported. That left 48 checks, 44 within tolerance. The four misses are measured values
the model does not reach, each pinned as a strict xfail. The scores got worse because
the targets became true, not because the model changed. A separate mechanistic back-test
predicts measured drug-interaction factors from enzyme fractions and inhibitor potency:
**4 of 4 addressable pairs within 2×, GMFE 1.11** (5 of 5 at 1.17 before one pair's
inhibitor regimen was corrected and the pair no longer matched).

Whole curves: **13 of 23 published concentration-time curves within 1.5-fold AAFE,
geometric mean AAFE 1.50** (§10).

Reproduce with `python -m core.validation`.

---

## 2. What the model can and cannot predict

Worth stating plainly, because the numbers above flatter it.

**It predicts well:** the *average* patient's exposure and peak from a documented dose,
accumulation on repeated dosing, the direction and rough size of a renal-impairment or
body-weight effect, and the size of a documented enzyme-inhibition interaction.

**It predicts a range, not a value, for an individual.** Inter-individual variability is
real and large — risperidone's active moiety spans ~20-fold across patients. The band is
modelled honestly, which means it is often wide enough to be unsatisfying. That is the
data's fault, not the model's, and narrowing it without evidence would be a lie.

**It cannot predict** absorption in a fed state, hepatic impairment, altered protein
binding, transporter-mediated interactions, or anything about a drug whose label does
not quantify its elimination. Curve *shape* is now tested against 23 published curves, and
it is the weakest part: absorption onset and, for three drugs, the post-peak decline
are wrong in ways summary metrics did not show (§10).

---

## 3. What was built

### 3.1 Measurement (the foundation)

- **`data/reference_pk.json`** — 50 published values across 15 presets, each with a
  source URL, a `verified` flag and a per-value fold tolerance.
- **`core/validation.py`** — scores predictions against them, grouped into exposure,
  Cmax, Tmax and redosing. Printed by `python -m core.validation` and surfaced in the
  Streamlit **Model** tab, so the app states its own accuracy to the user.
- **`validate_interactions()`** — the mechanistic gate described in §3.6.
- **`validate_curves()`** — whole-curve scorer, §6.

### 3.2 Engine (`core/pk_engine.py`)

- **Two-compartment option** with a closed-form solution — a real distribution phase
  instead of one lumped half-life. Used by paracetamol, lithium and mirtazapine.
- **Absorption lag** (`tlag_h`), which fixed several Tmax predictions.
- **Saturable (Michaelis–Menten) elimination** and **mechanism-based enzyme
  inactivation** (`kinact`, `KI`, `kdeg`) — MDMA's disproportionate rise with dose is a
  consequence of its own CYP2D6 inactivation, not a fitted curve.
- **Active-metabolite compartment** — a fifth ODE state. `Metabolite` is a frozen
  dataclass of scalars only, which keeps `PKModel` hashable for its three `lru_cache`s.
- `needs_integration` was separated from `is_nonlinear`: a metabolite breaks
  superposition without being nonlinear, so the profile must be integrated either way.

### 3.3 Patient covariates (`core/population.py`)

Clearance is one chokepoint, split across elimination pathways:

```
CL / CL_ref = Σ_p f_p · a_p  +  (1 − Σ_p f_p) · size
    size     = (WT / 70) ** 0.75          allometric
    a_renal  = CrCl / 100, or size when CrCl is unknown
    a_enzyme = size × enzyme activity     (genotype; 1.0 = normal metaboliser)
```

Renal clearance deliberately gets no allometric term, because CrCl already carries body
size. The unassigned remainder is size-scaled and genotype-blind. Creatinine clearance
comes from Cockcroft–Gault in SI units. The reference patient (70 kg, CrCl 100 mL/min,
normal metaboliser) returns exactly 1.0, and the pathway form was verified algebraically
identical to the renal-only form it replaced to 2e-16 over a weight × CrCl × fraction
grid — so no validated prediction moved when it landed.

Nine presets declare pathways; the rest behave exactly as before. Encoded only where
the data is solid:

| drug | pathways |
|---|---|
| paracetamol | UGT 0.50, sulfation 0.40, CYP2E1 0.05, renal 0.05 |
| levetiracetam | renal 0.66, hydrolysis 0.24 |
| pregabalin | renal 0.90 |
| lithium | renal 0.95 |
| amoxicillin | renal 0.70 |
| methylphenidate | CES1 0.80, renal 0.01 |
| diazepam | CYP2C19 0.50, CYP3A4 0.40 |
| mirtazapine | CYP2D6 0.36, CYP3A4 0.35, CYP1A2 0.19 |
| risperidone (parent) | CYP2D6 0.78 |

### 3.4 Variability band

Rebuilt on log-normal clearance and bioavailability rather than a half-life range. The
z-shift of total log exposure is split between CL and F in proportion to their
variances, so each edge of the band is the most likely (CL, F) pair for a patient at
that exposure percentile — not an impossible corner case. Volume is held constant, so
half-life moves with clearance.

### 3.5 Genotype (`CYP2D6`, `CYP2C19`, `CYP2C9`)

Phenotype acts as a patient covariate on its own pathway, with activity multipliers
back-derived from AUC ratios of sensitive substrates. CYP3A4 and CYP1A2 get inhibition
but **no** genotype selector, because neither has a validated phenotype system —
offering one would be fiction.

When the phenotype is *unknown*, the band widens in quadrature by the spread of the
population phenotype mixture, without moving its centre. Stating a phenotype narrows
it. Only drugs that actually carry a polymorphic pathway are affected.

### 3.6 Interactions (`core/interactions.py`, `core/enzymes.py`, `data/enzymes.json`)

19 measured interactions, plus derived ones for perpetrator × victim pairs with no
measured row, labelled `estimated` so the UI can distinguish them. 5 enzymes and 13
perpetrators in the registry.

Two design points that matter:

- **Residual activity is derived, never stored.** Reading it off the FDA class band
  ("strong ⇒ 0.2") is wrong, because the band describes the effect on an index
  substrate that is not wholly cleared by that enzyme. It is normalised instead:
  `residual = (1/AUC_ratio − (1 − fm_index)) / fm_index`. Paroxetine on desipramine
  (fm ≈ 0.9) gives 0.07, not 0.2. The loader asserts each derived residual falls inside
  its declared FDA class band, which makes the classification enforceable rather than
  decorative.
- **A pre-registered falsifiable gate.** Before building any of it: *if the mechanism
  cannot reproduce the measured factors within 2×, do not ship derived factors.* The
  expectation was recorded in advance (3–5 of 19 pairs addressable; most of our
  measured pairs are renal, first-pass or lack a usable fm). Outcome: 5 pairs
  addressable, all within 2×, GMFE 1.17.

### 3.7 Bayesian individualisation (`core/individualize.py`)

Maximum *a posteriori* fit of two random effects (CL/F and V) to measured blood levels,
by Laplace approximation. Measured levels win; the mechanistic model fills the gaps.

### 3.8 Interface

- FastAPI endpoints for simulation, harm reduction, co-medication and validation.
- Streamlit app: genotype selector, TDM fit tab, co-medication tab, metabolite curve,
  validation table, bilingual (English/Czech, key parity enforced by test).
- The Next.js portal was **routed through the API** rather than mirroring the model in
  TypeScript (§5).

---

## 4. Findings, and the places I was wrong

The useful part of the log. Each of these was caught by the measurement, not by review.

**Claimed a number before verifying it.** I asserted that a diazepam→nordazepam
conversion fraction of 0.5 reproduced the measured steady-state ratio of 1.26. It gives
0.83. Solved numerically instead: 0.761 — and the data note now says the value is
*calibrated, not observed*, with a test pinning it.

**Cut a variability parameter on an assumption.** During the risperidone parent/metabolite
split I dropped `cv_cl` from 0.6 to 0.35, assuming genotype would absorb the difference.
It does not. Re-derived from the published 20-fold spread: total log-SD 0.91, of which
CYP2D6 explains 0.23 — about 6% of the variance. Within-phenotype `cv_cl` ≈ 1.08.

**The band moved only the parent.** For a drug reported as an active moiety (parent +
metabolite), shifting only the parent's clearance left risperidone's band at 1.3-fold
against a published ~20. The band now moves the metabolite too when it is part of the
reported quantity → 19.8-fold.

**The static interaction formula is invalid for sum-reported victims.** Risperidone pairs
failed at 2.28× and 2.61×. Those pairs are now *simulated* rather than computed from the
static formula → 0.90× and 1.03×, and the gate went from 3 pairs at 1.25 to 5 at 1.17.

**Research contradicted the plan, twice, after it was approved.**
- *Risperidone's active moiety is CYP2D6-independent.* Poor metabolisers form less
  9-OH but clear parent far slower, and the sum barely moves. The plan's headline
  justification — "the ~20-fold spread is mostly CYP2D6" — was wrong. The spread is real
  but comes from elsewhere. The plan was rewritten and the phase order changed.
- *Propranolol must not get genotype.* Its label gives clean CYP2D6/1A2/UGT fractions,
  but poor metabolisers show no difference in oral clearance, because first-pass
  extraction dominates. Encoding the fractions would predict an effect documented not to
  exist, so propranolol has no pathways.

**Pushed back on the premise of a request.** Asked to model two drugs sharing a metabolic
path as interacting, the answer is that they generally do not: substrate competition is
negligible at therapeutic concentrations, which sit far below Ki. Only inhibition,
induction and genotype move enzyme activity. Encoding "shared pathway ⇒ interaction"
would have invented interactions that do not exist.

**Then nearly made the structural error the harness warns about.** Having built the curve
harness, I proposed correcting diazepam's `ka` from 4.0/h to 1.3/h, from a published
absorption half-life of 32 min. Checked first: it puts Tmax at 3.37 h against the
label's measured 1.25 h, a 2.7-fold miss, where 4.0/h gives 1.37 h. That study fitted
*two* compartments, where the distribution phase carries part of the early rise. **A
published parameter is only valid inside the structure it was fitted in.** No change
made; the worked example is recorded in the harness README.

**Two silent data errors surfaced by encoding pathways.** Levetiracetam's renal fraction
was 0.75 against a label saying 66%; pregabalin's was 0.98, which was % of *urinary
radioactivity* — the correct fraction of dose is 90%. Both change renal-impairment
predictions and no reference check covered them.

**One pre-existing bug in the TypeScript engine:** it served a wrong pregabalin renal
fraction. Fixed by deletion — see §5.

**A source inconsistency, reported rather than fitted around.** The risperidone label says
paroxetine raises the active moiety 1.3–1.8×, and also that poor metabolisers match
extensive metabolisers. Inhibition is a phenocopy of poor metabolism, so both cannot be
true. The mechanism predicts ≈1.0. This is recorded, not tuned away.

---

## 5. Two structural decisions

**The portal no longer has its own model.** The Next.js app previously reimplemented the
kinetics in TypeScript, which drifted and served at least one wrong number. It now calls
the Python API, with the local engine kept only as a visible, labelled fallback while
loading or when the backend is unreachable. Metabolites, enzyme pathways, genotype and
enzyme inactivation cannot be represented in the TypeScript engine at all, so mirroring
them was never going to work.

**Live-network tests now skip instead of failing.** Four tests call the live SÚKL
open-data API. They are genuine integration tests and were kept, not mocked — but they
cannot pass behind this container's proxy, and a suite with four permanently red tests is
not a usable signal. A session-scoped probe asks once whether SÚKL answers; if not, they
skip with the reason stated. If it answers, they run exactly as before, so a *wrong*
answer still fails.

---

## 6. The curve harness

Summary metrics have a blind spot: a model can match Cmax, Tmax and AUC and still draw
the wrong curve. Scale `ka` down 8× and AUC is **identically** unchanged, Cmax stays
close, and no existing check notices.

`validate_curves()` closes it. Curves are JSON files in `data/reference_curves/`, found
by filename with no code change; the README is the format spec. Each is scored point by
point for **AFE** (bias), **AAFE** (precision — the GMFE of the PBPK literature, and the
same quantity the summary rows report, so the two are comparable) and the fraction within
2-fold.

Three deliberate choices: units are converted from the paper's to the drug's by an
explicit table and an unlisted pair *raises*, because a silent 1000-fold slip would look
exactly like a model failure; measured curves (`observed`) are scored separately from
curves generated from someone else's fitted parameters (`published_model`) and never
averaged together; and where a preset's absorption lag makes the model predict exactly
zero, scoring stops and names the times rather than averaging away an infinite fold
error.

The harness was built with no network, so it was first **proven rather than
demonstrated**: its tests build curves from deliberately wrong models — absorption 10×
slow, half-life 4× long, volume 3× small, and the AUC-neutral `ka` error above — and
assert each is rejected. The verification pass then seeded it with 23 published curves
(§10). Its very first run stopped, as designed, on the two presets that carried a lag.

---

## 7. Open items

| item | status |
|---|---|
| **Curve-shape findings** | Seen in §10, not acted on. No preset models the 0.2-0.5 h absorption delay seen for clonazepam, methylphenidate, mirtazapine and d-amphetamine. Paracetamol's ka was fitted with a lag that is now gone (Tmax 3.1× early). Amoxicillin's ka 0.55/h is too slow. Lithium and mirtazapine decline too slowly late on. Ibuprofen runs about 30 % low at 200 mg. Each needs a structural decision, not a re-fit. |
| **28 values still unverified** | 11 have primaries that conflict, 12 have no reachable primary, and 5 were corrected from a primary but keep a field no source states. Each note says which (§10). |
| Risperidone curve | Risperidone was a priority, but no published mean concentration-time table was found in about 30 sources. |
| Portal fallback | `frontend/lib/sukl-data.ts` still carries diazepam's wrong default code (`0002479`, an antihistamine) and the dead `download.php` SPC links. Reported, not changed. |
| **CI does not run** | `.github/workflows/ci.yml` triggers on `push: branches: [main]`, but the remote default branch is `master` and no `main` exists. The one recorded run predates the rename. No commit on this branch has been through CI. |
| Known data gaps | Closed as far as the sources allow (§10). Nordazepam potency has only rat in-vitro binding data (not adopted). 9-OH-risperidone's volume has no primary, and its first-pass formation is not modelled. |
| Not modelled | protein binding (recorded per drug, no UI — every TDM target here is a *total* concentration, so free fraction would change no prediction), hepatic impairment, food effects, transporters. |
| Deliberately not encoded | pathways for clonazepam, ketamine, MDMA, ibuprofen, d-amphetamine, propranolol — each for a stated reason: unresolved sources, concentration-dependent shares, urine-pH dependence, or a mechanism already modelled dynamically. |

---

## 8. Verification

```bash
pytest -q                       # 571 passed, 6 xfailed, 0 skipped (live SÚKL included)
python -m core.validation       # 48 reference checks + interaction gate + 23 curves
cd frontend && npx tsc --noEmit -p . && npx next build
```

During the build, the backward-compatibility contract held at every step: the reference
patient stayed numerically identical, and no fold error in the 50 checks worsened. The
verification pass deliberately broke that contract where a target or the engine was wrong.
Each worsening is listed in §10.

---

## 9. Commits

| commit | subject |
|---|---|
| `440892c` | Validate PK presets against published data and extend the model |
| `a81719a` | Add interactions, two-compartment, CYP2D6 inactivation and formulations |
| `4b27523` | Generalise clearance to elimination pathways |
| `4d64621` | Add enzyme registry and a mechanistic back-test of interactions |
| `4e17e54` | Add an active-metabolite compartment to the engine |
| `264a023` | Simulate nordazepam alongside diazepam |
| `99ef6b0` | Split risperidone into parent and 9-hydroxyrisperidone |
| `db1329e` | Serve the web portal's kinetics from the backend |
| `b980179` | Add metaboliser phenotype as a patient covariate |
| `a7300b7` | Skip the live SÚKL tests where the network is blocked |
| `74e809f` | Add a whole-curve validation harness |
| `d8d2e9b` | Record why diazepam's ka must not be taken from the published fit |
| `bd0b1e7` | Document what this branch changed and what it is measured to be worth |
| `22b5915` | Stop showing ibuprofen's product data under paracetamol |
| `ee4c276` | Point SPC links at the live SÚKL document service and fix two wrong presets |
| `4c8590e` | Form metabolites in their own mass units, not an inverted molar factor |
| `69188c0` | Verify 104 of 111 data values against their primary documents |
| `ae7d10f` | Verify the CYP1A2 and CYP3A4 index-inhibitor entries |
| `c7bda28` | Let curve files convert mass to molar through a stated molar mass |
| `7d7b21c` | Remove the absorption lags the measured curves contradict |
| `6f5f7b5` | Seed 23 published concentration-time curves for 12 presets |
| `736101e` | Make the live SÚKL tests assert, and retry a dropped connection |

---

## 10. Verification pass

Three jobs, in order: verify every value that came from a search excerpt, seed the
curve harness with published tables, and run the live SÚKL tests. The rules were that
measured data wins over the model and no parameter is tuned to pass a check. A
correction that makes a score worse is reported, not reverted. No new modelling layer
was added.

### 10.1 Data: 111 values

One agent checked each value against its primary document: the SÚKL SPC, a regulator
label, the study's full text, or registry results. A second agent was then told to
refute it, and disagreements were settled by reading the source. An abstract counts only
when it states the number with its definition, and the note then says "abstract only".
Aggregators and label mirrors never count; `tests/test_data_provenance.py` enforces this.

| outcome | n |
|---|---|
| verified | **81** (44 matched directly, 17 by recomputed derivation, 20 after correction) |
| not verified: primaries conflict | 11 |
| not verified: no reachable primary | 12 |
| corrected from a primary, but a field no source states | 5 |
| removed | 2 (Peiró 2013 MDMA redosing: it gave 50 then 100 mg, and its tables do not reproduce its own +16 %) |

The message of `ae7d10f` gives the not-verified split as 8 / 11; that is the split for
the first 104 values only. The table above covers all 111.

Corrected reference values, pred/obs before → after (nothing re-fitted):

| value | correction | pred/obs |
|---|---|---|
| lithium Cmax 300 mg | 0.28 was the **fed** value; fasted is 0.40 mmol/L | 0.72 → **0.50** (xfail) |
| lithium Tmax | 2.5 → 1.21 h, same table | 0.77 → **1.60** (xfail, just outside 1.6) |
| lithium cavg_ss | 0.85 was the midpoint of an ER dosing target; now the measured 300 mg AUC0-inf 4.97 | 0.80 → 1.10 |
| mirtazapine 0.0368 | a TDM median **trough**, so cmin_ss, not cavg_ss | 0.87 → **0.52** (xfail) |
| diazepam AUC | 5.8 applied F to an oral clearance; Hogan 2020 measured 8.27 | 1.10 → 0.77 |
| risperidone | cavg_ss 0.028 was from long-acting injectables; now Health Canada AUC 0.132 per 1 mg | 1.07 → 1.39 |
| livizux | restated at 70 mg from NCT01890785 (n=30) | AUC 0.99 → 1.12 |
| MDMA | 75/125 mg Cmax were secondary restatements; 100 mg was 2 pilot subjects, now n=8 | 100 mg Cmax and AUC 0.80 |
| clonazepam, methylphenidate Tmax | 3.0 → 3.3 h; 1.75 → 1.8 h (test and reference arms had been averaged) | within tolerance |

Model data corrected:
- paracetamol pathways from Critchley 1986 (n=111);
- mirtazapine CYP1A2 0.19 → 0.245;
- risperidone CYP2D6 0.78 → 0.80;
- diazepam's nordazepam formation fraction 0.761 (calibrated) → 0.53 (observed, Greenblatt 1988).

Interaction corrections:
- risperidone + carbamazepine 0.51 (US label);
- amoxicillin + probenecid 2.0 in healthy adults, not 1.37 from children with meningitis;
- cimetidine acts on propranolol and mirtazapine through bioavailability, not clearance;
- MDMA + paroxetine 1.27;
- rifampicin's index ratio 0.18 → 0.04. The 0.18 was a PBPK simulation; Backman 1996 and
  the midazolam label measured 0.04. The derived diazepam estimate becomes ×0.09 against
  the SPC's roughly ×0.25. This is the oral-index gut-wall confound the file already
  warns about, now visible.

Fluvoxamine's CYP2C19 potency was measured at 50 mg/day, not the 100-150 mg/day its name
claimed. It was renamed, so the diazepam pair dropped out of the back-test: 5 pairs at
1.17 → 4 pairs at 1.11. The dropped pair was the worst one, at 0.68.

**The five gaps.**
- **Diazepam → nordazepam fraction:** observed, 0.53. It gives a steady-state ratio of
  0.79 against the measured 1.26: a strict xfail, not re-solved.
- **Nordazepam potency:** only rat in-vitro binding data (0.85) exist; not adopted.
- **9-OH:risperidone AUC ratio:** the primaries give about 6, not "about 5"; the note is
  corrected.
- **9-OH volume:** no primary reachable. First-pass formation of 9-OH is not modelled, a
  structural limit.
- **AGNP ranges:** total concentrations, and troughs unless a row says otherwise
  (methylphenidate's are Cmax-based). The `drugs.json` header now says so.

### 10.2 Engine and SÚKL defects found on the way

**The metabolite was formed in the wrong mass units** (`4c8590e`). The formation term
multiplied parent mass by MW_parent / MW_metabolite. That is neither the metabolite's mass
nor parent-equivalents, and the mass-balance test encoded the same factor. Nordazepam was
11 % too high and 9-OH-risperidone 7 % too low. After the fix, two values calibrated under
the bug miss and are strict xfails: diazepam's nordazepam:diazepam ratio, and
risperidone's solved active moiety (+6.8 %).

**SÚKL** (`ee4c276`, `736101e`):
- Every SPC/PIL link was dead: `download.php` now returns an HTML notice. Documents are now
  fetched by id from `prehledy.sukl.cz`.
- diazepam's default code was an antihistamine (DITHIADEN, bisulepin), and amoxicillin's
  was co-amoxiclav. risperidone's had no SPC on file. paracetamol is re-added as PARALEN
  500 MG. Each default code now declares the ATC it must resolve to, and a live test
  checks it.
- Nine live tests wrapped their assertions in `if results:`, so they passed when SÚKL
  returned nothing. The `736101e` message says seven; the correct count is nine. They now
  assert.
- Run repeatedly, they exposed about one search in 40 hitting a TLS handshake timeout or
  a reset. The client reads any failure as "no data", so the app showed one Livizux
  package instead of two. Each read-only request now retries once.

All 27 live tests pass.

### 10.3 Curves: 23 published tables, 12 presets

Each table was transcribed twice, independently, from the primary document, and every
value agreed. Each curve uses the reference arm and the paper's own unit. Lithium is
printed as mass of Li⁺, so `molar_mass_g_per_mol` was added; methylphenidate in pg/mL
got a conversion factor. Risperidone, a priority, has no published mean table in about
30 sources checked.

**The first run stopped on the two absorption lags** (`7d7b21c`). Paracetamol (0.44 h) and
ibuprofen (0.5 h) are measured well before their lags. Paracetamol is 2.6 ug/mL at 15 min
(1.6 ug/mL at 5 min in a second study), and ibuprofen 3.8-4.4 ug/mL at 15 min. Both lags were removed and nothing was re-fitted. Paracetamol's ka was fitted
together with its lag, so its Tmax is now 3.1× early, a strict xfail. Ibuprofen's Tmax
went 0.96 → 0.63, still inside its 1.6 tolerance. The Tmax group GMFE moved 1.13 → 1.28.

**Scores: 13 of 23 within 1.5-fold AAFE, geometric mean 1.50.**
`data/reference_curves/README.md` has the per-curve table. The findings:

- **Match:** levetiracetam (1.05) and d-amphetamine from 1 h on (1.21). Clonazepam and
  methylphenidate are within 1.5× from 1 h on. Propranolol matches Salem 1984
  (1.15, 1.31), while Parsons 1976 disagrees with both Salem and the model.
- **Absorption onset is wrong for five presets.** Clonazepam, methylphenidate, mirtazapine
  and d-amphetamine are 3.5-29× over-predicted in the first 20-30 min, a 0.2-0.5 h delay
  no preset models. Paracetamol shows the opposite: it is measurable from 3 min.
- **Post-peak shape is wrong for three, consistently across independent studies.**
  - Amoxicillin (3 studies): ka 0.55/h is below k, so the modelled decline is
    absorption-limited, and the curve is 1.8-3.8× high from 5 h.
  - Lithium (2 studies): low early, 1.4-1.7× high from 36 h.
  - Mirtazapine (2 studies): 1.6-4.8× high from 48 h.

  The summary metrics saw none of this.
- **Exposure:** ibuprofen at 200 mg runs about 30 % low in both FDA studies.
  Paracetamol's GSK 204779 study is 1.5-3.9× over-predicted, but that study's AUC is about
  60 % of the literature's.
- **Suspect source:** the pregabalin table (0.43×) prints about twice Pfizer's own
  concentrations and non-zero "pre-dose" samples. It is seeded with `verified: false`.

None of these was tuned away. Each is a decision for a later modelling pass (§7).

---

## 11. Co-medication: several drugs, interactions over time, effect clashes

Built on `develop` from `453b042`. The question: what happens to each drug's level when
several are taken together, and when do two drugs add up in effect without either level
changing? Until then the co-medication tab drew each drug in its own chart and applied
every interaction as a constant factor from the first hour.

### 11.1 What the model now does

All of it sits in one function, `core.regimen.simulate_regimen`. The Streamlit tab and
`POST /api/v1/pk/multi-drug` both call it, so the lab and the API give the same answer.

- **One scale for all drugs.** Each curve is shown as its position in its own
  therapeutic window, (C − MEC) / (MTC − MEC). Drugs in mg/L, ng/mL and mmol/L then share
  one axis, and the zones line up.
- **A shared enzyme pool.** An enzyme a drug in the regimen destroys (MDMA → CYP2D6) is
  one ODE state for every drug cleared by it. Its activity follows the destroyer's own
  level and recovers at the enzyme's turnover rate (`core.pk_engine.coupled_profiles`). A
  single-drug run reproduces the existing engine to about 1e-7.
- **Interaction timing.** An added perpetrator keeps its magnitude: the measured AUC
  ratio where a study exists, the mechanistic estimate otherwise. That magnitude is now
  reached over time, factor(t) = 1 + (F − 1)·g(t) (`core.perpetrator_timecourse`).
  - For a reversible inhibitor, g(t) builds at the perpetrator's own half-life.
  - For an inducer, g(t) follows an enzyme pool, dE/dt = kdeg·(S(t) − E), so induction
    lags the inducer by days.
  - Fading uses a separate half-life where the effect outlasts the drug.
  - The default is "already being taken", which is what a static factor assumed.
- **Effect clashes.** Four effect classes: serotonergic, CNS-depressant,
  QT-prolonging and sympathomimetic. They come from primary labels (or, for MDMA, the
  original studies). A clash is flagged while two tagged participants are active
  together:
  - a simulated drug counts as active above its MEC;
  - an added substance counts as active from half its steady-state exposure.

  A participant that was never assessed is listed as such, never left out.

### 11.2 Decisions where the answer could have gone the other way

- **Exposure, not one clearance.** The factor scales the victim's reported exposure
  curve.
  - Scaling the parent's clearance turned risperidone + paroxetine into about 1.1×,
    because the parent-plus-9-OH sum is buffered. The Risperdal label's measured
    active-moiety ratio is 1.6×. Applied to exposure, any victim reproduces its measured
    factor exactly.
  - Cost: the victim's own re-equilibration half-life is not modelled. The API's
    `t_half_h` is therefore the drug's own.
- **Timing calibrated to the verified magnitude, not to in vitro Ki.** At steady state
  the residual activity is exactly the one `data/enzymes.json` already derives and
  verifies. Published Ki values span about 100-fold between test systems, and would have
  contradicted those residuals.
- **Two mechanisms on one enzyme.** An added inhibitor on an enzyme that a regimen drug
  destroys goes into that pool as an activity multiplier, not as a factor on top.
  - Stacking the two gave risperidone up to about ×8 with MDMA + paroxetine. The ceiling
    for its CYP2D6 share of 0.8 is 1/(1 − 0.8) = ×5.
  - The cost is that this pair then uses the mechanistic size, not the measured 1.6×.
- **Bioavailability.** A bioavailability interaction does nothing to an IV dose, and is
  capped so F stays ≤ 1.

### 11.3 Data added

- **Perpetrator PK.** Each of the 13 perpetrators now has its elimination half-life and
  dosing interval, from its FDA/DailyMed label.
- **Enzyme turnover.** CYP3A4 turnover is 36 h (kdeg 0.0193 /h, Fan et al.,
  PMC6310869). This sets how fast induction builds.
- **Fade half-lives:**
  - paroxetine 70 h and fluoxetine 168 h: the CYP2D6 inhibition half-life measured after
    stopping (Liston 2002, abstract only);
  - terbinafine 300 h: the label's 200-400 h tissue half-life, with the interaction still
    present at 4 weeks;
  - clarithromycin 36 h: a mechanism-based CYP3A4 inhibitor, so it fades as CYP3A4 is
    resynthesised.
  - Bupropion has no sourced value and uses its own half-life.
- **Aliases.** "Carbamazepine" and "Fluvoxamine 100-150 mg/day" are measured pairs
  recorded under other names. Both are now linked to their PK entry.
- **PD tags.** `data/pd_effects.json` covers all 15 presets (11 tagged, 4 none) and the
  18 added perpetrators (9 tagged, 9 none). Every tag is `verified: true` against the
  primary document and covered by `tests/test_data_provenance.py`.
  - A symptom that appears only in an adverse-event table is not tagged. That is why
    pregabalin has no QT tag and ibuprofen's NSAID hypertension is not counted as
    sympathomimetic.
- **MDMA's CYP2D6 recovery** stays at 46.6 h (O'Mathúna 2008), now with its source. An
  independent study found 27.6-36.6 h (Yubero-Lahoz 2011), so the model sits at the slow
  end and errs toward overstating how long a substrate stays raised.

### 11.4 Checked numbers

| case | result |
|---|---|
| lithium + ibuprofen (measured 1.15×) | 1.15× |
| risperidone + paroxetine, already taken (measured 1.6×) | 1.60× |
| risperidone + paroxetine started day 5 | 1.00× on day 3, 1.41× on day 6, 1.60× by day 20 |
| diazepam + rifampicin started day 0 | 0.60× day 1, 0.29× day 3, 0.12× day 7, 0.09× day 20 |
| risperidone + MDMA, 72 h | 1.33× (CYP2D6 to about 5 % after one 100 mg dose) |
| ketamine IV + clarithromycin (bioavailability 2.6×) | 1.00× (no first pass to inhibit) |

### 11.5 Known gaps

- **Not built:** per-indication target windows.
- **Fluoxetine's build-up** uses the parent's 120 h, though norfluoxetine accumulates
  more slowly still.
- **Terbinafine's fade** comes from its tissue half-life. The label gives persistence
  (at least 4 weeks), not a half-life for the inhibition.
- **Clash presence is a proxy.** QT and sedation risk rise with concentration, and some
  effects outlast the level.
- **No timing for some perpetrators.** Perpetrators without a PK entry (phenytoin,
  cimetidine, omeprazole, probenecid, hydrochlorothiazide, ACE inhibitors) still act
  from the first hour. All of them are PD-assessed, and none carries one of the four
  effect classes, so this affects only their level effect, not the clash check.

### 11.6 Portal page

`/kombinace` in the Next.js portal calls the same `simulate_regimen` through
`POST /api/v1/pk/multi-drug`. Two things were added to the API for it:
`GET /api/v1/pk/perpetrators`, which lists the substances that can be added, and
`output_points`. The simulation still runs on the fine grid, but the curves come back
resampled: a full regimen is about 100 kB instead of about 600 kB.

The page has no offline fallback, because the portal's TypeScript engine cannot run the
coupled enzyme model or the interaction time course. When the service does not answer, the
page says so. The regimen is stored in the link (`?lek=…&latka=…`), with substance names as
repeated parameters because one of them contains a comma. Each drug page links there with
that drug preselected.

The six series colours were validated against the portal's own surfaces. Three light-mode
steps fall below 3:1, so the chart always has a value list in each drug's own unit.

## 12. Several substances on one enzyme

Built on `develop` from `599b045`, on `claude/ddi-model-on-develop`. It carries over the
part of the dynamic interaction plan (`docs/multi-drug-ddi-plan.md`) that develop did not
already have. The plan was written on a branch cut before §11; that branch is retired.

### 12.1 The defect

The regimen multiplied every added substance's factor. Two inhibitors of one enzyme
could then push a drug past what losing the whole enzyme can do, which is 1 / (1 − fm):

| regimen (all already taken) | multiplied | whole-enzyme ceiling | now |
|---|---|---|---|
| mirtazapine + paroxetine + fluoxetine + terbinafine (CYP2D6, fm 0.36) | 3.42× | 1.56× | 1.54× |
| diazepam + fluconazole + fluvoxamine 50 mg (CYP2C19, fm 0.5) | 3.83× | 2.0× | 1.98× |
| diazepam + ketoconazole + clarithromycin (CYP3A4, fm 0.4) | 2.59× | 1.67× | 1.65× |
| risperidone + paroxetine + fluoxetine (both measured) | 2.24× | — | 1.65× |
| diazepam + ketoconazole + carbamazepine (inhibitor + inducer) | 0.80× | — | 1.58× |

### 12.2 What the model now does

- **One activity per enzyme** (`core.perpetrator_timecourse.combined_activity`).
  - Competing inhibitors add their occupancy: A = E / (1 + Σ βⱼ gⱼ), βⱼ = (1 − rⱼ) / rⱼ.
    The product of their activities would overstate the combination.
  - Inducers add to one synthesis term in one pool, dE/dt = kdeg·(1 + Σ (rⱼ − 1)·hⱼ − E).
  - Inhibition acts on the induced pool. That is the net-effect model, where the
    inhibition and induction terms multiply on activity (Fahmi 2008, DMD 36:1698).
  - With one substance this is exactly the existing curve, to machine precision.
- **The drug's factor on that enzyme.** Each substance keeps its own factor and time
  course, measured where a study exists. When two or more act on one enzyme of the drug,
  their product is multiplied by the mechanism's ratio:
  F_mech(A_combined) / Π F_mech(Aⱼ).
  - F_mech is the existing `mechanistic_exposure_factor`, now also over an array
    (`mechanistic_exposure_response`). For a parent-plus-metabolite sum it is tabulated at
    11 activities from 0.001 to 100 and interpolated in log-log. Worst error is 0.9 %,
    and the tabulation costs about 2 s once per drug and enzyme.
  - The enzyme chart and an enzyme a regimen drug destroys (MDMA) use the same combined
    activity.
- **Shown, not hidden.** The lab, `POST /api/v1/pk/multi-drug` (`applied[].pooled`,
  `series[].pooled`) and `/kombinace` state each pooled enzyme as "×combined, not
  ×multiplied".

### 12.3 Decisions where the answer could have gone the other way

- **Correct the measured product, don't replace it.** Replacing combined pairs with the
  mechanistic size would be simpler. But adding a weak second substance would then
  silently swap risperidone's measured 1.6× for a different number. As built, a single
  substance gives exactly what it gave before, so every §11.4 number holds, and the
  mechanism only supplies what no study measured: how two substances interact.
- **Inhibitor plus inducer: the inhibitor wins.** Ketoconazole leaves 2.4 % of CYP3A4. On
  an enzyme induced to 365 %, that is still only 8.8 %, so diazepam rises (1.58×). The
  multiplied factors said it falls (0.80×). Neither has a clinical study in this
  repository; the net-effect model is the standard basis for such predictions, but it is
  a prediction.
- **Only substances with known PK pool.** All 13 registry perpetrators have it. A
  measured perpetrator outside the registry (phenytoin, cimetidine, omeprazole) has no
  enzyme assigned, so it stays a factor.

### 12.4 The interaction gate lost a pair

`validate_interactions` looked perpetrators up by name only. When "Fluvoxamine
100-150 mg/day" became an alias (§11.3), diazepam + fluvoxamine silently left the gate.
It now resolves aliases and a test pins the pair in. The gate is back to the
5 pairs / geometric mean fold error 1.17 its own comment recorded (5/5 within 2×;
diazepam + fluvoxamine predicted 1.94× against 2.86× measured).

### 12.5 Known gaps

- A substance listed twice under its name and an alias ("Carbamazepine" and
  "Carbamazepine 400 mg b.i.d.") counts twice. It did before as well.
- Risperidone + carbamazepine (0.51× measured) is outside the gate: the model gives
  risperidone no CYP3A4 pathway, so the mechanism cannot address it.
- Gut-wall first pass and per-enzyme turnover beyond CYP3A4 are not built (plan §7,
  Phases 3 and 0).
