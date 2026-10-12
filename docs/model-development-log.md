# PK model development log

What was built on `claude/pk-model-validation`, why, and what it is measured to be worth.

**Base:** `cfb2bb8` (tip of `master`) · **Head:** `d8d2e9b` · 12 commits · 43 files,
+5 864 / −256 lines · tests 97 → 235 functions (441 passed, 4 skipped)

The starting question was *"how good can predictions from public documentation be, and
can we improve the model?"* It could not be answered, because nothing measured the
model against anything. So the first thing built was the measurement, and every change
after it was justified or rejected by that measurement.

---

## 1. Headline result

50 published reference values (Cmax, Tmax, AUC, steady-state averages), scored as
geometric mean fold error — 1.00 is perfect, 1.50 is the tolerance:

| group | n | at session start | now |
|---|---|---|---|
| exposure (AUC, Cavg,ss, Cmin,ss) | 16 | 1.56 | **1.08** |
| Cmax | 16 | 1.45 | **1.14** |
| Tmax | 12 | 1.41 | **1.10** |
| redosing (accumulation vs superposition) | 6 | — | **1.06** |

All 50 are now within tolerance; at session start several were not. A separate
mechanistic back-test predicts measured drug-interaction factors from enzyme fractions
and inhibitor potency: **5 of 5 addressable pairs within 2×, GMFE 1.17**.

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
not quantify its elimination. Curve *shape* is now testable but untested, because no
published concentration-time table was reachable (§6).

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

## 6. The curve harness, and why it is empty

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

**No curve is seeded.** Published per-timepoint tables live in journal and label PDFs.
This container has no outbound network beyond web search — `curl` and page fetches fail
at DNS — and search returns snippets, which never carry tables. Five targeted extended
searches (paracetamol, ibuprofen, levetiracetam, methylphenidate, amoxicillin) returned
summary parameters only. Generating a curve from a partial parameter set would have
tested nothing but my own structural assumption (§4), so the directory holds only its
README and a list of where to get real data with full-text access.

The harness is therefore **proven rather than demonstrated**: its tests build curves from
deliberately wrong models — absorption 10× slow, half-life 4× long, volume 3× small, and
the AUC-neutral `ka` error above — and assert each is rejected. A real curve added later
is scored by machinery known to have teeth.

---

## 7. Open items

| item | status |
|---|---|
| **Seed a real curve** | Needs one published mean concentration-time table. Everything else is in place. |
| **`verified: false` on all new data** | Every new value came from a search excerpt with a source URL. `scripts/spc_check.py` verifies against SPC section 5.2 but cannot run behind this proxy. |
| **CI does not run** | `.github/workflows/ci.yml` triggers on `push: branches: [main]`, but the remote default branch is `master` and no `main` exists. The one recorded run predates the rename. No commit on this branch has been through CI. |
| Known data gaps | diazepam conversion fraction, nordazepam potency, the 9-OH:risperidone AUC ratio, 9-OH volume of distribution, AGNP total-vs-free wording. |
| Not modelled | protein binding (recorded per drug, no UI — every TDM target here is a *total* concentration, so free fraction would change no prediction), hepatic impairment, food effects, transporters. |
| Deliberately not encoded | pathways for clonazepam, ketamine, MDMA, ibuprofen, d-amphetamine, propranolol — each for a stated reason: unresolved sources, concentration-dependent shares, urine-pH dependence, or a mechanism already modelled dynamically. |

---

## 8. Verification

```bash
pytest -q                       # 441 passed, 4 skipped (the 4 are live-SÚKL, §5)
python -m core.validation       # 50 reference checks + interaction gate + curve table
cd frontend && npx tsc --noEmit -p . && npx next build
```

The backward-compatibility contract held at every step: the reference patient stayed
numerically identical, and no fold error in the 50 checks worsened.

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
