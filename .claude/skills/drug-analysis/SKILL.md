---
name: drug-analysis
description: Analyse a drug for this PK/PD modelling app - mechanism of action, receptor profile and effects in the body, pharmacokinetics and pharmacodynamics from primary sources - and turn it into a sourced preset the model can simulate. Use when asked to analyse, research, profile or add a drug (e.g. "analyse quetiapine", "how does X work and what are its PK", "add X as a preset", "model X in the app"), or to audit an existing preset against its sources.
---

# Drug analysis

One run produces two things for one drug:

1. **A monograph**, `docs/drugs/<id>.md`: what the drug does in the body and why,
   its PK and PD with a source on every number, and how those numbers map onto the model.
2. **A preset the app can simulate**: entries in `data/*.json` and `core/sukl_client.py`,
   tests, and the README tables. Every number carries `source`, `verified` and `note`.

The monograph explains; the data files are what the model runs on. They must agree.
Write the monograph from `references/monograph-template.md` and the data from
`references/app-data.md`. Read both before starting.

## Modes

- **Add** (default): the drug is not in `data/drugs.json`. Do every step below.
- **Audit**: the drug is already a preset. Do steps 1-4 *blind*. Until your own values
  are written down, don't open any file that names the drug: its entries in `data/*.json`,
  tests in any file (several pin its numbers), `docs/model-development-log.md`, the
  README tables, `core/sukl_client.py`, and code comments. If `references/app-data.md`
  names the audited drug as the template to copy, use the other template it names. Use
  the raw SÚKL route in step 2, not `core.sukl_client`. Then diff yours against the
  stored values in a table: value, yours, stored, source of each, verdict (agree /
  stored is wrong / mine is wrong / sources conflict / not comparable). Fitted values
  (V, ka) are "not comparable": compare their predictions instead, scoring both
  parameter sets against every primary reference you found. Change stored data only
  where you can name the primary document that settles it, and only once the person
  agrees. Report the diff even when everything agrees: it measures this skill's
  accuracy.

## Ground rules

These come from mistakes recorded in `docs/model-development-log.md`. Each one has cost
this project a wrong number before.

- **Primary sources only.** In order of preference: the SÚKL SPC (Czech product
  information) for the default product; the FDA/DailyMed label; the EMA EPAR or Health
  Canada monograph; the original study (full text, or its abstract if the note says
  "abstract only"); registry results. Aggregators and label mirrors never count, even
  when they are right. The list is in `tests/test_data_provenance.py` (drugs.com,
  DrugBank, UpToDate, medlibrary, rxreasoner, ...). Use them only to find the primary
  source.
- **`verified: true` only after you've read the number in the primary document.** The
  note says where it was read (section, table, page) and quotes or restates it. A value
  recomputed by plain arithmetic from numbers read in a primary (e.g. AUCτ/τ) may be
  verified if the note shows the arithmetic. It's `verified: false` when the derivation
  adds an assumption (no residual activity in poor metabolisers, equal volumes, an oral
  index substrate), uses an unverified input, or is calibrated or solved; also when
  primaries conflict or none was reachable. The note says which applies.
- **Find the definition behind each number.** Fed or fasted, single dose or steady
  state, trough or average, total or unbound, parent or active moiety (and for a sum:
  plain mass sum, as AGNP and the engine use, or a molar-corrected total active moiety),
  test or reference arm, dose as salt or base, which population. The development log's corrections were
  almost all definition errors (a fed Cmax stored as fasted, a trough stored as an
  average, a long-acting-injectable average used for tablets).
- **A published parameter is only valid inside the structure it was fitted in.** A ka
  or volume from a two-compartment fit doesn't transfer to a one-compartment preset.
  Fit the preset to observed summary values (AUC, Cmax, Tmax, troughs) instead.
- **Never tune to pass.** If a sourced value is missed after an honest fit, keep the
  reference and mark the miss as a strict xfail in `tests/test_reference_pk.py`
  (`KNOWN_MISSES`) with the reason. Don't widen `max_fold` to hide it.
- **Encode elimination pathways only where a label mass balance or in-vivo data
  supports them** (an inhibitor or phenotype study). In-vitro shares alone are
  `verified: false`, and the note says so. Leave pathways out when first pass dominates
  and the label shows no phenotype effect (propranolol), or when the share depends on
  concentration.
- **Deriving a pathway share from an inhibitor study counts as in-vivo evidence**, but
  it's circular for that pair. fm = (1 − 1/AUC ratio) / (1 − r), where r is the
  inhibitor's residual activity derived in `data/enzymes.json`. Name that measured
  interaction row for the dose actually given (e.g. `Ketoconazole 200 mg/day`), not the
  registry's name: `validate_interactions()` matches registry names exactly, and a pair
  that produced its own fm would otherwise score as a perfect, meaningless back-test.
  Say so in both notes.
- **A shared enzyme doesn't mean an interaction.** Only inhibition, induction and
  genotype change enzyme activity. Record measured interaction pairs, not inferred ones.
- **Healthy volunteers and patients can differ by 2× in clearance.** Quetiapine: 138
  against 65 L/h in the same paper. Fit the population the app simulates (patients at
  steady state, where therapeutic ranges apply). Cross-check the result against an
  independent value (the label's Vd), and keep the other population's reference as a
  strict known miss.
- **Therapeutic ranges are AGNP 2017 consensus values.** They are total plasma
  concentrations and, unless stated otherwise, troughs. The preset's MEC is the lower
  limit of the range. MTC is the AGNP laboratory alert level, unless a primary source
  names a different adverse threshold for the drug. `mec` and `mtc` have no
  `source`/`verified` fields, so no test catches an unsourced window: the `notes` must
  cite the document, give the full range and the alert level, and say which value MTC
  is and whether it was read in the consensus itself or second-hand.
- **Units:** `drugs.json` concentrations in µg/mL (= mg/L). Divide ng/mL by 1000. AUC is
  in µg·h/mL.

## Workflow

### 1. Identify

Work out the INN, Czech name (SÚKL spelling, e.g. KVETIAPIN), ATC code, class, and the
products on the Czech market (salt, strengths, IR/PR forms). Choose a default SÚKL
product (`kod_sukl`) whose ATC matches the substance, which is currently supplied
(`jeDodavka`), and which has an SPC on file. A wrong default code has twice pointed at
another substance.

To list products, POST `{"filtr": "<term>", "pocet": 50, "stranka": 1}` to
`https://prehledy.sukl.cz/prehledy/v1/dlprc`. The filter matches product names, so
search the product-name root ("quetiapin", or a brand), not the Czech INN: "KVETIAPIN"
returned nothing. `GET https://prehledy.sukl.cz/prehledy/v1/dlp/<kodSUKL>` gives the
product, `jeDodavka` and `dokumenty[]`. The entry with `typDokumentu` "SPC" has the id
for `GET .../v1/dokumenty/<id>` (a PDF). One SPC usually covers every strength of a
product. Prefer a supplied product whose SPC section 5.2 is complete; the originator's
often is, so check it when a generic's is thin.

Before changing anything, run `pytest -q` once and note what already fails. Live SÚKL
tests can fail on data drift that has nothing to do with your change.

### 2. Gather primary documents

Collect the URLs before reading anything in detail:
- the SÚKL SPC for the default product (`python scripts/spc_check.py` only covers
  existing presets, so fetch the document from `prehledy.sukl.cz` directly);
- the US label from DailyMed, plus any separate XR/PR label;
- the AGNP 2017 consensus (Hiemke et al., Pharmacopsychiatry 2018;51:9-62) and its
  2026 update (doi 10.1055/a-2860-7861) for psychotropic drugs. Thieme blocks
  automated access with a bot check; don't try to get around it. Unless the user
  supplies the text, record the range as unverified, say which secondary sources
  disagree, and ask the user to confirm it in the report;
- original PK studies with mean concentration-time tables (for the curve harness),
  mass balance, interaction studies, and receptor binding (Ki) data.

For open-access articles whose PMC or Europe PMC page shows a bot check, the official
API returns the full text:
`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=<PMCID>&rettype=xml`.
If it returns a stub saying the publisher doesn't allow downloading, record the value as
"abstract only".

Keep a scratch list of URL, document type and what it gives. Every value cited later
must come from this list.

### 3. Mechanism and actions in the body (monograph sections 2-3)

- Receptor/target profile: a table of Ki (nM) per target, with the source and assay
  (human cloned receptors preferred). Say which targets are reached at therapeutic free
  concentrations. Compare Ki with the unbound Cmax: a target with Ki far above it isn't
  engaged. Use PET receptor occupancy where it has been measured.
- For an active metabolite, give its own profile. It often explains effects the parent
  can't, as with norquetiapine's NET inhibition.
- Effects by organ system: what each receptor action does (sedation, antipsychotic,
  antidepressant, orthostatic hypotension, weight gain, anticholinergic, QT,
  prolactin...), and at which dose or occupancy it appears. Link every effect to a
  mechanism, and mark where the link is inferred rather than shown.
- Dose-dependent pharmacology: say so when the drug acts differently at low and high
  doses (for example, antihistaminergic only at low doses).

### 4. Pharmacokinetics and pharmacodynamics (monograph sections 4-6)

Extract with full definitions: F (absolute vs relative), Tmax, Cmax and AUC per dose
(IR and PR separately), t½ (parent and metabolites), Vd, CL/F, protein binding,
metabolism (enzymes with fractions and their evidence), excretion (mass balance),
linearity range, food effect, special populations (renal, hepatic, age, sex, smoking),
CYP genotype effects, measured DDIs (AUC ratios, with perpetrator dose), therapeutic
range and alert level, exposure-response or occupancy-response data.

Steps 3-4 run in **Audit** mode too. The diff happens after step 4.

### 5. Map onto the model

Decide the structure first, and write the reasoning down in monograph section 7:
- one or two compartments (two only if published curves show a distribution phase that
  matters at the time scale of dosing);
- which quantity the therapeutic window refers to (`reported`: parent, metabolite or
  sum). A metabolite gets its own compartment only when it's active and its exposure
  isn't proportional to the parent's on the time scale of dosing, or when the range
  refers to it;
- `formulations` for IR versus PR (the first entry is the default and must match the
  references);
- `elimination.pathways` under the evidence rule above;
- `cv_cl` / `cv_f` from a published population spread, never guessed. Prefer a
  population-PK between-subject CV. A percentile range of dose-normalised levels also
  works if its percentiles are stated; a min–max range isn't a percentile range. For a
  sum-reported drug the CV must describe the sum, which is usually dominated by the
  metabolite's clearance. The band already adds the spread of unstated metaboliser
  phenotypes, so a CV that includes genotype effects counts them twice.

For a metabolite reported as `parent`, its formation fraction and volume can't be
separated from plasma data. Set its volume equal to the parent's (an assumption, stated
in the note), solve the formation fraction from the steady-state metabolite/parent AUC
ratio, and check the Cmax ratio as the one independent test. Mass and molar ratios
differ by `molar_ratio`, and labels mix them.

The engine forms a metabolite only from systemic parent, so first-pass formation (as
with risperidone) can't be represented. When the reported quantity's Cmax and the
metabolite/parent AUC ratio can't both be met, keep the ratio (it's structural), and
record the Cmax miss as a strict xfail.

If the summary data show a misfit that only a more complex structure would fix (for
example, steady-state peaks at half the measured value), don't fit that structure to
summary numbers alone. Pin the misfit as a strict xfail test with the numbers in its
reason, and explain it in the monograph.

Then fit. Fix t½ and F from the sources. Choose `vd_l_per_kg` so that F·dose/CL matches
the observed AUC (CL = ln2·V/t½), and `ka_per_h` (plus `tlag_h` only if early
concentrations are truly zero) so that Tmax and Cmax match. Check the result with:

```bash
python .claude/skills/drug-analysis/scripts/preset_report.py <id>
```

It prints the 70 kg model, the single-dose and steady-state numbers against MEC/MTC,
and every reference row and curve for the drug with its fold error.

### 6. Write the data

Follow `references/app-data.md` file by file. At minimum: `data/drugs.json`, one or
more exposure references in `data/reference_pk.json` (a test requires one per preset),
`data/pd_effects.json`, the SÚKL maps in `core/sukl_client.py`, and a test in
`tests/test_drug_database.py`. Add `data/interactions.json` rows for measured pairs
whose perpetrator is already in `data/enzymes.json`, and a curve in
`data/reference_curves/` when a primary prints a mean table.

### 7. Verify

```bash
pytest -q                     # whole suite; live SÚKL tests skip without network
python -m core.validation     # the drug's rows must appear, and misses must be explained
```

Re-read every `verified: true` and confirm the URL is primary and the note says where
the number is. Every Cmax, Tmax and AUC printed in a cited document for the default
formulation should be a reference row, or the monograph should say why it isn't. Then check the monograph against the data. Every number in section 7
must equal the stored one.

### 8. Report

Summarise in chat: the structure chosen and the one decision that could have gone the
other way, the fit (pred/obs per reference), what's verified and what isn't and why,
known misses, and what the model can't represent for this drug. In Audit mode, lead
with the diff table.
