# Reference concentration-time curves

Whole published curves, scored by `core.validation.validate_curves()` and printed by
`python -m core.validation`. Drop a `.json` file in here and it is picked up; no code
changes, and `tests/test_curve_validation.py` will validate its shape.

## Why this exists alongside `data/reference_pk.json`

`reference_pk.json` checks summary numbers — Cmax, Tmax, AUC, steady-state averages.
A model can get every one of those right and still have the wrong curve: too fast an
absorption rate with a compensating volume, or a missing distribution phase hidden
inside a lumped half-life. Those errors show up only point by point, which is what
this directory tests.

## File format

One curve per file; the filename stem becomes the `curve_id`.

```json
{
  "drug_id": "ibuprofen",
  "label": "ibuprofen 400 mg oral, fasted, n=20",
  "kind": "observed",
  "doses": [{"time_h": 0.0, "dose_mg": 400.0}],
  "route": "oral",
  "body_weight_kg": 70.0,
  "unit": "ug/mL",
  "lloq": 0.1,
  "max_aafe": 1.5,
  "n_subjects": 20,
  "statistic": "arithmetic mean",
  "source": "https://...",
  "verified": false,
  "note": "What was measured, and anything that limits the comparison.",
  "points": [[0.25, 4.1], [0.5, 12.0], [1.0, 24.3], [2.0, 26.1]]
}
```

| field | required | meaning |
|---|---|---|
| `drug_id` | yes | preset id in `data/drugs.json` |
| `kind` | no | `observed` (default) or `published_model` — see below |
| `doses` | yes | list of `{time_h, dose_mg}`; a single dose is `[{"time_h": 0, "dose_mg": X}]` |
| `unit` | yes | the paper's unit; converted to the drug's own `concentration_unit` |
| `points` | yes | `[[time_h, concentration], ...]`, sorted on load |
| `route` | no | `oral` (default), `iv_bolus`, `insufflated`, `smoked` |
| `body_weight_kg` | no | 70 by default |
| `lloq` | no | in `unit`; points at or below it are dropped before scoring |
| `molar_mass_g_per_mol` | only for a mass unit against a molar preset | molar mass of the species the paper weighed: 6.941 for lithium printed as Li+, never lithium carbonate's 73.89 |
| `max_aafe` | no | fold tolerance on AAFE, 1.5 by default |
| `source` | yes | URL or citation |
| `verified` | no | `true` only once checked against the primary document, not a search excerpt |
| `n_subjects`, `statistic`, `label`, `note` | no | metadata; recorded, not used in scoring |

Supported units are in `core.validation.UNIT_FACTORS`. An unlisted pair raises rather
than guessing — a silent 1000-fold slip would look exactly like a model failure. Mass
to molar (ug/mL of Li+ against lithium's mmol/L) also needs `molar_mass_g_per_mol`,
because the unit never says which species was weighed; without it, it raises too.

## `observed` versus `published_model`

`observed` is the real test: concentrations measured in subjects.

`published_model` is a curve generated from someone else's fitted compartmental
parameters. It is not a measurement, but it is still an external check, because their
model was fitted to data this one has never seen. It is scored and reported
separately so the two are never averaged together. Only use it when the source
reports a **complete** parameter set for the structure it fitted — generating a
one-compartment curve from a paper that fitted two compartments invents a curve the
paper does not support, and tests nothing but the assumption you just made.

The same trap applies to lifting a single parameter. A multiple-dose oral diazepam
study in healthy volunteers reports an absorption half-life of 32 min, i.e. ka ≈
1.3/h, against this preset's 4.0/h — an apparent 3-fold error. But that study fitted
two compartments, where the distribution phase carries part of the early rise. Moving
its ka into this one-compartment preset makes it do both jobs and predicts Tmax at
3.37 h against the label's measured 1.25 h, a 2.7-fold miss; the preset's own 4.0/h
gives 1.37 h. The preset is right and the "correction" is an artefact of structure.

## Early samples and absorption lag

Several presets carry an absorption lag (`tlag_h`), and before it the model predicts
exactly zero — an infinite fold error. Rather than averaging that away, scoring stops
and names the offending times. Either the early points belong below `lloq`, or the
preset's lag is too long, and which it is matters: the second is a real finding.

## Scores

* **AFE** — average fold error, `10^mean(log10(pred/obs))`. Bias: >1 the model runs
  high, <1 low. Over- and under-predictions cancel.
* **AAFE** — `10^mean|log10(pred/obs)|`, ≥1. Precision, direction-blind. Called GMFE
  in PBPK papers, and the same quantity `ValidationRow.fold_error` reports for a
  single metric, so curve and summary scores are directly comparable.
* **fraction within 2-fold** — the conventional acceptance measure.

## Current contents: 23 observed curves, 12 presets

Each table was transcribed twice, independently, from the primary document (FDA review
PDFs, GSK study reports, full-text papers and their scanned pages, one patent); every
value agreed. Each file uses the reference or brand arm and the paper's own unit. Pre-dose
samples are left out. Points the source flags as non-representative are also dropped, and
the file's `note` says so. No preset parameter was changed to fit any of them.

`python -m core.validation` scores 13 of 23 within the default 1.5-fold AAFE, with a
geometric mean AAFE of 1.50. A curve outside tolerance is reported, not hidden: the test
suite checks only that every file loads and scores. What the curves show:

| preset | curves (AAFE) | what the curve shows |
|---|---|---|
| levetiracetam | Ruijs 2022 (1.05) | Matches throughout. |
| livizux | NDA 21-977 (1.21) | Matches from 1 h on. Only the 0.5 h point is 5.9x high: d-amphetamine forms from the prodrug later than the model. |
| rivotril | ANDA 74-925 (1.21), 74-940 (1.54) | Within 1.0-1.5x from 1 h on. The first samples (0.33-0.5 h) are 4-29x high. |
| ritalin | ANDA 40-300 (1.42), US 7,247,318 (1.30) | Within 1.0-1.5x from 1 h on, after both doses. The 20 min point is 14x high, and 30 min after the second dose is 2x high. |
| propranolol | Salem 1984 40 mg (1.31) and 80 mg (1.15); Parsons 1976 (2.50) | Salem matches (1.0-1.6x). Against Parsons, at the same 40 mg, the model is 1.7-3.3x high from 1 h on. Parsons used an old fluorimetric assay, so the two sources disagree with each other. |
| lithium | ANDA 76-243 (1.30), 76-121 (1.30) | Both studies agree: 0.6-0.9x at 1-4 h, rising to 1.4-1.7x at 36-96 h. The model peaks too low and declines too slowly late on. |
| diazepam | Swift 1985 (1.24), Bittencourt 1981 (1.54) | Both studies agree: 0.45-0.8x at 1-6.5 h, then near or above 1 from 6-10 h. This is consistent with a distribution phase that a one-compartment preset cannot show. |
| mirtazapine | ANDA 76-119 (1.44), 76-241 (1.63) | 0.5 h is 6-10x high, 0.6-0.9x at 4-12 h, then 1.6-4.8x high from 48 h. The fitted two-compartment tail declines too slowly. Part of this is below-limit zeros in the 15 mg study's 72-96 h means. |
| ibuprofen | ANDA 76-117 (1.45), 75-010 (1.52) | Both 200 mg studies: 0.5-0.9x from 1 h on (AFE 0.72-0.75). The model's exposure at 200 mg runs about 30 % low. |
| amoxicillin | Arancibia 1980 (1.32), ANDA 65-080 (1.82), Gordon 1972 (1.91) | All three: below 1x around the peak, 1.8-3.8x high from 5 h. Absorption is too slow: ka 0.55/h is below k (t1/2 1.0 h), so the decline is absorption-limited. |
| paracetamol | GSK 204779 (1.95), GSK E3870651 (2.10) | 204779 (fasted 500 mg) is 1.5-3.9x high everywhere. That study's AUC is about 60 % of what the preset and other studies give. E3870651 (semi-fed, combination product) matches at 1-3 h but is 2-9x high in the first 20 min. |
| pregabalin | Tujii 2014 (2.30, `verified: false`) | Runs 0.43x throughout. The source is suspect: its values are about 2x Pfizer's own data, and its 'pre-dose' samples are non-zero. |

Two patterns run across presets.

First, the earliest sample is often far above the model's prediction: clonazepam,
methylphenidate, mirtazapine and d-amphetamine are 3.5-29x high in the first 20-30 min.
That is a delay of about 0.2-0.5 h that none of these presets models. Paracetamol shows
the opposite: drug was measurable at 3-15 min, which is why its former 0.44 h lag
(and ibuprofen's 0.5 h lag) were removed. No new lag was added, because that would fit
the presets to these curves.

Second, three presets get the shape after the peak wrong, consistently across
independent studies: amoxicillin (absorption too slow), lithium (late phase too high)
and mirtazapine (tail too long). These are findings for the next modelling pass. They
were not tuned here.

Not seeded: risperidone, a priority. About 30 sources were checked (FDA and EMA reviews,
bioequivalence papers, registries), and none prints a mean concentration-time table for
risperidone or 9-hydroxyrisperidone. The others lack a published table this search
reached.
