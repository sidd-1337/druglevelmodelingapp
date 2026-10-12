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
| `max_aafe` | no | fold tolerance on AAFE, 1.5 by default |
| `source` | yes | URL or citation |
| `verified` | no | `true` only once checked against the primary document, not a search excerpt |
| `n_subjects`, `statistic`, `label`, `note` | no | metadata; recorded, not used in scoring |

Supported units are in `core.validation.UNIT_FACTORS`. An unlisted pair raises rather
than guessing — a silent 1000-fold slip would look exactly like a model failure.

## `observed` versus `published_model`

`observed` is the real test: concentrations measured in subjects.

`published_model` is a curve generated from someone else's fitted compartmental
parameters. It is not a measurement, but it is still an external check, because their
model was fitted to data this one has never seen. It is scored and reported
separately so the two are never averaged together. Only use it when the source
reports a **complete** parameter set for the structure it fitted — generating a
one-compartment curve from a paper that fitted two compartments invents a curve the
paper does not support, and tests nothing but the assumption you just made.

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

## Current contents: none

No curve has been seeded yet. Published per-timepoint tables live in journal PDFs,
label PDFs and full-text articles, and the development container reaches none of
them — `curl` and page fetches have no network at all, and search returns only
snippets, which never carry tables. Five targeted searches for paracetamol,
ibuprofen, levetiracetam, methylphenidate and amoxicillin tables returned summary
parameters only.

The harness is therefore complete and tested but unseeded. The tests prove it works
by generating curves from the engine and checking that a deliberately wrong
absorption rate or half-life is caught, so a real curve added later is scored by
machinery that is known to have teeth.

Known-good sources to pull from with full-text access:

* FDA clinical pharmacology reviews on `accessdata.fda.gov` — usually print mean
  concentration-time tables.
* ClinicalTrials.gov posted results for bioequivalence studies.
* Open-access bioequivalence papers (ibuprofen 400 mg and paracetamol 1000 mg both
  have several).
* Amoxicillin: *Antimicrob Agents Chemother* 1980;17(2):199-202, Table 3 — a single
  500 mg oral dose in healthy fasted subjects.
