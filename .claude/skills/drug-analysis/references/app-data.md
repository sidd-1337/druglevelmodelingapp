# Writing a preset into the app

Each data file starts with a `_comment` that is its spec, so read it first. The
provenance rule is the same for every block that carries `verified`, and
`tests/test_data_provenance.py` enforces it: `source` is an `http` URL to a primary
document, `note` says where the value was read, and aggregator domains fail the test.

## `data/drugs.json`

Copy the shape of the closest existing preset. A metabolite plus CYP pathways is
`risperidone`, `diazepam` or `quetiapine`; IR plus PR is `propranolol` or `quetiapine`. Required fields: `id`, `name`
(Czech, as the app shows it), `category`, `t_half_h`, `vd_l_per_kg`, `bioavailability`,
`ka_per_h`, `mec`, `mtc`. Also give `t_half_range_h`, `concentration_unit`,
`dose_conversion`, `default_dose_mg`, `default_interval_h`, `dose_basis` and `notes`.
Most presets write `notes` in Czech, since the app shows them; a few older ones are in
English.

| block | fields | rule |
|---|---|---|
| `elimination` | `pathways` {enzyme: fraction}, source, verified, note | Sum ≤ 1. Enzyme names must exist in `data/enzymes.json` or be a non-enzyme route (`renal`, `UGT`, `sulfation`, `hydrolysis`, `CES1`, `CYP2E1`). `renal` must equal `fe_renal`. |
| `metabolite` | id, name, formation_fraction, t_half_h, vd_l_per_kg, molar_ratio, potency_ratio, reported, source, verified, note | `molar_ratio` = MW_parent / MW_metabolite. `reported`: `parent`, `metabolite` or `sum`, whichever the therapeutic window refers to. potency_ratio is used only when the reported quantity is the sum. |
| `protein_binding` | fu, source, verified, note | Recorded only; no prediction uses it. |
| `formulations` | {id: {name, name_cs, overrides}} | The first entry is the default and has no overrides. A PR form's id is `pr` (the SÚKL dosage-form matcher relies on it). It may override `bioavailability`, `ka_per_h`, `tlag_h`, `default_dose_mg` and `default_interval_h`. |
| `two_compartment` | v1_fraction, q_l_per_h | Only when curves show a distribution phase. |
| `cv_cl`, `cv_f` | floats | Log-normal CVs from a published spread. The note in `notes` or the monograph says how they were derived. |

`t_half_h` is the **parent's** half-life when a metabolite is simulated. The reported
quantity's apparent half-life comes out of the model.

## `data/reference_pk.json`

Add at least one exposure row (`auc_inf`, `auc_0_24`, `cavg_ss` or `cmin_ss`); a test
fails without it. Add `cmax` and `tmax_h` where a primary gives them. Rows are for the
**default formulation**, oral, 70 kg adult. Validation can't select a formulation, so
PR values go in the monograph and a dedicated test instead.

```json
{"drug_id": "<id>", "metric": "auc_inf", "dose_mg": 100, "value": 1.23,
 "source": "https://...", "verified": true,
 "note": "<Doc> <table/section>: <what, n, fed/fasted, statistic>; converted ng*h/mL -> ug*h/mL."}
```

Steady-state rows need `interval_h`. `max_fold` defaults to 1.5; don't raise it to pass.
When the reported quantity is a sum, the reference must be the sum too.

## `data/pd_effects.json`

Add the preset under `drugs`. Each key in `effects` is one of `serotonergic`,
`cns_depressant`, `qt_prolonging` or `sympathomimetic`, with source, verified and note.
A class needs a label **warning, stated interaction or described mechanism**; an
adverse-event table row isn't enough. An empty `effects` means "assessed, none". Leaving
the drug out of the file means "not assessed", which is worse.

## `data/interactions.json`

Add a row only for a **measured** pair, under `interactions`, with fields `victim_id`,
`perpetrator` (the exact name of an entry in `data/enzymes.json` `perpetrators` or one
of its `aliases`, so the time course can be simulated), `perpetrator_cs`,
`exposure_factor` (the AUC ratio), `site` (`clearance` or `bioavailability`),
`onset_h`, `qualitative`, `source`, `verified` and `note`. When the measured perpetrator
dose differs from the registry name's dose, name the row for the dose actually given
(`Phenytoin 300 mg/day`) rather than reusing a registry name with a different dose; the
time course then isn't simulated, but the label is true. Add every new perpetrator name
to its substance's `names` list in `data/pd_effects.json`, so the clash check doesn't
report it as unassessed. Don't add a new perpetrator to `enzymes.json` as part of a
preset unless asked; list the pair in the monograph instead.

How the app applies a factor: the co-medication path (`core/regimen.py`) scales the
victim's **reported** concentration curve by `exposure_factor`, so a measured ratio,
including an active-moiety ratio, is reproduced exactly and the half-life isn't changed.
`site` decides only two things: a `bioavailability` factor is skipped for IV dosing and
capped so F stays ≤ 1. With `bioavailability: 1.0` (unknown absolute F), use `clearance`,
or the cap removes the effect. The `core/interactions.py` docstring and
`apply_interactions` (only tests call it) describe an older clearance-scaling behaviour.
The one exception: when another drug in the regimen inactivates the enzyme (MDMA and
CYP2D6), a timed registry perpetrator on that enzyme acts through the shared enzyme
pool, at its mechanistic size.

## `core/sukl_client.py`

- the name map near the top: the Czech and English INN in upper case -> preset id;
- `PRESET_DEFAULT_SUKL_CODES[<id>]` = the 7-digit code, with a comment giving the
  product name;
- `PRESET_ATC_CODES[<id>]` = the ATC that the code must resolve to (a live test checks
  it).

## Tests

In `tests/test_drug_database.py`:
- add the id to `test_all_presets_load_and_are_valid`;
- add a drug-specific test that pins the structural facts the monograph argues for
  (reported quantity, the metabolite's share of the window, the PR form flattening the
  profile, the pathway response to an inhibitor). Pin numbers with a tolerance and say
  in the docstring where each target comes from.

If a sourced reference is missed after an honest fit, add it to `KNOWN_MISSES` in
`tests/test_reference_pk.py` with the reason. The xfail is strict, so a later fix
shows up.

## README

Add a row to the preset table in `README.md` and `README.cs.md` (same columns as the
existing rows).

## `data/reference_curves/` (optional, high value)

When a primary prints a mean concentration-time **table**, transcribe it twice,
independently, into `<id>_<source>_<dose>mg.json`. The format is in that directory's
README. Don't digitise figures. Set `verified: true` only if both transcriptions agree
with each other and with the source.
