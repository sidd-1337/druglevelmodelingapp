# Monograph template: `docs/drugs/<id>.md`

Write in English. Every number gets a citation in the form `[n]`, pointing to the
numbered source list at the end. Mark values not read in a primary document with
`(unverified: <reason>)`. Use tables where a list of values would otherwise be prose.
Educational use only: the first line after the title says so.

```markdown
# <INN> (<Czech name>)

> Educational summary for the drug-level modelling app. Not clinical guidance.

## 1. Identity
INN, Czech name, ATC, class, salt and base molecular weights, products on the Czech
market (forms and strengths), the preset id and default SÚKL code.

## 2. Mechanism of action
- One paragraph: what the drug is for and the core mechanism in one sentence.
- Table: target | Ki parent (nM) | Ki metabolite (nM) | action (antagonist, partial
  agonist, inhibitor...) | engaged at therapeutic levels? | source.
- Occupancy: PET data (target, dose, % occupancy, time after dose) where measured.
- Dose dependence: which targets come in as the dose rises.
- Active metabolites: what they add, and their share of the effect.

## 3. Actions in the body
Table: system | effect | mechanism (receptor) | dose/exposure where it appears | source.
Cover CNS (therapeutic and sedative), cardiovascular (BP, HR, QT), metabolic and
endocrine (weight, glucose, lipids, prolactin), autonomic (anticholinergic),
haematological, and anything on the label's boxed warning or warnings section.

## 4. Pharmacokinetics
Table: parameter | value | definition (dose, fed/fasted, single/steady state,
population) | source | verified.
Rows: F, Tmax, Cmax/dose, AUC/dose, t½ (parent, metabolites), Vd, CL/F, protein
binding, linearity, food effect, then the same for each PR form.
Then: metabolism (enzymes, fractions, evidence type: label mass balance / inhibitor
study / genotype study / in vitro), excretion (mass balance), special populations,
genotype.

## 5. Interactions
- PK: table of measured pairs: perpetrator (dose) | AUC ratio | Cmax ratio | source.
  Note which ones the app can simulate (perpetrator present in data/enzymes.json).
- PD: which of the app's four clash classes apply (serotonergic, CNS depressant,
  QT-prolonging, sympathomimetic), each with the label sentence that justifies it,
  plus important PD interactions outside those classes.

## 6. Pharmacodynamics and therapeutic drug monitoring
Therapeutic range and alert level (AGNP: total concentration, trough unless stated),
sampling time, exposure-response and occupancy-response data, dose-response, time
course of effect versus concentration (hysteresis, tolerance).

## 7. How the app models it
- Structure chosen (compartments, metabolite, reported quantity, formulations) and the
  one decision that could have gone the other way.
- Table: preset field | value | how it was obtained (sourced / derived: formula /
  fitted to: reference) | source.
- Fit: reference | observed | predicted | ratio (output of preset_report.py).
- What the model can't represent for this drug (food, first-pass formation, nonlinear
  absorption, CYP3A5, ...), and the known misses.

## 8. Open items
Values left unverified and why; conflicts between primaries; data that would close them.

## Sources
1. <Document title>, <issuer/journal, year>, <section or table used>. <URL>
```
