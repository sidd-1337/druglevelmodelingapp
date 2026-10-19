"""A co-medication regimen simulated as a whole.

Levels of every drug, each documented or derived interaction applied over time,
the activity of the enzymes involved, and pharmacodynamic clashes - one function,
used by the Streamlit lab and the API alike so both give the same answer.

How an interaction is applied:

* An enzyme that a drug of the regimen destroys (MDMA and CYP2D6) is simulated as a
  shared pool (:func:`core.pk_engine.coupled_profiles`).
* An added substance with known PK applies its exposure factor over time
  (:mod:`core.perpetrator_timecourse`); without PK, or when the perpetrator is a drug
  of the regimen, the factor is constant.
* An added substance acting on an enzyme that is already such a pool is fed into the
  pool as an activity multiplier instead of as a factor, so the two cannot add up to
  more than losing the whole enzyme would.
* A bioavailability interaction does nothing to an IV dose, and never raises
  bioavailability above 1.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from core.drug_database import Drug
from core.enzymes import load_enzymes, perpetrator_pk
from core.interactions import (
    BIOAVAILABILITY,
    CLEARANCE,
    Interaction,
    derived_interactions,
    find,
    interactions_for,
    load_interactions,
)
from core.pd_effects import Clash, EffectTag, find_clashes, load_pd_effects, load_substance_effects
from core.perpetrator_timecourse import (
    ActingFactor,
    PerpetratorCourse,
    activity_from_fraction,
    apply_factor_course,
    effect_fraction,
    exposure_fraction,
)
from core.pk_engine import CoDrug, DoseEvent, PKModel, Route, auc_trapezoid, coupled_profiles, concentration_profile

# An added substance counts as present for the clash check once at half its steady exposure.
ACTIVE_FRACTION = 0.5
SUBSTANCE = "substance:"  # key prefix of an added substance in clash results


@dataclass(frozen=True)
class RegimenDrug:
    """A simulated drug: this patient's model before interactions, and its doses."""

    drug: Drug
    model: PKModel
    events: tuple[DoseEvent, ...]
    interval_h: float
    enzyme_shares: dict[str, float]

    @property
    def key(self) -> str:
        return self.drug.id


@dataclass(frozen=True)
class AddedSubstance:
    """A perpetrator that is not simulated itself (a name from data/interactions.json)."""

    name: str
    start_h: float | None = None  # None: already being taken at steady state
    stop_h: float | None = None


@dataclass(frozen=True)
class AppliedInteraction:
    interaction: Interaction
    factor: float | None  # steady-state exposure factor applied; None if not applied as a factor
    timed: bool  # follows the perpetrator's PK
    via_enzyme: bool = False  # applied through the shared enzyme pool instead
    skipped: str = ""  # "qualitative" or "iv" when it changes nothing


@dataclass
class DrugResult:
    entry: RegimenDrug
    concentrations: np.ndarray
    alone: np.ndarray
    applied: list[AppliedInteraction]
    drivers: list[str]  # names whose time course shapes this drug's interaction effect
    exposure_factor: float  # over the last dosing interval, with co-medication / alone


@dataclass
class RegimenResult:
    t: np.ndarray
    drugs: list[DrugResult]
    enzyme_activity: dict[str, np.ndarray]
    enzyme_drivers: dict[str, list[str]]  # what moves each enzyme ("destroys" or a substance name)
    destroyers: dict[str, str]  # enzyme -> regimen drug that destroys it
    clashes: list[Clash]
    tags: dict[str, tuple[EffectTag, ...]]  # clash key -> PD tags of every assessed participant
    labels: dict[str, tuple[str, str]]  # clash key -> (English, Czech) name
    unassessed: list[str]  # clash keys never assessed for PD effects
    unmatched: list[str]  # added substances that act on no drug of the regimen


def _acting(drug: Drug, others: set[str], added: list[str]) -> list[Interaction]:
    acting = [i for i in interactions_for(drug.id) if i.perpetrator in others]
    acting += [i for i in find(drug.id, added, drug=drug) if i not in acting]
    acting += [i for i in derived_interactions(drug, others) if i.perpetrator not in {a.perpetrator for a in acting}]
    return list({(i.perpetrator, i.site): i for i in acting}.values())


def _cs_name(name: str) -> str:
    hit = next((i.perpetrator_cs for i in load_interactions() if i.perpetrator == name and i.perpetrator_cs), "")
    if not hit and (p := perpetrator_pk(name)) is not None:
        hit = p.name_cs if name == p.name else ""
    return hit or name


def simulate_regimen(t: np.ndarray, entries: list[RegimenDrug], added: Iterable[AddedSubstance] = ()
                     ) -> RegimenResult:
    t = np.asarray(t, dtype=float)
    added = list(added)
    added_names = [a.name for a in added]
    enzymes = load_enzymes()
    names = {e.drug.name for e in entries} | {e.drug.id for e in entries} | set(added_names)

    co_drugs = [CoDrug(e.model, e.events, e.enzyme_shares, e.drug.inactivated_enzyme) for e in entries]
    destroyers = {e.drug.inactivated_enzyme: e.drug.name for e in entries if e.drug.inactivated_enzyme}
    used_enzymes = {enz for c in co_drugs for enz, s in c.shares.items() if s > 0}

    # Each added substance's course and its development g(t), computed once.
    courses: dict[str, tuple[PerpetratorCourse, np.ndarray]] = {}
    for a in added:
        p = perpetrator_pk(a.name)
        course = p.course(a.start_h, a.stop_h) if p is not None else None
        if course is not None:
            courses[a.name] = (course, effect_fraction(course, t, enzymes[course.enzyme].turnover_per_h))

    # Substances on an enzyme a regimen drug destroys go into that pool.
    multipliers: dict[str, np.ndarray] = {}
    enzyme_drivers: dict[str, list[str]] = {e: [d] for e, d in destroyers.items()}
    perp_activity: dict[str, np.ndarray] = {}
    for name, (course, g) in courses.items():
        if course.enzyme not in used_enzymes:
            continue
        activity = activity_from_fraction(course, g)
        target = multipliers if course.enzyme in destroyers else perp_activity
        target[course.enzyme] = target.get(course.enzyme, np.ones_like(t)) * activity
        enzyme_drivers.setdefault(course.enzyme, []).append(name)

    coupled = coupled_profiles(t, co_drugs, multipliers)

    results: list[DrugResult] = []
    acted: set[str] = set()
    for e, co, c_base in zip(entries, co_drugs, coupled.concentrations):
        acting = _acting(e.drug, names - {e.drug.name, e.drug.id}, added_names)
        acted |= {i.perpetrator for i in acting}
        applied, factors, drivers = [], [], []
        for i in acting:
            timed = i.perpetrator in courses
            if i.exposure_factor is None:
                applied.append(AppliedInteraction(i, None, timed, skipped="qualitative"))
                continue
            course_g = courses.get(i.perpetrator)
            if (course_g and course_g[0].enzyme in multipliers and i.site == CLEARANCE
                    and co.shares.get(course_g[0].enzyme, 0) > 0):
                applied.append(AppliedInteraction(i, None, timed, via_enzyme=True))
                drivers.append(i.perpetrator)
                continue
            factor = i.exposure_factor
            if i.site == BIOAVAILABILITY:
                if e.model.route is Route.IV_BOLUS:
                    applied.append(AppliedInteraction(i, None, timed, skipped="iv"))
                    continue
                factor = min(factor, 1.0 / e.model.f_eff)
            factors.append(ActingFactor(factor, course_g[1] if course_g else np.ones_like(t)))
            applied.append(AppliedInteraction(i, factor, timed))
            if timed:
                drivers.append(i.perpetrator)
        drivers = [destroyers[enz] for enz in coupled.enzyme_activity
                   if co.shares.get(enz, 0) > 0 and destroyers[enz] != e.drug.name] + drivers
        c_with = apply_factor_course(c_base, factors)
        alone = concentration_profile(t, e.events, e.model)
        lo = max(float(t[0]), float(t[-1]) - e.interval_h)
        base_auc = auc_trapezoid(t, alone, lo, float(t[-1]))
        exposure = auc_trapezoid(t, c_with, lo, float(t[-1])) / base_auc if base_auc > 0 else 1.0
        results.append(DrugResult(e, c_with, alone, applied, list(dict.fromkeys(drivers)), exposure))

    # Clashes over regimen drugs and added substances alike.
    drug_tags, sub_tags = load_pd_effects(), load_substance_effects()
    levels: dict[str, tuple[np.ndarray, float]] = {r.entry.key: (r.concentrations, r.entry.drug.mec) for r in results}
    tags = {k: drug_tags[k] for k in levels if k in drug_tags}
    labels = {r.entry.key: (r.entry.drug.name, r.entry.drug.name) for r in results}
    for a in added:
        key = SUBSTANCE + a.name
        presence = exposure_fraction(courses[a.name][0], t) if a.name in courses else np.ones_like(t)
        levels[key] = (presence, ACTIVE_FRACTION)
        labels[key] = (a.name, _cs_name(a.name))
        if a.name in sub_tags:
            tags[key] = sub_tags[a.name]
    clashes = find_clashes(t, levels, tags)

    enzyme_activity = {**perp_activity, **coupled.enzyme_activity}
    return RegimenResult(
        t=t, drugs=results, enzyme_activity=enzyme_activity,
        enzyme_drivers={e: enzyme_drivers.get(e, []) for e in enzyme_activity}, destroyers=destroyers,
        clashes=clashes, tags=tags, labels=labels, unassessed=[k for k in levels if k not in tags],
        unmatched=[n for n in added_names if n not in acted])
