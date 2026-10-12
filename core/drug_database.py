"""Drug presets: loading, validation and conversion to PK models."""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Mapping

from core.pk_engine import Metabolite, PKModel, Route
from core.population import (
    REFERENCE_WEIGHT_KG,
    RENAL,
    clearance_factor,
    exposure_band_models,
    pathway_clearance_factor,
)

DEFAULT_DB_PATH = Path(__file__).resolve().parent.parent / "data" / "drugs.json"
CUSTOM_ID = "custom"


@dataclass(frozen=True)
class Drug:
    id: str
    name: str
    t_half_h: float
    vd_l_per_kg: float
    bioavailability: float
    ka_per_h: float
    mec: float
    mtc: float
    category: str = ""
    t_half_range_h: tuple[float, float] | None = None
    concentration_unit: str = "ug/mL"
    dose_conversion: float = 1.0
    default_dose_mg: float = 100.0
    default_interval_h: float = 12.0
    dose_basis: str = "mg"
    notes: str = ""
    routes: dict[str, dict[str, float]] | None = None
    formulations: dict[str, dict] | None = None
    cv_cl: float | None = None  # inter-individual CV of clearance (log-normal)
    cv_f: float = 0.0  # inter-individual CV of bioavailability
    fe_renal: float = 0.0  # fraction of clearance that scales with CrCl
    # Fraction of clearance per elimination pathway (enzyme name or "renal"); the
    # unassigned remainder is size-scaled and genotype-free. Kept consistent with
    # fe_renal, which is the "renal" entry.
    pathways: dict[str, float] | None = None
    # Unbound fraction. Recorded for reference only: every therapeutic range here is a
    # total plasma concentration, so no prediction depends on it.
    fu: float = 1.0
    # One active metabolite, as raw fields; vd is per kg so it scales with the patient.
    metabolite: dict | None = None
    km: float | None = None  # Michaelis-Menten constant, concentration unit
    saturable_fraction: float = 0.0  # share of low-concentration CL that saturates
    v1_fraction: float = 1.0  # central / steady-state volume (< 1: two compartments)
    q_l_per_h: float = 0.0  # inter-compartmental clearance of a 70 kg adult
    tlag_h: float = 0.0  # oral absorption lag
    # Mechanism-based enzyme inactivation; t_half_h is then the inactivated-enzyme half-life.
    inactivation_fraction: float = 0.0
    kinact_per_h: float = 0.0
    ki: float = 1.0
    kdeg_per_h: float = 0.0
    inactivated_enzyme: str = ""

    def __post_init__(self) -> None:
        if self.t_half_h <= 0 or self.vd_l_per_kg <= 0 or self.ka_per_h <= 0:
            raise ValueError(f"{self.id}: t_half_h, vd_l_per_kg and ka_per_h must be > 0")
        if not 0 < self.bioavailability <= 1:
            raise ValueError(f"{self.id}: bioavailability must be in (0, 1]")
        if not 0 < self.mec < self.mtc:
            raise ValueError(f"{self.id}: require 0 < MEC < MTC")
        if self.dose_conversion <= 0:
            raise ValueError(f"{self.id}: dose_conversion must be > 0")
        if self.t_half_range_h is not None:
            lo, hi = self.t_half_range_h
            if not 0 < lo <= self.t_half_h <= hi:
                raise ValueError(f"{self.id}: t_half_range_h must bracket t_half_h")
        if (self.cv_cl is not None and self.cv_cl < 0) or self.cv_f < 0:
            raise ValueError(f"{self.id}: cv_cl and cv_f must be >= 0")
        if not 0 <= self.fe_renal <= 1:
            raise ValueError(f"{self.id}: fe_renal must be in [0, 1]")
        if self.pathways is not None:
            if any(not 0 <= v <= 1 for v in self.pathways.values()):
                raise ValueError(f"{self.id}: pathway fractions must be in [0, 1]")
            if sum(self.pathways.values()) > 1 + 1e-9:
                raise ValueError(f"{self.id}: pathway fractions must sum to <= 1")
            if abs(self.pathways.get(RENAL, 0.0) - self.fe_renal) > 1e-9:
                raise ValueError(f"{self.id}: pathways[{RENAL!r}] must equal fe_renal")
        if self.km is not None and (self.km <= 0 or not 0 < self.saturable_fraction <= 1):
            raise ValueError(f"{self.id}: km must be > 0 and saturable_fraction in (0, 1]")

    @property
    def half_life_band_h(self) -> tuple[float, float]:
        """Min/max half-life; defaults to +-30 % when no range is defined."""
        return self.t_half_range_h or (0.7 * self.t_half_h, 1.3 * self.t_half_h)

    def amount_from_dose_mg(self, dose_mg: float) -> float:
        """Dose in mg -> amount in the concentration's amount unit."""
        return dose_mg * self.dose_conversion

    @property
    def formulation_ids(self) -> list[str]:
        return list(self.formulations) if self.formulations else []

    def formulation_name(self, formulation_id: str, lang: str = "en") -> str:
        cfg = (self.formulations or {}).get(formulation_id, {})
        return (cfg.get("name_cs") if lang == "cs" else None) or cfg.get("name", formulation_id)

    def with_formulation(self, formulation_id: str | None) -> "Drug":
        """Preset with a product form's overrides applied (default form or None: unchanged)."""
        if not formulation_id or not self.formulations:
            return self
        if formulation_id not in self.formulations:
            raise KeyError(f"{self.id}: unknown formulation {formulation_id!r}")
        cfg = {k: v for k, v in self.formulations[formulation_id].items() if k not in ("name", "name_cs")}
        return replace(self, **cfg) if cfg else self

    def to_pk_model(
        self,
        body_weight_kg: float,
        route: Route | str = Route.ORAL,
        crcl_ml_min: float | None = None,
        genotypes: Mapping[str, str] | None = None,
        enzyme_activity: Mapping[str, float] | None = None,
    ) -> PKModel:
        """PK model for a patient; t1/2 reflects weight, renal function and genotype.

        ``genotypes`` maps an enzyme to a phenotype label (resolved through
        ``core.enzymes``); ``enzyme_activity`` sets activities directly and wins
        where both are given. Inhibitors do NOT belong here - they are applied as
        an :class:`~core.interactions.Interaction`, so their effect is never
        counted twice.
        """
        if body_weight_kg <= 0:
            raise ValueError("body_weight_kg must be > 0")
        if isinstance(route, str):
            r_str = route.strip().lower()
            if r_str in ("inhaled", "smoked"):
                route = Route.SMOKED
            elif r_str in ("iv", "iv_bolus", "intravenous"):
                route = Route.IV_BOLUS
            elif r_str in ("insufflated", "nasal", "snorted"):
                route = Route.INSUFFLATED
            elif r_str in ("oral", "po", "peroral"):
                route = Route.ORAL
            else:
                route = Route(r_str)

        bio = self.bioavailability
        ka = self.ka_per_h
        if self.routes and route.value in self.routes:
            r_cfg = self.routes[route.value]
            bio = float(r_cfg.get("bioavailability", bio))
            ka = float(r_cfg.get("ka_per_h", ka))
        elif self.routes and route == Route.SMOKED and "inhaled" in self.routes:
            r_cfg = self.routes["inhaled"]
            bio = float(r_cfg.get("bioavailability", bio))
            ka = float(r_cfg.get("ka_per_h", ka))
        elif route == Route.INSUFFLATED:
            ka = max(self.ka_per_h * 2.5, 3.5)
            bio = min(1.0, max(self.bioavailability, 0.70))
        elif route == Route.SMOKED:
            ka = max(self.ka_per_h * 5.0, 8.0)
            bio = min(1.0, max(self.bioavailability * 0.9, 0.60))
        elif route == Route.IV_BOLUS:
            bio = 1.0
            ka = 1.0

        # V scales with weight, CL (and Q) with clearance_factor: t1/2 = ln2 * V / CL.
        if genotypes:
            from core.enzymes import resolve_enzyme_activity  # local: avoids a cycle via data loading

            activity = dict(resolve_enzyme_activity(genotypes))
        else:
            activity = {}
        activity.update(enzyme_activity or {})
        cl_rel = pathway_clearance_factor(
            body_weight_kg, self.pathways or ({RENAL: self.fe_renal} if self.fe_renal else None),
            crcl_ml_min, activity)
        t_half = self.t_half_h * (body_weight_kg / REFERENCE_WEIGHT_KG) / cl_rel
        # The engine's t1/2 is the fresh-enzyme value: CL_fresh = CL_inactivated / (1 - fm).
        t_half *= 1.0 - self.inactivation_fraction
        return PKModel(
            t_half_h=t_half,
            vd_l=self.vd_l_per_kg * body_weight_kg,
            bioavailability=bio,
            ka_per_h=ka,
            route=route,
            km=self.km,
            saturable_fraction=self.saturable_fraction if self.km is not None else 0.0,
            v1_fraction=self.v1_fraction,
            # Q is organ blood flow, not metabolism: pure allometric scaling.
            q_l_per_h=self.q_l_per_h * clearance_factor(body_weight_kg),
            tlag_h=self.tlag_h if route is Route.ORAL else 0.0,
            metabolite=(Metabolite(
                name=self.metabolite["name"],
                t_half_h=float(self.metabolite["t_half_h"]),
                vd_l=float(self.metabolite["vd_l_per_kg"]) * body_weight_kg,
                formation_fraction=float(self.metabolite.get("formation_fraction", 1.0)),
                molar_ratio=float(self.metabolite.get("molar_ratio", 1.0)),
                potency_ratio=float(self.metabolite.get("potency_ratio", 1.0)),
                reported=self.metabolite.get("reported", "parent"),
            ) if self.metabolite else None),
            inactivation_fraction=self.inactivation_fraction,
            kinact_per_h=self.kinact_per_h,
            ki=self.ki,
            kdeg_per_h=self.kdeg_per_h,
        )

    def exposure_band(self, model: PKModel) -> tuple[PKModel, PKModel]:
        """(low, high) 90 % inter-individual exposure models around ``model``.

        Uses the CL/F variability when known, else the t1/2 range relative to
        the preset t1/2 (so covariate effects on ``model`` are kept).
        """
        if self.cv_cl is not None:
            return exposure_band_models(model, self.cv_cl, self.cv_f)
        lo, hi = self.half_life_band_h
        scale = model.t_half_h / self.t_half_h
        return model.with_half_life(lo * scale), model.with_half_life(hi * scale)

    def with_overrides(self, **changes) -> "Drug":
        return replace(self, **changes)


# SÚKL/EU dosage-form markers of a modified-release product (CPS PRO, TBL RET, ...).
_PROLONGED_RELEASE = ("PRO", "RET", "PROLONG", "RETARD", "MODIFIK", "MODIFIED", "EXTENDED")


def formulation_from_dosage_form(drug: Drug, dosage_form: str) -> str | None:
    """Formulation id matching a dosage-form description, or None if nothing fits."""
    if not drug.formulations or not dosage_form:
        return None
    text = dosage_form.upper()
    prolonged = any(marker in text for marker in _PROLONGED_RELEASE)
    for fid in drug.formulation_ids:
        if (fid == "pr") == prolonged:
            return fid
    return None


def custom_drug() -> Drug:
    """Neutral template for user-defined parameters."""
    return Drug(
        id=CUSTOM_ID,
        name="Custom drug",
        category="User defined",
        t_half_h=6.0,
        t_half_range_h=None,
        vd_l_per_kg=0.5,
        bioavailability=0.8,
        ka_per_h=1.0,
        mec=5.0,
        mtc=25.0,
        default_dose_mg=100,
        default_interval_h=12,
        cv_cl=0.3,
    )


def _parse_pathways(entry: dict) -> tuple[dict[str, float] | None, float]:
    """(pathways, fe_renal) from an ``elimination`` block, or from a bare ``fe_renal``.

    Either spelling works: the renal fraction is always both ``fe_renal`` and the
    ``renal`` pathway, so readers of either stay correct.
    """
    elimination = entry.get("elimination")
    if elimination:
        pathways = {str(k): float(v) for k, v in elimination["pathways"].items()}
        return pathways, pathways.get(RENAL, 0.0)
    fe_renal = float(entry.get("fe_renal", 0.0))
    return ({RENAL: fe_renal} if fe_renal else None), fe_renal


def _parse(entry: dict) -> Drug:
    rng = entry.get("t_half_range_h")
    pathways, fe_renal = _parse_pathways(entry)
    return Drug(
        id=entry["id"],
        name=entry["name"],
        category=entry.get("category", ""),
        t_half_h=float(entry["t_half_h"]),
        t_half_range_h=tuple(map(float, rng)) if rng else None,
        vd_l_per_kg=float(entry["vd_l_per_kg"]),
        bioavailability=float(entry["bioavailability"]),
        ka_per_h=float(entry["ka_per_h"]),
        mec=float(entry["mec"]),
        mtc=float(entry["mtc"]),
        concentration_unit=entry.get("concentration_unit", "ug/mL"),
        dose_conversion=float(entry.get("dose_conversion", 1.0)),
        default_dose_mg=float(entry.get("default_dose_mg", 100)),
        default_interval_h=float(entry.get("default_interval_h", 12)),
        dose_basis=entry.get("dose_basis", "mg"),
        notes=entry.get("notes", ""),
        routes=entry.get("routes"),
        formulations=entry.get("formulations"),
        cv_cl=float(entry["cv_cl"]) if "cv_cl" in entry else None,
        cv_f=float(entry.get("cv_f", 0.0)),
        fe_renal=fe_renal,
        pathways=pathways,
        fu=float(binding["fu"]) if (binding := entry.get("protein_binding")) else 1.0,
        metabolite=entry.get("metabolite"),
        km=float(sat["km"]) if (sat := entry.get("saturable")) else None,
        saturable_fraction=float(sat["fraction"]) if sat else 0.0,
        v1_fraction=float(two["v1_fraction"]) if (two := entry.get("two_compartment")) else 1.0,
        q_l_per_h=float(two["q_l_per_h"]) if two else 0.0,
        tlag_h=float(entry.get("tlag_h", 0.0)),
        inactivation_fraction=float(inact["fraction"]) if (inact := entry.get("enzyme_inactivation")) else 0.0,
        kinact_per_h=float(inact["kinact_per_h"]) if inact else 0.0,
        ki=float(inact["ki"]) if inact else 1.0,
        kdeg_per_h=float(inact["kdeg_per_h"]) if inact else 0.0,
        inactivated_enzyme=inact.get("enzyme", "") if inact else "",
    )


@lru_cache(maxsize=4)
def _load(path: str) -> tuple[Drug, ...]:
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    drugs = tuple(_parse(e) for e in raw["drugs"])
    ids = [d.id for d in drugs]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate drug ids in database")
    return drugs


def load_drugs(path: str | Path = DEFAULT_DB_PATH) -> list[Drug]:
    return list(_load(str(path)))


def get_drug(drug_id: str, path: str | Path = DEFAULT_DB_PATH) -> Drug:
    if drug_id == CUSTOM_ID:
        return custom_drug()
    for d in load_drugs(path):
        if d.id == drug_id:
            return d
    raise KeyError(f"unknown drug id: {drug_id!r}")
