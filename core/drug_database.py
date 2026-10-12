"""Drug presets: loading, validation and conversion to PK models."""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path

from core.pk_engine import PKModel, Route
from core.population import REFERENCE_WEIGHT_KG, clearance_factor, exposure_band_models

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
    cv_cl: float | None = None  # inter-individual CV of clearance (log-normal)
    cv_f: float = 0.0  # inter-individual CV of bioavailability
    fe_renal: float = 0.0  # fraction of clearance that scales with CrCl
    km: float | None = None  # Michaelis-Menten constant, concentration unit
    saturable_fraction: float = 0.0  # share of low-concentration CL that saturates

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
        if self.km is not None and (self.km <= 0 or not 0 < self.saturable_fraction <= 1):
            raise ValueError(f"{self.id}: km must be > 0 and saturable_fraction in (0, 1]")

    @property
    def half_life_band_h(self) -> tuple[float, float]:
        """Min/max half-life; defaults to +-30 % when no range is defined."""
        return self.t_half_range_h or (0.7 * self.t_half_h, 1.3 * self.t_half_h)

    def amount_from_dose_mg(self, dose_mg: float) -> float:
        """Dose in mg -> amount in the concentration's amount unit."""
        return dose_mg * self.dose_conversion

    def to_pk_model(
        self,
        body_weight_kg: float,
        route: Route | str = Route.ORAL,
        crcl_ml_min: float | None = None,
    ) -> PKModel:
        """PK model for a patient; t1/2 reflects allometric weight and renal function."""
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

        # V scales with weight, CL with clearance_factor: t1/2 = ln2 * V / CL.
        cl_rel = clearance_factor(body_weight_kg, crcl_ml_min, self.fe_renal)
        t_half = self.t_half_h * (body_weight_kg / REFERENCE_WEIGHT_KG) / cl_rel
        return PKModel(
            t_half_h=t_half,
            vd_l=self.vd_l_per_kg * body_weight_kg,
            bioavailability=bio,
            ka_per_h=ka,
            route=route,
            km=self.km,
            saturable_fraction=self.saturable_fraction if self.km is not None else 0.0,
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


def _parse(entry: dict) -> Drug:
    rng = entry.get("t_half_range_h")
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
        cv_cl=float(entry["cv_cl"]) if "cv_cl" in entry else None,
        cv_f=float(entry.get("cv_f", 0.0)),
        fe_renal=float(entry.get("fe_renal", 0.0)),
        km=float(sat["km"]) if (sat := entry.get("saturable")) else None,
        saturable_fraction=float(sat["fraction"]) if sat else 0.0,
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
