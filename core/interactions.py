"""Drug-drug interactions as exposure factors (data/interactions.json).

An interaction multiplies the victim drug's exposure (AUC, steady-state level) by
``exposure_factor``, applied either to clearance (half-life changes with it) or to
bioavailability (half-life unchanged). Effects of several perpetrators multiply.

The factors describe the full effect of the stated perpetrator regimen; the onset
(days to weeks for enzyme induction) is documented but not simulated.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

from core.pk_engine import PKModel

DEFAULT_INTERACTIONS_PATH = Path(__file__).resolve().parent.parent / "data" / "interactions.json"
CLEARANCE, BIOAVAILABILITY = "clearance", "bioavailability"


@dataclass(frozen=True)
class Interaction:
    victim_id: str
    perpetrator: str
    perpetrator_cs: str
    site: str
    source: str
    exposure_factor: float | None = None
    onset_h: float = 0.0
    qualitative: bool = False
    verified: bool = False
    note: str = ""

    def __post_init__(self) -> None:
        if self.site not in (CLEARANCE, BIOAVAILABILITY):
            raise ValueError(f"{self.perpetrator}: site must be {CLEARANCE!r} or {BIOAVAILABILITY!r}")
        if self.qualitative:
            if self.exposure_factor is not None:
                raise ValueError(f"{self.perpetrator}: qualitative interactions carry no exposure_factor")
        elif not self.exposure_factor or self.exposure_factor <= 0:
            raise ValueError(f"{self.perpetrator}: exposure_factor must be > 0")

    def name(self, lang: str = "en") -> str:
        return self.perpetrator_cs if lang == "cs" and self.perpetrator_cs else self.perpetrator

    @property
    def direction(self) -> str:
        """'increase', 'decrease' or 'unknown' (qualitative)."""
        if self.exposure_factor is None:
            return "unknown"
        return "increase" if self.exposure_factor > 1.0 else "decrease"


@lru_cache(maxsize=4)
def _load(path: str) -> tuple[Interaction, ...]:
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    return tuple(Interaction(**entry) for entry in raw["interactions"])


def load_interactions(path: str | Path = DEFAULT_INTERACTIONS_PATH) -> list[Interaction]:
    return list(_load(str(path)))


def interactions_for(drug_id: str, path: str | Path = DEFAULT_INTERACTIONS_PATH) -> list[Interaction]:
    """Known interactions where ``drug_id`` is the victim."""
    return [i for i in _load(str(path)) if i.victim_id == drug_id]


def find(drug_id: str, perpetrators: Iterable[str],
         path: str | Path = DEFAULT_INTERACTIONS_PATH) -> list[Interaction]:
    """Interactions on ``drug_id`` from the named perpetrators (exact ``perpetrator`` strings)."""
    wanted = set(perpetrators)
    return [i for i in interactions_for(drug_id, path) if i.perpetrator in wanted]


def combined_factor(interactions: Sequence[Interaction], site: str) -> float:
    """Product of the quantified exposure factors acting on ``site``."""
    factor = 1.0
    for i in interactions:
        if i.site == site and i.exposure_factor is not None:
            factor *= i.exposure_factor
    return factor


def apply_interactions(model: PKModel, interactions: Sequence[Interaction]) -> PKModel:
    """Model with the quantified interaction effects applied.

    Clearance effects scale the half-life (volume is unchanged), bioavailability
    effects scale F and are capped at 1. Qualitative entries change nothing.
    """
    cl_factor = combined_factor(interactions, CLEARANCE)
    f_factor = combined_factor(interactions, BIOAVAILABILITY)
    if cl_factor == 1.0 and f_factor == 1.0:
        return model
    return replace(
        model,
        t_half_h=model.t_half_h * cl_factor,
        bioavailability=min(1.0, model.bioavailability * f_factor),
    )
