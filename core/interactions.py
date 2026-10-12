"""Drug-drug interactions as exposure factors (data/interactions.json).

An interaction multiplies the victim drug's exposure (AUC, steady-state level) by
``exposure_factor``, applied either to clearance (half-life changes with it) or to
bioavailability (half-life unchanged). Effects of several perpetrators multiply.

The factors describe the full effect of the stated perpetrator regimen; the onset
(days to weeks for enzyme induction) is documented but not simulated.

Measured factors always win. Where no study exists for a pair, one can be derived
from the victim's fraction metabolised and the perpetrator's residual enzyme
activity (:mod:`core.enzymes`); such entries are marked ``estimated`` and carry
the inputs they came from in ``basis``. Derived entries only ever act on
clearance - a first-pass effect needs a gut model this does not have.
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
    estimated: bool = False  # derived from fm x residual activity, not measured
    basis: str = ""  # how a derived factor was obtained

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
         path: str | Path = DEFAULT_INTERACTIONS_PATH, drug=None) -> list[Interaction]:
    """Interactions on ``drug_id`` from the named perpetrators (exact ``perpetrator`` strings).

    With a ``drug`` (a :class:`~core.drug_database.Drug`), pairs that have no
    measured study also get a derived estimate, where the perpetrator's enzyme is
    one the drug is actually cleared by.
    """
    wanted = set(perpetrators)
    measured = [i for i in interactions_for(drug_id, path) if i.perpetrator in wanted]
    if drug is None:
        return measured
    covered = {i.perpetrator for i in measured}
    return measured + derived_interactions(drug, wanted - covered)


def mechanistic_exposure_factor(drug, enzyme: str, residual: float,
                                body_weight_kg: float = 70.0) -> tuple[float, str]:
    """(exposure factor, how it was obtained) for inhibiting ``enzyme`` to ``residual``.

    For an ordinary drug this is the static relation 1 / (fm * residual + 1 - fm).
    That relation does NOT hold when the reported quantity is a parent plus its
    metabolite: inhibiting the enzyme that forms the metabolite raises the parent
    but lowers the metabolite, so the sum is buffered and the static formula
    overstates the effect. In that case the factor is simulated instead, which the
    engine can do because parent and metabolite are separate compartments.
    """
    from core.enzymes import derived_exposure_factor  # avoids an import cycle
    from core.pk_engine import steady_state_summary

    fm = (drug.pathways or {}).get(enzyme, 0.0)
    static = derived_exposure_factor(fm, residual)
    if drug.metabolite is None or drug.to_pk_model(body_weight_kg).reported_quantity == "parent":
        return static, (f"{enzyme} carries {100 * fm:.0f} % of this drug's clearance and "
                        f"{residual:.2f} of that activity is left")
    amount = drug.amount_from_dose_mg(drug.default_dose_mg)
    tau = drug.default_interval_h
    base = steady_state_summary(amount, tau, drug.to_pk_model(body_weight_kg))
    hit = steady_state_summary(amount, tau,
                               drug.to_pk_model(body_weight_kg, enzyme_activity={enzyme: residual}))
    return (hit.cavg_ss / base.cavg_ss,
            f"simulated: inhibiting {enzyme} to {residual:.2f} of normal raises the parent but forms less "
            f"metabolite, so the reported sum moves far less than the parent alone")


def derived_interactions(drug, perpetrator_names: Iterable[str]) -> list[Interaction]:
    """Mechanistic estimates for perpetrators acting on a pathway the drug uses."""
    from core.enzymes import find_perpetrators  # avoids an import cycle

    out = []
    for p in find_perpetrators(perpetrator_names):
        if (drug.pathways or {}).get(p.enzyme, 0.0) <= 0:
            continue
        factor, how = mechanistic_exposure_factor(drug, p.enzyme, p.residual_activity)
        out.append(Interaction(
            victim_id=drug.id, perpetrator=p.name, perpetrator_cs=p.name_cs, site=CLEARANCE,
            source=p.source, exposure_factor=factor, verified=False, estimated=True,
            basis=(f"{how}; the perpetrator's effect is derived from a {p.index_auc_ratio:g}x AUC change "
                   f"in {p.index_substrate}"),
            note="Estimated from the mechanism, not measured for this pair."))
    return out


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
