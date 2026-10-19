"""Enzyme registry: phenotype activities and perpetrator inhibition/induction.

Both act on the same quantity - the activity of one elimination pathway, relative
to a normal metaboliser with no co-medication - so both feed the pathway clearance
model in :mod:`core.population`.

An inhibitor's residual activity is derived from the AUC increase it causes in a
sensitive index substrate, corrected for how much of that substrate's clearance
the enzyme actually carries::

    AUC_ratio = 1 / (fm * residual + (1 - fm))
    residual  = (1 / AUC_ratio - (1 - fm)) / fm

Reading the residual straight off the AUC ratio (i.e. assuming fm = 1) would
overstate the remaining activity, so the index substrate's own fm is required.
The same relation run forwards predicts the effect on any other substrate of that
enzyme, which is what :func:`derived_exposure_factor` does.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping

DEFAULT_ENZYMES_PATH = Path(__file__).resolve().parent.parent / "data" / "enzymes.json"
NORMAL_PHENOTYPE = "normal"
# FDA classification bands, as residual activity of the enzyme.
CLASS_BANDS: dict[str, tuple[float, float]] = {
    "strong_inhibitor": (0.0, 0.2),
    "moderate_inhibitor": (0.2, 0.5),
    "weak_inhibitor": (0.5, 0.8),
    "inducer": (1.0, 100.0),
}


@dataclass(frozen=True)
class Enzyme:
    name: str
    genotype_activity: dict[str, float]
    phenotype_frequency: dict[str, float]
    source: str
    verified: bool = False
    note: str = ""
    turnover_h: float | None = None  # enzyme degradation half-life, for the onset of induction

    @property
    def turnover_per_h(self) -> float | None:
        return math.log(2.0) / self.turnover_h if self.turnover_h else None

    def __post_init__(self) -> None:
        if self.turnover_h is not None and self.turnover_h <= 0:
            raise ValueError(f"{self.name}: turnover_h must be > 0")
        if self.genotype_activity.get(NORMAL_PHENOTYPE) != 1.0:
            raise ValueError(f"{self.name}: the normal phenotype must have activity 1.0")
        if any(v < 0 for v in self.genotype_activity.values()):
            raise ValueError(f"{self.name}: activities must be >= 0")
        if self.phenotype_frequency:
            if set(self.phenotype_frequency) - set(self.genotype_activity):
                raise ValueError(f"{self.name}: a frequency names an unknown phenotype")
            total = sum(self.phenotype_frequency.values())
            if abs(total - 1.0) > 0.05:
                raise ValueError(f"{self.name}: phenotype frequencies sum to {total:.3f}, not ~1")


@dataclass(frozen=True)
class Perpetrator:
    """A drug that changes one enzyme's activity, with the evidence it rests on."""

    name: str
    name_cs: str
    enzyme: str
    drug_class: str
    index_substrate: str
    index_auc_ratio: float
    index_fm: float
    source: str
    verified: bool = False
    note: str = ""
    # Own pharmacokinetics, for the onset/offset of its effect (core.perpetrator_timecourse).
    # t_half_h None: no time course available, so the effect is treated as static.
    t_half_h: float | None = None
    interval_h: float = 24.0
    pk_source: str = ""
    # Slower fade after stopping (active metabolite, tissue reservoir, destroyed enzyme).
    offset_t_half_h: float | None = None
    offset_source: str = ""
    # Other perpetrator strings in data/interactions.json naming this same drug, so a
    # measured pair recorded under a different regimen label still gets this PK.
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.drug_class not in CLASS_BANDS:
            raise ValueError(f"{self.name}: unknown class {self.drug_class!r}")
        if self.index_auc_ratio <= 0:
            raise ValueError(f"{self.name}: index_auc_ratio must be > 0")
        if not 0 < self.index_fm <= 1:
            raise ValueError(f"{self.name}: index_fm must be in (0, 1]")
        if self.t_half_h is not None and (self.t_half_h <= 0 or self.interval_h <= 0):
            raise ValueError(f"{self.name}: t_half_h and interval_h must be > 0")
        if self.offset_t_half_h is not None and self.offset_t_half_h <= 0:
            raise ValueError(f"{self.name}: offset_t_half_h must be > 0")
        object.__setattr__(self, "aliases", tuple(self.aliases))
        lo, hi = CLASS_BANDS[self.drug_class]
        if not lo - 1e-9 <= self.residual_activity <= hi + 1e-9:
            raise ValueError(
                f"{self.name}: residual {self.residual_activity:.3f} derived from "
                f"{self.index_auc_ratio}x on {self.index_substrate} falls outside the "
                f"{self.drug_class} band {lo}-{hi}")

    @property
    def residual_activity(self) -> float:
        """Enzyme activity left under this perpetrator (1.0 = unaffected, > 1 = induced)."""
        return max(0.0, (1.0 / self.index_auc_ratio - (1.0 - self.index_fm)) / self.index_fm)

    @property
    def induces(self) -> bool:
        return self.drug_class == "inducer"

    def course(self, start_h: float | None = None, stop_h: float | None = None):
        """A :class:`~core.perpetrator_timecourse.PerpetratorCourse`, or None without known PK.

        ``start_h`` None (the default) means the patient was already taking it at
        steady state, which is what a static interaction factor assumes.
        """
        if self.t_half_h is None:
            return None
        from core.perpetrator_timecourse import PerpetratorCourse  # avoids an import cycle

        return PerpetratorCourse(self.enzyme, self.residual_activity, self.t_half_h, self.induces,
                                 interval_h=self.interval_h, start_h=start_h, stop_h=stop_h,
                                 offset_t_half_h=self.offset_t_half_h)

    def display_name(self, lang: str = "en") -> str:
        return self.name_cs if lang == "cs" and self.name_cs else self.name


@lru_cache(maxsize=4)
def _load(path: str) -> tuple[dict[str, Enzyme], tuple[Perpetrator, ...]]:
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    enzymes = {name: Enzyme(name=name, **cfg) for name, cfg in raw["enzymes"].items()}
    perpetrators = tuple(Perpetrator(**entry) for entry in raw["perpetrators"])
    unknown = {p.enzyme for p in perpetrators} - set(enzymes)
    if unknown:
        raise ValueError(f"perpetrators name unknown enzymes: {sorted(unknown)}")
    # Induction onset needs the enzyme's turnover; fail at load, not in the UI.
    missing = sorted({p.name for p in perpetrators
                      if p.induces and p.t_half_h is not None and enzymes[p.enzyme].turnover_h is None})
    if missing:
        raise ValueError(f"inducers on an enzyme without turnover_h: {missing}")
    names = [n for p in perpetrators for n in (p.name, *p.aliases)]
    if len(names) != len(set(names)):
        raise ValueError("a perpetrator name or alias is used twice")
    return enzymes, perpetrators


def load_enzymes(path: str | Path = DEFAULT_ENZYMES_PATH) -> dict[str, Enzyme]:
    return dict(_load(str(path))[0])


def load_perpetrators(path: str | Path = DEFAULT_ENZYMES_PATH) -> list[Perpetrator]:
    return list(_load(str(path))[1])


def polymorphic_enzymes(path: str | Path = DEFAULT_ENZYMES_PATH) -> set[str]:
    """Enzymes with a usable phenotype system (others get inhibition but no genotype)."""
    return {name for name, e in _load(str(path))[0].items() if len(e.genotype_activity) > 1}


def genotype_activity(enzyme: str, phenotype: str | None,
                      path: str | Path = DEFAULT_ENZYMES_PATH) -> float:
    """Activity multiplier for a phenotype; 1.0 for a normal or unstated one."""
    if not phenotype or phenotype == NORMAL_PHENOTYPE:
        return 1.0
    enzymes = _load(str(path))[0]
    if enzyme not in enzymes:
        raise KeyError(f"unknown enzyme {enzyme!r}")
    try:
        return enzymes[enzyme].genotype_activity[phenotype]
    except KeyError:
        raise KeyError(f"{enzyme}: unknown phenotype {phenotype!r}")


def resolve_enzyme_activity(genotypes: Mapping[str, str],
                            path: str | Path = DEFAULT_ENZYMES_PATH) -> dict[str, float]:
    """Phenotype labels per enzyme -> activity multipliers."""
    return {enzyme: genotype_activity(enzyme, phenotype, path)
            for enzyme, phenotype in genotypes.items()}


def perpetrator_pk(name: str, path: str | Path = DEFAULT_ENZYMES_PATH) -> Perpetrator | None:
    """The perpetrator whose name or alias is ``name`` (for its PK), or None."""
    for p in _load(str(path))[1]:
        if name == p.name or name in p.aliases:
            return p
    return None


def find_perpetrators(names: Iterable[str],
                      path: str | Path = DEFAULT_ENZYMES_PATH) -> list[Perpetrator]:
    wanted = set(names)
    return [p for p in _load(str(path))[1] if p.name in wanted]


def residual_activity(enzyme: str, perpetrators: Iterable[Perpetrator]) -> float:
    """Combined residual activity of one enzyme under several perpetrators."""
    residual = 1.0
    for p in perpetrators:
        if p.enzyme == enzyme:
            residual *= p.residual_activity
    return residual


def derived_exposure_factor(fm: float, residual: float) -> float:
    """Victim exposure factor from the static mechanistic model.

    ``fm`` is the victim's fraction of clearance through the affected enzyme and
    ``residual`` that enzyme's remaining activity. Induction (residual > 1) uses
    the same relation and returns a factor below 1.
    """
    if not 0 <= fm <= 1:
        raise ValueError("fm must be in [0, 1]")
    if residual < 0:
        raise ValueError("residual must be >= 0")
    denominator = fm * residual + (1.0 - fm)
    if denominator <= 0:
        raise ValueError("a fully inhibited, wholly enzyme-cleared drug has no finite exposure")
    return 1.0 / denominator
