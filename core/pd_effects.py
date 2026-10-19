"""Pharmacodynamic effect classes, and clashes between drugs taken together.

Some combinations are dangerous with no change in blood levels at all: two
serotonergic drugs (serotonin syndrome), two CNS depressants (sedation, respiratory
depression), two QT-prolonging drugs (torsades de pointes), two sympathomimetics
(blood pressure, heart rate). ``data/pd_effects.json`` tags each preset with the
classes its primary document states - the label, or for a drug without one the
original studies. A drug missing from the file has not been assessed, which is not
the same as having none of these effects.

A clash is flagged where at least two drugs of one class are above their MEC at the
same time, the simulated level standing in for "pharmacologically active". That is
a proxy: the risk of QT prolongation and of sedation rises with concentration, and
some effects outlast the level.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Mapping

import numpy as np

DEFAULT_PD_PATH = Path(__file__).resolve().parent.parent / "data" / "pd_effects.json"
EFFECTS = ("serotonergic", "cns_depressant", "qt_prolonging", "sympathomimetic")


@dataclass(frozen=True)
class EffectTag:
    effect: str
    source: str
    verified: bool
    note: str

    def __post_init__(self) -> None:
        if self.effect not in EFFECTS:
            raise ValueError(f"unknown effect class {self.effect!r}")


@dataclass(frozen=True)
class Clash:
    effect: str
    drug_ids: tuple[str, ...]  # drugs of the regimen tagged with the effect, in regimen order
    overlap_h: float  # time at least two of them are above their MEC together
    first_h: float | None  # when that first happens


@lru_cache(maxsize=4)
def _load(path: str) -> dict[str, tuple[EffectTag, ...]]:
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    return {drug_id: tuple(EffectTag(effect, **tag) for effect, tag in entry["effects"].items())
            for drug_id, entry in raw["drugs"].items()}


def load_pd_effects(path: str | Path = DEFAULT_PD_PATH) -> dict[str, tuple[EffectTag, ...]]:
    """Drug id -> its effect tags, for every assessed drug (an empty tuple: assessed, none stated)."""
    return dict(_load(str(path)))


def find_clashes(
    t: np.ndarray,
    levels: Mapping[str, tuple[np.ndarray, float]],
    tags: Mapping[str, tuple[EffectTag, ...]] | None = None,
) -> list[Clash]:
    """Effect classes shared by at least two drugs of a regimen.

    ``levels`` maps a drug id to (concentration over ``t``, its MEC). A class shared
    by drugs that are never active together is still returned, with ``overlap_h`` 0:
    a different schedule could bring them together.
    """
    tags = load_pd_effects() if tags is None else tags
    t = np.asarray(t, dtype=float)
    out = []
    for effect in EFFECTS:
        ids = tuple(i for i in levels if any(tag.effect == effect for tag in tags.get(i, ())))
        if len(ids) < 2:
            continue
        active = np.sum([np.asarray(levels[i][0]) >= levels[i][1] for i in ids], axis=0) >= 2
        overlap = float(np.sum(np.diff(t)[active[:-1]]))  # left-point rule, as in metrics._time_where
        first = float(t[np.argmax(active)]) if active.any() else None
        out.append(Clash(effect, ids, overlap, first))
    return out
