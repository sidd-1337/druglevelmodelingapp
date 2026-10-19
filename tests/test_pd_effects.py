"""Pharmacodynamic effect tags and the clash detector."""
import json
from pathlib import Path

import numpy as np
import pytest

from core.pd_effects import EFFECTS, EffectTag, find_clashes, load_pd_effects

DATA = Path(__file__).resolve().parent.parent / "data" / "pd_effects.json"


def _tag(effect):
    return EffectTag(effect, "https://dailymed.nlm.nih.gov/x", True, "quote")


def test_clash_needs_two_drugs_above_mec_at_once():
    t = np.linspace(0.0, 10.0, 1001)
    a = np.where(t < 6, 2.0, 0.0)  # above its MEC (1.0) for the first 6 h
    b = np.where(t > 4, 2.0, 0.0)  # above its MEC from 4 h on
    tags = {"a": (_tag("serotonergic"),), "b": (_tag("serotonergic"),)}
    (clash,) = find_clashes(t, {"a": (a, 1.0), "b": (b, 1.0)}, tags)
    assert clash.effect == "serotonergic"
    assert clash.drug_ids == ("a", "b")
    assert clash.overlap_h == pytest.approx(2.0, abs=0.05)  # overlap window 4-6 h
    assert clash.first_h == pytest.approx(4.0, abs=0.05)


def test_shared_effect_but_never_concurrent_is_reported_with_zero_overlap():
    t = np.linspace(0.0, 10.0, 1001)
    a = np.where(t < 4, 2.0, 0.0)
    b = np.where(t > 6, 2.0, 0.0)  # disjoint from a
    tags = {"a": (_tag("qt_prolonging"),), "b": (_tag("qt_prolonging"),)}
    (clash,) = find_clashes(t, {"a": (a, 1.0), "b": (b, 1.0)}, tags)
    assert clash.overlap_h == 0.0 and clash.first_h is None


def test_no_clash_when_effect_not_shared():
    t = np.linspace(0.0, 10.0, 101)
    c = np.full_like(t, 2.0)
    tags = {"a": (_tag("serotonergic"),), "b": (_tag("cns_depressant"),)}
    assert find_clashes(t, {"a": (c, 1.0), "b": (c, 1.0)}, tags) == []


def test_three_drugs_one_class_reported_once():
    t = np.linspace(0.0, 10.0, 101)
    c = np.full_like(t, 2.0)
    tags = {k: (_tag("cns_depressant"),) for k in ("a", "b", "c")}
    (clash,) = find_clashes(t, {k: (c, 1.0) for k in ("a", "b", "c")}, tags)
    assert clash.drug_ids == ("a", "b", "c")


def test_unknown_effect_rejected():
    with pytest.raises(ValueError):
        EffectTag("hepatotoxic", "x", True, "y")


def test_data_file_parses_and_covers_only_known_effects_and_drugs():
    tags = load_pd_effects()
    ids = {d["id"] for d in json.loads(
        (DATA.parent / "drugs.json").read_text(encoding="utf-8"))["drugs"]}
    for drug_id, drug_tags in tags.items():
        assert drug_id in ids, drug_id
        assert all(tag.effect in EFFECTS for tag in drug_tags)


def test_substance_tags_parse_and_name_known_perpetrators():
    from core.enzymes import load_perpetrators
    from core.interactions import load_interactions
    from core.pd_effects import EFFECTS, load_substance_effects
    subs = load_substance_effects()
    known = ({p.name for p in load_perpetrators()} | {i.perpetrator for i in load_interactions()}
             | {n for p in load_perpetrators() for n in p.aliases})
    for name, tags in subs.items():
        assert name in known, name
        assert all(tag.effect in EFFECTS for tag in tags)


def test_every_timed_perpetrator_is_pd_assessed():
    # A perpetrator we can place on the timeline should also be checked for PD clashes.
    from core.enzymes import load_perpetrators
    from core.pd_effects import load_substance_effects
    subs = load_substance_effects()
    for p in load_perpetrators():
        if p.t_half_h:
            assert p.name in subs, p.name
