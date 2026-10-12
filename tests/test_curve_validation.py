"""The curve harness, tested against curves the engine itself generates.

No published curve has been seeded (see data/reference_curves/README.md), so these
tests prove the scorer works rather than that the model fits. The one that matters is
``test_a_wrong_absorption_rate_is_caught``: it builds a curve from a deliberately
wrong model and checks the harness rejects it. Without that, a passing score on a
real curve added later would mean nothing.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from core.drug_database import get_drug
from core.pk_engine import DoseEvent, Route, concentration_profile
from core.validation import (
    CURVE_KINDS,
    REFERENCE_CURVE_DIR,
    CurveFit,
    convert_concentration,
    load_reference_curves,
    predict_curve,
    score_curve,
    validate_curve,
    validate_curves,
)

# Ibuprofen's preset carries a 0.5 h absorption lag, so sampling starts after it;
# test_a_sample_before_the_absorption_lag_is_reported covers what happens if it does not.
TIMES = [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0]


def _curve_from_model(drug_id: str, dose_mg: float, times=TIMES, **overrides) -> dict:
    """A curve file whose points come from the engine, optionally with a wrong preset."""
    drug = get_drug(drug_id)
    if overrides:
        drug = drug.with_overrides(**overrides)
    model = drug.to_pk_model(70.0)
    c = concentration_profile(np.array(times, dtype=float),
                             (DoseEvent(0.0, drug.amount_from_dose_mg(dose_mg)),), model)
    return {
        "curve_id": f"synthetic_{drug_id}",
        "drug_id": drug_id,
        "kind": "observed",
        "doses": [{"time_h": 0.0, "dose_mg": dose_mg}],
        "unit": drug.concentration_unit,
        "source": "generated from the engine inside the test suite",
        "points": [[t, float(v)] for t, v in zip(times, c)],
    }


# ------------------------------------------------------------------ the scorer
def test_a_curve_generated_from_the_preset_scores_perfectly():
    row = validate_curve(_curve_from_model("ibuprofen", 400.0))
    assert row.fit.aafe == pytest.approx(1.0, abs=1e-9)
    assert row.fit.afe == pytest.approx(1.0, abs=1e-9)
    assert row.fit.fraction_within_2fold == 1.0
    assert row.within and row.n_points == len(TIMES)


@pytest.mark.parametrize("overrides,label", [
    ({"ka_per_h": 0.4}, "absorption 10x too slow"),
    ({"t_half_h": 32.0, "t_half_range_h": None}, "half-life 4x too long"),
    ({"vd_l_per_kg": 0.04}, "volume 3x too small"),
])
def test_a_wrong_model_is_caught(overrides, label):
    """The harness must reject a curve the preset cannot reproduce."""
    row = validate_curve(_curve_from_model("ibuprofen", 400.0, **overrides))
    assert not row.within, (label, row.fit.aafe)
    assert row.fit.aafe > 1.5, (label, row.fit.aafe)


def test_a_wrong_absorption_rate_is_caught_where_cmax_and_auc_would_not_be():
    """The point of curve validation: shape error the summary metrics miss.

    Halving ka and halving the volume keeps AUC identical and Cmax close, so the
    existing reference checks stay happy while the curve is visibly wrong.
    """
    from core.pk_engine import cmax_single_dose

    drug = get_drug("ibuprofen")
    wrong = drug.with_overrides(ka_per_h=drug.ka_per_h / 8.0)
    amount = drug.amount_from_dose_mg(400.0)
    right_m, wrong_m = drug.to_pk_model(70.0), wrong.to_pk_model(70.0)
    # AUC is ka-independent, so that check cannot see this at all.
    assert (right_m.f_eff * amount / right_m.clearance_l_per_h
            == pytest.approx(wrong_m.f_eff * amount / wrong_m.clearance_l_per_h))
    assert cmax_single_dose(amount, wrong_m) < cmax_single_dose(amount, right_m)
    # ... but the curve shows it.
    row = validate_curve(_curve_from_model("ibuprofen", 400.0, ka_per_h=drug.ka_per_h / 8.0))
    assert row.fit.aafe > 1.5, row.fit.aafe


def test_afe_shows_direction_and_aafe_does_not():
    fit = score_curve([1.0, 2.0], [10.0, 10.0], [20.0, 5.0])  # one 2x high, one 2x low
    assert fit.afe == pytest.approx(1.0)      # the bias cancels
    assert fit.aafe == pytest.approx(2.0)     # the error does not
    assert fit.fraction_within_2fold == 1.0   # exactly 2x counts as within
    high = score_curve([1.0], [10.0], [30.0])
    assert high.afe == pytest.approx(3.0) and high.aafe == pytest.approx(3.0)
    assert high.fraction_within_2fold == 0.0


@pytest.mark.parametrize("bad", [
    ([1.0], [10.0], [1.0, 2.0]),   # mismatched shapes
    ([], [], []),                  # nothing to score
    ([1.0], [0.0], [1.0]),         # zero observation
    ([1.0], [10.0], [-1.0]),       # negative prediction
])
def test_unscoreable_input_is_rejected(bad):
    with pytest.raises(ValueError):
        score_curve(*bad)


# ------------------------------------------------------------------ units
def test_the_same_curve_in_different_units_scores_the_same():
    ug = _curve_from_model("ibuprofen", 400.0)
    ng = {**ug, "unit": "ng/mL",
          "points": [[t, c * 1000.0] for t, c in ug["points"]]}
    mg_per_l = {**ug, "unit": "mg/L"}  # identical to ug/mL
    base = validate_curve(ug).fit.aafe
    assert validate_curve(ng).fit.aafe == pytest.approx(base)
    assert validate_curve(mg_per_l).fit.aafe == pytest.approx(base)


def test_molar_units_work_for_lithium():
    mmol = _curve_from_model("lithium", 300.0)
    assert mmol["unit"] == "mmol/L"
    umol = {**mmol, "unit": "umol/L",
            "points": [[t, c * 1000.0] for t, c in mmol["points"]]}
    meq = {**mmol, "unit": "mEq/L"}  # lithium is monovalent
    assert validate_curve(umol).fit.aafe == pytest.approx(1.0, abs=1e-9)
    assert validate_curve(meq).fit.aafe == pytest.approx(1.0, abs=1e-9)


def test_lithium_printed_as_mass_of_li_converts_only_with_its_molar_mass():
    """FDA reviews print lithium as ug/mL of Li+; mmol/L is that divided by 6.941."""
    mmol = _curve_from_model("lithium", 300.0)
    ug_li = {**mmol, "unit": "ug/mL", "molar_mass_g_per_mol": 6.941,
             "points": [[t, c * 6.941] for t, c in mmol["points"]]}
    assert validate_curve(ug_li).fit.aafe == pytest.approx(1.0, abs=1e-9)
    ng_li = {**ug_li, "unit": "ng/mL", "points": [[t, c * 1000.0] for t, c in ug_li["points"]]}
    assert validate_curve(ng_li).fit.aafe == pytest.approx(1.0, abs=1e-9)
    # Read as lithium carbonate the same numbers would be tenfold low - so no default.
    with pytest.raises(ValueError, match="molar_mass_g_per_mol"):
        validate_curve({k: v for k, v in ug_li.items() if k != "molar_mass_g_per_mol"})


def test_an_unknown_unit_pair_raises_instead_of_guessing():
    with pytest.raises(ValueError, match="no conversion"):
        convert_concentration([1.0], "nmol/L", "ug/mL")
    with pytest.raises(ValueError, match="no conversion"):
        validate_curve({**_curve_from_model("ibuprofen", 400.0), "unit": "parts per gill"})


# ------------------------------------------------------------------ file handling
def test_points_below_the_lloq_are_dropped():
    curve = _curve_from_model("ibuprofen", 400.0)
    tail = min(c for _, c in curve["points"])
    row = validate_curve({**curve, "lloq": tail * 1.5})
    assert row.n_points < len(TIMES)
    with pytest.raises(ValueError, match="LLOQ"):
        validate_curve({**curve, "lloq": 1e6})


def test_a_sample_before_the_absorption_lag_is_reported():
    """A zero prediction is an infinite fold error, so it must be named, not averaged."""
    drug = get_drug("ibuprofen")
    assert drug.tlag_h > 0, "this test needs a preset with an absorption lag"
    curve = _curve_from_model("ibuprofen", 400.0, times=[0.25, *TIMES])
    curve["points"][0][1] = 0.4  # a real assay would have measured something
    with pytest.raises(ValueError, match="predicts zero at t = 0.25"):
        validate_curve(curve)
    # Declaring it below the assay's limit resolves it.
    assert validate_curve({**curve, "lloq": 0.5}).n_points == len(TIMES)


def test_points_are_sorted_by_time():
    curve = _curve_from_model("ibuprofen", 400.0)
    shuffled = {**curve, "points": list(reversed(curve["points"]))}
    row = validate_curve(shuffled)
    assert np.all(np.diff(row.fit.times_h) > 0)
    assert row.fit.aafe == pytest.approx(1.0, abs=1e-9)


def test_a_multiple_dose_schedule_accumulates():
    drug = get_drug("levetiracetam")
    times = [12.5, 13.0, 14.0, 18.0, 24.0]
    single = predict_curve("levetiracetam", times, [(0.0, 500.0)])
    twice = predict_curve("levetiracetam", times, [(0.0, 500.0), (12.0, 500.0)])
    assert np.all(twice > single)
    assert drug.concentration_unit == "ug/mL"


def test_route_is_honoured():
    oral = predict_curve("propranolol", TIMES, [(0.0, 40.0)], route=Route.ORAL)
    iv = predict_curve("propranolol", TIMES, [(0.0, 40.0)], route=Route.IV_BOLUS)
    assert iv[0] > oral[0]  # no absorption delay and no first-pass loss


def test_an_unknown_kind_is_rejected():
    with pytest.raises(ValueError, match="kind"):
        validate_curve({**_curve_from_model("ibuprofen", 400.0), "kind": "vibes"})
    assert CURVE_KINDS == ("observed", "published_model")


def test_published_model_curves_are_scored_but_kept_separate():
    curve = {**_curve_from_model("ibuprofen", 400.0), "kind": "published_model"}
    row = validate_curve(curve)
    assert row.kind == "published_model" and row.within


def test_a_missing_directory_yields_no_curves(tmp_path):
    assert load_reference_curves(tmp_path / "nope") == []
    assert validate_curves(tmp_path / "nope") == []


def test_curve_files_are_discovered_and_named_after_the_file(tmp_path):
    curve = _curve_from_model("ibuprofen", 400.0)
    curve.pop("curve_id")
    (tmp_path / "ibuprofen_400mg.json").write_text(json.dumps(curve), encoding="utf-8")
    rows = validate_curves(tmp_path)
    assert [r.curve_id for r in rows] == ["ibuprofen_400mg"]
    assert rows[0].within


# ------------------------------------------------------------------ the real directory
def test_the_reference_curve_directory_documents_itself():
    assert (REFERENCE_CURVE_DIR / "README.md").is_file()


def test_every_seeded_curve_loads_and_is_scoreable():
    """A malformed seed file must fail here, not be silently skipped."""
    rows = validate_curves()
    assert len(rows) == len(load_reference_curves())
    for row in rows:
        assert row.n_points > 0 and row.source
        assert math.isfinite(row.fit.aafe) and row.fit.aafe >= 1.0
        assert row.kind in CURVE_KINDS


def test_curve_fit_requires_matching_shapes_on_construction():
    with pytest.raises(ValueError, match="same shape"):
        CurveFit(np.array([1.0]), np.array([1.0, 2.0]), np.array([1.0, 2.0]))
