"""Presets must reproduce published PK data (data/reference_pk.json) within tolerance."""
import pytest

import numpy as np

from core.drug_database import load_drugs
from core.validation import load_references, validate, validate_interactions

ROWS = validate()

# Published values the model misses once checked against the primary document. Each was
# hidden by a data error until the 2026 verification pass; none is tuned away. Strict, so a
# model change that fixes one fails here until the entry is removed.
KNOWN_MISSES = {
    ("lithium", "cmax", 300): "fasted Cmax is 0.40 mmol/L; the stored 0.28 was the fed value. The "
                              "two-compartment preset predicts 0.20",
    ("lithium", "tmax_h", 300): "fasted Tmax is 1.21 h; the preset peaks at 1.94 h (1.60-fold, limit 1.6)",
    ("mirtazapine", "cmin_ss", 30): "the TDM median 36.8 ng/mL is a trough, stored as cavg_ss until verified; "
                                    "the preset's trough is 19 ng/mL",
    ("paracetamol", "tmax_h", 1000): "measured curves show paracetamol from 5 minutes, so the 0.44 h lag was removed; "
                                     "ka 3.8/h was fitted together with that lag and now peaks at 0.21 h against 0.65",
    ("quetiapine", "auc_inf", 25): "healthy young men after one 25 mg dose clear quetiapine at 138 L/h; the preset "
                                   "is fitted to patients at steady state (65 L/h, same study), so it predicts 2.1x",
}


def _row_param(row):
    reason = KNOWN_MISSES.get((row.drug_id, row.metric, row.dose_mg))
    marks = [pytest.mark.xfail(strict=True, reason=reason)] if reason else []
    return pytest.param(row, marks=marks, id=f"{row.drug_id}-{row.metric}-{row.dose_mg:g}mg")


def test_every_known_miss_names_a_real_reference():
    keys = {(r.drug_id, r.metric, r.dose_mg) for r in ROWS}
    assert set(KNOWN_MISSES) <= keys


@pytest.mark.parametrize("row", [_row_param(r) for r in ROWS])
def test_prediction_within_published_range(row):
    assert row.within, (
        f"{row.drug_id} {row.metric}: predicted {row.predicted:.4g}, published {row.observed:.4g} "
        f"(ratio {row.ratio:.2f}, allowed {row.max_fold:g}x) - {row.source}"
    )


def test_every_preset_has_an_exposure_reference():
    exposure = {"auc_inf", "auc_0_24", "cavg_ss"}
    covered = {r["drug_id"] for r in load_references() if r["metric"] in exposure}
    # ketamine: no oral AUC found; amoxicillin: published AUCs conflict (17-29 ug*h/mL)
    missing = {d.id for d in load_drugs()} - covered - {"ketamine", "amoxicillin"}
    assert not missing


def test_references_cite_a_source():
    assert all(r["source"].startswith("http") for r in load_references())


def test_saturable_exposure_rises_more_than_dose():
    auc = {r.dose_mg: r.predicted for r in ROWS if r.drug_id == "mdma" and r.metric == "auc_0_24"}
    assert auc[125] / auc[75] > 125 / 75


# --------------------------------------------- mechanistic interaction back-test
DDI = validate_interactions()
# Pinned from the observed run: 5 pairs, geometric mean fold error 1.17. It was 1.25
# over 3 pairs before parent and metabolite were split; simulating the sum-reported
# victims both widened coverage and improved accuracy. Raise the ceiling only with a
# reason; never by tuning fm to fit.
MAX_DDI_GMFE = 1.35


@pytest.mark.parametrize("row", DDI, ids=lambda r: f"{r.victim_id}-{r.perpetrator.split()[0]}")
def test_mechanism_reproduces_each_measured_factor(row):
    assert row.within, (
        f"{row.victim_id} + {row.perpetrator}: measured {row.observed:.2f}, "
        f"mechanism predicts {row.predicted:.2f} from fm={row.fm:.2f} x residual={row.residual:.3f}")


def test_mechanistic_back_test_stays_within_its_pinned_error():
    assert DDI, "the mechanism can address no measured pair at all"
    gmfe = float(np.exp(np.mean([np.log(r.fold_error) for r in DDI])))
    assert gmfe <= MAX_DDI_GMFE, f"geometric mean fold error {gmfe:.2f} exceeds the pinned {MAX_DDI_GMFE}"


def test_back_test_only_covers_pairs_the_mechanism_can_address():
    # Renal handling, transporter and first-pass interactions are out of scope by
    # construction - they must not silently count as mechanistic successes.
    covered = {(r.victim_id, r.perpetrator) for r in DDI}
    assert ("lithium", "Ibuprofen") not in covered  # renal clearance, not an enzyme
    assert ("amoxicillin", "Probenecid") not in covered  # tubular secretion
    assert ("ketamine", "Rifampicin") not in covered  # first-pass bioavailability
