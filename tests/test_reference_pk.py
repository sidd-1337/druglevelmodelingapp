"""Presets must reproduce published PK data (data/reference_pk.json) within tolerance."""
import pytest

from core.drug_database import load_drugs
from core.validation import load_references, validate

ROWS = validate()


@pytest.mark.parametrize("row", ROWS, ids=lambda r: f"{r.drug_id}-{r.metric}-{r.dose_mg:g}mg")
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
