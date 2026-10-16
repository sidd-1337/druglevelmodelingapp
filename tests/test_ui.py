import re
import tomllib
from pathlib import Path

import numpy as np
import pytest

from core.drug_database import get_drug
from core.metrics import TherapeuticWindow
from core.pk_engine import Route, concentration_profile, regular_schedule, time_grid
from ui.charts import build_figure
from ui.i18n import _TEXT, make_translator
from ui.reporting import build_dataframe

ROOT = Path(__file__).resolve().parent.parent


def test_translation_tables_have_same_keys():
    assert set(_TEXT["en"]) == set(_TEXT["cs"])


def test_all_keys_used_in_code_exist():
    used = set()
    for rel in ("app.py", "ui/charts.py", "ui/reporting.py"):
        src = (ROOT / rel).read_text(encoding="utf-8")
        used |= set(re.findall(r"\bT\(\s*\"([a-z_0-9]+)\"", src))
        used |= set(re.findall(r"T\(\s*'([a-z_0-9]+)'", src))
    used |= {
        "route_oral", "route_iv_bolus", "mode_skip", "mode_delay",
        "sukl_status_rx", "sukl_status_otc", "sukl_status_blue_stripe",
        "sukl_reimb_full", "sukl_reimb_partial", "sukl_reimb_none",
        "sukl_shortage_active", "sukl_shortage_interrupted", "sukl_shortage_ended", "sukl_shortage_unknown",
    }
    used = {k for k in used if not k.endswith("_")}  # dynamic prefixes covered above
    assert not (used - set(_TEXT["en"]))


def test_czech_toxic_alert_matches_spec_wording():
    T = make_translator("cs")
    msg = T("alert_toxic", t=6.5, mtc=50, unit="ug/mL", hours=1.0, cmax=60)
    assert "Pozor" in msg and "T=6.5 h" in msg


@pytest.mark.parametrize("log", [False, True])
def test_figure_builds_with_zones_band_and_reference(log):
    T = make_translator("en")
    drug = get_drug("ibuprofen")
    model = drug.to_pk_model(70)
    events = regular_schedule(400.0, 8, 72)
    t = time_grid(72)
    c = concentration_profile(t, events, model)
    w = TherapeuticWindow(drug.mec, drug.mtc)
    fig = build_figure(t, c, w, "ug/mL", T, log_scale=log, band=(0.8 * c, 1.2 * c), reference=c,
                       dose_times=[e.time_h for e in events], missed_time=8.0)
    names = {tr.name for tr in fig.data}
    assert {"Toxic (> MTC)", "Therapeutic window", "Sub-therapeutic (< MEC)", "Concentration"} <= names
    assert fig.layout.yaxis.type == ("log" if log else "linear")


def test_dataframe_columns_and_zones():
    w = TherapeuticWindow(2.0, 8.0)
    t = np.array([0.0, 1.0, 2.0])
    df = build_dataframe(t, np.array([1.0, 5.0, 9.0]), w, "mg/L")
    assert list(df["zone"]) == ["subtherapeutic", "therapeutic", "toxic"]
    assert "concentration_mg/L" in df.columns
    assert Route.ORAL.value == "oral"


def test_streamlit_config_defines_light_and_dark_themes():
    cfg = tomllib.loads((ROOT / ".streamlit" / "config.toml").read_text(encoding="utf-8"))
    light, dark = cfg["theme"]["light"], cfg["theme"]["dark"]
    for theme in (light, dark):
        assert {"primaryColor", "backgroundColor", "secondaryBackgroundColor", "textColor"} <= set(theme)
    assert light["backgroundColor"] != dark["backgroundColor"]


def test_dose_markers_stay_visible_on_dark_backgrounds():
    # The viewer can switch theme without a rerun, so the marker carries its own light halo.
    T = make_translator("en")
    drug = get_drug("ibuprofen")
    model = drug.to_pk_model(70)
    events = regular_schedule(400.0, 8, 24)
    t = time_grid(24)
    w = TherapeuticWindow(drug.mec, drug.mtc)
    fig = build_figure(t, concentration_profile(t, events, model), w, "ug/mL", T,
                       dose_times=[e.time_h for e in events])
    markers = next(tr for tr in fig.data if tr.name == T("doses")).marker
    assert markers.line.width >= 1 and markers.line.color.startswith("rgba(248")
