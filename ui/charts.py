"""Plotly figure builder for concentration-time curves with therapeutic zones."""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import plotly.graph_objects as go

from core.metrics import TherapeuticWindow

RED = "220,38,38"
GREEN = "22,163,74"
GREY = "100,116,139"
BLUE = "37,99,235"


def _rgba(rgb: str, a: float) -> str:
    return f"rgba({rgb},{a})"


def build_figure(
    t: np.ndarray,
    c: np.ndarray,
    window: TherapeuticWindow,
    unit: str,
    T: Callable[..., str],
    *,
    log_scale: bool = False,
    band: tuple[np.ndarray, np.ndarray] | None = None,
    reference: np.ndarray | None = None,
    dose_times: Sequence[float] = (),
    missed_time: float | None = None,
    x_range: tuple[float, float] | None = None,
    show_zones: bool = True,
    height: int = 520,
    observations: Sequence[tuple[float, float]] = (),
    conc_label: str | None = None,
    extra: tuple[np.ndarray, str] | None = None,
    reference_label: str | None = None,
    band_label: str | None = None,
) -> go.Figure:
    t = np.asarray(t, dtype=float)
    x0, x1 = x_range if x_range else (float(t[0]), float(t[-1]))

    peaks = [float(np.max(c))]
    if reference is not None:
        peaks.append(float(np.max(reference)))
    if band is not None:
        peaks.append(float(np.max(band[1])))
    peaks += [obs_c for _, obs_c in observations]
    if extra is not None:
        peaks.append(float(np.max(extra[0])))
    y_top = max(window.mtc * 1.3, max(peaks) * 1.1)
    y_bot = window.mec / 100.0 if log_scale else 0.0

    def prep(arr: np.ndarray) -> np.ndarray:
        arr = np.asarray(arr, dtype=float)
        return np.where(arr > 0, arr, np.nan) if log_scale else arr

    fig = go.Figure()

    if show_zones:
        lo = y_bot if log_scale else 0.0
        fig.add_trace(go.Scatter(
            x=[x0, x1, x1, x0], y=[lo, lo, window.mec, window.mec], fill="toself", mode="lines",
            line=dict(width=0), fillcolor=_rgba(GREY, 0.22), name=T("zone_sub"), hoverinfo="skip"))
        fig.add_trace(go.Scatter(
            x=[x0, x1, x1, x0], y=[window.mec, window.mec, window.mtc, window.mtc], fill="toself",
            mode="lines", line=dict(width=0), fillcolor=_rgba(GREEN, 0.16), name=T("zone_ther"),
            hoverinfo="skip"))
        fig.add_trace(go.Scatter(
            x=[x0, x1, x1, x0], y=[window.mtc, window.mtc, y_top, y_top], fill="toself", mode="lines",
            line=dict(width=0), name=T("zone_tox"), hoverinfo="skip",
            fillpattern=dict(shape="/", fgcolor=_rgba(RED, 0.55), bgcolor=_rgba(RED, 0.07), size=8)))
        fig.add_hline(y=window.mec, line=dict(color=_rgba(GREEN, 0.9), width=1, dash="dash"),
                      annotation_text=f"MEC {window.mec:g}", annotation_position="bottom right")
        fig.add_hline(y=window.mtc, line=dict(color=_rgba(RED, 0.9), width=1, dash="dash"),
                      annotation_text=f"MTC {window.mtc:g}", annotation_position="top right")

    if band is not None:
        fig.add_trace(go.Scatter(x=t, y=prep(band[1]), mode="lines", line=dict(width=0),
                                 hoverinfo="skip", showlegend=False))
        fig.add_trace(go.Scatter(x=t, y=prep(band[0]), mode="lines", line=dict(width=0), fill="tonexty",
                                 fillcolor=_rgba(BLUE, 0.14), name=band_label or T("band"), hoverinfo="skip"))

    hover = "t = %{x:.2f} h<br>C = %{y:.4g} " + unit + "<extra>%{fullData.name}</extra>"
    if reference is not None:
        fig.add_trace(go.Scatter(x=t, y=prep(reference), mode="lines", name=reference_label or T("reference"),
                                 line=dict(color=_rgba(GREY, 0.9), width=2, dash="dash"),
                                 hovertemplate=hover))
    fig.add_trace(go.Scatter(x=t, y=prep(c), mode="lines", name=conc_label or T("conc"),
                             line=dict(color=f"rgb({BLUE})", width=3), hovertemplate=hover))
    if observations:
        fig.add_trace(go.Scatter(
            x=[o[0] for o in observations], y=[o[1] for o in observations], mode="markers",
            name=T("measured"), marker=dict(symbol="diamond", size=12, color=f"rgb({RED})",
                                            line=dict(width=1.5, color="#0f172a")),
            hovertemplate=hover))

    if extra is not None:
        fig.add_trace(go.Scatter(x=t, y=prep(extra[0]), mode="lines", name=extra[1],
                                 line=dict(color=_rgba(GREY, 0.95), width=2, dash="dot"),
                                 hovertemplate=hover))

    marker_y = y_bot * 1.6 if log_scale else y_top * 0.015
    shown = [d for d in dose_times if x0 <= d <= x1]
    if shown:
        fig.add_trace(go.Scatter(
            x=shown, y=[marker_y] * len(shown), mode="markers", name=T("doses"),
            # Dark fill with a light halo: readable on both themes (the viewer can switch
            # theme without a rerun, so the colour cannot depend on the current one).
            marker=dict(symbol="triangle-up", size=9, color="#0f172a",
                        line=dict(width=1, color="rgba(248,250,252,0.9)")),
            hovertemplate="t = %{x:.2f} h<extra>" + T("dose_marker") + "</extra>"))
    if missed_time is not None and x0 <= missed_time <= x1:
        fig.add_trace(go.Scatter(
            x=[missed_time], y=[marker_y], mode="markers", name=T("missed_marker"),
            marker=dict(symbol="x", size=13, color=f"rgb({RED})", line=dict(width=2)),
            hovertemplate="t = %{x:.2f} h<extra>" + T("missed_marker") + "</extra>"))

    fig.update_layout(
        template="plotly_white", height=height, hovermode="x unified", dragmode="zoom",
        margin=dict(l=10, r=10, t=30, b=10),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        xaxis=dict(title=T("axis_time"), range=[x0, x1], showspikes=True, spikemode="across"),
        yaxis=dict(title=f"{T('axis_conc')} [{unit}]",
                   type="log" if log_scale else "linear",
                   range=[np.log10(y_bot), np.log10(y_top)] if log_scale else [0, y_top]),
    )
    return fig
