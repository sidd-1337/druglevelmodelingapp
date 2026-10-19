"use client";

import React, { useMemo, useRef } from "react";
import { num, sig } from "@/lib/format";
import Status, { Tone } from "@/components/ui/Status";
import { HALO, decimalsOf, hourStep, niceStep, scrubHandlers, useWidth } from "@/components/chart-utils";

export interface OverlaySeries {
  key: string;
  name: string;
  unit: string;
  mec: number;
  mtc: number;
  /** Concentration over the shared time grid, in the drug's own unit. */
  values: number[];
}

interface OverlayChartProps {
  timeGrid: number[];
  series: OverlaySeries[];
  /** Grid index shown in the readout; null = end of the simulation. Shared with other charts. */
  selected: number | null;
  onSelect: (index: number | null) => void;
}

const MARGIN = { top: 22, right: 40, bottom: 28, left: 2 };

/** Fixed colour per position in the regimen (never by rank); tokens in globals.css. */
export const seriesColor = (i: number) => `rgb(var(--series-${(i % 6) + 1}))`;

/** Where a level sits in its own window: 0 at MEC, 1 at MTC. */
const position = (c: number, s: OverlaySeries) => (c - s.mec) / (s.mtc - s.mec);

function zone(c: number, s: OverlaySeries): [Tone, string] {
  if (c > s.mtc) return ["bad", "nad MTC"];
  if (c >= s.mec) return ["ok", "v terapeutickém okně"];
  return ["neutral", "pod MEC"];
}

/**
 * Several drugs on one axis. Each is drawn as its position in its own therapeutic
 * window, the only scale drugs in different units share; the readout gives every
 * level in its own unit.
 */
export default function OverlayChart({ timeGrid, series, selected, onSelect }: OverlayChartProps) {
  const [wrapRef, width] = useWidth<HTMLDivElement>();
  const svgRef = useRef<SVGSVGElement>(null);

  const n = timeGrid.length;
  const durationH = n ? timeGrid[n - 1] : 1;
  const height = Math.round(Math.min(400, Math.max(260, width * 0.45)));
  const plotW = Math.max(0, width - MARGIN.left - MARGIN.right);
  const plotH = height - MARGIN.top - MARGIN.bottom;

  const positions = useMemo(() => series.map((s) => s.values.map((c) => position(c, s))), [series]);
  const [yMin, yMax] = useMemo(() => {
    let lo = 0;
    let hi = 1;
    positions.forEach((p) => p.forEach((v) => ((lo = Math.min(lo, v)), (hi = Math.max(hi, v)))));
    return [Math.min(-0.15, lo * 1.05), Math.max(1.3, hi * 1.08)];
  }, [positions]);

  const getX = (t: number) => MARGIN.left + (Math.max(0, Math.min(durationH, t)) / durationH) * plotW;
  const getY = (v: number) =>
    MARGIN.top + plotH - ((Math.max(yMin, Math.min(yMax, v)) - yMin) / (yMax - yMin)) * plotH;

  const paths = useMemo(
    () =>
      positions.map((p) =>
        p.map((v, i) => `${i === 0 ? "M" : "L"}${getX(timeGrid[i]).toFixed(1)} ${getY(v).toFixed(1)}`).join(" ")
      ),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [positions, timeGrid, width, height, yMin, yMax]
  );

  if (!n || series.length === 0) return null;

  const at = Math.min(selected ?? n - 1, n - 1);
  const handlers = scrubHandlers({
    svgRef, left: MARGIN.left, plotW, durationH, n, timeAt: (i) => timeGrid[i], at, select: onSelect,
  });

  const yStep = niceStep(yMax - yMin, 5);
  const yDecimals = decimalsOf(yStep);
  const yTicks: number[] = [];
  for (let v = Math.ceil(yMin / yStep) * yStep; v <= yMax + 1e-9; v += yStep) yTicks.push(v);
  const xStep = hourStep(durationH, Math.max(3, Math.floor(plotW / 72)));
  const xTicks: number[] = [];
  for (let t = 0; t <= durationH + 1e-9; t += xStep) xTicks.push(t);

  const yTop = MARGIN.top;
  const y0 = getY(0);
  const y1 = getY(1);
  const yAxis = MARGIN.top + plotH;
  const label =
    `Graf: ${series.length} léčiva v čase, každé jako poloha ve vlastním terapeutickém okně ` +
    `(0 = MEC, 1 = MTC). Šipkami vlevo a vpravo procházíte jednotlivé časy.`;

  return (
    <div>
      <div
        ref={wrapRef}
        tabIndex={0}
        role="group"
        aria-label={label}
        onKeyDown={handlers.onKeyDown}
        className="select-none rounded-ctl"
        style={{ height: height || 260 }}
      >
        {width > 0 && (
          <svg
            ref={svgRef}
            width={width}
            height={height}
            viewBox={`0 0 ${width} ${height}`}
            className="block overflow-visible"
            style={{ touchAction: "pan-y", cursor: "crosshair" }}
            onPointerMove={handlers.onPointerMove}
            onPointerDown={handlers.onPointerDown}
            onPointerLeave={handlers.onPointerLeave}
            aria-hidden
          >
            {/* Shared zones: every drug's window maps to 0..1 */}
            <rect x={MARGIN.left} y={yTop} width={plotW} height={Math.max(0, y1 - yTop)}
                  fill="rgb(var(--bad-fill) / 0.08)" />
            <rect x={MARGIN.left} y={y1} width={plotW} height={Math.max(0, y0 - y1)}
                  fill="rgb(var(--ok-fill) / 0.13)" />

            {yTicks.map((v) => (
              <g key={v}>
                <line x1={MARGIN.left} x2={MARGIN.left + plotW} y1={getY(v)} y2={getY(v)}
                      stroke="rgb(var(--line-soft))" />
                <text x={width - 2} y={getY(v) + 4} textAnchor="end" fontSize={12} fill="rgb(var(--ink-2))"
                      style={{ fontVariantNumeric: "tabular-nums" }}>
                  {num(v, yDecimals)}
                </text>
              </g>
            ))}

            <line x1={MARGIN.left} x2={MARGIN.left + plotW} y1={y1} y2={y1}
                  stroke="rgb(var(--bad-fill) / 0.7)" strokeDasharray="4 4" />
            <text x={MARGIN.left + 6} y={y1 - yTop < 16 ? y1 + 15 : y1 - 6} fontSize={12} fontWeight={500}
                  fill="rgb(var(--bad))" {...HALO}>
              MTC
            </text>
            <line x1={MARGIN.left} x2={MARGIN.left + plotW} y1={y0} y2={y0}
                  stroke="rgb(var(--ok-fill) / 0.8)" strokeDasharray="4 4" />
            <text x={MARGIN.left + 6} y={y0 - 6} fontSize={12} fontWeight={500} fill="rgb(var(--ok))" {...HALO}>
              MEC
            </text>

            <line x1={MARGIN.left} x2={MARGIN.left + plotW} y1={yAxis} y2={yAxis} stroke="rgb(var(--line))" />
            {xTicks.map((t) => (
              <text key={t} x={Math.min(getX(t), width - MARGIN.right - 2)} y={yAxis + 19}
                    textAnchor={t === 0 ? "start" : t >= durationH ? "end" : "middle"} fontSize={12}
                    fill="rgb(var(--ink-2))" style={{ fontVariantNumeric: "tabular-nums" }}>
                {t === 0 ? "0 h" : t}
              </text>
            ))}

            {/* Curves: a surface-coloured underlay keeps crossings readable */}
            {paths.map((d, i) => (
              <g key={series[i].key}>
                <path d={d} fill="none" stroke="rgb(var(--bg))" strokeWidth={5} strokeLinejoin="round" />
                <path d={d} fill="none" stroke={seriesColor(i)} strokeWidth={2.25} strokeLinejoin="round"
                      strokeLinecap="round" />
              </g>
            ))}

            <line x1={getX(timeGrid[at])} x2={getX(timeGrid[at])} y1={yTop} y2={yAxis}
                  stroke="rgb(var(--ink) / 0.35)" />
            {positions.map((p, i) => (
              <circle key={series[i].key} cx={getX(timeGrid[at])} cy={getY(p[at])} r={4.5} fill={seriesColor(i)}
                      stroke="rgb(var(--bg))" strokeWidth={2} />
            ))}
          </svg>
        )}
      </div>

      {/* Readout and legend: every drug's level in its own unit at the chosen time */}
      <p className="mt-6 text-ui text-ink-2">
        {selected === null ? "Na konci simulace" : "V čase"} {"t = "}
        <span className="figures">{num(timeGrid[at], 1)}</span>
        {" h"}
      </p>
      <ul className="mt-3 divide-y divide-line-soft border-y border-line-soft">
        {series.map((s, i) => {
          const c = s.values[at];
          const [tone, text] = zone(c, s);
          return (
            <li key={s.key} className="flex flex-wrap items-baseline justify-between gap-x-6 gap-y-1 py-3">
              <span className="flex items-center gap-3 text-ui">
                <span aria-hidden className="h-[3px] w-5 rounded-full" style={{ background: seriesColor(i) }} />
                {s.name}
              </span>
              <span className="flex items-baseline gap-5 text-ui">
                <span className="figures font-semibold">
                  {sig(c)} <span className="font-normal text-ink-2">{s.unit}</span>
                </span>
                <Status tone={tone}>{text}</Status>
              </span>
            </li>
          );
        })}
      </ul>
    </div>
  );
}
