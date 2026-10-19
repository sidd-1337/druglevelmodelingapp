"use client";

import React, { useMemo, useRef } from "react";
import { num } from "@/lib/format";
import { HALO, hourStep, scrubHandlers, useWidth } from "@/components/chart-utils";

interface EnzymeChartProps {
  timeGrid: number[];
  /** Active enzyme as a fraction of baseline (1 = normal) over the time grid. */
  activity: Record<string, number[]>;
  /** What changes each enzyme, already as display names. */
  drivers: Record<string, string[]>;
  selected: number | null;
  onSelect: (index: number | null) => void;
}

const MARGIN = { top: 16, right: 48, bottom: 28, left: 2 };
// Monochrome on purpose: the drug colours belong to the drugs in the overlay above.
const STROKES = [{ dash: undefined, ink: "--ink" }, { dash: "6 4", ink: "--ink-2" }, { dash: "2 3", ink: "--ink-2" }];

/**
 * Activity of the enzymes the regimen changes, in % of this patient's baseline.
 * Once an inducer takes an enzyme above baseline (rifampicin: about 26-fold) the axis
 * becomes logarithmic, so a few % of inhibition and a 2600 % induction both stay readable.
 */
export default function EnzymeChart({ timeGrid, activity, drivers, selected, onSelect }: EnzymeChartProps) {
  const [wrapRef, width] = useWidth<HTMLDivElement>();
  const svgRef = useRef<SVGSVGElement>(null);
  const enzymes = Object.keys(activity);
  const n = timeGrid.length;
  const durationH = n ? timeGrid[n - 1] : 1;
  const height = 200;
  const plotW = Math.max(0, width - MARGIN.left - MARGIN.right);
  const plotH = height - MARGIN.top - MARGIN.bottom;

  const pct = useMemo(
    () => Object.fromEntries(enzymes.map((e) => [e, activity[e].map((a) => 100 * a)])),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [activity]
  );
  const all = enzymes.flatMap((e) => pct[e]);
  const top = Math.max(100, ...all);
  const bottom = Math.min(100, ...all);
  const log = top > 105;
  const [lo, hi] = log ? [Math.max(bottom, 1) / 1.5, top * 1.3] : [0, 105];
  const scale = (v: number) => (log ? Math.log10(Math.max(v, lo)) : v);
  const getX = (t: number) => MARGIN.left + (t / durationH) * plotW;
  const getY = (v: number) =>
    MARGIN.top + plotH - ((scale(Math.min(Math.max(v, lo), hi)) - scale(lo)) / (scale(hi) - scale(lo))) * plotH;

  if (!n || enzymes.length === 0) return null;

  const at = Math.min(selected ?? n - 1, n - 1);
  const handlers = scrubHandlers({
    svgRef, left: MARGIN.left, plotW, durationH, n, timeAt: (i) => timeGrid[i], at, select: onSelect,
  });
  const yTicks = log
    ? [1, 3, 10, 30, 100, 300, 1000, 3000].filter((v) => v >= lo && v <= hi)
    : [0, 25, 50, 75, 100];
  const xStep = hourStep(durationH, Math.max(3, Math.floor(plotW / 72)));
  const xTicks: number[] = [];
  for (let t = 0; t <= durationH + 1e-9; t += xStep) xTicks.push(t);
  const yAxis = MARGIN.top + plotH;

  return (
    <div>
      <div
        ref={wrapRef}
        tabIndex={0}
        role="group"
        aria-label={`Graf: aktivita enzymů ${enzymes.join(", ")} v procentech výchozí hodnoty. Šipkami procházíte časy.`}
        onKeyDown={handlers.onKeyDown}
        className="select-none rounded-ctl"
        style={{ height }}
      >
        {width > 0 && (
          <svg ref={svgRef} width={width} height={height} viewBox={`0 0 ${width} ${height}`}
               className="block overflow-visible" style={{ touchAction: "pan-y", cursor: "crosshair" }}
               onPointerMove={handlers.onPointerMove} onPointerDown={handlers.onPointerDown}
               onPointerLeave={handlers.onPointerLeave} aria-hidden>
            {yTicks.map((v) => (
              <g key={v}>
                <line x1={MARGIN.left} x2={MARGIN.left + plotW} y1={getY(v)} y2={getY(v)} stroke="rgb(var(--line-soft))" />
                <text x={width - 2} y={getY(v) + 4} textAnchor="end" fontSize={12} fill="rgb(var(--ink-2))"
                      style={{ fontVariantNumeric: "tabular-nums" }}>
                  {num(v)} %
                </text>
              </g>
            ))}
            <line x1={MARGIN.left} x2={MARGIN.left + plotW} y1={getY(100)} y2={getY(100)}
                  stroke="rgb(var(--ink-2))" strokeDasharray="2 3" />
            <text x={MARGIN.left + 6} y={getY(100) - 6} fontSize={12} fill="rgb(var(--ink-2))" {...HALO}>
              výchozí aktivita
            </text>
            <line x1={MARGIN.left} x2={MARGIN.left + plotW} y1={yAxis} y2={yAxis} stroke="rgb(var(--line))" />
            {xTicks.map((t) => (
              <text key={t} x={Math.min(getX(t), width - MARGIN.right - 2)} y={yAxis + 19}
                    textAnchor={t === 0 ? "start" : t >= durationH ? "end" : "middle"} fontSize={12}
                    fill="rgb(var(--ink-2))" style={{ fontVariantNumeric: "tabular-nums" }}>
                {t === 0 ? "0 h" : t}
              </text>
            ))}
            {enzymes.map((e, i) => {
              const s = STROKES[i % STROKES.length];
              const d = pct[e].map((v, k) => `${k ? "L" : "M"}${getX(timeGrid[k]).toFixed(1)} ${getY(v).toFixed(1)}`).join(" ");
              return (
                <path key={e} d={d} fill="none" stroke={`rgb(var(${s.ink}))`} strokeWidth={2}
                      strokeDasharray={s.dash} strokeLinejoin="round" />
              );
            })}
            <line x1={getX(timeGrid[at])} x2={getX(timeGrid[at])} y1={MARGIN.top} y2={yAxis}
                  stroke="rgb(var(--ink) / 0.35)" />
          </svg>
        )}
      </div>
      <ul className="mt-4 space-y-1.5 text-ui">
        {enzymes.map((e, i) => (
          <li key={e} className="flex flex-wrap items-baseline gap-x-3">
            <span aria-hidden className="w-5 border-t-2"
                  style={{ borderColor: `rgb(var(${STROKES[i % STROKES.length].ink}))`,
                           borderStyle: i % STROKES.length ? "dashed" : "solid" }} />
            <span className="font-medium">{e}</span>
            <span className="figures">{num(pct[e][at])} %</span>
            {drivers[e]?.length > 0 && <span className="text-ink-2">mění: {drivers[e].join(", ")}</span>}
          </li>
        ))}
      </ul>
    </div>
  );
}
