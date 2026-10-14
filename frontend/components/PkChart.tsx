"use client";

import React, { KeyboardEvent, PointerEvent, useEffect, useMemo, useRef, useState } from "react";
import { PKCalculationPoint } from "@/lib/types";
import { num, sig } from "@/lib/format";
import Status from "@/components/ui/Status";

interface PkChartProps {
  data: PKCalculationPoint[];
  mec: number;
  mtc: number;
  unit: string;
  durationH: number;
  doseTimes?: number[];
  secondaryData?: PKCalculationPoint[];
  secondaryLabel?: string;
  primaryLabel?: string;
}

// Keeps a label legible where the curve passes behind it.
const HALO = { stroke: "rgb(var(--bg))", strokeWidth: 4, paintOrder: "stroke" } as const;

const MARGIN = { top: 22, right: 46, bottom: 28, left: 2 };

/** Width of the element, tracked so the SVG is drawn at real pixel size and text stays legible. */
function useWidth<T extends HTMLElement>() {
  const ref = useRef<T>(null);
  const [width, setWidth] = useState(0);
  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    setWidth(Math.round(el.getBoundingClientRect().width));
    const observer = new ResizeObserver(([entry]) => setWidth(Math.round(entry.contentRect.width)));
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  return [ref, width] as const;
}

function niceStep(max: number, count: number): number {
  const raw = max / count;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / mag;
  const m = n <= 1 ? 1 : n <= 2 ? 2 : n <= 2.5 ? 2.5 : n <= 5 ? 5 : 10;
  return m * mag;
}

function decimalsOf(step: number): number {
  let d = 0;
  while (Math.abs(step * 10 ** d - Math.round(step * 10 ** d)) > 1e-9 && d < 6) d++;
  return d;
}

function hourStep(durationH: number, maxTicks: number): number {
  const candidates = [1, 2, 3, 4, 6, 8, 12, 24, 48, 96];
  return candidates.find((s) => durationH / s <= maxTicks) ?? 96;
}

function nearestIndex(data: PKCalculationPoint[], t: number): number {
  let lo = 0;
  let hi = data.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (data[mid].timeH < t) lo = mid;
    else hi = mid;
  }
  return Math.abs(data[lo].timeH - t) <= Math.abs(data[hi].timeH - t) ? lo : hi;
}

export default function PkChart({
  data,
  mec,
  mtc,
  unit,
  durationH,
  doseTimes,
  secondaryData,
  secondaryLabel = "2. dávka",
  primaryLabel = "Hladina v krvi",
}: PkChartProps) {
  const [wrapRef, width] = useWidth<HTMLDivElement>();
  const svgRef = useRef<SVGSVGElement>(null);
  const [selected, setSelected] = useState<number | null>(null);

  const height = Math.round(Math.min(380, Math.max(240, width * 0.42)));
  const plotW = Math.max(0, width - MARGIN.left - MARGIN.right);
  const plotH = height - MARGIN.top - MARGIN.bottom;

  const maxC = useMemo(() => {
    const all = secondaryData ? [...data, ...secondaryData] : data;
    const maxData = all.reduce((m, p) => Math.max(m, p.concentration), 0);
    return Math.max(mtc * 1.15, maxData * 1.1) || 1;
  }, [data, secondaryData, mtc]);

  const peakIndex = useMemo(() => {
    let best = 0;
    data.forEach((p, i) => {
      if (p.concentration > data[best].concentration) best = i;
    });
    return best;
  }, [data]);

  const getX = (t: number) => MARGIN.left + (Math.max(0, Math.min(durationH, t)) / durationH) * plotW;
  const getY = (c: number) => MARGIN.top + plotH - (Math.max(0, Math.min(maxC, c)) / maxC) * plotH;

  const toPath = (points: PKCalculationPoint[]) =>
    points
      .map((p, i) => `${i === 0 ? "M" : "L"}${getX(p.timeH).toFixed(1)} ${getY(p.concentration).toFixed(1)}`)
      .join(" ");

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const primaryPath = useMemo(() => toPath(data), [data, width, height, maxC, durationH]);
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const secondaryPath = useMemo(
    () => (secondaryData ? toPath(secondaryData) : null),
    [secondaryData, width, height, maxC, durationH]
  );

  if (!data || data.length === 0) return null;

  const cur = data[Math.min(selected ?? peakIndex, data.length - 1)];
  const curSecondary =
    secondaryData && secondaryData.length
      ? secondaryData[nearestIndex(secondaryData, cur.timeH)]
      : null;

  const scrub = (clientX: number) => {
    const rect = svgRef.current?.getBoundingClientRect();
    if (!rect || plotW <= 0) return;
    const t = ((clientX - rect.left - MARGIN.left) / plotW) * durationH;
    setSelected(nearestIndex(data, Math.max(0, Math.min(durationH, t))));
  };

  const onPointerMove = (e: PointerEvent<SVGSVGElement>) => {
    // Touch scrubs only while a finger is down; a mouse follows the pointer.
    if (e.pointerType === "mouse" || e.buttons > 0) scrub(e.clientX);
  };
  const onPointerDown = (e: PointerEvent<SVGSVGElement>) => scrub(e.clientX);
  const onPointerLeave = (e: PointerEvent<SVGSVGElement>) => {
    if (e.pointerType === "mouse") setSelected(null);
  };

  const onKeyDown = (e: KeyboardEvent<HTMLDivElement>) => {
    const step = Math.max(1, Math.floor(data.length / 120));
    const at = selected ?? peakIndex;
    if (e.key === "ArrowRight") setSelected(Math.min(data.length - 1, at + step));
    else if (e.key === "ArrowLeft") setSelected(Math.max(0, at - step));
    else if (e.key === "Home") setSelected(0);
    else if (e.key === "End") setSelected(data.length - 1);
    else if (e.key === "Escape") setSelected(null);
    else return;
    e.preventDefault();
  };

  // Axes
  const yStep = niceStep(maxC, 4);
  const yDecimals = decimalsOf(yStep);
  const yTicks: number[] = [];
  for (let v = 0; v <= maxC + 1e-9; v += yStep) yTicks.push(v);

  const xStep = hourStep(durationH, Math.max(3, Math.floor(plotW / 72)));
  const xTicks: number[] = [];
  for (let t = 0; t <= durationH + 1e-9; t += xStep) xTicks.push(t);

  const yMtc = getY(mtc);
  const yMec = getY(mec);
  const yZero = getY(0);
  const yTop = MARGIN.top;

  const zone =
    cur.concentration > mtc ? "bad" : cur.concentration >= mec ? "ok" : "neutral";
  const zoneText =
    zone === "bad"
      ? "Toxická hladina, nad MTC"
      : zone === "ok"
      ? "Terapeutické rozmezí"
      : "Pod terapeutickou hladinou";

  const peak = data[peakIndex];
  const label = `Graf: ${primaryLabel.toLowerCase()} v čase. Vrchol ${sig(peak.concentration)} ${unit} v ${num(peak.timeH, 1)} h. Terapeutické okno ${sig(mec)} až ${sig(mtc)} ${unit}.`;

  return (
    <div>
      {/* Readout: the value under the finger or pointer; the peak until you scrub. */}
      <div className="flex flex-wrap items-end justify-between gap-x-8 gap-y-3">
        <div>
          <p className="text-ui text-ink-2">
            {selected === null ? "Vrchol" : "V čase"} {"t = "}
            <span className="figures">{num(cur.timeH, 1)}</span>
            {" h"}
          </p>
          <p className="figures text-figure">
            {sig(cur.concentration)}
            <span className="ml-2 text-h3 font-normal text-ink-2">{unit}</span>
          </p>
        </div>
        <div className="text-ui sm:text-right">
          <Status tone={zone}>{zoneText}</Status>
          {curSecondary && (
            <p className="figures mt-0.5 text-ink-2">
              {secondaryLabel}: {sig(curSecondary.concentration)}
              {" "}
              {unit}
            </p>
          )}
        </div>
      </div>

      <div
        ref={wrapRef}
        tabIndex={0}
        role="group"
        aria-label={`${label} Šipkami vlevo a vpravo procházíte jednotlivé časy.`}
        onKeyDown={onKeyDown}
        className="mt-5 select-none rounded-ctl"
        style={{ height: height || 240 }}
      >
        {width > 0 && (
          <svg
            ref={svgRef}
            width={width}
            height={height}
            viewBox={`0 0 ${width} ${height}`}
            className="block overflow-visible"
            style={{ touchAction: "pan-y", cursor: "crosshair" }}
            onPointerMove={onPointerMove}
            onPointerDown={onPointerDown}
            onPointerLeave={onPointerLeave}
            aria-hidden
          >
            {/* Zones */}
            <rect
              x={MARGIN.left}
              y={yTop}
              width={plotW}
              height={Math.max(0, yMtc - yTop)}
              fill="rgb(var(--bad-fill) / 0.08)"
            />
            <rect
              x={MARGIN.left}
              y={Math.min(yMtc, yMec)}
              width={plotW}
              height={Math.abs(yMec - yMtc)}
              fill="rgb(var(--ok-fill) / 0.13)"
            />

            {/* Grid and y axis, labels on the right */}
            {yTicks.map((v) => (
              <g key={v}>
                <line
                  x1={MARGIN.left}
                  x2={MARGIN.left + plotW}
                  y1={getY(v)}
                  y2={getY(v)}
                  stroke="rgb(var(--line-soft))"
                />
                <text
                  x={width - 2}
                  y={getY(v) + 4}
                  textAnchor="end"
                  fontSize={12}
                  fill="rgb(var(--ink-2))"
                  style={{ fontVariantNumeric: "tabular-nums" }}
                >
                  {num(v, yDecimals)}
                </text>
              </g>
            ))}
            <text
              x={width - 2}
              y={11}
              textAnchor="end"
              fontSize={12}
              fill="rgb(var(--ink-2))"
            >
              {unit}
            </text>

            {/* Therapeutic limits */}
            <line
              x1={MARGIN.left}
              x2={MARGIN.left + plotW}
              y1={yMtc}
              y2={yMtc}
              stroke="rgb(var(--bad-fill) / 0.7)"
              strokeDasharray="4 4"
            />
            <text
              x={MARGIN.left + 6}
              y={yMtc - yTop < 16 ? yMtc + 15 : yMtc - 6}
              fontSize={12}
              fontWeight={500}
              fill="rgb(var(--bad))"
              {...HALO}
            >
              MTC {sig(mtc)}
            </text>
            <line
              x1={MARGIN.left}
              x2={MARGIN.left + plotW}
              y1={yMec}
              y2={yMec}
              stroke="rgb(var(--ok-fill) / 0.8)"
              strokeDasharray="4 4"
            />
            <text
              x={MARGIN.left + 6}
              y={yMec - 6}
              fontSize={12}
              fontWeight={500}
              fill="rgb(var(--ok))"
              {...HALO}
            >
              MEC {sig(mec)}
            </text>

            {/* x axis */}
            <line
              x1={MARGIN.left}
              x2={MARGIN.left + plotW}
              y1={yZero}
              y2={yZero}
              stroke="rgb(var(--line))"
            />
            {xTicks.map((t) => (
              <text
                key={t}
                x={Math.min(getX(t), width - MARGIN.right - 2)}
                y={yZero + 19}
                textAnchor={t === 0 ? "start" : t >= durationH ? "end" : "middle"}
                fontSize={12}
                fill="rgb(var(--ink-2))"
                style={{ fontVariantNumeric: "tabular-nums" }}
              >
                {t === 0 ? "0 h" : t}
              </text>
            ))}

            {/* Doses */}
            {doseTimes?.map((td, i) =>
              td > durationH ? null : (
                <g key={`${td}-${i}`}>
                  <line
                    x1={getX(td)}
                    x2={getX(td)}
                    y1={yTop}
                    y2={yZero}
                    stroke="rgb(var(--ink) / 0.1)"
                  />
                  <circle cx={getX(td)} cy={yZero} r={3} fill="rgb(var(--ink-2))" />
                </g>
              )
            )}

            {/* Curves */}
            {secondaryPath && (
              <path
                d={secondaryPath}
                fill="none"
                stroke="rgb(var(--ink-2))"
                strokeWidth={2}
                strokeDasharray="5 4"
                strokeLinejoin="round"
              />
            )}
            <path
              d={primaryPath}
              fill="none"
              stroke="rgb(var(--ink))"
              strokeWidth={2.25}
              strokeLinejoin="round"
              strokeLinecap="round"
            />

            {/* Marker */}
            <line
              x1={getX(cur.timeH)}
              x2={getX(cur.timeH)}
              y1={yTop}
              y2={yZero}
              stroke="rgb(var(--ink) / 0.35)"
            />
            <circle
              cx={getX(cur.timeH)}
              cy={getY(cur.concentration)}
              r={5.5}
              fill="rgb(var(--ink))"
              stroke="rgb(var(--bg))"
              strokeWidth={2.5}
            />
          </svg>
        )}
      </div>

      <ul className="mt-4 flex flex-wrap gap-x-6 gap-y-1.5 text-small text-ink-2">
        <li className="flex items-center gap-2">
          <span aria-hidden className="h-0.5 w-4 rounded-full bg-ink" />
          {primaryLabel}
        </li>
        {secondaryData && (
          <li className="flex items-center gap-2">
            <span
              aria-hidden
              className="h-0 w-4 border-t-2 border-dashed border-ink-2"
            />
            {secondaryLabel}
          </li>
        )}
        {doseTimes && (
          <li className="flex items-center gap-2">
            <span aria-hidden className="h-1.5 w-1.5 rounded-full bg-ink-2" />
            Dávka
          </li>
        )}
        <li className="flex items-center gap-2">
          <span aria-hidden className="h-2.5 w-4 rounded-sm bg-ok-fill/25" />
          Terapeutické okno (MEC až MTC)
        </li>
        <li className="flex items-center gap-2">
          <span aria-hidden className="h-2.5 w-4 rounded-sm bg-bad-fill/20" />
          Toxická zóna (nad MTC)
        </li>
      </ul>
    </div>
  );
}
