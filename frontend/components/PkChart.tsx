"use client";

import React, { useState } from "react";
import { PKCalculationPoint } from "@/lib/types";

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

export default function PkChart({
  data,
  mec,
  mtc,
  unit,
  durationH,
  doseTimes,
  secondaryData,
  secondaryLabel = "2. dávka (redose)",
  primaryLabel = "Hladina v krvi C(t)",
}: PkChartProps) {
  const [hoverIndex, setHoverIndex] = useState<number | null>(null);

  if (!data || data.length === 0) return null;

  const width = 800;
  const height = 360;
  const padding = { top: 25, right: 30, bottom: 45, left: 60 };
  const chartW = width - padding.left - padding.right;
  const chartH = height - padding.top - padding.bottom;

  const allPoints = secondaryData ? [...data, ...secondaryData] : data;
  const maxDataC = Math.max(...allPoints.map((d) => d.concentration));
  const maxC = Math.max(mtc * 1.3, maxDataC * 1.15, 0.1);

  const getX = (t: number) => padding.left + (Math.max(0, Math.min(durationH, t)) / durationH) * chartW;
  const getY = (c: number) => padding.top + chartH - (Math.max(0, c) / maxC) * chartH;

  // Build SVG path
  const pathD = data
    .map((pt, i) => `${i === 0 ? "M" : "L"} ${getX(pt.timeH)} ${getY(pt.concentration)}`)
    .join(" ");

  const secondaryPathD = secondaryData
    ? secondaryData
        .map((pt, i) => `${i === 0 ? "M" : "L"} ${getX(pt.timeH)} ${getY(pt.concentration)}`)
        .join(" ")
    : null;

  // Find Cmax
  let peakPoint = data[0];
  for (const pt of data) {
    if (pt.concentration > peakPoint.concentration) {
      peakPoint = pt;
    }
  }

  // Zone coordinates
  const yMtc = getY(mtc);
  const yMec = getY(mec);
  const yZero = getY(0);

  const handleMouseMove = (e: React.MouseEvent<SVGSVGElement>) => {
    const rect = e.currentTarget.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const relX = (mouseX / rect.width) * width;
    const tHover = ((relX - padding.left) / chartW) * durationH;

    if (tHover >= 0 && tHover <= durationH) {
      const idx = Math.round((tHover / durationH) * (data.length - 1));
      if (idx >= 0 && idx < data.length) {
        setHoverIndex(idx);
      }
    } else {
      setHoverIndex(null);
    }
  };

  const hovered = hoverIndex !== null ? data[hoverIndex] : null;
  const secondaryHovered = hoverIndex !== null && secondaryData ? secondaryData[hoverIndex] : null;

  return (
    <div className="apple-card p-5 select-none">
      <div className="flex flex-wrap items-center justify-between mb-4 text-xs gap-3">
        <div className="flex flex-wrap items-center gap-4">
          <div className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-full bg-blue-600 dark:bg-blue-400"></span>
            <span className="text-zinc-800 dark:text-zinc-200 font-medium">{primaryLabel}</span>
          </div>
          {secondaryData && (
            <div className="flex items-center gap-1.5">
              <span className="w-2.5 h-2.5 rounded-full bg-amber-500"></span>
              <span className="text-zinc-800 dark:text-zinc-200 font-medium">{secondaryLabel}</span>
            </div>
          )}
          <div className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm bg-emerald-500/20 border border-emerald-500/50"></span>
            <span className="text-zinc-500 dark:text-zinc-400">Terapeutické okno (MEC–MTC)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm bg-rose-500/15 border border-rose-500/40"></span>
            <span className="text-zinc-500 dark:text-zinc-400">Toxická zóna (&gt; MTC)</span>
          </div>
        </div>

        {peakPoint && (
          <div className="text-xs px-2.5 py-1 rounded-lg bg-zinc-100 dark:bg-zinc-800 text-zinc-800 dark:text-zinc-200 border border-zinc-200/80 dark:border-zinc-700 font-medium">
            Cmax: <strong>{peakPoint.concentration.toFixed(2)} {unit}</strong> (t = {peakPoint.timeH.toFixed(1)} h)
          </div>
        )}
      </div>

      <div className="relative w-full aspect-[800/360]">
        <svg
          viewBox={`0 0 ${width} ${height}`}
          className="w-full h-full cursor-crosshair overflow-visible"
          onMouseMove={handleMouseMove}
          onMouseLeave={() => setHoverIndex(null)}
        >
          {/* Toxic Zone (> MTC) */}
          <rect
            x={padding.left}
            y={padding.top}
            width={chartW}
            height={Math.max(0, yMtc - padding.top)}
            fill="rgba(239, 68, 68, 0.05)"
          />

          {/* Therapeutic Window (MEC - MTC) */}
          <rect
            x={padding.left}
            y={Math.min(yMtc, yMec)}
            width={chartW}
            height={Math.abs(yMec - yMtc)}
            fill="rgba(16, 185, 129, 0.06)"
          />

          {/* Sub-therapeutic Zone (< MEC) */}
          <rect
            x={padding.left}
            y={yMec}
            width={chartW}
            height={Math.max(0, yZero - yMec)}
            fill="rgba(113, 113, 122, 0.03)"
          />

          {/* Grid lines */}
          {[0, 0.25, 0.5, 0.75, 1].map((pct) => {
            const tVal = pct * durationH;
            const x = getX(tVal);
            return (
              <g key={pct}>
                <line
                  x1={x}
                  y1={padding.top}
                  x2={x}
                  y2={padding.top + chartH}
                  stroke="rgba(113, 113, 122, 0.15)"
                  strokeDasharray="4 4"
                />
                <text
                  x={x}
                  y={height - 15}
                  fill="#71717a"
                  fontSize="10"
                  textAnchor="middle"
                  fontFamily="system-ui, sans-serif"
                >
                  {tVal.toFixed(0)} h
                </text>
              </g>
            );
          })}

          {/* Dose event markers */}
          {doseTimes &&
            doseTimes.map((td, idx) => {
              if (td > durationH) return null;
              const x = getX(td);
              return (
                <g key={`dose-${idx}`}>
                  <line
                    x1={x}
                    y1={padding.top}
                    x2={x}
                    y2={padding.top + chartH}
                    stroke="rgba(0, 113, 227, 0.25)"
                    strokeWidth="1.5"
                    strokeDasharray="3 3"
                  />
                  <circle cx={x} cy={padding.top + chartH + 8} r="3" fill="#0071e3" />
                  <text
                    x={x}
                    y={padding.top + chartH + 22}
                    fill="#0071e3"
                    fontSize="9"
                    textAnchor="middle"
                    fontFamily="monospace"
                  >
                    D{idx + 1}
                  </text>
                </g>
              );
            })}

          {/* MTC line */}
          <line
            x1={padding.left}
            y1={yMtc}
            x2={padding.left + chartW}
            y2={yMtc}
            stroke="#ef4444"
            strokeWidth="1.5"
            strokeDasharray="5 3"
          />
          <text
            x={padding.left + chartW - 5}
            y={yMtc - 6}
            fill="#ef4444"
            fontSize="10"
            textAnchor="end"
            fontFamily="system-ui, sans-serif"
            fontWeight="500"
          >
            MTC: {mtc.toFixed(2)} {unit}
          </text>

          {/* MEC line */}
          <line
            x1={padding.left}
            y1={yMec}
            x2={padding.left + chartW}
            y2={yMec}
            stroke="#10b981"
            strokeWidth="1.5"
            strokeDasharray="5 3"
          />
          <text
            x={padding.left + chartW - 5}
            y={yMec - 6}
            fill="#10b981"
            fontSize="10"
            textAnchor="end"
            fontFamily="system-ui, sans-serif"
            fontWeight="500"
          >
            MEC: {mec.toFixed(2)} {unit}
          </text>

          {/* Y-axis ticks */}
          {[0, maxC * 0.25, maxC * 0.5, maxC * 0.75, maxC].map((cVal, i) => {
            const y = getY(cVal);
            return (
              <g key={i}>
                <line
                  x1={padding.left - 4}
                  y1={y}
                  x2={padding.left + chartW}
                  y2={y}
                  stroke="rgba(113, 113, 122, 0.12)"
                />
                <text
                  x={padding.left - 8}
                  y={y + 3}
                  fill="#71717a"
                  fontSize="10"
                  textAnchor="end"
                  fontFamily="monospace"
                >
                  {cVal.toFixed(1)}
                </text>
              </g>
            );
          })}

          {/* Secondary curve */}
          {secondaryPathD && (
            <path
              d={secondaryPathD}
              fill="none"
              stroke="#f59e0b"
              strokeWidth="2.5"
              strokeDasharray="4 2"
              className="transition-all"
            />
          )}

          {/* Primary curve */}
          <path
            d={pathD}
            fill="none"
            stroke="#0071e3"
            strokeWidth="2.5"
            className="transition-all"
          />

          {/* Hover indicator */}
          {hovered && (
            <g>
              <line
                x1={getX(hovered.timeH)}
                y1={padding.top}
                x2={getX(hovered.timeH)}
                y2={padding.top + chartH}
                stroke="#0071e3"
                strokeWidth="1"
                strokeDasharray="2 2"
              />
              <circle
                cx={getX(hovered.timeH)}
                cy={getY(hovered.concentration)}
                r="4.5"
                fill="#0071e3"
                stroke="#ffffff"
                strokeWidth="2"
              />
              {secondaryHovered && (
                <circle
                  cx={getX(secondaryHovered.timeH)}
                  cy={getY(secondaryHovered.concentration)}
                  r="4.5"
                  fill="#f59e0b"
                  stroke="#ffffff"
                  strokeWidth="2"
                />
              )}
            </g>
          )}
        </svg>

        {/* Hover Tooltip Popup */}
        {hovered && (
          <div
            className="absolute pointer-events-none z-20 backdrop-blur-md bg-white/95 dark:bg-zinc-800/95 border border-zinc-200 dark:border-zinc-700 shadow-xl rounded-xl p-3 text-xs space-y-1 transition-all"
            style={{
              left: `${Math.min(80, Math.max(10, (hovered.timeH / durationH) * 100))}%`,
              top: "15%",
            }}
          >
            <div className="font-semibold text-zinc-900 dark:text-zinc-100 flex items-center justify-between gap-4">
              <span>Čas t:</span>
              <span className="font-semibold text-blue-600 dark:text-blue-400">{hovered.timeH.toFixed(2)} h</span>
            </div>
            <div className="flex items-center justify-between gap-4">
              <span className="text-zinc-500 dark:text-zinc-400">Koncentrace:</span>
              <span className="font-bold text-zinc-900 dark:text-zinc-100">
                {hovered.concentration.toFixed(3)} {unit}
              </span>
            </div>
            {secondaryHovered && (
              <div className="flex items-center justify-between gap-4 pt-1 border-t border-zinc-100 dark:border-zinc-700/60">
                <span className="text-amber-600 dark:text-amber-400">Po redosu:</span>
                <span className="font-bold text-amber-600 dark:text-amber-400">
                  {secondaryHovered.concentration.toFixed(3)} {unit}
                </span>
              </div>
            )}
            <div className="pt-1 text-[10px] border-t border-zinc-100 dark:border-zinc-700/60 font-medium">
              {hovered.concentration > mtc ? (
                <span className="text-rose-600 dark:text-rose-400">⚠️ Toxická hladina (&gt; MTC)</span>
              ) : hovered.concentration >= mec ? (
                <span className="text-emerald-600 dark:text-emerald-400">✅ Terapeutické rozmezí</span>
              ) : (
                <span className="text-zinc-400 dark:text-zinc-500">⚪ Subterapeutická hladina</span>
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
