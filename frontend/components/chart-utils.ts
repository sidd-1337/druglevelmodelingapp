"use client";

import React, { useEffect, useRef, useState } from "react";

// Keeps a label legible where a curve passes behind it.
export const HALO = { stroke: "rgb(var(--bg))", strokeWidth: 4, paintOrder: "stroke" } as const;

/** Width of the element, tracked so the SVG is drawn at real pixel size and text stays legible. */
export function useWidth<T extends HTMLElement>() {
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

export function niceStep(max: number, count: number): number {
  const raw = max / count;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / mag;
  const m = n <= 1 ? 1 : n <= 2 ? 2 : n <= 2.5 ? 2.5 : n <= 5 ? 5 : 10;
  return m * mag;
}

export function decimalsOf(step: number): number {
  let d = 0;
  while (Math.abs(step * 10 ** d - Math.round(step * 10 ** d)) > 1e-9 && d < 6) d++;
  return d;
}

export function hourStep(durationH: number, maxTicks: number): number {
  const candidates = [1, 2, 3, 4, 6, 8, 12, 24, 48, 96];
  return candidates.find((s) => durationH / s <= maxTicks) ?? 96;
}

/** Index of the point nearest to ``t`` among ``n`` sorted times read through ``timeAt``. */
export function nearestTime(n: number, timeAt: (i: number) => number, t: number): number {
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (timeAt(mid) < t) lo = mid;
    else hi = mid;
  }
  return Math.abs(timeAt(lo) - t) <= Math.abs(timeAt(hi) - t) ? lo : hi;
}

/**
 * Pointer and keyboard handlers that pick a time on a chart by scrubbing.
 * A mouse follows the pointer; touch scrubs only while a finger is down.
 */
export function scrubHandlers(opts: {
  svgRef: React.RefObject<SVGSVGElement | null>;
  left: number;
  plotW: number;
  durationH: number;
  n: number;
  timeAt: (i: number) => number;
  at: number;
  select: (i: number | null) => void;
}) {
  const { svgRef, left, plotW, durationH, n, timeAt, at, select } = opts;
  const scrub = (clientX: number) => {
    const rect = svgRef.current?.getBoundingClientRect(); // read at event time, not render time
    if (!rect || plotW <= 0 || n === 0) return;
    const t = ((clientX - rect.left - left) / plotW) * durationH;
    select(nearestTime(n, timeAt, Math.max(0, Math.min(durationH, t))));
  };
  return {
    onPointerMove: (e: React.PointerEvent<SVGSVGElement>) => {
      if (e.pointerType === "mouse" || e.buttons > 0) scrub(e.clientX);
    },
    onPointerDown: (e: React.PointerEvent<SVGSVGElement>) => scrub(e.clientX),
    onPointerLeave: (e: React.PointerEvent<SVGSVGElement>) => {
      if (e.pointerType === "mouse") select(null);
    },
    onKeyDown: (e: React.KeyboardEvent<HTMLElement>) => {
      const step = Math.max(1, Math.floor(n / 120));
      if (e.key === "ArrowRight") select(Math.min(n - 1, at + step));
      else if (e.key === "ArrowLeft") select(Math.max(0, at - step));
      else if (e.key === "Home") select(0);
      else if (e.key === "End") select(n - 1);
      else if (e.key === "Escape") select(null);
      else return;
      e.preventDefault();
    },
  };
}
