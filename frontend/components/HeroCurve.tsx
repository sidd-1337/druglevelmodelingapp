import React from "react";
import { cz } from "@/lib/typo";

// A real one-compartment oral model (Bateman function, four doses), normalised so the
// highest peak is 1. It is an illustration, not a drug: no axes, no numbers.
const W = 1440;
const H = 250;
const TOP = 52;
const BASE = 214;

const KA = 1.3; // absorption rate, 1/h
const T_HALF = 5; // elimination half-life, h
const TAU = 9; // dosing interval, h
const DOSES = 4;
const HORIZON = 38; // h
const STEP = 0.1; // h

const WINDOW_LOW = 0.3; // illustrative MEC, relative to the highest peak
const WINDOW_HIGH = 1.18; // illustrative MTC

function buildCurve() {
  const ke = Math.LN2 / T_HALF;
  const single = (t: number) =>
    t <= 0 ? 0 : (KA / (KA - ke)) * (Math.exp(-ke * t) - Math.exp(-KA * t));
  const total = (t: number) => {
    let c = 0;
    for (let i = 0; i < DOSES; i++) c += single(t - i * TAU);
    return c;
  };

  const points: Array<[number, number]> = [];
  let max = 0;
  for (let t = 0; t <= HORIZON + 1e-9; t += STEP) {
    const c = total(t);
    points.push([t, c]);
    if (c > max) max = c;
  }

  const x = (t: number) => (t / HORIZON) * W;
  const y = (c: number) => BASE - (c / max) * (BASE - TOP);
  const d = points
    .map(([t, c], i) => `${i === 0 ? "M" : "L"}${x(t).toFixed(1)} ${y(c).toFixed(1)}`)
    .join(" ");

  return { d, x, y };
}

export default function HeroCurve() {
  const { d, x, y } = buildCurve();
  const yHigh = y(WINDOW_HIGH);
  const yLow = y(WINDOW_LOW);

  return (
    <figure className="mx-auto w-full max-w-[1680px]">
      <div
        className="aspect-[1440/250] max-h-[300px] min-h-[190px] w-full"
        style={{
          WebkitMaskImage:
            "linear-gradient(to right, transparent, #000 5%, #000 92%, transparent)",
          maskImage: "linear-gradient(to right, transparent, #000 5%, #000 92%, transparent)",
        }}
      >
        <svg
          viewBox={`0 0 ${W} ${H}`}
          preserveAspectRatio="xMinYMid slice"
          className="h-full w-full"
          aria-hidden
        >
          <rect
            className="curve-band"
            x={0}
            y={yHigh}
            width={W}
            height={yLow - yHigh}
            fill="rgb(var(--ok-fill) / 0.1)"
          />
          <line x1={0} x2={W} y1={BASE} y2={BASE} stroke="rgb(var(--line))" strokeWidth={1} />
          {Array.from({ length: DOSES }, (_, i) => {
            const p = (i * TAU) / HORIZON;
            const eased = p < 0.5 ? 2 * p * p : 1 - 2 * (1 - p) * (1 - p);
            return (
              <line
                key={i}
                className="curve-tick"
                x1={x(i * TAU)}
                x2={x(i * TAU)}
                y1={BASE}
                y2={BASE + 12}
                stroke="rgb(var(--ink-2))"
                strokeWidth={1.5}
                style={{ animationDelay: `${(0.35 + 2.6 * eased).toFixed(2)}s` }}
              />
            );
          })}
          <path
            className="curve-line"
            d={d}
            pathLength={1}
            fill="none"
            stroke="rgb(var(--ink))"
            strokeWidth={3}
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
      </div>
      <figcaption className="page mt-3 text-center text-small text-ink-2">
        {cz("Ilustrace: hladina léku při opakovaném dávkování. Zeleně je terapeutické okno, čárky označují dávky.")}
      </figcaption>
    </figure>
  );
}
