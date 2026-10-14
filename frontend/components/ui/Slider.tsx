"use client";

import React, { CSSProperties, useId } from "react";

interface SliderProps {
  label: string;
  /** Formatted current value shown beside the label. */
  display: string;
  value: number;
  min: number;
  max: number;
  step: number;
  onChange: (value: number) => void;
}

export default function Slider({ label, display, value, min, max, step, onChange }: SliderProps) {
  const id = useId();
  const pct = max > min ? ((value - min) / (max - min)) * 100 : 0;

  return (
    <div>
      <div className="flex items-baseline justify-between gap-4">
        <label htmlFor={id} className="text-ui text-ink-2">
          {label}
        </label>
        <output htmlFor={id} className="figures text-ui font-semibold text-ink">
          {display}
        </output>
      </div>
      <input
        id={id}
        type="range"
        className="slider mt-1"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(Number(e.target.value))}
        style={{ "--pct": `${pct}%` } as CSSProperties}
      />
    </div>
  );
}
