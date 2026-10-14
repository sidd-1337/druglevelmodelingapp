"use client";

import React, { KeyboardEvent, useRef } from "react";

export interface SegmentOption<T extends string> {
  value: T;
  label: string;
  /** Optional second line, e.g. the limit for a tier. */
  hint?: string;
}

interface SegmentedProps<T extends string> {
  label: string;
  value: T;
  options: SegmentOption<T>[];
  onChange: (value: T) => void;
}

/** Equal-width segmented control. The thumb slides to the chosen option. */
export default function Segmented<T extends string>({
  label,
  value,
  options,
  onChange,
}: SegmentedProps<T>) {
  const index = Math.max(0, options.findIndex((o) => o.value === value));
  const refs = useRef<(HTMLButtonElement | null)[]>([]);

  const onKeyDown = (e: KeyboardEvent<HTMLDivElement>) => {
    const dir =
      e.key === "ArrowRight" || e.key === "ArrowDown"
        ? 1
        : e.key === "ArrowLeft" || e.key === "ArrowUp"
        ? -1
        : 0;
    if (!dir) return;
    e.preventDefault();
    const next = (index + dir + options.length) % options.length;
    onChange(options[next].value);
    refs.current[next]?.focus();
  };

  return (
    <div
      role="radiogroup"
      aria-label={label}
      onKeyDown={onKeyDown}
      className="relative grid rounded-ctl bg-ink/[0.07] p-0.5"
      style={{ gridTemplateColumns: `repeat(${options.length}, minmax(0, 1fr))` }}
    >
      <span
        aria-hidden
        className="absolute inset-y-0.5 left-0.5 rounded-[10px] bg-thumb shadow-thumb transition-transform duration-200 ease-out motion-reduce:transition-none"
        style={{
          width: `calc((100% - 4px) / ${options.length})`,
          transform: `translateX(${index * 100}%)`,
        }}
      />
      {options.map((option, i) => {
        const active = i === index;
        return (
          <button
            key={option.value}
            ref={(el) => {
              refs.current[i] = el;
            }}
            type="button"
            role="radio"
            aria-checked={active}
            tabIndex={active ? 0 : -1}
            onClick={() => onChange(option.value)}
            className={`relative z-10 rounded-[10px] px-2 py-2 text-center text-ui font-medium transition-colors ${
              active ? "text-ink" : "text-ink-2 hover:text-ink"
            }`}
          >
            <span className="block">{option.label}</span>
            {option.hint && (
              <span className="figures mt-0.5 block text-small font-normal text-ink-2">
                {option.hint}
              </span>
            )}
          </button>
        );
      })}
    </div>
  );
}
