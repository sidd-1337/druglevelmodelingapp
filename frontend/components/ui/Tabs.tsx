"use client";

import React, { KeyboardEvent, useRef } from "react";

export interface TabItem<T extends string> {
  id: T;
  label: string;
}

interface TabsProps<T extends string> {
  label: string;
  /** Prefix that links each tab to its panel through aria-controls. */
  idPrefix: string;
  tabs: readonly TabItem<T>[];
  value: T;
  onChange: (id: T) => void;
}

export const tabId = (prefix: string, id: string) => `${prefix}-tab-${id}`;
export const panelId = (prefix: string, id: string) => `${prefix}-panel-${id}`;

/** Text tabs with an ink underline on the current one, set on a full-width hairline. */
export default function Tabs<T extends string>({
  label,
  idPrefix,
  tabs,
  value,
  onChange,
}: TabsProps<T>) {
  const index = Math.max(0, tabs.findIndex((t) => t.id === value));
  const refs = useRef<(HTMLButtonElement | null)[]>([]);

  const onKeyDown = (e: KeyboardEvent<HTMLDivElement>) => {
    const dir = e.key === "ArrowRight" ? 1 : e.key === "ArrowLeft" ? -1 : 0;
    if (!dir) return;
    e.preventDefault();
    const next = (index + dir + tabs.length) % tabs.length;
    onChange(tabs[next].id);
    refs.current[next]?.focus();
  };

  return (
    <div className="border-b border-line-soft">
      <div className="page">
        <div
          role="tablist"
          aria-label={label}
          onKeyDown={onKeyDown}
          className="no-scrollbar -mb-px flex gap-7 overflow-x-auto"
        >
          {tabs.map((tab, i) => {
            const active = tab.id === value;
            return (
              <button
                key={tab.id}
                ref={(el) => {
                  refs.current[i] = el;
                }}
                id={tabId(idPrefix, tab.id)}
                type="button"
                role="tab"
                aria-selected={active}
                aria-controls={panelId(idPrefix, tab.id)}
                tabIndex={active ? 0 : -1}
                onClick={() => onChange(tab.id)}
                className={`shrink-0 border-b-2 py-4 text-ui font-medium transition-colors ${
                  active
                    ? "border-ink text-ink"
                    : "border-transparent text-ink-2 hover:text-ink"
                }`}
              >
                {tab.label}
              </button>
            );
          })}
        </div>
      </div>
    </div>
  );
}
