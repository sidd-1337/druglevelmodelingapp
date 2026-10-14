import React from "react";

export interface SpecItem {
  label: string;
  value: React.ReactNode;
  /** Quieter second line under the value. */
  note?: React.ReactNode;
}

/** Label on the left, value on the right, a hairline between rows. No boxes. */
export default function SpecList({ items }: { items: SpecItem[] }) {
  return (
    <dl className="border-b border-line-soft">
      {items.map((item) => (
        <div
          key={item.label}
          className="grid gap-1 border-t border-line-soft py-4 sm:grid-cols-[13rem_minmax(0,1fr)] sm:gap-8"
        >
          <dt className="text-ui text-ink-2">{item.label}</dt>
          <dd className="text-body text-ink">
            {item.value}
            {item.note && <span className="figures mt-0.5 block text-ui text-ink-2">{item.note}</span>}
          </dd>
        </div>
      ))}
    </dl>
  );
}
