"use client";

import React, { useState } from "react";
import { CopayAgeTier, CopayTierInfo } from "@/lib/types";
import { kc } from "@/lib/format";
import Segmented, { SegmentOption } from "./ui/Segmented";
import SectionHead from "./ui/SectionHead";
import Slider from "./ui/Slider";

interface CopayCalculatorProps {
  maxDoplatek: number;
  zapocitatelnyDoplatek: number;
}

const TIERS: CopayTierInfo[] = [
  {
    id: "adult",
    label: "Dospělý, 19–64 let",
    limitKč: 5000,
    description: "Standardní zákonný limit podle § 16b zákona č. 48/1997 Sb.",
  },
  {
    id: "child_senior65",
    label: "Dítě do 18 let, senior 65–69 let",
    limitKč: 1000,
    description: "Zvýhodněný limit pro děti do 18 let a seniory od 65 do 69 let.",
  },
  {
    id: "senior70_disabled",
    label: "Senior od 70 let, invalidní",
    limitKč: 500,
    description: "Nejnižší limit pro pojištěnce od 70 let a invalidní osoby.",
  },
];

const OPTIONS: SegmentOption<CopayAgeTier>[] = TIERS.map((t) => ({
  value: t.id,
  label: t.label,
  hint: `${kc(t.limitKč)} ročně`,
}));

export default function CopayCalculator({
  maxDoplatek,
  zapocitatelnyDoplatek,
}: CopayCalculatorProps) {
  const [selectedTier, setSelectedTier] = useState<CopayAgeTier>("adult");
  const [packsPerYear, setPacksPerYear] = useState<number>(12);

  const activeTier = TIERS.find((t) => t.id === selectedTier) || TIERS[0];
  const totalPatientPaid = packsPerYear * maxDoplatek;
  const totalCountedToLimit = packsPerYear * zapocitatelnyDoplatek;

  const limitReached = totalCountedToLimit >= activeTier.limitKč;
  const refundedByInsurance = Math.max(0, totalCountedToLimit - activeTier.limitKč);
  const netPatientExpense = totalPatientPaid - refundedByInsurance;
  const percentOfLimit = Math.min(100, Math.round((totalCountedToLimit / activeTier.limitKč) * 100));

  const results = [
    { label: "Zaplatíte v lékárně", value: kc(totalPatientPaid) },
    { label: "Započte se do limitu", value: kc(totalCountedToLimit) },
    { label: "Pojišťovna vrátí", value: kc(refundedByInsurance) },
  ];

  return (
    <section aria-labelledby="limit">
      <SectionHead id="limit" title="Roční ochranný limit na doplatky">
        Částku započitatelných doplatků nad limit vám pojišťovna vrátí (§ 16b zákona č. 48/1997 Sb.).
      </SectionHead>

      <div className="panel mt-8 space-y-8 p-6 sm:p-8">
        <div>
          <Segmented
            label="Věková skupina"
            value={selectedTier}
            options={OPTIONS}
            onChange={setSelectedTier}
          />
          <p className="mt-3 text-ui text-ink-2">{activeTier.description}</p>
        </div>

        <Slider
          label="Předpokládaná spotřeba za rok"
          display={`${packsPerYear} balení`}
          value={packsPerYear}
          min={1}
          max={36}
          step={1}
          onChange={setPacksPerYear}
        />

        <div>
          <div className="flex flex-wrap items-baseline justify-between gap-x-6 text-ui">
            <span className="text-ink-2">Vyčerpání ročního limitu</span>
            <span className="figures font-semibold">
              {kc(totalCountedToLimit)} z {kc(activeTier.limitKč)} ({percentOfLimit}&nbsp;%)
            </span>
          </div>
          <div
            role="progressbar"
            aria-label="Vyčerpání ročního limitu"
            aria-valuemin={0}
            aria-valuemax={100}
            aria-valuenow={percentOfLimit}
            className="mt-2 h-1 overflow-hidden rounded-full bg-line"
          >
            <div
              className={`h-full rounded-full transition-[width] duration-300 ${
                limitReached ? "bg-ok-fill" : "bg-accent"
              }`}
              style={{ width: `${percentOfLimit}%` }}
            />
          </div>
        </div>

        <dl className="border-b border-line">
          {results.map((r) => (
            <div key={r.label} className="flex items-baseline justify-between gap-6 border-t border-line py-3">
              <dt className="text-ui text-ink-2">{r.label}</dt>
              <dd className="figures text-body">{r.value}</dd>
            </div>
          ))}
          <div className="flex items-baseline justify-between gap-6 border-t border-line py-4">
            <dt className="text-body font-semibold">Čistý náklad za rok</dt>
            <dd className="figures text-sub">{kc(netPatientExpense)}</dd>
          </div>
        </dl>
      </div>
    </section>
  );
}
