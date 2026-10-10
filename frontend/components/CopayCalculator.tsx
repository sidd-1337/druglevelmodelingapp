"use client";

import React, { useState } from "react";
import { CopayAgeTier, CopayTierInfo } from "@/lib/types";

interface CopayCalculatorProps {
  maxDoplatek: number;
  zapocitatelnyDoplatek: number;
}

const TIERS: CopayTierInfo[] = [
  {
    id: "adult",
    label: "Dospělý (19–64 let)",
    limitKč: 5000,
    description: "Standardní zákonný limit dle § 16b zákona č. 48/1997 Sb.",
  },
  {
    id: "child_senior65",
    label: "Dítě < 18 / Senior 65–69 let",
    limitKč: 1000,
    description: "Zvýhodněný limit pro děti do 18 let a seniory od 65 do 69 let.",
  },
  {
    id: "senior70_disabled",
    label: "Senior 70+ / Invalida",
    limitKč: 500,
    description: "Maximálně chráněná kategorie pojištěnců od 70 let a invalidních osob.",
  },
];

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

  return (
    <div className="apple-card p-6 space-y-6">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <div>
          <h4 className="font-semibold text-base text-zinc-950 dark:text-zinc-50">
            Kalkulačka ročního ochranného limitu na doplatky
          </h4>
          <p className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5">
            Zákon č. 48/1997 Sb. (§ 16b). Částku započitatelných doplatků nad limit pojišťovna vrací.
          </p>
        </div>
        <span className="text-xs px-3 py-1 rounded-full bg-zinc-100 dark:bg-zinc-800 text-zinc-700 dark:text-zinc-300 border border-zinc-200 dark:border-zinc-700 font-medium">
          Limit: {activeTier.limitKč.toLocaleString("cs-CZ")} Kč / rok
        </span>
      </div>

      {/* Tier selector */}
      <div className="grid grid-cols-1 sm:grid-cols-3 gap-2.5">
        {TIERS.map((tier) => {
          const isSelected = selectedTier === tier.id;
          return (
            <button
              key={tier.id}
              onClick={() => setSelectedTier(tier.id)}
              className={`text-left p-3.5 rounded-xl border transition-all ${
                isSelected
                  ? "bg-blue-50/70 border-blue-600 dark:bg-blue-950/30 dark:border-blue-500 text-blue-950 dark:text-blue-100 ring-1 ring-blue-600 dark:ring-blue-500"
                  : "bg-zinc-50 hover:bg-zinc-100 dark:bg-zinc-800/60 dark:hover:bg-zinc-800 border-zinc-200/80 dark:border-zinc-700/60 text-zinc-600 dark:text-zinc-400"
              }`}
            >
              <div className="text-xs font-semibold">{tier.label}</div>
              <div className="text-sm font-bold tracking-tight text-zinc-900 dark:text-zinc-100 mt-1">
                {tier.limitKč.toLocaleString("cs-CZ")} Kč
              </div>
            </button>
          );
        })}
      </div>

      {/* Pack quantity slider */}
      <div className="bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 rounded-xl p-4">
        <div className="flex items-center justify-between text-xs mb-2">
          <span className="text-zinc-600 dark:text-zinc-300 font-medium">Předpokládaná spotřeba balení za rok:</span>
          <span className="font-bold text-blue-600 dark:text-blue-400 text-sm">
            {packsPerYear} balení / rok
          </span>
        </div>
        <input
          type="range"
          min="1"
          max="36"
          step="1"
          value={packsPerYear}
          onChange={(e) => setPacksPerYear(Number(e.target.value))}
          className="w-full accent-blue-600 dark:accent-blue-500 cursor-pointer"
        />
        <div className="flex justify-between text-[10px] text-zinc-400 dark:text-zinc-500 mt-1">
          <span>1 balení</span>
          <span>12 balení (1× měsíčně)</span>
          <span>24 balení</span>
          <span>36 balení</span>
        </div>
      </div>

      {/* Progress towards annual limit */}
      <div>
        <div className="flex items-center justify-between text-xs mb-1.5">
          <span className="text-zinc-600 dark:text-zinc-400">Vyčerpání ročního ochranného limitu:</span>
          <span className="font-semibold text-zinc-900 dark:text-zinc-100">
            {totalCountedToLimit.toFixed(0)} / {activeTier.limitKč} Kč ({percentOfLimit} %)
          </span>
        </div>
        <div className="h-2 w-full bg-zinc-100 dark:bg-zinc-800 rounded-full overflow-hidden">
          <div
            className={`h-full transition-all duration-300 ${
              limitReached ? "bg-emerald-500" : "bg-blue-600 dark:bg-blue-500"
            }`}
            style={{ width: `${percentOfLimit}%` }}
          />
        </div>
      </div>

      {/* Results grid */}
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 pt-2">
        <div className="p-3 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
          <div className="text-[11px] text-zinc-500 dark:text-zinc-400">Zaplaceno v lékárně</div>
          <div className="text-base font-bold tracking-tight text-zinc-900 dark:text-zinc-100 mt-0.5">
            {totalPatientPaid.toFixed(0)} Kč
          </div>
        </div>

        <div className="p-3 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
          <div className="text-[11px] text-zinc-500 dark:text-zinc-400">Započteno do limitu</div>
          <div className="text-base font-bold tracking-tight text-zinc-900 dark:text-zinc-100 mt-0.5">
            {totalCountedToLimit.toFixed(0)} Kč
          </div>
        </div>

        <div className="p-3 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
          <div className="text-[11px] text-zinc-500 dark:text-zinc-400">Pojišťovna vrátí</div>
          <div className={`text-base font-bold tracking-tight mt-0.5 ${refundedByInsurance > 0 ? "text-emerald-600 dark:text-emerald-400" : "text-zinc-400"}`}>
            {refundedByInsurance.toFixed(0)} Kč
          </div>
        </div>

        <div className="p-3 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
          <div className="text-[11px] text-zinc-500 dark:text-zinc-400">Čistý náklad pacienta</div>
          <div className="text-base font-bold tracking-tight text-zinc-900 dark:text-zinc-100 mt-0.5">
            {netPatientExpense.toFixed(0)} Kč
          </div>
        </div>
      </div>
    </div>
  );
}
