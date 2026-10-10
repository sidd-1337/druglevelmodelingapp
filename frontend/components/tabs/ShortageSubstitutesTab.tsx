import React from "react";
import Link from "next/link";
import { SuklPackageDetail } from "@/lib/types";

interface ShortageSubstitutesTabProps {
  suklDetail: SuklPackageDetail;
}

export default function ShortageSubstitutesTab({ suklDetail }: ShortageSubstitutesTabProps) {
  const isInterrupted = suklDetail.dostupnostStatus === "interrupted";
  const isEnded = suklDetail.dostupnostStatus === "ended";
  const isActive = suklDetail.dostupnostStatus === "active";

  return (
    <div className="space-y-6">
      {/* Shortage Status Banner */}
      <div className="apple-card p-6">
        <div className="flex flex-wrap items-center justify-between gap-4 mb-4">
          <div>
            <h4 className="font-semibold text-base text-zinc-900 dark:text-zinc-100 flex items-center gap-2">
              Hlášení o dostupnosti na trhu ČR
            </h4>
            <p className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5">
              Evidováno dle § 33 zákona č. 378/2007 Sb., o léčivech
            </p>
          </div>

          <span
            className={`inline-flex items-center gap-1.5 px-3 py-1 rounded-full text-xs font-medium ${
              isActive
                ? "bg-emerald-50 text-emerald-700 border border-emerald-200/80 dark:bg-emerald-950/40 dark:text-emerald-300 dark:border-emerald-800/60"
                : isInterrupted
                ? "bg-amber-50 text-amber-800 border border-amber-200/80 dark:bg-amber-950/40 dark:text-amber-300 dark:border-amber-800/60"
                : isEnded
                ? "bg-rose-50 text-rose-700 border border-rose-200/80 dark:bg-rose-950/40 dark:text-rose-300 dark:border-rose-800/60"
                : "bg-zinc-100 text-zinc-700 border border-zinc-200 dark:bg-zinc-800 dark:text-zinc-300 dark:border-zinc-700"
            }`}
          >
            <span
              className={`w-1.5 h-1.5 rounded-full ${
                isActive
                  ? "bg-emerald-500"
                  : isInterrupted
                  ? "bg-amber-500"
                  : isEnded
                  ? "bg-rose-500"
                  : "bg-zinc-400"
              }`}
            />
            {isActive
              ? "Dodávky aktivní"
              : isInterrupted
              ? "Výpadek / Přerušení"
              : isEnded
              ? "Ukončeno"
              : "Bez hlášení"}
          </span>
        </div>

        <p className="text-sm text-zinc-600 dark:text-zinc-300 mb-4">{suklDetail.dostupnostText}</p>

        {isInterrupted && (
          <div className="p-4 rounded-xl bg-amber-50/70 border border-amber-200/80 dark:bg-amber-950/30 dark:border-amber-800/50 text-xs space-y-2">
            {suklDetail.terminObnoveni && (
              <div className="flex flex-wrap items-center gap-2 text-amber-900 dark:text-amber-200">
                <span className="font-medium">Předpokládaný termín obnovení dodávek:</span>
                <span className="font-mono bg-white dark:bg-amber-900/40 px-2 py-0.5 rounded border border-amber-300 dark:border-amber-700 font-semibold">
                  {suklDetail.terminObnoveni}
                </span>
              </div>
            )}
            {suklDetail.duvodPreruseni && (
              <div className="text-zinc-700 dark:text-zinc-300">
                <span className="font-medium text-amber-900 dark:text-amber-200">Důvod přerušení:</span>{" "}
                {suklDetail.duvodPreruseni}
              </div>
            )}
          </div>
        )}
      </div>

      {/* Available Generic Substitutes Table */}
      <div className="apple-card p-6">
        <div className="flex items-center justify-between mb-4">
          <div>
            <h4 className="font-semibold text-base text-zinc-900 dark:text-zinc-100 flex items-center gap-2">
              Dostupné generické náhrady bez výpadku
            </h4>
            <p className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5">
              Přípravky se shodnou nebo terapeuticky ekvivalentní účinnou látkou, které jsou aktuálně dodávány do lékáren.
            </p>
          </div>
        </div>

        {suklDetail.substitutes && suklDetail.substitutes.length > 0 ? (
          <div className="overflow-x-auto rounded-xl border border-zinc-200/80 dark:border-zinc-800">
            <table className="w-full text-left text-xs divide-y divide-zinc-200/80 dark:divide-zinc-800">
              <thead>
                <tr className="text-zinc-500 dark:text-zinc-400 font-medium bg-zinc-50 dark:bg-zinc-800/60">
                  <th className="py-3 px-4">Přípravek</th>
                  <th className="py-3 px-4">Balení & Síla</th>
                  <th className="py-3 px-4">Stav dodávek</th>
                  <th className="py-3 px-4 text-right">Doplatek</th>
                  <th className="py-3 px-4 text-right">Akce</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-zinc-200/60 dark:divide-zinc-800/60">
                {suklDetail.substitutes.map((sub) => (
                  <tr key={sub.kodSukl} className="hover:bg-zinc-50 dark:hover:bg-zinc-800/40 transition-colors">
                    <td className="py-3 px-4 font-semibold text-zinc-900 dark:text-zinc-100">
                      {sub.nazev}
                    </td>
                    <td className="py-3 px-4 text-zinc-500 dark:text-zinc-400">{sub.doplnek}</td>
                    <td className="py-3 px-4">
                      <span className="inline-flex items-center gap-1.5 px-2.5 py-0.5 rounded-full text-[11px] font-medium bg-emerald-50 text-emerald-700 border border-emerald-200/80 dark:bg-emerald-950/40 dark:text-emerald-300 dark:border-emerald-800/60">
                        <span className="w-1.5 h-1.5 rounded-full bg-emerald-500" />
                        Bez výpadku
                      </span>
                    </td>
                    <td className="py-3 px-4 text-right font-mono font-medium text-zinc-900 dark:text-zinc-100">
                      {sub.doplatek.toFixed(2)} Kč
                    </td>
                    <td className="py-3 px-4 text-right">
                      <Link
                        href={`/leky/${sub.kodSukl}`}
                        className="px-2.5 py-1 rounded-lg bg-zinc-100 hover:bg-zinc-200 dark:bg-zinc-800 dark:hover:bg-zinc-700 text-zinc-800 dark:text-zinc-200 font-medium text-xs border border-zinc-200/80 dark:border-zinc-700/80 transition-colors inline-block"
                      >
                        Detail ↗
                      </Link>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        ) : (
          <div className="text-center py-8 text-xs text-zinc-500 dark:text-zinc-400 border border-dashed border-zinc-200 dark:border-zinc-800 rounded-xl">
            Pro tento přípravek nejsou evidovány žádné hlášené výpadky nebo je jedinečným monopolem na trhu ČR.
          </div>
        )}
      </div>
    </div>
  );
}
