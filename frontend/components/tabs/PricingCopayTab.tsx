import React from "react";
import { SuklPackageDetail } from "@/lib/types";
import CopayCalculator from "../CopayCalculator";
import { Stethoscope } from "lucide-react";

interface PricingCopayTabProps {
  suklDetail: SuklPackageDetail;
}

export default function PricingCopayTab({ suklDetail }: PricingCopayTabProps) {
  const uhrada = suklDetail.uhradaZp ?? 0;
  const maxCena = suklDetail.maxCenaLek ?? 0;
  const maxDoplatek = suklDetail.maxDoplatek ?? 0;
  const zapocitatelny = suklDetail.zapocitatelnyDoplatek ?? 0;

  return (
    <div className="space-y-6">
      {/* 4 Pricing KPI Cards */}
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-3.5">
        <div className="apple-card p-5">
          <div className="text-xs text-zinc-500 dark:text-zinc-400 font-normal">Úhrada pojišťovny (SCAU)</div>
          <div className="text-2xl font-semibold tracking-tight text-emerald-600 dark:text-emerald-400 mt-1">
            {uhrada > 0 ? `${uhrada.toFixed(2)} Kč` : "0.00 Kč"}
          </div>
          <div className="text-[11px] text-zinc-400 dark:text-zinc-500 mt-1">
            {suklDetail.uhradaRezim === "full"
              ? "Plně hrazeno ze ZP"
              : suklDetail.uhradaRezim === "partial"
              ? "Částečně hrazeno"
              : "Nehrazeno (samoplátce)"}
          </div>
        </div>

        <div className="apple-card p-5">
          <div className="text-xs text-zinc-500 dark:text-zinc-400 font-normal">Maximální cena v lékárně</div>
          <div className="text-2xl font-semibold tracking-tight text-zinc-950 dark:text-zinc-50 mt-1">
            {maxCena > 0 ? `${maxCena.toFixed(2)} Kč` : "Není stanovena"}
          </div>
          <div className="text-[11px] text-zinc-400 dark:text-zinc-500 mt-1">Úředně stanovený strop SÚKL</div>
        </div>

        <div className="apple-card p-5">
          <div className="text-xs text-zinc-500 dark:text-zinc-400 font-normal">Odhadovaný doplatek pacienta</div>
          <div className="text-2xl font-semibold tracking-tight text-amber-600 dark:text-amber-400 mt-1">
            {maxDoplatek > 0 ? `${maxDoplatek.toFixed(2)} Kč` : "0.00 Kč"}
          </div>
          <div className="text-[11px] text-zinc-400 dark:text-zinc-500 mt-1">Orientační doplatek za 1 balení</div>
        </div>

        <div className="apple-card p-5">
          <div className="text-xs text-zinc-500 dark:text-zinc-400 font-normal">Započitatelný doplatek do limitu</div>
          <div className="text-2xl font-semibold tracking-tight text-blue-600 dark:text-blue-400 mt-1">
            {zapocitatelny > 0 ? `${zapocitatelny.toFixed(2)} Kč` : "0.00 Kč"}
          </div>
          <div className="text-[11px] text-zinc-400 dark:text-zinc-500 mt-1">Částka počítaná do ročního limitu</div>
        </div>
      </div>

      {/* Interactive Annual Copay Protection Limit Calculator */}
      <CopayCalculator
        maxDoplatek={maxDoplatek}
        zapocitatelnyDoplatek={zapocitatelny}
      />

      {/* Prescriptions & Conditions */}
      {(suklDetail.preskripcniOmezeni || suklDetail.indikacniOmezeni) && (
        <div className="apple-card p-6 space-y-4">
          <div className="flex items-center gap-2">
            <Stethoscope className="w-4 h-4 text-blue-600 dark:text-blue-400" />
            <h4 className="font-semibold text-sm text-zinc-950 dark:text-zinc-50">
              Podmínky úhrady a preskripce ze zdravotního pojištění
            </h4>
          </div>

          {suklDetail.preskripcniOmezeni && (
            <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 text-xs">
              <span className="font-medium text-zinc-800 dark:text-zinc-200 block mb-1">
                Preskripční omezení (oprávněná odbornost lékaře):
              </span>
              <span className="text-zinc-500 dark:text-zinc-400 font-mono">{suklDetail.preskripcniOmezeni}</span>
            </div>
          )}

          {suklDetail.indikacniOmezeni && (
            <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 text-xs leading-relaxed">
              <span className="font-medium text-zinc-800 dark:text-zinc-200 block mb-1">
                Indikační omezení úhrady:
              </span>
              <p className="text-zinc-600 dark:text-zinc-300 whitespace-pre-line font-normal">{suklDetail.indikacniOmezeni}</p>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
