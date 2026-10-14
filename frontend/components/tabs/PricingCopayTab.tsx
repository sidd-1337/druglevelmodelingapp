import React from "react";
import { SuklPackageDetail } from "@/lib/types";
import { kc } from "@/lib/format";
import CopayCalculator from "../CopayCalculator";
import SectionHead from "../ui/SectionHead";

interface PricingCopayTabProps {
  suklDetail: SuklPackageDetail;
}

export default function PricingCopayTab({ suklDetail }: PricingCopayTabProps) {
  const uhrada = suklDetail.uhradaZp ?? 0;
  const maxCena = suklDetail.maxCenaLek ?? 0;
  const maxDoplatek = suklDetail.maxDoplatek ?? 0;
  const zapocitatelny = suklDetail.zapocitatelnyDoplatek ?? 0;

  const rows = [
    {
      label: "Nejvyšší cena v lékárně",
      hint: "Úředně stanovený strop SÚKL",
      value: maxCena > 0 ? kc(maxCena, 2) : "Není stanovena",
      lead: false,
    },
    {
      label: "Úhrada pojišťovny",
      hint:
        suklDetail.uhradaRezim === "full"
          ? "Plně hrazeno ze zdravotního pojištění"
          : suklDetail.uhradaRezim === "partial"
          ? "Částečně hrazeno"
          : "Nehrazeno, platí samoplátce",
      value: kc(uhrada, 2),
      lead: false,
    },
    {
      label: "Odhadovaný doplatek",
      hint: "Orientačně za jedno balení",
      value: kc(maxDoplatek, 2),
      lead: true,
    },
    {
      label: "Započitatelný doplatek",
      hint: "Částka, která se počítá do ročního limitu",
      value: kc(zapocitatelny, 2),
      lead: false,
    },
  ];

  return (
    <div className="space-y-20">
      <section>
        <SectionHead title="Cena a úhrada">Částky za jedno balení.</SectionHead>
        <dl className="mt-6 border-b border-line-soft">
          {rows.map((row) => (
            <div
              key={row.label}
              className="grid gap-1 border-t border-line-soft py-5 sm:grid-cols-[1fr_auto] sm:items-baseline sm:gap-8"
            >
              <dt>
                <span className="block text-body">{row.label}</span>
                <span className="block text-ui text-ink-2">{row.hint}</span>
              </dt>
              <dd className={`figures sm:text-right ${row.lead ? "text-figure" : "text-sub"}`}>
                {row.value}
              </dd>
            </div>
          ))}
        </dl>
      </section>

      <CopayCalculator maxDoplatek={maxDoplatek} zapocitatelnyDoplatek={zapocitatelny} />

      {(suklDetail.preskripcniOmezeni || suklDetail.indikacniOmezeni) && (
        <section>
          <SectionHead title="Podmínky úhrady a preskripce" />
          <div className="mt-6 space-y-8">
            {suklDetail.preskripcniOmezeni && (
              <div>
                <h3 className="text-h3">Preskripční omezení</h3>
                <p className="mt-1 text-ui text-ink-2">Odbornost lékaře, který smí lék předepsat.</p>
                <p className="mt-3 max-w-prose text-body">{suklDetail.preskripcniOmezeni}</p>
              </div>
            )}
            {suklDetail.indikacniOmezeni && (
              <div>
                <h3 className="text-h3">Indikační omezení úhrady</h3>
                <p className="mt-3 max-w-prose whitespace-pre-line text-body text-ink-2">
                  {suklDetail.indikacniOmezeni}
                </p>
              </div>
            )}
          </div>
        </section>
      )}
    </div>
  );
}
