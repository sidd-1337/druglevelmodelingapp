import React from "react";
import Link from "next/link";
import { SuklPackageDetail } from "@/lib/types";
import { kc, packageText, sentenceCase } from "@/lib/format";
import SectionHead from "../ui/SectionHead";
import SpecList from "../ui/SpecList";
import Status, { Tone } from "../ui/Status";

interface ShortageSubstitutesTabProps {
  suklDetail: SuklPackageDetail;
}

export default function ShortageSubstitutesTab({ suklDetail }: ShortageSubstitutesTabProps) {
  const status = suklDetail.dostupnostStatus;
  const isInterrupted = status === "interrupted";

  const tone: Tone =
    status === "active" ? "ok" : status === "interrupted" ? "warn" : status === "ended" ? "bad" : "neutral";
  const label =
    status === "active"
      ? "Dodávky jsou aktivní"
      : status === "interrupted"
      ? "Dodávky jsou přerušené"
      : status === "ended"
      ? "Dodávky byly ukončeny"
      : "Bez hlášení o dostupnosti";

  const details = [
    suklDetail.terminObnoveni && {
      label: "Předpokládané obnovení dodávek",
      value: suklDetail.terminObnoveni,
    },
    suklDetail.duvodPreruseni && { label: "Důvod přerušení", value: suklDetail.duvodPreruseni },
  ].filter(Boolean) as Array<{ label: string; value: string }>;

  const substitutes = suklDetail.substitutes ?? [];

  return (
    <div className="space-y-20">
      <section>
        <SectionHead title="Dostupnost na trhu v ČR">
          Evidováno podle § 33 zákona č. 378/2007 Sb., o léčivech.
        </SectionHead>
        <p className="mt-8 text-sub">
          <Status tone={tone}>{label}</Status>
        </p>
        <p className="mt-2 max-w-prose text-body text-ink-2">{suklDetail.dostupnostText}</p>
        {isInterrupted && details.length > 0 && (
          <div className="mt-8">
            <SpecList items={details} />
          </div>
        )}
      </section>

      <section>
        <SectionHead title="Dostupné náhrady">
          Přípravky se stejnou nebo terapeuticky rovnocennou látkou, které se teď dodávají do lékáren.
        </SectionHead>

        {substitutes.length > 0 ? (
          <div className="-mx-5 mt-6 overflow-x-auto px-5 sm:mx-0 sm:px-0">
            <table className="w-full min-w-[480px] border-collapse text-left">
              <thead>
                <tr className="text-ui text-ink-2">
                  <th scope="col" className="border-b border-line py-3 pr-6 font-normal">
                    Přípravek
                  </th>
                  <th scope="col" className="border-b border-line py-3 pr-6 font-normal">
                    Balení a síla
                  </th>
                  <th scope="col" className="border-b border-line py-3 text-right font-normal">
                    Doplatek
                  </th>
                </tr>
              </thead>
              <tbody>
                {substitutes.map((sub) => (
                  <tr key={sub.kodSukl}>
                    <td className="border-b border-line-soft py-4 pr-6 text-body">
                      <Link href={`/leky/${sub.kodSukl}`} className="link">
                        {sentenceCase(sub.nazev)}
                      </Link>
                    </td>
                    <td className="border-b border-line-soft py-4 pr-6 text-ui text-ink-2">{packageText(sub.doplnek)}</td>
                    <td className="figures border-b border-line-soft py-4 text-right text-body">
                      {kc(sub.doplatek, 2)}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        ) : (
          <p className="mt-6 max-w-prose text-body text-ink-2">
            Pro tento přípravek nejsou evidovány žádné náhrady. Buď nebyl hlášen výpadek, nebo jde o jediný
            přípravek s touto látkou na českém trhu.
          </p>
        )}
      </section>
    </div>
  );
}
