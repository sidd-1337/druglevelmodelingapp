import React from "react";
import { DrugPreset, SuklPackageDetail } from "@/lib/types";
import { packageText, substanceCase, upperFirst } from "@/lib/format";
import SectionHead from "../ui/SectionHead";
import SpecList from "../ui/SpecList";

interface AboutDrugTabProps {
  drug: DrugPreset | null;
  suklDetail: SuklPackageDetail;
}

export default function AboutDrugTab({ drug, suklDetail }: AboutDrugTabProps) {
  const substance = substanceCase(suklDetail.lecivaLatka || drug?.substance || suklDetail.atcNazev);
  const atc = [
    suklDetail.atcNazev && suklDetail.atcNazev !== (suklDetail.lecivaLatka || "")
      ? substanceCase(suklDetail.atcNazev)
      : "",
    suklDetail.atcKod ? `ATC ${suklDetail.atcKod}` : "",
  ]
    .filter(Boolean)
    .join(", ");

  const documents = [
    suklDetail.spcUrl && {
      title: "Souhrn údajů o přípravku (SPC)",
      text: "Odborná informace pro lékaře a lékárníky",
      action: "Otevřít PDF",
      href: suklDetail.spcUrl,
    },
    suklDetail.pilUrl && {
      title: "Příbalová informace (PIL)",
      text: "Srozumitelný návod pro pacienty",
      action: "Otevřít PDF",
      href: suklDetail.pilUrl,
    },
    {
      title: "Karta v databázi SÚKL",
      text: "prehledy.sukl.cz",
      action: "Otevřít kartu",
      href: suklDetail.suklUrl,
    },
  ].filter(Boolean) as Array<{ title: string; text: string; action: string; href: string }>;

  return (
    <div className="space-y-20">
      <section>
        <SectionHead title="Údaje o přípravku" />
        <div className="mt-6">
          <SpecList
            items={[
              { label: "Účinná látka", value: substance, note: atc || undefined },
              { label: "Léková forma a síla", value: packageText(suklDetail.sila), note: suklDetail.forma },
              { label: "Cesta podání", value: suklDetail.cesta },
              { label: "Způsob výdeje", value: upperFirst(suklDetail.vydejText) },
              { label: "Registrační číslo", value: suklDetail.registracniCislo },
              {
                label: "Držitel rozhodnutí o registraci",
                value: suklDetail.drzitel || "Neuvedeno v databázi",
              },
            ]}
          />
        </div>
      </section>

      <section>
        <SectionHead title="Oficiální dokumentace">
          Dokumenty schválené Státním ústavem pro kontrolu léčiv.
        </SectionHead>
        <ul className="mt-6 border-b border-line-soft">
          {documents.map((doc) => (
            <li
              key={doc.title}
              className="flex flex-wrap items-baseline justify-between gap-x-8 gap-y-1 border-t border-line-soft py-4"
            >
              <div>
                <p className="text-body">{doc.title}</p>
                <p className="text-ui text-ink-2">{doc.text}</p>
              </div>
              <a href={doc.href} target="_blank" rel="noopener noreferrer" className="link text-body">
                {doc.action}
                <span className="sr-only"> (otevře se v novém okně)</span>
              </a>
            </li>
          ))}
        </ul>
      </section>

      {drug && drug.notes && (
        <section>
          <SectionHead title="Farmakologický souhrn" />
          <p className="mt-4 max-w-prose text-body text-ink-2">{drug.notes}</p>
        </section>
      )}
    </div>
  );
}
