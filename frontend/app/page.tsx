import React from "react";
import Link from "next/link";
import HeroCurve from "@/components/HeroCurve";
import SuklSearch from "@/components/SuklSearch";
import { cz } from "@/lib/typo";

const FEATURED = [
  { name: "Ritalin", code: "0268766", note: "10 mg tablety" },
  { name: "Livizux", code: "0265888", note: "70 mg tobolky" },
  { name: "Rivotril", code: "0014957", note: "0,5 mg tablety" },
  { name: "Pragiola", code: "0211488", note: "75 mg tobolky" },
  { name: "Ibalgin", code: "0254292", note: "200 mg tablety" },
  { name: "Propranolol", code: "0284404", note: "80 mg tobolky" },
];

const PARTS = [
  {
    title: "Přehled a dokumentace",
    text: "Účinná látka, léková forma, síla, držitel registrace a způsob výdeje. Odkazy na příbalovou informaci (PIL) a souhrn údajů o přípravku (SPC) ve formátu PDF.",
  },
  {
    title: "Cena a úhrada",
    text: "Úhrada pojišťovny, nejvyšší cena v lékárně a odhad doplatku za balení. Kalkulačka ukáže, kdy dosáhnete ročního ochranného limitu: 5 000 Kč pro dospělé, 1 000 Kč pro děti a seniory 65–69 let, 500 Kč od 70 let.",
  },
  {
    title: "Dostupnost a náhrady",
    text: "Hlášená přerušení a ukončení dodávek podle § 33 zákona o léčivech, předpokládaný termín obnovení a seznam dostupných přípravků se stejnou látkou.",
  },
  {
    title: "Hladina v těle",
    text: "Jednokompartmentový model ukáže, jak rychle lék nastupuje, kdy je v krvi nejvíc a jak dlouho působí. Čtyři cesty podání, bezpečný odstup mezi dávkami a varování před vrstvením dávek.",
  },
];

export default function HomePage() {
  return (
    <>
      <section className="overflow-hidden pt-14 sm:pt-20">
        <div className="page text-center">
          <h1 className="mx-auto max-w-[16ch] text-display">Co se s lékem děje v těle.</h1>
          <p className="mx-auto mt-6 max-w-[38rem] text-lead text-ink-2">
            {cz(
              "Vyhledejte lék v registru SÚKL. Uvidíte jeho cenu a úhradu, dostupnost v lékárnách a modelovanou hladinu v krvi."
            )}
          </p>

          <div className="mt-10">
            <SuklSearch />
          </div>

          <p className="mt-5 text-ui text-ink-2">
            Zkuste{" "}
            {FEATURED.map((d, i) => (
              <React.Fragment key={d.code}>
                <Link href={`/leky/${d.code}`} title={d.note} className="link">
                  {d.name}
                </Link>
                {i < FEATURED.length - 2 ? ", " : i === FEATURED.length - 2 ? " nebo " : "."}
              </React.Fragment>
            ))}
          </p>
        </div>

        <div className="mt-12 sm:mt-14">
          <HeroCurve />
        </div>
      </section>

      <section className="mt-20 border-t border-line-soft sm:mt-28">
        <div className="page grid gap-10 py-20 sm:py-28 lg:grid-cols-[5fr_7fr] lg:gap-20">
          <h2 className="text-headline lg:sticky lg:top-24 lg:self-start">
            Všechno o léku na jedné stránce.
          </h2>
          <dl className="border-b border-line-soft">
            {PARTS.map((part) => (
              <div key={part.title} className="border-t border-line-soft py-7">
                <dt className="text-h3">{part.title}</dt>
                <dd className="mt-2 max-w-prose text-body text-ink-2">{cz(part.text)}</dd>
              </div>
            ))}
          </dl>
        </div>
      </section>
    </>
  );
}
