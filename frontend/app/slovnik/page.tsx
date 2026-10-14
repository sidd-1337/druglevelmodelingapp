import React from "react";
import { cz } from "@/lib/typo";

interface Term {
  /** Mathematical symbol as [base, subscript?]. */
  symbol?: [string, string?];
  name: string;
  en: string;
  desc: string;
  note?: string;
}

const GROUPS: Array<{ id: string; title: string; terms: Term[] }> = [
  {
    id: "kinetika",
    title: "Kinetika",
    terms: [
      {
        symbol: ["C", "max"],
        name: "Vrcholová koncentrace v krvi",
        en: "Peak concentration",
        desc: "Nejvyšší hladina léčiva v krevní plazmě po jednorázovém nebo opakovaném podání. Pokud Cmax překročí toxický limit MTC, hrozí akutní nežádoucí účinky nebo toxicita.",
        note: "U stimulancií (například methylfenidát) určuje intenzitu žádoucího i vedlejšího účinku.",
      },
      {
        symbol: ["T", "max"],
        name: "Čas do dosažení maxima",
        en: "Time to peak",
        desc: "Doba od podání léčiva do nejvyšší plazmatické koncentrace (Cmax). Závisí na rychlosti absorpce (ka) a cestě podání. Inhalace dosahuje maxima za 1–3 minuty, nosní sliznice za 10–20 minut, perorální tablety obvykle za 1–3 hodiny.",
        note: "Klíčový údaj pro bezpečný rozestup: další dávku neužívejte dříve, než nastane Tmax.",
      },
      {
        symbol: ["t½"],
        name: "Eliminační poločas",
        en: "Half-life",
        desc: "Doba, za kterou organismus (játra a ledviny) vyloučí nebo přemění přesně polovinu léčiva z krevního oběhu. Určuje, jak často je třeba dávkovat a jak dlouho látka tělo ovlivňuje.",
        note: "Po 1× t½ zbývá 50 %, po 2× t½ 25 %, po 3× t½ 12,5 %. Po 4–5× t½ je látka z 95–97 % odbourána.",
      },
    ],
  },
  {
    id: "distribuce",
    title: "Distribuce a eliminace",
    terms: [
      {
        symbol: ["V", "d"],
        name: "Distribuční objem",
        en: "Volume of distribution",
        desc: "Zdánlivý objem tělesné tekutiny, ve kterém by se musela celá dávka rozpustit, aby vznikla naměřená plazmatická koncentrace. Vyjadřuje se v litrech nebo v litrech na kilogram hmotnosti.",
        note: "Vysoký Vd (nad 2–4 l/kg) znamená, že léčivo proniká hluboko do tukové tkáně, svalů a mozku (například propranolol, lisdexamfetamin).",
      },
      {
        symbol: ["CL"],
        name: "Celková clearance",
        en: "Clearance",
        desc: "Objem krve, který organismus za jednotku času (obvykle l/h) úplně a nevratně očistí od léčiva. Rovná se součinu eliminační konstanty a distribučního objemu: CL = ke × Vd.",
        note: "Snížená funkce ledvin nebo jater clearance výrazně snižuje a vede k nebezpečnému hromadění léčiva.",
      },
    ],
  },
  {
    id: "akumulace",
    title: "Akumulace a terapeutické okno",
    terms: [
      {
        symbol: ["C", "ss"],
        name: "Ustálený stav",
        en: "Steady state",
        desc: "Rovnováha při pravidelném opakovaném dávkování, kdy se rychlost přívodu látky vyrovná rychlosti jejího odbourávání. Nastává obvykle po 4 až 5 poločasech.",
        note: "V ustáleném stavu hladina pravidelně kolísá mezi vrcholem (Css,max) a údolím (Css,min těsně před další dávkou).",
      },
      {
        symbol: ["R"],
        name: "Index akumulace",
        en: "Accumulation ratio",
        desc: "Poměr hladiny v ustáleném stavu k hladině po první dávce: R = 1 / (1 − exp(−ke × τ)). Je-li interval dávkování τ kratší než poločas, léčivo se v těle výrazně hromadí (R > 1).",
        note: "U látek s dlouhým poločasem (Rivotril 35 h, diazepam 40 h) je hladina v ustáleném stavu několikanásobně vyšší než po první tabletě.",
      },
      {
        name: "Terapeutické okno (MEC a MTC)",
        en: "Therapeutic window",
        desc: "Rozmezí mezi minimální efektivní koncentrací (MEC), pod kterou léčivo nezabírá, a maximální tolerovanou koncentrací (MTC), nad kterou hrozí toxické a nežádoucí účinky.",
        note: "Mimořádně úzké okno má například lithium (0,6–1,2 mmol/l). Toxicita nastává už od 1,5 mmol/l.",
      },
    ],
  },
  {
    id: "harm-reduction",
    title: "Harm reduction",
    terms: [
      {
        name: "Vrstvení dávek",
        en: "Dose stacking",
        desc: "Další dávka užitá dříve, než se předchozí stihla vstřebat nebo odbourat. Vrcholy se nečekaně sečtou, což je nejčastější příčina akutního předávkování psychoaktivními látkami.",
        note: "Pokud látka po 45 minutách „nezabírá“, nebere se další dávka. Počkejte alespoň 2–3 hodiny na plný nástup a zhodnoťte účinek.",
      },
      {
        name: "Očišťovací interval",
        en: "Washout period",
        desc: "Čas potřebný k téměř úplnému vymizení léčiva z těla před nasazením jiné látky nebo dalším cyklem, obvykle 5 poločasů. Brání nebezpečným interakcím, například serotoninovému syndromu při kombinaci MAO inhibitorů se stimulancii nebo SSRI.",
        note: "U fluoxetinu s poločasem metabolitu až 16 dnů musí očišťovací interval trvat až 5 týdnů.",
      },
    ],
  },
];

function Sym({ symbol }: { symbol: [string, string?] }) {
  return (
    <>
      {symbol[0]}
      {symbol[1] && <sub className="text-[0.6em] font-medium">{symbol[1]}</sub>}
    </>
  );
}

export default function GlossaryPage() {
  return (
    <div className="page pb-8 pt-14 sm:pt-20">
      <h1 className="text-title">Slovník</h1>
      <p className="mt-4 max-w-prose text-lead text-ink-2">
        Veličiny, pojmy a bezpečnostní pravidla, na které narazíte u každého léku v portálu.
      </p>

      <div className="mt-16 space-y-20 sm:mt-20">
        {GROUPS.map((group) => (
          <section key={group.id} aria-labelledby={group.id}>
            <h2 id={group.id} className="text-sub">
              {group.title}
            </h2>
            <dl className="mt-6 border-b border-line-soft">
              {group.terms.map((term) => (
                <div
                  key={term.en}
                  className="grid gap-3 border-t border-line-soft py-8 sm:grid-cols-[15rem_minmax(0,1fr)] sm:gap-12"
                >
                  <dt>
                    {term.symbol ? (
                      <>
                        <span className="figures block text-sub">
                          <Sym symbol={term.symbol} />
                        </span>
                        <span className="mt-1 block text-ui">{term.name}</span>
                      </>
                    ) : (
                      <span className="block text-h3">{term.name}</span>
                    )}
                    <span className="mt-0.5 block text-small text-ink-2">{term.en}</span>
                  </dt>
                  <dd className="max-w-prose">
                    <p className="text-body">{cz(term.desc)}</p>
                    {term.note && (
                      <p className="mt-3 text-ui text-ink-2">
                        <span className="font-semibold text-ink">V praxi: </span>
                        {cz(term.note)}
                      </p>
                    )}
                  </dd>
                </div>
              ))}
            </dl>
          </section>
        ))}
      </div>
    </div>
  );
}
