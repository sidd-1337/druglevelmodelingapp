import React from "react";
import Link from "next/link";

export default function GlossaryPage() {
  const GLOSSARY_ITEMS = [
    {
      term: "C_max – Vrcholová koncentrace v krvi",
      category: "Kinetika",
      badge: "Peak concentration",
      desc: "Nejvyšší dosažená hladina léčiva v krevní plazmě po jednorázovém nebo opakovaném podání dávky. Pokud Cmax překročí toxický limit MTC, hrozí akutní nežádoucí účinky nebo toxicita.",
      clinicalNote: "U stimulancií (např. methylfenidát) determinuje intenzitu žádoucího i vedlejšího účinku.",
    },
    {
      term: "T_max – Čas do dosažení maxima",
      category: "Kinetika",
      badge: "Time to peak",
      desc: "Doba, která uplyne od podání léčiva do dosažení nejvyšší plazmatické koncentrace (Cmax). Odvíjí se od rychlosti absorpce (ka) a cesty podání. Inhalace dosahuje maxima za 1–3 minuty, nosní sliznice za 10–20 minut, perorální tablety obvykle za 1–3 hodiny.",
      clinicalNote: "Klíčový údaj pro bezpečný rozestup: nikdy neužívejte další dávku dříve, než nastane Tmax!",
    },
    {
      term: "t½ – Eliminační poločas rozpadu",
      category: "Kinetika",
      badge: "Half-life",
      desc: "Doba nutná k tomu, aby organismus (játra a ledviny) vyloučil nebo metabolizoval přesně 50 % léčiva z krevního oběhu. Přímo určuje potřebnou frekvenci dávkování a dobu, po kterou látka ovlivňuje tělo.",
      clinicalNote: "Po 1× t½ zbývá 50 %, po 2× t½ 25 %, po 3× t½ 12.5 %, po 4–5× t½ je látka z 95–97 % odbourána.",
    },
    {
      term: "V_d – Distribuční objem",
      category: "Distribuce",
      badge: "Volume of distribution",
      desc: "Zdánlivý teoretický objem tělesné tekutiny, ve kterém by se celková dávka musela rozpustit, aby vytvořila naměřenou plazmatickou koncentraci. Vyjadřuje se v L nebo L/kg tělesné hmotnosti.",
      clinicalNote: "Vysoký Vd (> 2–4 L/kg) značí, že léčivo proniká hluboko do tukové tkáně, svalů a mozku (např. propranolol, lisdexamfetamin).",
    },
    {
      term: "CL – Celková clearance",
      category: "Eliminace",
      badge: "Clearance",
      desc: "Objem krve, který organismus kompletně a nevratně očistí od léčiva za jednotku času (obvykle L/h). Matematicky se rovná součinu eliminační konstanty a distribučního objemu: CL = ke × Vd.",
      clinicalNote: "Snížená funkce ledvin nebo jater clearance dramaticky snižuje a vede k nebezpečnému hromadění.",
    },
    {
      term: "C_ss – Ustálený stav (Steady-State)",
      category: "Akumulace",
      badge: "Steady state",
      desc: "Rovnovážný stav dosažený při pravidelném opakovaném dávkování, kdy se rychlost přívodu látky přesně vyrovná rychlosti její eliminace. Nastává typicky po 4 až 5 biologických poločasech (4–5 × t½).",
      clinicalNote: "V ustáleném stavu koncentrace pravidelně osciluje mezi vrcholem (Css,max) a údolím (Css,min těsně před další dávkou).",
    },
    {
      term: "R – Index akumulace",
      category: "Akumulace",
      badge: "Accumulation ratio",
      desc: "Poměr mezi koncentrací v ustáleném stavu a koncentrací po první samostatné dávce: R = 1 / (1 - exp(-ke × τ)). Pokud je interval podávání (τ) kratší než poločas t½, lék se v těle výrazně hromadí (R > 1).",
      clinicalNote: "U látek s dlouhým poločasem (Rivotril t½ = 35 h, Diazepam t½ = 40 h) je hladina v ustáleném stavu několikanásobně vyšší než po první tabletě.",
    },
    {
      term: "Terapeutické okno (MEC & MTC)",
      category: "TDM",
      badge: "Therapeutic window",
      desc: "Rozmezí mezi Minimální efektivní koncentrací (MEC – pod kterou léčivo nezabírá) a Maximální tolerovanou koncentrací (MTC – nad kterou hrozí toxické a nežádoucí účinky).",
      clinicalNote: "Extrémně úzké okno má např. Lithium (0.6–1.2 mmol/L; toxicita nastává již od 1.5 mmol/L).",
    },
    {
      term: "Dose Stacking (Vrstvení dávek)",
      category: "Harm Reduction",
      badge: "Kumulace",
      desc: "Situace, kdy uživatel vezme další dávku dříve, než se předchozí dávka stihla vstřebat nebo odbourat. Dochází k nečekanému sčítání vrcholů, což je nejčastější příčina akutního předávkování u psychoaktivních látek.",
      clinicalNote: "Pravidlo harm reduction: pokud látka po 45 minutách 'nezabírá', neberte další dávku – počkejte alespoň 2–3 hodiny na plný nástup a zhodnocení.",
    },
    {
      term: "Washout Period (Očišťovací interval)",
      category: "Harm Reduction",
      badge: "Vyloučení",
      desc: "Čas potřebný k téměř úplnému vymizení léčiva z těla před nasazením jiné látky nebo dalšího cyklu (obvykle 5 biologických poločasů). Zabraňuje fatálním lékovým interakcím (např. serotoninový syndrom při kombinaci MAOI a stimulancií/SSRI).",
      clinicalNote: "Např. u fluoxetinu s poločasem metabolitu až 16 dnů musí washout trvat až 5 týdnů.",
    },
  ];

  return (
    <div className="space-y-8 max-w-4xl mx-auto py-2">
      {/* Header */}
      <div>
        <Link
          href="/"
          className="text-xs text-zinc-500 hover:text-zinc-900 dark:hover:text-zinc-100 transition-colors mb-4 inline-flex items-center gap-1 font-medium"
        >
          ← Zpět na přehled
        </Link>
        <h1 className="text-3xl font-bold tracking-tight text-zinc-900 dark:text-zinc-100">
          Farmakokinetický & Harm Reduction Slovník
        </h1>
        <p className="text-sm text-zinc-500 dark:text-zinc-400 mt-2 leading-relaxed">
          Srozumitelný průvodce klíčovými proměnnými, matematickými koncepty a bezpečnostními pravidly pro pacienty i odbornou veřejnost.
        </p>
      </div>

      {/* Glossary Items List */}
      <div className="space-y-4">
        {GLOSSARY_ITEMS.map((item, idx) => (
          <div key={idx} className="apple-card p-6 space-y-3">
            <div className="flex flex-wrap items-center justify-between gap-2 border-b border-zinc-200/80 dark:border-zinc-800 pb-3">
              <h3 className="font-semibold text-base text-zinc-900 dark:text-zinc-100">{item.term}</h3>
              <div className="flex items-center gap-2">
                <span className="text-[11px] font-mono px-2 py-0.5 rounded-md bg-zinc-100 dark:bg-zinc-800 text-zinc-600 dark:text-zinc-400 border border-zinc-200/80 dark:border-zinc-700/80">
                  {item.badge}
                </span>
                <span className="text-xs px-2.5 py-0.5 rounded-full bg-blue-50 text-blue-700 border border-blue-200/80 dark:bg-blue-950/40 dark:text-blue-300 dark:border-blue-800/60 font-medium">
                  {item.category}
                </span>
              </div>
            </div>

            <p className="text-sm text-zinc-600 dark:text-zinc-300 leading-relaxed pt-0.5">
              {item.desc}
            </p>

            {item.clinicalNote && (
              <div className="p-3.5 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 text-xs text-zinc-600 dark:text-zinc-400 leading-relaxed">
                <strong className="font-semibold text-zinc-900 dark:text-zinc-200">Klinický kontext:</strong>{" "}
                {item.clinicalNote}
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
}
