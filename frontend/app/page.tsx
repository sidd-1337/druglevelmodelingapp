import React from "react";
import Link from "next/link";
import SuklSearch from "@/components/SuklSearch";
import { CLINICAL_DRUGS } from "@/lib/sukl-data";
import { FileText, ShieldAlert, Activity, CreditCard } from "lucide-react";

export default function HomePage() {
  const featuredDrugs = [
    { id: "ritalin", name: "Ritalin", code: "0268766", subtitle: "10 mg tablety" },
    { id: "livizux", name: "Livizux", code: "0265888", subtitle: "70 mg tobolky" },
    { id: "rivotril", name: "Rivotril", code: "0014957", subtitle: "0.5 mg tablety" },
    { id: "pregabalin", name: "Pragiola", code: "0211488", subtitle: "75 mg tobolky" },
    { id: "ibuprofen", name: "Ibalgin", code: "0254292", subtitle: "200 mg tablety" },
    { id: "propranolol", name: "Propranolol", code: "0284404", subtitle: "80 mg tobolky" },
  ];

  return (
    <div className="space-y-16 sm:space-y-20">
      {/* Apple Style Hero */}
      <section className="text-center max-w-3xl mx-auto pt-6 sm:pt-10">
        <h1 className="text-4xl sm:text-5xl lg:text-6xl font-semibold tracking-tight text-zinc-950 dark:text-zinc-50 leading-tight">
          Lékový portál ČR
        </h1>
        <p className="mt-4 text-zinc-500 dark:text-zinc-400 text-base sm:text-lg max-w-xl mx-auto font-normal">
          Otevřená data SÚKL, úhrady pojišťoven, monitoring výpadků a klinická farmakokinetika s harm reduction.
        </p>

        {/* Spotlight Search Bar */}
        <div className="mt-8">
          <SuklSearch />
        </div>

        {/* Quick jump pills */}
        <div className="mt-6 flex flex-wrap items-center justify-center gap-2 text-xs">
          <span className="text-zinc-400 dark:text-zinc-500 font-medium">Častá léčiva:</span>
          {featuredDrugs.map((d) => (
            <Link
              key={d.id}
              href={`/leky/${d.code}`}
              className="px-3 py-1.5 rounded-full bg-zinc-100 hover:bg-zinc-200/80 dark:bg-zinc-800 dark:hover:bg-zinc-700 text-zinc-800 dark:text-zinc-200 border border-zinc-200/60 dark:border-zinc-700 transition-all font-medium inline-flex items-center gap-1.5"
            >
              <span>{d.name}</span>
              <span className="text-zinc-400 dark:text-zinc-500 text-[11px]">({d.subtitle})</span>
            </Link>
          ))}
        </div>
      </section>

      {/* Clean Feature Overview */}
      <section className="space-y-6">
        <div className="flex items-center justify-between">
          <h2 className="text-lg font-semibold tracking-tight text-zinc-900 dark:text-zinc-100">
            Klíčové moduly portálu
          </h2>
          <span className="text-xs text-zinc-400 dark:text-zinc-500">
            Pro pacienty, lékaře i studenty
          </span>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          {/* Card 1: Registr SÚKL */}
          <div className="apple-card p-6 flex flex-col justify-between">
            <div>
              <div className="w-10 h-10 rounded-xl bg-zinc-100 dark:bg-zinc-800 text-zinc-900 dark:text-zinc-100 flex items-center justify-center mb-4">
                <FileText className="w-5 h-5 text-blue-600 dark:text-blue-400" />
              </div>
              <h3 className="text-base font-semibold text-zinc-900 dark:text-zinc-100 mb-1">
                Registr léčiv & Dokumentace
              </h3>
              <p className="text-sm text-zinc-500 dark:text-zinc-400 leading-relaxed font-normal">
                Přesná data z databáze DLP SÚKL ČR. Léková forma, síla, držitel registrace,
                výdejní režim (volně prodejné, na recept, návykové látky) a přímé odkazy na oficiální PDF příbalového letáku (PIL) i SPC pro lékaře.
              </p>
            </div>
            <div className="mt-5 pt-4 border-t border-zinc-100 dark:border-zinc-800/80 text-xs text-zinc-400">
              Data synchronizována s přehledy SÚKL
            </div>
          </div>

          {/* Card 2: Ceny a úhrady */}
          <div className="apple-card p-6 flex flex-col justify-between">
            <div>
              <div className="w-10 h-10 rounded-xl bg-zinc-100 dark:bg-zinc-800 text-zinc-900 dark:text-zinc-100 flex items-center justify-center mb-4">
                <CreditCard className="w-5 h-5 text-emerald-600 dark:text-emerald-400" />
              </div>
              <h3 className="text-base font-semibold text-zinc-900 dark:text-zinc-100 mb-1">
                Ceny, úhrady & Ochranný limit
              </h3>
              <p className="text-sm text-zinc-500 dark:text-zinc-400 leading-relaxed font-normal">
                Kalkulace úhrad zdravotní pojišťovny a odhadovaného doplatku v lékárně.
                Interaktivní kalkulačka započitatelného doplatku do ročního ochranného limitu
                (5 000 Kč dospělí, 1 000 Kč děti a senioři, 500 Kč senioři 70+).
              </p>
            </div>
            <div className="mt-5 pt-4 border-t border-zinc-100 dark:border-zinc-800/80 text-xs text-zinc-400">
              Dle zákona č. 48/1997 Sb. o veřejném zdravotním pojištění
            </div>
          </div>

          {/* Card 3: Výpadky & Náhrady */}
          <div className="apple-card p-6 flex flex-col justify-between">
            <div>
              <div className="w-10 h-10 rounded-xl bg-zinc-100 dark:bg-zinc-800 text-zinc-900 dark:text-zinc-100 flex items-center justify-center mb-4">
                <ShieldAlert className="w-5 h-5 text-amber-600 dark:text-amber-400" />
              </div>
              <h3 className="text-base font-semibold text-zinc-900 dark:text-zinc-100 mb-1">
                Výpadky dodávek & Generika
              </h3>
              <p className="text-sm text-zinc-500 dark:text-zinc-400 leading-relaxed font-normal">
                Průběžný monitoring hlášení přerušení a ukončení dodávek na český trh dle § 33 zákona o léčivech.
                Předpokládaný termín obnovení a automatické vyhledání dostupných generických alternativ bez výpadku.
              </p>
            </div>
            <div className="mt-5 pt-4 border-t border-zinc-100 dark:border-zinc-800/80 text-xs text-zinc-400">
              Hlášení držitelů rozhodnutí o registraci
            </div>
          </div>

          {/* Card 4: Kinetika & Harm Reduction */}
          <div className="apple-card p-6 flex flex-col justify-between">
            <div>
              <div className="w-10 h-10 rounded-xl bg-zinc-100 dark:bg-zinc-800 text-zinc-900 dark:text-zinc-100 flex items-center justify-center mb-4">
                <Activity className="w-5 h-5 text-indigo-600 dark:text-indigo-400" />
              </div>
              <h3 className="text-base font-semibold text-zinc-900 dark:text-zinc-100 mb-1">
                Kinetika v těle & Harm Reduction
              </h3>
              <p className="text-sm text-zinc-500 dark:text-zinc-400 leading-relaxed font-normal">
                Jednokompartmentové modelování koncentrace v krvi se 4 cestami podání [Perorální, IV bolus, Nosní, Inhalace].
                Výpočet bezpečného odstupu před další dávkou, prevence toxického vrstvení a analýza fází účinku.
              </p>
            </div>
            <div className="mt-5 pt-4 border-t border-zinc-100 dark:border-zinc-800/80 text-xs text-zinc-400">
              Validovaný Batemanův model & TDM směrnice
            </div>
          </div>
        </div>
      </section>
    </div>
  );
}
