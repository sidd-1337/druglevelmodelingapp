"use client";

import React, { use, useEffect, useState } from "react";
import Link from "next/link";
import { getDrugPresetByCode, getDrugPresetById, getSuklDetailByCode } from "@/lib/sukl-data";
import { fetchDrugDetailApi, routesFromApi } from "@/lib/api";
import { SuklPackageDetail } from "@/lib/types";
import AboutDrugTab from "@/components/tabs/AboutDrugTab";
import PricingCopayTab from "@/components/tabs/PricingCopayTab";
import ShortageSubstitutesTab from "@/components/tabs/ShortageSubstitutesTab";
import KineticsHarmReductionTab from "@/components/tabs/KineticsHarmReductionTab";
import { ArrowLeft, ExternalLink } from "lucide-react";

interface PageProps {
  params: Promise<{ kod: string }>;
}

export default function DrugProfilePage({ params }: PageProps) {
  const resolvedParams = use(params);
  const kod = resolvedParams.kod;

  const [activeTab, setActiveTab] = useState<"about" | "pricing" | "shortage" | "kinetics">("about");

  const [suklDetail, setSuklDetail] = useState<SuklPackageDetail>(() => {
    return (
      getSuklDetailByCode(kod) || {
        kodSukl: kod,
        nazev: `Léčivý přípravek [${kod}]`,
        doplnek: "Balení z registru SÚKL",
        sila: "Neuvedena",
        forma: "Nespecifikováno",
        cesta: "Perorální podání",
        atcKod: "N06BA",
        atcNazev: "Léčivo",
        registracniCislo: "SÚKL/REG",
        drzitel: "Držitel registrace",
        vydejKod: "rx",
        vydejText: "Na lékařský předpis (Rx)",
        dostupnostStatus: "active",
        dostupnostText: "Běžně dodáváno na trh ČR",
        cenaLek: 250,
        maxCenaLek: 290,
        uhradaZp: 180,
        maxDoplatek: 110,
        zapocitatelnyDoplatek: 45,
        uhradaRezim: "partial",
        suklUrl: `https://prehledy.sukl.cz/prehled_leciv.html#/detail/dlp/${kod}`,
      }
    );
  });

  useEffect(() => {
    let active = true;
    fetchDrugDetailApi(kod).then((live) => {
      if (active && live) {
        setSuklDetail(live);
      }
    });
    return () => {
      active = false;
    };
  }, [kod]);

  const drugPreset =
    getDrugPresetByCode(kod) ||
    (suklDetail.matchedPresetId ? getDrugPresetById(suklDetail.matchedPresetId) : null) ||
    (suklDetail.matchedPreset
      ? {
          id: suklDetail.matchedPreset.id,
          name: suklDetail.matchedPreset.name,
          substance: suklDetail.matchedPreset.substance || suklDetail.matchedPreset.name,
          category: suklDetail.matchedPreset.category,
          defaultSuklCode: kod,
          tHalfHours: suklDetail.matchedPreset.t_half_h,
          tHalfRangeHours: suklDetail.matchedPreset.t_half_range_h,
          vdLPerKg: suklDetail.matchedPreset.vd_l_per_kg,
          bioavailability: suklDetail.matchedPreset.bioavailability,
          kaPerH: suklDetail.matchedPreset.ka_per_h,
          cvCl: suklDetail.matchedPreset.cv_cl ?? undefined,
          cvF: suklDetail.matchedPreset.cv_f,
          feRenal: suklDetail.matchedPreset.fe_renal,
          mec: suklDetail.matchedPreset.mec,
          mtc: suklDetail.matchedPreset.mtc,
          unit: suklDetail.matchedPreset.concentration_unit,
          defaultDoseMg: suklDetail.matchedPreset.default_dose_mg,
          defaultIntervalH: suklDetail.matchedPreset.default_interval_h,
          doseConversion: suklDetail.matchedPreset.dose_conversion,
          doseBasis: suklDetail.matchedPreset.dose_basis,
          notes: suklDetail.matchedPreset.notes,
          routes: routesFromApi(suklDetail.matchedPreset.routes),
        }
      : null);

  const TABS = [
    { id: "about", label: "O léku & SÚKL" },
    { id: "pricing", label: "Ceny & Úhrady" },
    { id: "shortage", label: "Výpadky & Náhrady" },
    { id: "kinetics", label: "Kinetika & Harm Reduction" },
  ] as const;

  return (
    <div className="space-y-6">
      {/* Breadcrumb Navigation */}
      <div className="flex items-center justify-between text-xs text-zinc-500 dark:text-zinc-400">
        <div className="flex items-center gap-2">
          <Link
            href="/"
            className="inline-flex items-center gap-1 hover:text-zinc-900 dark:hover:text-zinc-100 transition-colors"
          >
            <ArrowLeft className="w-3.5 h-3.5" />
            <span>Zpět na vyhledávání</span>
          </Link>
          <span>/</span>
          <span className="text-zinc-700 dark:text-zinc-300 font-medium">{suklDetail.nazev}</span>
          <span className="text-zinc-400 font-medium">[{kod}]</span>
        </div>

        <a
          href={suklDetail.suklUrl}
          target="_blank"
          rel="noopener noreferrer"
          className="inline-flex items-center gap-1 hover:text-zinc-900 dark:hover:text-zinc-100 transition-colors"
        >
          <span>Oficiální registr SÚKL</span>
          <ExternalLink className="w-3 h-3" />
        </a>
      </div>

      {/* Apple Style Drug Header Card */}
      <div className="apple-card p-6 flex flex-wrap items-center justify-between gap-4">
        <div>
          <div className="flex items-center gap-3">
            <h1 className="text-2xl sm:text-3xl font-semibold tracking-tight text-zinc-950 dark:text-zinc-50">
              {suklDetail.nazev}
            </h1>
            <span className="text-xs px-2.5 py-1 rounded-full bg-zinc-100 dark:bg-zinc-800 text-zinc-600 dark:text-zinc-300 border border-zinc-200/80 dark:border-zinc-700 font-medium">
              SÚKL {kod}
            </span>
          </div>
          <p className="text-sm text-zinc-500 dark:text-zinc-400 mt-1 font-normal">
            {suklDetail.doplnek || "Registrováno v databázi SÚKL"}
            {suklDetail.drzitel ? ` • ${suklDetail.drzitel}` : ""}
          </p>
        </div>

        {/* Status Badges */}
        <div className="flex flex-wrap items-center gap-2">
          {/* Dispensing badge */}
          <span
            className={`text-xs px-3 py-1 rounded-full font-medium border ${
              suklDetail.vydejKod === "blue_stripe"
                ? "bg-blue-50 text-blue-700 border-blue-200/80 dark:bg-blue-950/40 dark:text-blue-300 dark:border-blue-800/60"
                : suklDetail.vydejKod === "otc"
                ? "bg-emerald-50 text-emerald-700 border-emerald-200/80 dark:bg-emerald-950/40 dark:text-emerald-300 dark:border-emerald-800/60"
                : "bg-zinc-100 text-zinc-700 border-zinc-200 dark:bg-zinc-800 dark:text-zinc-300 dark:border-zinc-700"
            }`}
          >
            {suklDetail.vydejKod === "blue_stripe"
              ? "Recept s modrým pruhem"
              : suklDetail.vydejKod === "otc"
              ? "Volně prodejné (OTC)"
              : "Na lékařský předpis (Rx)"}
          </span>

          {/* Supply status badge */}
          <span
            className={`text-xs px-3 py-1 rounded-full font-medium border ${
              suklDetail.dostupnostStatus === "active"
                ? "bg-emerald-50 text-emerald-700 border-emerald-200/80 dark:bg-emerald-950/40 dark:text-emerald-300 dark:border-emerald-800/60"
                : suklDetail.dostupnostStatus === "interrupted"
                ? "bg-amber-50 text-amber-700 border-amber-200/80 dark:bg-amber-950/40 dark:text-amber-300 dark:border-amber-800/60"
                : "bg-rose-50 text-rose-700 border-rose-200/80 dark:bg-rose-950/40 dark:text-rose-300 dark:border-rose-800/60"
            }`}
          >
            {suklDetail.dostupnostStatus === "active"
              ? "Dodáváno v ČR"
              : suklDetail.dostupnostStatus === "interrupted"
              ? "Výpadek dodávek"
              : "Ukončeno"}
          </span>
        </div>
      </div>

      {/* Apple Style Segmented Pill Control */}
      <div className="overflow-x-auto pb-1">
        <div className="inline-flex p-1 bg-zinc-100 dark:bg-zinc-900 rounded-2xl border border-zinc-200/80 dark:border-zinc-800">
          {TABS.map((tab) => {
            const isActive = activeTab === tab.id;
            return (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`px-4 sm:px-5 py-2 rounded-xl text-xs sm:text-sm transition-all shrink-0 ${
                  isActive
                    ? "bg-white dark:bg-zinc-800 text-zinc-950 dark:text-zinc-50 shadow-sm border border-black/5 dark:border-white/10 font-semibold"
                    : "text-zinc-600 dark:text-zinc-400 hover:text-zinc-900 dark:hover:text-zinc-100 font-medium"
                }`}
              >
                {tab.label}
              </button>
            );
          })}
        </div>
      </div>

      {/* Tab Panels */}
      <div className="pt-1">
        {activeTab === "about" && (
          <AboutDrugTab drug={drugPreset} suklDetail={suklDetail} />
        )}
        {activeTab === "pricing" && (
          <PricingCopayTab suklDetail={suklDetail} />
        )}
        {activeTab === "shortage" && (
          <ShortageSubstitutesTab suklDetail={suklDetail} />
        )}
        {activeTab === "kinetics" && (
          <KineticsHarmReductionTab drug={drugPreset} suklDetail={suklDetail} />
        )}
      </div>
    </div>
  );
}
