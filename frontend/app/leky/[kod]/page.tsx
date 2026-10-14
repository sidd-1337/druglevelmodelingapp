"use client";

import React, { use, useEffect, useState } from "react";
import Link from "next/link";
import { getDrugPresetByCode, getDrugPresetById, getSuklDetailByCode } from "@/lib/sukl-data";
import { fetchDrugDetailApi } from "@/lib/api";
import { SuklPackageDetail } from "@/lib/types";
import AboutDrugTab from "@/components/tabs/AboutDrugTab";
import PricingCopayTab from "@/components/tabs/PricingCopayTab";
import ShortageSubstitutesTab from "@/components/tabs/ShortageSubstitutesTab";
import KineticsHarmReductionTab from "@/components/tabs/KineticsHarmReductionTab";
import { packageText, sentenceCase } from "@/lib/format";
import Status, { Tone } from "@/components/ui/Status";
import Tabs, { panelId, tabId } from "@/components/ui/Tabs";

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
          mec: suklDetail.matchedPreset.mec,
          mtc: suklDetail.matchedPreset.mtc,
          unit: suklDetail.matchedPreset.concentration_unit,
          defaultDoseMg: suklDetail.matchedPreset.default_dose_mg,
          defaultIntervalH: suklDetail.matchedPreset.default_interval_h,
          doseConversion: suklDetail.matchedPreset.dose_conversion,
          doseBasis: suklDetail.matchedPreset.dose_basis,
          notes: suklDetail.matchedPreset.notes,
          routes: suklDetail.matchedPreset.routes,
        }
      : null);

  const TABS = [
    { id: "about", label: "Přehled" },
    { id: "pricing", label: "Cena a úhrada" },
    { id: "shortage", label: "Dostupnost" },
    { id: "kinetics", label: "Hladina v těle" },
  ] as const;

  const supplyTone: Tone =
    suklDetail.dostupnostStatus === "active"
      ? "ok"
      : suklDetail.dostupnostStatus === "interrupted"
      ? "warn"
      : suklDetail.dostupnostStatus === "ended"
      ? "bad"
      : "neutral";
  const supplyText =
    suklDetail.dostupnostStatus === "active"
      ? "Dodáváno v ČR"
      : suklDetail.dostupnostStatus === "interrupted"
      ? "Výpadek dodávek"
      : suklDetail.dostupnostStatus === "ended"
      ? "Dodávky ukončeny"
      : "Dostupnost neznámá";
  const dispensingText =
    suklDetail.vydejKod === "blue_stripe"
      ? "Recept s modrým pruhem"
      : suklDetail.vydejKod === "otc"
      ? "Volně prodejné"
      : "Na lékařský předpis";

  return (
    <div>
      <div className="page flex flex-wrap items-baseline justify-between gap-x-8 gap-y-1 pt-8 text-ui text-ink-2">
        <nav aria-label="Drobečková navigace">
          <Link href="/" className="link">
            Hledání
          </Link>
          <span aria-hidden className="mx-2">
            /
          </span>
          <span>{sentenceCase(suklDetail.nazev)}</span>
        </nav>
        <a href={suklDetail.suklUrl} target="_blank" rel="noopener noreferrer" className="link">
          Registr SÚKL
          <span className="sr-only"> (otevře se v novém okně)</span>
        </a>
      </div>

      <header className="page pb-12 pt-6 sm:pb-16">
        <h1 className="text-title">{sentenceCase(suklDetail.nazev)}</h1>
        <p className="mt-3 max-w-prose text-lead text-ink-2">
          {packageText(suklDetail.doplnek) || "Registrováno v databázi SÚKL"}
        </p>
        <div className="mt-6 flex flex-wrap items-baseline gap-x-8 gap-y-2 text-ui">
          <Status tone={supplyTone}>{supplyText}</Status>
          <span>{dispensingText}</span>
          {suklDetail.drzitel && <span className="text-ink-2">{suklDetail.drzitel}</span>}
          <span className="text-ink-2">
            Kód SÚKL <span className="figures text-ink">{kod}</span>
          </span>
        </div>
      </header>

      <Tabs
        label="Části stránky léku"
        idPrefix="lek"
        tabs={TABS}
        value={activeTab}
        onChange={setActiveTab}
      />

      <div
        key={activeTab}
        role="tabpanel"
        id={panelId("lek", activeTab)}
        aria-labelledby={tabId("lek", activeTab)}
        className="page py-14 sm:py-20"
      >
        {activeTab === "about" && <AboutDrugTab drug={drugPreset} suklDetail={suklDetail} />}
        {activeTab === "pricing" && <PricingCopayTab suklDetail={suklDetail} />}
        {activeTab === "shortage" && <ShortageSubstitutesTab suklDetail={suklDetail} />}
        {activeTab === "kinetics" && (
          <KineticsHarmReductionTab drug={drugPreset} suklDetail={suklDetail} />
        )}
      </div>
    </div>
  );
}
