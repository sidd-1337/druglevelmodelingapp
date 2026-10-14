"use client";

import React, { useMemo, useState } from "react";
import {
  assessDoseStacking,
  computeSafeRedosingInterval,
  getRouteTimelineInfo,
} from "@/lib/harm-reduction";
import {
  PKModelParams,
  calculateMultipleDoseCurve,
  singleDosePoint,
} from "@/lib/pk-engine";
import { DrugPreset, PKCalculationPoint, RouteType, SuklPackageDetail } from "@/lib/types";
import { czDecimals, num, sig, stripEmoji } from "@/lib/format";
import PkChart from "../PkChart";
import RouteSelector from "../RouteSelector";
import SectionHead from "../ui/SectionHead";
import Slider from "../ui/Slider";
import SpecList from "../ui/SpecList";
import Status, { Tone } from "../ui/Status";

interface KineticsHarmReductionTabProps {
  drug: DrugPreset | null;
  suklDetail: SuklPackageDetail;
}

const NBSP = " ";
const hours = (h: number) => `${num(h, 1)}${NBSP}h`;

const RISK_TONE: Record<string, Tone> = {
  safe: "ok",
  moderate: "neutral",
  high: "warn",
  critical_toxic: "bad",
};

const RISK_SURFACE: Record<string, string> = {
  safe: "bg-band",
  moderate: "bg-band",
  high: "bg-warn-fill/15",
  critical_toxic: "bg-bad-fill/12",
};

const ROUTE_LABEL: Record<RouteType, string> = {
  oral: "Perorálně (spolknutí)",
  iv_bolus: "Nitrožilně (injekce)",
  insufflated: "Nazálně (šňupání)",
  smoked: "Inhalačně (kouření)",
};

const PRINCIPLES = [
  {
    title: "Otestujte látku",
    text: "Před užitím ji ověřte reagenciemi (Marquis, Mecke) a testovacími proužky na fentanyl, aby se vyloučily nebezpečné příměsi.",
  },
  {
    title: "Začněte nízko, postupujte pomalu",
    text: "Zkuste čtvrtinu obvyklé dávky. Než zvážíte další, počkejte alespoň do času Tmax.",
  },
  {
    title: "Nekombinujte tlumivé látky",
    text: "Benzodiazepiny, alkohol a opioidy dohromady mohou fatálně utlumit dech. Mějte po ruce naloxon.",
  },
];

export default function KineticsHarmReductionTab({
  drug,
  suklDetail,
}: KineticsHarmReductionTabProps) {
  // Model state
  const baseT12 = drug ? drug.tHalfHours : 4.0;
  const baseVdKg = drug ? drug.vdLPerKg : 0.5;
  const baseBio = drug ? drug.bioavailability : 0.8;
  const baseKa = drug ? drug.kaPerH : 1.5;
  const mec = drug ? drug.mec : 5.0;
  const mtc = drug ? drug.mtc : 25.0;
  const unit = drug ? drug.unit : "µg/mL";
  const defaultDose = drug ? drug.defaultDoseMg : 100;
  const defaultInterval = drug ? drug.defaultIntervalH : 12;

  // Extract dose from SÚKL strength text if available (e.g. 70MG -> 70)
  const parsedDetailDose = useMemo(() => {
    const raw = suklDetail?.sila || suklDetail?.doplnek || "";
    const m = raw.match(/(\d+(?:[.,]\d+)?)\s*MG/i);
    return m ? parseFloat(m[1].replace(",", ".")) : null;
  }, [suklDetail]);

  const initialDose = parsedDetailDose || defaultDose;

  const [route, setRoute] = useState<RouteType>("oral");
  const [doseMg, setDoseMg] = useState<number>(initialDose);
  const [intervalH, setIntervalH] = useState<number>(defaultInterval);
  const [weightKg, setWeightKg] = useState<number>(90);
  const [elapsedH, setElapsedH] = useState<number>(defaultInterval);

  const modelParams: PKModelParams = useMemo(() => {
    return {
      tHalfHours: baseT12,
      vdLiters: baseVdKg * weightKg,
      bioavailability: baseBio,
      kaPerH: baseKa,
      route,
      doseConversion: drug ? drug.doseConversion : 1.0,
      routes: drug?.routes,
    };
  }, [baseT12, baseVdKg, weightKg, baseBio, baseKa, route, drug]);

  const effectiveAmount = doseMg * (drug ? drug.doseConversion : 1.0);
  const durationH = Math.max(48, Math.ceil(intervalH * 4));

  // Compute repeated doses schedule
  const doseTimes = useMemo(() => {
    const numDoses = Math.max(1, Math.ceil(durationH / intervalH));
    const times: number[] = [];
    for (let i = 0; i < numDoses; i++) {
      const td = i * intervalH;
      if (td < durationH) times.push(td);
    }
    return times;
  }, [durationH, intervalH]);

  // Compute repeated curve
  const curveData = useMemo(() => {
    return calculateMultipleDoseCurve(durationH, effectiveAmount, intervalH, modelParams);
  }, [durationH, effectiveAmount, intervalH, modelParams]);

  // Harm reduction calculations
  const safeIntervalH = useMemo(() => {
    return computeSafeRedosingInterval(effectiveAmount, modelParams, { mec, mtc });
  }, [effectiveAmount, modelParams, mec, mtc]);

  const stackingAssessment = useMemo(() => {
    return assessDoseStacking(
      effectiveAmount,
      elapsedH,
      modelParams,
      { mec, mtc },
      unit
    );
  }, [effectiveAmount, elapsedH, modelParams, mec, mtc, unit]);

  // Route comparison timelines (each calculated with its own route model)
  const routeTimelines = useMemo(() => {
    return (["oral", "iv_bolus", "insufflated", "smoked"] as RouteType[]).map((r) =>
      getRouteTimelineInfo(r, effectiveAmount, { ...modelParams, route: r })
    );
  }, [effectiveAmount, modelParams]);

  // Stacking curve simulation for Harm Reduction section
  const { stackedCurveData, dose2CurveData, stackDurationH } = useMemo(() => {
    const horizon = Math.max(36, Math.ceil(elapsedH + baseT12 * 4));
    const steps = 240;
    const dt = horizon / steps;
    const sCurve: PKCalculationPoint[] = [];
    const d2Curve: PKCalculationPoint[] = [];

    for (let i = 0; i <= steps; i++) {
      const t = i * dt;
      const c1 = singleDosePoint(t, effectiveAmount, modelParams);
      const c2 = t >= elapsedH ? singleDosePoint(t - elapsedH, effectiveAmount, modelParams) : 0;
      sCurve.push({ timeH: Math.round(t * 100) / 100, concentration: c1 + c2 });
      d2Curve.push({ timeH: Math.round(t * 100) / 100, concentration: c2 });
    }

    return { stackedCurveData: sCurve, dose2CurveData: d2Curve, stackDurationH: horizon };
  }, [elapsedH, baseT12, effectiveAmount, modelParams]);

  const remainingToSafe = Math.max(0, safeIntervalH - elapsedH);
  const safeProgressPct = Math.min(100, Math.round((elapsedH / safeIntervalH) * 100));

  const risk = stackingAssessment.riskLevel;
  const doseDisplay = `${num(doseMg, Number.isInteger(doseMg) ? 0 : 1)}${NBSP}mg`;

  return (
    <div className="space-y-24">
      {/* ------------------------------------------------------------- curve */}
      <section aria-labelledby="krivka">
        <SectionHead id="krivka" title="Hladina léku v krvi">
          Simulace koncentrace v plazmě při zvoleném schématu dávkování.
        </SectionHead>

        <div className="mt-8 max-w-lg">
          <RouteSelector value={route} onChange={setRoute} />
        </div>

        <div className="mt-10 grid gap-10 lg:grid-cols-[minmax(0,1fr)_19rem] lg:gap-14">
          <PkChart
            data={curveData}
            mec={mec}
            mtc={mtc}
            unit={unit}
            durationH={durationH}
            doseTimes={doseTimes}
          />

          <div className="panel space-y-5 p-6 lg:sticky lg:top-24 lg:self-start">
            <Slider
              label="Dávka"
              display={doseDisplay}
              value={doseMg}
              min={Math.max(1, Math.round(defaultDose * 0.2))}
              max={Math.round(defaultDose * 3)}
              step={defaultDose > 50 ? 5 : 0.5}
              onChange={setDoseMg}
            />
            <Slider
              label="Dávkovací interval"
              display={`každých${NBSP}${intervalH}${NBSP}h`}
              value={intervalH}
              min={4}
              max={48}
              step={2}
              onChange={setIntervalH}
            />
            <Slider
              label="Tělesná hmotnost"
              display={`${weightKg}${NBSP}kg`}
              value={weightKg}
              min={40}
              max={140}
              step={1}
              onChange={setWeightKg}
            />
          </div>
        </div>
      </section>

      {/* ---------------------------------------------------- safe interval */}
      <section aria-labelledby="odstup">
        <SectionHead id="odstup" title="Bezpečný odstup mezi dávkami">
          Kolik času má uplynout, než užijete další dávku, aby se hladiny nesečetly nad toxickou mez (MTC).
        </SectionHead>

        <div className="mt-10 grid gap-8 lg:grid-cols-2 lg:gap-14">
          <div className="panel space-y-6 p-6 sm:p-8 lg:self-start">
            <div>
              <p className="text-ui text-ink-2">Doporučený odstup mezi dávkami</p>
              <p className="figures text-figure">
                {"≥"}
                {NBSP}
                {hours(safeIntervalH)}
              </p>
            </div>

            <div>
              <Slider
                label="Od poslední dávky uplynulo"
                display={hours(elapsedH)}
                value={elapsedH}
                min={0}
                max={Math.max(24, Math.ceil(baseT12 * 4))}
                step={0.5}
                onChange={setElapsedH}
              />
              <p className="mt-1 text-small text-ink-2">Poločas léčiva je {hours(baseT12)}.</p>
            </div>

            <div>
              <div className="flex items-baseline justify-between text-ui">
                <span className="text-ink-2">Odbourání k bezpečnému odstupu</span>
                <span className="figures font-semibold">{safeProgressPct}{NBSP}%</span>
              </div>
              <div
                role="progressbar"
                aria-label="Odbourání k bezpečnému odstupu"
                aria-valuemin={0}
                aria-valuemax={100}
                aria-valuenow={safeProgressPct}
                className="mt-2 h-1 overflow-hidden rounded-full bg-line"
              >
                <div
                  className={`h-full rounded-full transition-[width] duration-300 ${
                    safeProgressPct >= 100 ? "bg-ok-fill" : "bg-accent"
                  }`}
                  style={{ width: `${safeProgressPct}%` }}
                />
              </div>
              <p className="mt-3 text-ui">
                {remainingToSafe > 0
                  ? `Vyčkejte ještě ${hours(remainingToSafe)}.`
                  : "Další dávku už lze bezpečně podat."}
              </p>
            </div>
          </div>

          <div className="space-y-8">
            <div className={`rounded-panel p-6 sm:p-8 ${RISK_SURFACE[risk] ?? "bg-band"}`}>
              <h3 className="text-h3">
                <Status tone={RISK_TONE[risk] ?? "neutral"}>{czDecimals(stripEmoji(stackingAssessment.riskHeadline))}</Status>
              </h3>
              <p className="mt-2 text-ui text-ink-2">{czDecimals(stripEmoji(stackingAssessment.riskExplanation))}</p>
            </div>

            <SpecList
              items={[
                {
                  label: "Zbytková hladina",
                  value: `${sig(stackingAssessment.currentLevel)}${NBSP}${unit}`,
                  note: `v čase t = ${hours(elapsedH)}`,
                },
                {
                  label: "Vrchol po další dávce",
                  value: `${sig(stackingAssessment.stackedPeak)}${NBSP}${unit}`,
                  note: `limit MTC ${sig(mtc)}${NBSP}${unit}`,
                },
              ]}
            />
          </div>
        </div>

        <div className="mt-16">
          <h3 className="text-h3">Vrstvení dávek</h3>
          <p className="mb-8 mt-1 max-w-prose text-body text-ink-2">
            První dávka v čase 0 a druhá po {hours(elapsedH)}. Plná čára je součet obou.
          </p>
          <PkChart
            data={stackedCurveData}
            secondaryData={dose2CurveData}
            secondaryLabel="Samotná 2. dávka"
            primaryLabel="Součet obou dávek"
            mec={mec}
            mtc={mtc}
            unit={unit}
            durationH={stackDurationH}
            doseTimes={[0, elapsedH]}
          />
        </div>
      </section>

      {/* ------------------------------------------------------------ routes */}
      <section aria-labelledby="cesty">
        <SectionHead id="cesty" title="Porovnání cest podání">
          Stejná dávka, jiná cesta: liší se nástup, vrchol i délka účinku.
        </SectionHead>

        <div className="-mx-5 mt-8 overflow-x-auto px-5 sm:mx-0 sm:px-0">
          <table className="w-full min-w-[640px] border-collapse text-left">
            <thead>
              <tr>
                <td className="w-40 border-b border-line-soft" />
                {routeTimelines.map((tl) => (
                  <th
                    key={tl.route}
                    scope="col"
                    aria-current={route === tl.route ? "true" : undefined}
                    className={`border-b pb-3 pr-6 align-bottom text-ui font-semibold ${
                      route === tl.route ? "border-b-2 border-ink text-ink" : "border-line-soft text-ink-2"
                    }`}
                  >
                    {ROUTE_LABEL[tl.route] ?? tl.labelCs}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody className="text-ui">
              {[
                { label: "Nástup", cell: (tl: (typeof routeTimelines)[number]) => czDecimals(tl.onsetStr) },
                {
                  label: "Vrchol (Tmax)",
                  cell: (tl: (typeof routeTimelines)[number]) => `${num(tl.tmaxH, tl.tmaxH < 1 ? 2 : 1)}${NBSP}h`,
                },
                { label: "Trvání", cell: (tl: (typeof routeTimelines)[number]) => czDecimals(tl.durationStr) },
                {
                  label: "Biologická dostupnost",
                  cell: (tl: (typeof routeTimelines)[number]) => `${num(tl.bioavailabilityPct)}${NBSP}%`,
                },
                {
                  label: "Poznámka",
                  cell: (tl: (typeof routeTimelines)[number]) => czDecimals(tl.descriptionCs),
                },
              ].map((row) => (
                <tr key={row.label} className="align-top">
                  <th scope="row" className="border-b border-line-soft py-3.5 pr-6 font-normal text-ink-2">
                    {row.label}
                  </th>
                  {routeTimelines.map((tl) => (
                    <td
                      key={tl.route}
                      className={`figures border-b border-line-soft py-3.5 pr-6 ${
                        row.label === "Poznámka" ? "text-small text-ink-2" : "text-ink"
                      }`}
                    >
                      {row.cell(tl)}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      {/* -------------------------------------------------------- principles */}
      <section aria-labelledby="zasady">
        <SectionHead id="zasady" title="Zásady bezpečnějšího užití" />
        <div className="mt-8 grid gap-x-12 gap-y-8 md:grid-cols-3">
          {PRINCIPLES.map((p) => (
            <div key={p.title} className="border-t border-line pt-5">
              <h3 className="text-h3">{p.title}</h3>
              <p className="mt-2 text-ui text-ink-2">{p.text}</p>
            </div>
          ))}
        </div>
      </section>

      {/* ---------------------------------------------------------- advanced */}
      <section className="border-t border-line-soft pt-10">
        <h3 className="text-h3">Pokročilý simulátor</h3>
        <p className="mt-2 max-w-prose text-body text-ink-2">
          Vynechané nebo opožděné dávky, semilogaritmický graf, pásma variability mezi pacienty a podrobné
          TDM metriky.
        </p>
        <a
          href={`http://localhost:8501/?drug=${drug?.id || "livizux"}&sukl=${suklDetail?.kodSukl || ""}&dose=${doseMg}&interval=${intervalH}&weight=${weightKg}&route=${route}`}
          target="_blank"
          rel="noopener noreferrer"
          className="link mt-4 inline-block text-body"
        >
          Otevřít ve Streamlitu
          <span className="sr-only"> (otevře se v novém okně)</span>
        </a>
      </section>
    </div>
  );
}
