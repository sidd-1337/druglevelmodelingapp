"use client";

import React, { useMemo, useState } from "react";
import {
  assessDoseStacking,
  computeSafeRedosingInterval,
  getRouteTimelineInfo,
} from "@/lib/harm-reduction";
import {
  PKModelParams,
  calculateCmax,
  calculateMultipleDoseCurve,
  calculateTmax,
  patientHalfLife,
  singleDosePoint,
} from "@/lib/pk-engine";
import { DrugPreset, PKCalculationPoint, RouteType, SuklPackageDetail } from "@/lib/types";
import { ExternalLink } from "lucide-react";
import PkChart from "../PkChart";
import RouteSelector from "../RouteSelector";

interface KineticsHarmReductionTabProps {
  drug: DrugPreset | null;
  suklDetail: SuklPackageDetail;
}

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
      tHalfHours: patientHalfLife(baseT12, weightKg),
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

  const singleCmax = useMemo(() => {
    return calculateCmax(effectiveAmount, modelParams);
  }, [effectiveAmount, modelParams]);

  const singleTmax = useMemo(() => {
    return calculateTmax(modelParams);
  }, [modelParams]);

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

  // Route comparison timelines (each calculated with its own route model)
  const routeTimelines = useMemo(() => {
    return (["oral", "iv_bolus", "insufflated", "smoked"] as RouteType[]).map((r) =>
      getRouteTimelineInfo(r, effectiveAmount, { ...modelParams, route: r })
    );
  }, [effectiveAmount, modelParams]);

  const remainingToSafe = Math.max(0, safeIntervalH - elapsedH);
  const safeProgressPct = Math.min(100, Math.round((elapsedH / safeIntervalH) * 100));

  return (
    <div className="space-y-6">
      {/* Streamlit Deep Clinical Laboratory Banner */}
      <div className="apple-card p-5 bg-gradient-to-r from-blue-50/70 via-white to-indigo-50/40 dark:from-zinc-800/80 dark:via-zinc-850 dark:to-zinc-800/60 border border-blue-200/70 dark:border-blue-900/40 shadow-sm">
        <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
          <div className="space-y-1">
            <div className="flex items-center gap-2">
              <span className="text-base">🔬</span>
              <h4 className="font-semibold text-sm text-zinc-950 dark:text-zinc-50">
                Pokročilá klinická laboratoř (Streamlit)
              </h4>
              <span className="text-[10px] px-2 py-0.5 rounded-full font-medium bg-blue-100 text-blue-700 dark:bg-blue-950 dark:text-blue-300 border border-blue-200/60 dark:border-blue-800/60">
                Hloubková simulace
              </span>
            </div>
            <p className="text-xs text-zinc-600 dark:text-zinc-400 max-w-2xl leading-relaxed">
              Modelování <strong>vynechané či odložené dávky (Missed Dose)</strong>, semilogaritmický graf,
              populační variační pásma eliminace a detailní TDM metriky v expertním simulátoru.
            </p>
          </div>

          <a
            href={`http://localhost:8501/?drug=${drug?.id || "livizux"}&sukl=${suklDetail?.kodSukl || ""}&dose=${doseMg}&interval=${intervalH}&weight=${weightKg}&route=${route}`}
            target="_blank"
            rel="noopener noreferrer"
            className="inline-flex items-center justify-center gap-1.5 px-4 py-2 rounded-xl text-xs font-medium bg-zinc-900 hover:bg-zinc-800 text-white dark:bg-zinc-100 dark:hover:bg-white dark:text-zinc-900 transition-all shadow-sm shrink-0"
          >
            <span>Otevřít ve Streamlitu</span>
            <ExternalLink className="w-3.5 h-3.5" />
          </a>
        </div>
      </div>

      {/* Route Selector */}
      <div className="apple-card p-5">
        <div className="flex flex-wrap items-center justify-between gap-2 mb-3">
          <div>
            <label className="text-sm font-semibold text-zinc-900 dark:text-zinc-100">
              Cesta podání léčiva (Route of Administration)
            </label>
            <p className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5">
              Zvolte aplikační cestu pro okamžitý přepočet absorpční konstanty a biologické dostupnosti
            </p>
          </div>
        </div>
        <RouteSelector value={route} onChange={setRoute} />
      </div>

      {/* Regimen Sliders */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="apple-card p-4 space-y-2">
          <div className="flex justify-between items-center text-xs">
            <span className="text-zinc-600 dark:text-zinc-400 font-medium">Dávka léčiva:</span>
            <span className="font-semibold text-zinc-900 dark:text-zinc-100 text-sm">{doseMg} mg</span>
          </div>
          <input
            type="range"
            min={Math.max(1, Math.round(defaultDose * 0.2))}
            max={Math.round(defaultDose * 3)}
            step={defaultDose > 50 ? 5 : 0.5}
            value={doseMg}
            onChange={(e) => setDoseMg(parseFloat(e.target.value))}
            className="w-full accent-apple-blue dark:accent-apple-darkblue cursor-pointer"
          />
        </div>

        <div className="apple-card p-4 space-y-2">
          <div className="flex justify-between items-center text-xs">
            <span className="text-zinc-600 dark:text-zinc-400 font-medium">Dávkovací interval (τ):</span>
            <span className="font-semibold text-zinc-900 dark:text-zinc-100 text-sm">každých {intervalH} h</span>
          </div>
          <input
            type="range"
            min={4}
            max={48}
            step={2}
            value={intervalH}
            onChange={(e) => setIntervalH(parseInt(e.target.value, 10))}
            className="w-full accent-apple-blue dark:accent-apple-darkblue cursor-pointer"
          />
        </div>

        <div className="apple-card p-4 space-y-2">
          <div className="flex justify-between items-center text-xs">
            <span className="text-zinc-600 dark:text-zinc-400 font-medium">Tělesná hmotnost:</span>
            <span className="font-semibold text-zinc-900 dark:text-zinc-100 text-sm">{weightKg} kg</span>
          </div>
          <input
            type="range"
            min={40}
            max={140}
            step={1}
            value={weightKg}
            onChange={(e) => setWeightKg(parseInt(e.target.value, 10))}
            className="w-full accent-apple-blue dark:accent-apple-darkblue cursor-pointer"
          />
        </div>
      </div>

      {/* Interactive PK Curve Chart */}
      <div className="space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div>
            <h4 className="font-semibold text-base text-zinc-900 dark:text-zinc-100 flex items-center gap-2">
              Farmakokinetická křivka & Terapeutické okno (TDM)
            </h4>
            <p className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5">
              Simulace plazmatických koncentrací v závislosti na čase a zvoleném schématu
            </p>
          </div>
          <div className="flex items-center gap-3 text-xs text-zinc-500 dark:text-zinc-400">
            <span>Cmax (1. dávka): <strong className="text-zinc-900 dark:text-zinc-100">{singleCmax < 0.1 ? singleCmax.toFixed(3) : singleCmax.toFixed(2)} {unit}</strong></span>
            <span>•</span>
            <span>Tmax: <strong className="text-zinc-900 dark:text-zinc-100">{singleTmax.toFixed(2)} h</strong></span>
          </div>
        </div>

        <PkChart
          data={curveData}
          mec={mec}
          mtc={mtc}
          unit={unit}
          durationH={durationH}
          doseTimes={doseTimes}
        />
      </div>

      {/* Harm Reduction Calculator Module */}
      <div className="apple-card p-6 space-y-6">
        <div className="flex flex-wrap items-center justify-between border-b border-zinc-200/80 dark:border-zinc-800 pb-4 gap-3">
          <div>
            <h4 className="text-lg font-semibold text-zinc-900 dark:text-zinc-100 flex items-center gap-2">
              Harm Reduction Kalkulačka – Bezpečný odstup & Vrstvení dávek
            </h4>
            <p className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5">
              Predikce kumulace v krvi při opakovaném užití pro zamezení překročení toxického prahu (MTC).
            </p>
          </div>
          <span className="inline-flex items-center gap-1.5 px-3 py-1 rounded-full text-xs font-medium bg-blue-50 text-blue-700 border border-blue-200/80 dark:bg-blue-950/40 dark:text-blue-300 dark:border-blue-800/60">
            Doporučený odstup: ≥ {safeIntervalH.toFixed(1)} h
          </span>
        </div>

        {/* 4 Harm Reduction KPI Cards */}
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-3">
          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 text-center">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block mb-1">Doporučený bezpečný odstup</span>
            <span className="text-xl font-bold tracking-tight text-apple-blue dark:text-apple-darkblue block">
              {safeIntervalH.toFixed(1)} h
            </span>
            <span className="text-[11px] text-zinc-500 dark:text-zinc-400">Minimální doba před další dávkou</span>
          </div>

          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 text-center">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block mb-1">Aktuální zbytková hladina</span>
            <span className="text-xl font-bold tracking-tight text-zinc-900 dark:text-zinc-100 block">
              {stackingAssessment.currentLevel.toFixed(2)} {unit}
            </span>
            <span className="text-[11px] text-zinc-500 dark:text-zinc-400">V čase t = {elapsedH.toFixed(1)} h</span>
          </div>

          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 text-center">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block mb-1">Kumulativní vrchol (peak)</span>
            <span className="text-xl font-bold tracking-tight text-zinc-900 dark:text-zinc-100 block">
              {stackingAssessment.stackedPeak.toFixed(2)} {unit}
            </span>
            <span className="text-[11px] text-zinc-500 dark:text-zinc-400">MTC limit: {mtc} {unit}</span>
          </div>

          <div
            className={`p-4 rounded-xl border text-center transition-all ${
              remainingToSafe > 0
                ? "bg-amber-50/70 border-amber-200/80 dark:bg-amber-950/30 dark:border-amber-800/50"
                : "bg-emerald-50/70 border-emerald-200/80 dark:bg-emerald-950/30 dark:border-emerald-800/50"
            }`}
          >
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block mb-1">Odpočet do bezpečného užití</span>
            <span
              className={`text-xl font-bold tracking-tight block ${
                remainingToSafe > 0
                  ? "text-amber-700 dark:text-amber-400"
                  : "text-emerald-700 dark:text-emerald-400"
              }`}
            >
              {remainingToSafe > 0 ? `${remainingToSafe.toFixed(1)} h` : "0.0 h (Bezpečné)"}
            </span>
            <span className="text-[11px] text-zinc-500 dark:text-zinc-400">
              {remainingToSafe > 0 ? "Vyčkejte na odbourání" : "Lze bezpečně podat"}
            </span>
          </div>
        </div>

        {/* Elapsed Time Slider & Live Countdown Bar */}
        <div className="bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 rounded-xl p-4 space-y-3">
          <div className="flex items-center justify-between text-xs">
            <span className="text-zinc-700 dark:text-zinc-300 font-medium">
              Uplynulý čas od podání předchozí dávky:
            </span>
            <span className="font-mono font-semibold text-zinc-900 dark:text-zinc-100 text-sm">
              {elapsedH.toFixed(1)} hodin
            </span>
          </div>
          <input
            type="range"
            min={0}
            max={Math.max(24, Math.ceil(baseT12 * 4))}
            step={0.5}
            value={elapsedH}
            onChange={(e) => setElapsedH(parseFloat(e.target.value))}
            className="w-full accent-apple-blue dark:accent-apple-darkblue cursor-pointer"
          />
          <div className="flex justify-between text-[11px] text-zinc-500 dark:text-zinc-400">
            <span>0 h (okamžité redávkování)</span>
            <span>{baseT12.toFixed(1)} h (1× poločas t½)</span>
            <span>{(baseT12 * 2).toFixed(1)} h</span>
            <span>{(baseT12 * 4).toFixed(1)} h</span>
          </div>

          {/* Progress bar towards safe redose */}
          <div className="pt-2">
            <div className="flex justify-between text-[11px] text-zinc-500 dark:text-zinc-400 mb-1.5">
              <span>Postup odbourání k bezpečnému rozestupu:</span>
              <span className="font-mono font-medium text-zinc-800 dark:text-zinc-200">{safeProgressPct} %</span>
            </div>
            <div className="w-full h-2 bg-zinc-200 dark:bg-zinc-700 rounded-full overflow-hidden">
              <div
                className={`h-full transition-all duration-300 ${
                  safeProgressPct >= 100
                    ? "bg-emerald-500 dark:bg-emerald-400"
                    : "bg-apple-blue dark:bg-apple-darkblue"
                }`}
                style={{ width: `${safeProgressPct}%` }}
              />
            </div>
          </div>
        </div>

        {/* Stacking Evaluation Alert Box */}
        <div
          className={`p-4 rounded-xl border transition-all ${
            stackingAssessment.riskLevel === "critical_toxic"
              ? "bg-rose-50/80 border-rose-200 dark:bg-rose-950/40 dark:border-rose-800/60"
              : stackingAssessment.riskLevel === "high"
              ? "bg-amber-50/80 border-amber-200 dark:bg-amber-950/40 dark:border-amber-800/60"
              : stackingAssessment.riskLevel === "moderate"
              ? "bg-blue-50/80 border-blue-200 dark:bg-blue-950/40 dark:border-blue-800/60"
              : "bg-emerald-50/80 border-emerald-200 dark:bg-emerald-950/40 dark:border-emerald-800/60"
          }`}
        >
          <div className="flex flex-wrap items-start justify-between gap-4">
            <div className="flex-1 min-w-[240px]">
              <div
                className={`font-semibold text-sm mb-1 ${
                  stackingAssessment.riskLevel === "critical_toxic"
                    ? "text-rose-800 dark:text-rose-300"
                    : stackingAssessment.riskLevel === "high"
                    ? "text-amber-800 dark:text-amber-300"
                    : stackingAssessment.riskLevel === "moderate"
                    ? "text-blue-800 dark:text-blue-300"
                    : "text-emerald-800 dark:text-emerald-300"
                }`}
              >
                {stackingAssessment.riskHeadline}
              </div>
              <p className="text-xs text-zinc-600 dark:text-zinc-300 leading-relaxed">
                {stackingAssessment.riskExplanation}
              </p>
            </div>

            <div className="text-right shrink-0">
              <span className="text-xs block text-zinc-500 dark:text-zinc-400">Kumulativní vrchol:</span>
              <span className="font-mono text-lg font-bold text-zinc-900 dark:text-zinc-100">
                {stackingAssessment.stackedPeak.toFixed(2)} {unit}
              </span>
              <span className="text-[11px] block text-zinc-500 dark:text-zinc-400 mt-0.5">
                (MTC limit: {mtc} {unit})
              </span>
            </div>
          </div>
        </div>

        {/* Dose Stacking Simulation Chart */}
        <div className="space-y-2">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <h5 className="font-semibold text-xs uppercase tracking-wider text-zinc-600 dark:text-zinc-400 flex items-center gap-2">
              Simulace vrstvení hladin (Dávka 1 + Dávka 2 při t = {elapsedH.toFixed(1)} h)
            </h5>
            <span className="text-xs text-zinc-500 dark:text-zinc-400 font-mono">
              Čárkovaně: 2. dávka | Plně: Souhrnná hladina
            </span>
          </div>
          <PkChart
            data={stackedCurveData}
            secondaryData={dose2CurveData}
            secondaryLabel={`2. dávka (+${elapsedH.toFixed(1)} h)`}
            primaryLabel="Kumulativní hladina"
            mec={mec}
            mtc={mtc}
            unit={unit}
            durationH={stackDurationH}
            doseTimes={[0, elapsedH]}
          />
        </div>

        {/* Route Timeline Matrix */}
        <div className="space-y-3">
          <h5 className="font-semibold text-xs uppercase tracking-wider text-zinc-600 dark:text-zinc-400">
            Porovnání kinetiky a nástupu podle cest podání
          </h5>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-3">
            {routeTimelines.map((tl) => (
              <div
                key={tl.route}
                className={`p-3.5 rounded-xl border text-xs transition-all ${
                  route === tl.route
                    ? "bg-blue-50/60 border-apple-blue ring-1 ring-apple-blue dark:bg-blue-950/30 dark:border-apple-darkblue dark:ring-apple-darkblue"
                    : "bg-zinc-50 dark:bg-zinc-800/60 border-zinc-200/80 dark:border-zinc-700/60"
                }`}
              >
                <div className="font-semibold text-zinc-900 dark:text-zinc-100 mb-2 flex items-center justify-between">
                  <span>{tl.labelCs}</span>
                  {route === tl.route && (
                    <span className="text-[10px] text-apple-blue dark:text-apple-darkblue uppercase font-mono font-medium">Vybráno</span>
                  )}
                </div>
                <div className="space-y-1 text-zinc-600 dark:text-zinc-400">
                  <div><strong>Nástup:</strong> {tl.onsetStr}</div>
                  <div><strong>Vrchol (Tmax):</strong> {tl.tmaxH} h</div>
                  <div><strong>Trvání:</strong> {tl.durationStr}</div>
                  <div><strong>Dostupnost F:</strong> {tl.bioavailabilityPct} %</div>
                </div>
                <p className="mt-2 text-[11px] text-zinc-500 dark:text-zinc-400 border-t border-zinc-200/80 dark:border-zinc-800 pt-2">
                  {tl.descriptionCs}
                </p>
              </div>
            ))}
          </div>
        </div>

        {/* Harm Reduction Golden Rules */}
        <div className="space-y-3">
          <h5 className="font-semibold text-xs uppercase tracking-wider text-zinc-600 dark:text-zinc-400">
            Zásady Harm Reduction a prevence předávkování
          </h5>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-3 text-xs">
            <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
              <span className="font-semibold text-zinc-900 dark:text-zinc-100 block mb-1">1. Testujte látku (Drug Checking)</span>
              <p className="text-zinc-600 dark:text-zinc-400">
                Před užitím otestujte substanci reagenciemi (Marquis, Mecke) a proužky na fentanyl k vyloučení nebezpečných příměsí.
              </p>
            </div>
            <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
              <span className="font-semibold text-zinc-900 dark:text-zinc-100 block mb-1">2. Start low, go slow</span>
              <p className="text-zinc-600 dark:text-zinc-400">
                Začněte s malou testovací dávkou (1/4 obvyklého množství). Vyčkejte alespoň na Tmax před zvažováním další dávky.
              </p>
            </div>
            <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
              <span className="font-semibold text-rose-700 dark:text-rose-400 block mb-1">3. Nekombinujte tlumivé látky</span>
              <p className="text-zinc-600 dark:text-zinc-400">
                Kombinace benzodiazepinů, alkoholu a opioidů vede k fatálnímu útlumu dechu. Mějte po ruce Naloxone.
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
