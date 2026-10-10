import {
  PKModelParams,
  calculateCmax,
  calculateTmax,
  getEffectiveBioavailability,
  singleDosePoint,
} from "./pk-engine";
import { RouteTimelineInfo, RouteType, StackingAssessment } from "./types";

export function computeSafeRedosingInterval(
  amount: number,
  params: PKModelParams,
  window: { mec: number; mtc: number },
  safetyMarginFactor: number = 0.85
): number {
  const peak = calculateCmax(amount, params);
  const ceiling = safetyMarginFactor * window.mtc;
  const maxResidual = Math.max(0, ceiling - peak);

  if (maxResidual <= 0) {
    return Math.round(Math.max(5.0 * params.tHalfHours, 6.0) * 10) / 10;
  }

  // Consistent with assessDoseStacking: residual must satisfy stacked peak <= 0.85*MTC
  // AND be at or below MEC so previous dose is cleared (if single peak >= MEC)
  const targetResidual = peak >= window.mec ? Math.min(maxResidual, window.mec) : maxResidual;
  const tPeak = calculateTmax(params);
  const horizon = Math.max(48.0, 8.0 * params.tHalfHours);
  const step = 0.05;

  for (let t = tPeak; t <= horizon; t += step) {
    const c = singleDosePoint(t, amount, params);
    if (c <= targetResidual) {
      const minBound = params.route === "iv_bolus" ? 0.5 : tPeak + 0.5;
      return Math.round(Math.max(t, minBound) * 10) / 10;
    }
  }

  return Math.round(4.0 * params.tHalfHours * 10) / 10;
}

export function assessDoseStacking(
  amount: number,
  elapsedH: number,
  params: PKModelParams,
  window: { mec: number; mtc: number },
  unit: string = "µg/mL",
  redoseAmount?: number
): StackingAssessment {
  const nextDose = redoseAmount ?? amount;
  const currentLevel = singleDosePoint(elapsedH, amount, params);
  const singleCmax = calculateCmax(nextDose, params);

  // Peak of combined curves
  let maxStacked = 0;
  const tPeak = calculateTmax(params);
  const searchEnd = Math.max(24.0, tPeak + 4 * params.tHalfHours);
  const dt = 0.05;

  for (let dtOffset = 0; dtOffset <= searchEnd; dtOffset += dt) {
    const residual = singleDosePoint(elapsedH + dtOffset, amount, params);
    const added = singleDosePoint(dtOffset, nextDose, params);
    const sum = residual + added;
    if (sum > maxStacked) {
      maxStacked = sum;
    }
  }

  const willBreachMtc = maxStacked >= window.mtc;
  let riskLevel: "safe" | "moderate" | "high" | "critical_toxic";
  let headline: string;
  let explanation: string;
  let safeToRedose: boolean;

  if (willBreachMtc) {
    riskLevel = "critical_toxic";
    headline = "🚨 Kritické riziko toxicity (překročení MTC)";
    explanation = `Při podání další dávky po ${elapsedH.toFixed(1)} h kumulativní hladina dosáhne ${maxStacked.toFixed(2)} ${unit}, což překračuje toxický limit MTC (${window.mtc} ${unit})! Zbytková hladina z minulé dávky je stále ${currentLevel.toFixed(2)} ${unit}.`;
    safeToRedose = false;
  } else if (maxStacked > 0.85 * window.mtc) {
    riskLevel = "high";
    headline = "⚠️ Vysoké riziko kumulace (blízko toxické zóny)";
    explanation = `Předpokládaný kumulativní vrchol ${maxStacked.toFixed(2)} ${unit} se nachází v horních 15 % bezpečnostního okna (MTC: ${window.mtc} ${unit}). Biologická variabilita může vyvolat nežádoucí reakci.`;
    safeToRedose = false;
  } else if (currentLevel > window.mec) {
    riskLevel = "moderate";
    headline = "ℹ️ Střední riziko – látka je v těle stále plně aktivní";
    explanation = `Zbytková koncentrace (${currentLevel.toFixed(2)} ${unit}) je stále nad minimální účinnou hladinou MEC (${window.mec} ${unit}). Dávky se budou sčítat a prodlouží celkovou dobu působení.`;
    safeToRedose = false;
  } else {
    riskLevel = "safe";
    headline = "✅ Nízké riziko kumulace – bezpečná zóna pro podání";
    explanation = `Předchozí dávka je z větší části odbourána (${currentLevel.toFixed(2)} ${unit} pod MEC). Nový odhadovaný vrchol ${maxStacked.toFixed(2)} ${unit} nepřekročí terapeutické okno.`;
    safeToRedose = true;
  }

  return {
    elapsedH,
    currentLevel,
    singleCmax,
    stackedPeak: maxStacked,
    mtc: window.mtc,
    mec: window.mec,
    unit,
    riskLevel,
    riskHeadline: headline,
    riskExplanation: explanation,
    willBreachMtc,
    safeToRedose,
  };
}

export function getRouteTimelineInfo(
  route: RouteType,
  amount: number,
  baseParams: PKModelParams
): RouteTimelineInfo {
  const params: PKModelParams = { ...baseParams, route };
  const tmax = calculateTmax(params);
  const bioPct = Math.round(getEffectiveBioavailability(params) * 100);
  const tHalf = params.tHalfHours;

  switch (route) {
    case "iv_bolus":
      return {
        route,
        labelCs: "IV Bolus (injekce)",
        labelEn: "IV Bolus",
        onsetStr: "< 1 min",
        tmaxH: 0.0,
        durationStr: `${(2.0 * tHalf).toFixed(1)}–${(3.5 * tHalf).toFixed(1)} h`,
        bioavailabilityPct: 100,
        descriptionCs: "100% biologická dostupnost, okamžitý nástup do systémového krevního řečiště bez zpoždění vstřebáváním.",
      };
    case "insufflated":
      return {
        route,
        labelCs: "Insufflated (šňupání)",
        labelEn: "Insufflated (nasal)",
        onsetStr: "5–15 min",
        tmaxH: Math.round(tmax * 10) / 10,
        durationStr: `${(2.0 * tHalf).toFixed(1)}–${(3.8 * tHalf).toFixed(1)} h`,
        bioavailabilityPct: bioPct,
        descriptionCs: "Rychlý přestup nosní sliznicí do krevního oběhu, obchází primární metabolismus v trávicím traktu.",
      };
    case "smoked":
      return {
        route,
        labelCs: "Smoked (kouření/inhalace)",
        labelEn: "Smoked / Inhaled",
        onsetStr: "1–3 min",
        tmaxH: Math.round(tmax * 10) / 10,
        durationStr: `${(1.5 * tHalf).toFixed(1)}–${(3.2 * tHalf).toFixed(1)} h`,
        bioavailabilityPct: bioPct,
        descriptionCs: "Okamžitý přestup plicními alveoly do arteriálního oběhu a mozku. Velmi strmý nárůst a časné maximum.",
      };
    case "oral":
    default:
      return {
        route: "oral",
        labelCs: "Oral (spolknutí)",
        labelEn: "Oral",
        onsetStr: "20–45 min",
        tmaxH: Math.round(tmax * 10) / 10,
        durationStr: `${(2.5 * tHalf).toFixed(1)}–${(4.5 * tHalf).toFixed(1)} h`,
        bioavailabilityPct: bioPct,
        descriptionCs: "Průchod žaludkem a játry (first-pass metabolismus), pozvolnější náběh a stabilnější trvání hladiny.",
      };
  }
}
