import { PKCalculationPoint, RouteType } from "./types";

export interface PKModelParams {
  tHalfHours: number;
  vdLiters: number;
  bioavailability: number;
  kaPerH: number;
  route: RouteType;
  doseConversion?: number;
  routes?: Record<string, { bioavailability?: number; kaPerH?: number }>;
  /** Central volume as a fraction of vdLiters; < 1 adds a peripheral compartment. */
  v1Fraction?: number;
  /** Inter-compartmental clearance [L/h], required when v1Fraction < 1. */
  qLPerH?: number;
  /** Absorption lag [h] (extravascular routes only). */
  tlagH?: number;
}

/** Exponential terms of the single-dose curve: C(t) = sum(coefs[i] * exp(-rates[i] * t)). */
function expTerms(amount: number, params: PKModelParams): { coefs: number[]; rates: number[] } {
  const v1 = params.vdLiters * (params.v1Fraction ?? 1);
  const cl = getKe(params.tHalfHours) * params.vdLiters;
  const twoCompartment = (params.v1Fraction ?? 1) < 1;
  const k10 = cl / v1;
  const k12 = twoCompartment ? (params.qLPerH ?? 0) / v1 : 0;
  const k21 = twoCompartment ? (params.qLPerH ?? 0) / (params.vdLiters - v1) : 0;

  let rates: number[];
  if (twoCompartment) {
    const sum = k10 + k12 + k21;
    const root = Math.sqrt(sum * sum - 4 * k10 * k21);
    rates = [(sum + root) / 2, (sum - root) / 2];
  } else {
    rates = [k10];
  }

  if (params.route === "iv_bolus") {
    if (!twoCompartment) return { coefs: [amount / v1], rates };
    const [a, b] = rates;
    return { coefs: [((a - k21) / (a - b)) * (amount / v1), ((k21 - b) / (a - b)) * (amount / v1)], rates };
  }

  let ka = getEffectiveKa(params);
  for (const r of rates) {
    if (Math.abs(ka - r) <= REL_TOL * Math.max(ka, r)) ka = r * (1 + 1e-5);
  }
  const allRates = [...rates, ka];
  const scale = (getEffectiveBioavailability(params) * amount * ka) / v1;
  const coefs = allRates.map((r, i) => {
    const num = twoCompartment ? k21 - r : 1;
    const den = allRates.reduce((acc, o, j) => (j === i ? acc : acc * (o - r)), 1);
    return (scale * num) / den;
  });
  return { coefs, rates: allRates };
}

const LN2 = Math.LN2;
const REL_TOL = 1e-7;

const REFERENCE_WEIGHT_KG = 70;
const REFERENCE_CRCL_ML_MIN = 100;

/**
 * Half-life for a patient, mirroring core/population.py: V scales with weight,
 * CL with weight^0.75 (and its renal fraction with CrCl / 100 mL/min).
 * Preset half-lives describe a 70 kg adult with normal renal function.
 */
export function patientHalfLife(
  tHalfRefHours: number,
  weightKg: number,
  crclMlMin?: number,
  feRenal = 0
): number {
  const size = Math.pow(weightKg / REFERENCE_WEIGHT_KG, 0.75);
  const clRel =
    crclMlMin === undefined || feRenal === 0
      ? size
      : (1 - feRenal) * size + (feRenal * crclMlMin) / REFERENCE_CRCL_ML_MIN;
  return (tHalfRefHours * (weightKg / REFERENCE_WEIGHT_KG)) / clRel;
}

export function getKe(tHalfHours: number): number {
  return tHalfHours > 0 ? LN2 / tHalfHours : 0;
}

export function getEffectiveBioavailability(params: PKModelParams): number {
  if (params.routes && params.routes[params.route]?.bioavailability !== undefined) {
    return params.routes[params.route].bioavailability!;
  }
  if (params.route === "iv_bolus") return 1.0;
  if (params.route === "insufflated") return Math.min(1.0, Math.max(params.bioavailability, 0.70));
  if (params.route === "smoked") return Math.min(1.0, Math.max(params.bioavailability * 0.9, 0.60));
  return params.bioavailability;
}

export function getEffectiveKa(params: PKModelParams): number {
  if (params.routes && params.routes[params.route]?.kaPerH !== undefined) {
    return params.routes[params.route].kaPerH!;
  }
  if (params.route === "iv_bolus") return 1.0;
  if (params.route === "insufflated") return Math.max(params.kaPerH * 2.5, 3.5);
  if (params.route === "smoked") return Math.max(params.kaPerH * 5.0, 8.0);
  return params.kaPerH;
}

export function singleDosePoint(t: number, amount: number, params: PKModelParams): number {
  const lag = params.route === "iv_bolus" ? 0 : params.tlagH ?? 0;
  const tp = t - lag;
  if (tp < 0) return 0.0;
  const ke = getKe(params.tHalfHours);
  const vdL = params.vdLiters;
  if (vdL <= 0) return 0.0;

  if ((params.v1Fraction ?? 1) < 1) {
    const { coefs, rates } = expTerms(amount, params);
    return Math.max(0, coefs.reduce((acc, c, i) => acc + c * Math.exp(-rates[i] * tp), 0));
  }

  if (params.route === "iv_bolus") {
    return (amount / vdL) * Math.exp(-ke * tp);
  }

  const ka = getEffectiveKa(params);
  const fEff = getEffectiveBioavailability(params);
  const scale = (fEff * amount) / vdL;

  if (Math.abs(ka - ke) <= REL_TOL * Math.max(ka, ke)) {
    return scale * ka * tp * Math.exp(-ke * tp);
  }

  // Use expm1 to prevent numerical cancellation when ka ≈ ke
  return (scale * ka * Math.exp(-ke * tp) * (-Math.expm1(-(ka - ke) * tp))) / (ka - ke);
}

export function calculateTmax(params: PKModelParams): number {
  if (params.route === "iv_bolus") return 0.0;
  const lag = params.tlagH ?? 0;
  const ke = getKe(params.tHalfHours);
  const ka = getEffectiveKa(params);
  if ((params.v1Fraction ?? 1) < 1) {
    // Two compartments: locate the peak numerically.
    const horizon = lag + 10 / ka + 2 * terminalHalfLife(params);
    const steps = 2000;
    let best = 0;
    let bestC = -1;
    for (let i = 0; i <= steps; i++) {
      const t = (i * horizon) / steps;
      const c = singleDosePoint(t, 1, params);
      if (c > bestC) {
        bestC = c;
        best = t;
      }
    }
    return best;
  }
  if (Math.abs(ka - ke) <= REL_TOL * Math.max(ka, ke)) {
    return lag + 1.0 / ke;
  }
  return lag + Math.log(ka / ke) / (ka - ke);
}

/** Half-life of the elimination phase (longer than tHalfHours for two compartments). */
export function terminalHalfLife(params: PKModelParams): number {
  if ((params.v1Fraction ?? 1) >= 1) return params.tHalfHours;
  const { rates } = expTerms(1, { ...params, route: "iv_bolus" });
  return LN2 / Math.min(...rates);
}

export function calculateCmax(amount: number, params: PKModelParams): number {
  const tmax = calculateTmax(params);
  return singleDosePoint(tmax, amount, params);
}

export function calculateMultipleDoseCurve(
  durationH: number,
  amount: number,
  intervalH: number,
  params: PKModelParams,
  steps: number = 300
): PKCalculationPoint[] {
  const dt = durationH / steps;
  const numDoses = Math.max(1, Math.ceil(durationH / intervalH));
  const doseTimes: number[] = [];
  for (let i = 0; i < numDoses; i++) {
    const td = i * intervalH;
    if (td < durationH) doseTimes.push(td);
  }

  const points: PKCalculationPoint[] = [];
  for (let i = 0; i <= steps; i++) {
    const t = i * dt;
    let totalC = 0.0;
    for (const td of doseTimes) {
      if (t >= td) {
        totalC += singleDosePoint(t - td, amount, params);
      }
    }
    points.push({
      timeH: Math.round(t * 100) / 100,
      concentration: Math.max(0, totalC),
    });
  }

  return points;
}
