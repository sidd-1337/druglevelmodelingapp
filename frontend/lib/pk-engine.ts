import { PKCalculationPoint, RouteType } from "./types";

export interface PKModelParams {
  tHalfHours: number;
  vdLiters: number;
  bioavailability: number;
  kaPerH: number;
  route: RouteType;
  doseConversion?: number;
  routes?: Record<string, { bioavailability?: number; kaPerH?: number }>;
}

const LN2 = Math.LN2;
const REL_TOL = 1e-7;

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
  if (t < 0) return 0.0;
  const ke = getKe(params.tHalfHours);
  const vdL = params.vdLiters;
  if (vdL <= 0) return 0.0;

  if (params.route === "iv_bolus") {
    return (amount / vdL) * Math.exp(-ke * t);
  }

  const ka = getEffectiveKa(params);
  const fEff = getEffectiveBioavailability(params);
  const scale = (fEff * amount) / vdL;

  if (Math.abs(ka - ke) <= REL_TOL * Math.max(ka, ke)) {
    return scale * ka * t * Math.exp(-ke * t);
  }

  // Use expm1 to prevent numerical cancellation when ka ≈ ke
  return (scale * ka * Math.exp(-ke * t) * (-Math.expm1(-(ka - ke) * t))) / (ka - ke);
}

export function calculateTmax(params: PKModelParams): number {
  if (params.route === "iv_bolus") return 0.0;
  const ke = getKe(params.tHalfHours);
  const ka = getEffectiveKa(params);
  if (Math.abs(ka - ke) <= REL_TOL * Math.max(ka, ke)) {
    return 1.0 / ke;
  }
  return Math.log(ka / ke) / (ka - ke);
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
