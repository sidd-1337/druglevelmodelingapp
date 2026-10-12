export type RouteType = "oral" | "iv_bolus" | "insufflated" | "smoked";

export interface DrugPreset {
  id: string;
  name: string;
  substance: string;
  category: string;
  defaultSuklCode: string;
  tHalfHours: number;
  tHalfRangeHours?: [number, number];
  vdLPerKg: number;
  bioavailability: number;
  kaPerH: number;
  v1Fraction?: number; // central / steady-state volume (< 1: two compartments)
  qLPerH?: number; // inter-compartmental clearance [L/h]
  tlagH?: number; // oral absorption lag [h]
  cvCl?: number; // between-patient CV of clearance (log-normal)
  cvF?: number; // between-patient CV of bioavailability
  feRenal?: number; // fraction of clearance that scales with creatinine clearance
  mec: number;
  mtc: number;
  unit: string;
  defaultDoseMg: number;
  defaultIntervalH: number;
  doseConversion: number;
  doseBasis: string;
  notes: string;
  routes?: Record<string, { bioavailability?: number; kaPerH?: number }>;
}

export interface SuklPackageDetail {
  kodSukl: string;
  nazev: string;
  doplnek: string;
  sila: string;
  forma: string;
  cesta: string;
  lecivaLatka?: string;
  atcKod: string;
  atcNazev: string;
  registracniCislo: string;
  drzitel: string;
  matchedPresetId?: string;
  matchedPreset?: any;
  vydejKod: "rx" | "otc" | "blue_stripe";
  vydejText: string;
  dostupnostStatus: "active" | "interrupted" | "ended" | "unknown";
  dostupnostText: string;
  duvodPreruseni?: string;
  terminObnoveni?: string;
  cenaLek?: number;
  maxCenaLek?: number;
  uhradaZp?: number;
  maxDoplatek?: number;
  zapocitatelnyDoplatek?: number;
  uhradaRezim: "full" | "partial" | "none";
  indikacniOmezeni?: string;
  preskripcniOmezeni?: string;
  spcUrl?: string;
  pilUrl?: string;
  suklUrl: string;
  substitutes?: Array<{
    kodSukl: string;
    nazev: string;
    doplnek: string;
    doplatek: number;
    dostupnost: "active";
  }>;
}

export interface PKCalculationPoint {
  timeH: number;
  concentration: number;
}

export interface StackingAssessment {
  elapsedH: number;
  currentLevel: number;
  singleCmax: number;
  stackedPeak: number;
  mtc: number;
  mec: number;
  unit: string;
  riskLevel: "safe" | "moderate" | "high" | "critical_toxic";
  riskHeadline: string;
  riskExplanation: string;
  willBreachMtc: boolean;
  safeToRedose: boolean;
}

export interface RouteTimelineInfo {
  route: RouteType;
  labelCs: string;
  labelEn: string;
  onsetStr: string;
  tmaxH: number;
  durationStr: string;
  bioavailabilityPct: number;
  descriptionCs: string;
}

export type CopayAgeTier = "adult" | "child_senior65" | "senior70_disabled";

export interface CopayTierInfo {
  id: CopayAgeTier;
  label: string;
  limitKč: number;
  description: string;
}
