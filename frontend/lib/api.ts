import {
  DrugPreset,
  PKCalculationPoint,
  RouteType,
  StackingAssessment,
  SuklPackageDetail,
} from "./types";
import { getDrugPresetByCode, getSuklDetailByCode, searchSuklCatalog } from "./sukl-data";
import { PKModelParams, calculateMultipleDoseCurve } from "./pk-engine";
import { assessDoseStacking, computeSafeRedosingInterval, getRouteTimelineInfo } from "./harm-reduction";

const API_BASE =
  typeof window === "undefined"
    ? process.env.INTERNAL_API_URL || process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000"
    : process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";
const REQUEST_TIMEOUT_MS = 2500;
// The coupled multi-drug ODE takes about a second for a full regimen on a fast machine.
const MULTI_DRUG_TIMEOUT_MS = 10000;

async function fetchWithTimeout(
  url: string,
  options: RequestInit = {},
  timeoutMs: number = REQUEST_TIMEOUT_MS
): Promise<Response> {
  const controller = new AbortController();
  const id = setTimeout(() => controller.abort(), timeoutMs);
  try {
    const res = await fetch(url, {
      ...options,
      signal: controller.signal,
    });
    clearTimeout(id);
    return res;
  } catch (err) {
    clearTimeout(id);
    throw err;
  }
}

/** API route overrides use snake_case ({ka_per_h}); the PK engine expects {kaPerH}. */
export function routesFromApi(
  routes: Record<string, { bioavailability?: number; ka_per_h?: number }> | null | undefined
): DrugPreset["routes"] {
  if (!routes) return undefined;
  return Object.fromEntries(
    Object.entries(routes).map(([route, cfg]) => [
      route,
      { bioavailability: cfg.bioavailability, kaPerH: cfg.ka_per_h },
    ])
  );
}

export interface SearchResultItem {
  kodSukl: string;
  nazev: string;
  doplnek: string;
  matchedId?: string;
}

export async function searchDrugsApi(query: string): Promise<SearchResultItem[]> {
  const clean = query.trim();
  if (clean.length < 2) return [];

  try {
    const res = await fetchWithTimeout(`${API_BASE}/api/v1/drugs/search?q=${encodeURIComponent(clean)}&limit=20`);
    if (res.ok) {
      const data = await res.json();
      if (Array.isArray(data) && data.length > 0) {
        return data.map((item: any) => ({
          kodSukl: String(item.kod_sukl),
          nazev: item.nazev,
          doplnek: item.doplnek || "",
          matchedId: item.matched_preset_id,
        }));
      }
    }
  } catch {
    // Backend unavailable: fallback gracefully to local catalog
  }

  return searchSuklCatalog(clean);
}

export async function fetchDrugDetailApi(kod: string): Promise<SuklPackageDetail | null> {
  const clean = kod.trim();
  if (!clean) return null;

  try {
    const res = await fetchWithTimeout(`${API_BASE}/api/v1/drugs/${encodeURIComponent(clean)}`);
    if (res.ok) {
      const d = await res.json();
      let substitutes: SuklPackageDetail["substitutes"] = [];

      try {
        const altRes = await fetchWithTimeout(`${API_BASE}/api/v1/drugs/${encodeURIComponent(clean)}/alternatives?limit=5`);
        if (altRes.ok) {
          const alts = await altRes.json();
          if (Array.isArray(alts)) {
            substitutes = alts.map((a: any) => ({
              kodSukl: String(a.kod_sukl),
              nazev: a.nazev,
              doplnek: a.doplnek || "",
              doplatek: 60,
              dostupnost: "active",
            }));
          }
        }
      } catch {
        // Alternatives fetch non-critical
      }

      return {
        kodSukl: String(d.kod_sukl),
        nazev: d.nazev,
        doplnek: d.doplnek || "",
        sila: d.sila || d.doplnek || "",
        forma: d.lekova_forma || "Dle registrace",
        cesta: d.cesta_podani || "Perorální / systémové",
        lecivaLatka: d.leciva_latka || "",
        atcKod: d.atc_kod || "",
        atcNazev: d.atc_nazev || d.leciva_latka || d.nazev,
        registracniCislo: d.registracni_cislo || "SÚKL",
        drzitel: d.drzitel || "",
        matchedPresetId: d.matched_preset_id || (d.matched_preset ? d.matched_preset.id : undefined),
        matchedPreset: d.matched_preset || undefined,
        vydejKod: d.zpusob_vydeje_kod || "rx",
        vydejText: d.zpusob_vydeje_text || "Na lékařský předpis (Rx)",
        dostupnostStatus: d.dostupnost_status || "active",
        dostupnostText: d.dostupnost_text || "Běžně dodáváno v ČR",
        cenaLek: d.cena_lek,
        maxCenaLek: d.max_cena_lek,
        uhradaZp: d.uhrada_zp,
        maxDoplatek: d.max_doplatek,
        zapocitatelnyDoplatek: d.zapocitatelny_doplatek,
        uhradaRezim: d.uhrada_rezim || "partial",
        indikacniOmezeni: d.indikacni_omezeni,
        preskripcniOmezeni: d.preskripcni_omezeni,
        spcUrl: d.spc_url,
        pilUrl: d.pil_url,
        suklUrl: d.sukl_url || `https://prehledy.sukl.cz/prehled_leciv.html#/detail/dlp/${clean}`,
        substitutes: substitutes && substitutes.length > 0 ? substitutes : undefined,
      };
    }
  } catch {
    // Backend offline: fallback to internal static catalog
  }

  return getSuklDetailByCode(clean);
}

// ---------------------------------------------------------------- PK endpoints
// The Python core is the single source of truth for pharmacokinetics. It models
// things the local TypeScript engine cannot: active metabolites, enzyme
// inactivation, elimination pathways and genotype. The local engine is kept only
// as an offline fallback, and anything it cannot represent is flagged to the user.

export interface PkSimulation {
  timeGrid: number[];
  concentrations: number[];
  bandLow: number[];
  bandHigh: number[];
  doseTimes: number[];
  cmax: number;
  tmaxH: number;
  tHalfH: number;
  aucO24: number;
  cavgSs: number;
  cminSs: number;
  cmaxSs: number;
  fluctuationPct: number;
  tdmStatus: string;
}

export interface PkHarmReduction {
  safeIntervalH: number;
  stacking: StackingAssessment;
  singleCurve: PKCalculationPoint[];
  stackedCurve: PKCalculationPoint[];
  routeTimelines: Array<{
    route: string;
    nameCs: string;
    onsetMinutes: [number, number];
    tmaxH: number;
    durationHours: [number, number];
    cmax: number;
    bioavailabilityPct: number;
    descriptionCs: string;
  }>;
}

async function postJson<T>(path: string, body: unknown, timeoutMs?: number): Promise<T | null> {
  try {
    const res = await fetchWithTimeout(
      `${API_BASE}${path}`,
      { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) },
      timeoutMs
    );
    return res.ok ? ((await res.json()) as T) : null;
  } catch {
    return null; // backend unavailable: the caller falls back to the local engine
  }
}

export async function simulatePkApi(req: {
  drugId: string;
  doseMg: number;
  intervalH: number;
  durationH: number;
  bodyWeightKg: number;
  route: RouteType;
}): Promise<PkSimulation | null> {
  const d = await postJson<any>("/api/v1/pk/simulate", {
    drug_id: req.drugId, dose_mg: req.doseMg, interval_h: req.intervalH,
    duration_h: req.durationH, body_weight_kg: req.bodyWeightKg, route: req.route,
  });
  if (!d) return null;
  return {
    timeGrid: d.time_grid, concentrations: d.concentrations,
    bandLow: d.band_low, bandHigh: d.band_high, doseTimes: d.dose_times,
    cmax: d.metrics.cmax, tmaxH: d.metrics.tmax_h, tHalfH: d.metrics.t_half_h,
    aucO24: d.metrics.auc_0_24, cavgSs: d.steady_state.cavg_ss, cminSs: d.steady_state.cmin_ss,
    cmaxSs: d.steady_state.cmax_ss, fluctuationPct: d.steady_state.fluctuation_pct,
    tdmStatus: d.tdm.status,
  };
}

export async function harmReductionApi(req: {
  drugId: string;
  doseMg: number;
  bodyWeightKg: number;
  route: RouteType;
  elapsedH: number;
}): Promise<PkHarmReduction | null> {
  const d = await postJson<any>("/api/v1/pk/harm-reduction", {
    drug_id: req.drugId, dose_mg: req.doseMg, body_weight_kg: req.bodyWeightKg,
    route: req.route, elapsed_h: req.elapsedH,
  });
  if (!d) return null;
  const a = d.stacking_assessment;
  const points = (xs: any[]): PKCalculationPoint[] =>
    xs.map((p) => ({ timeH: p.time_h, concentration: p.concentration }));
  return {
    safeIntervalH: d.safe_redosing_interval_h,
    stacking: {
      elapsedH: a.elapsed_h, currentLevel: a.current_level, singleCmax: a.single_cmax,
      stackedPeak: a.stacked_peak, mtc: a.mtc, mec: a.mec, unit: a.unit,
      riskLevel: a.risk_level, riskHeadline: a.risk_headline,
      riskExplanation: a.risk_explanation, willBreachMtc: a.will_breach_mtc,
      safeToRedose: a.safe_to_redose, timeAboveMtcH: a.time_above_mtc_h,
    },
    singleCurve: points(d.single_curve),
    stackedCurve: points(d.stacked_curve),
    routeTimelines: d.route_timelines.map((t: any) => ({
      route: t.route, nameCs: t.route_name_cs, onsetMinutes: t.onset_minutes_range,
      tmaxH: t.tmax_h, durationHours: t.duration_hours_range, cmax: t.cmax,
      bioavailabilityPct: t.bioavailability_pct, descriptionCs: t.description_cs,
    })),
  };
}

// ------------------------------------------------------------ co-medication
// The coupled enzyme model and the interaction time course run only in the Python
// core, so the co-medication page has no local fallback: it says so when the
// simulation service cannot be reached.

export interface PresetSummary {
  id: string;
  name: string;
  unit: string;
  mec: number;
  mtc: number;
  defaultDoseMg: number;
  defaultIntervalH: number;
}

export interface Perpetrator {
  name: string;
  nameCs: string;
  timed: boolean;
  enzyme: string | null;
  induces: boolean;
  victims: string[];
}

export interface RegimenDrugIn {
  drugId: string;
  doseMg: number;
  intervalH: number;
}

export interface SubstanceIn {
  name: string;
  /** Hour first taken; null when it is already being taken at steady state. */
  startH: number | null;
  stopH: number | null;
}

export interface AppliedInteraction {
  perpetrator: string;
  site: string;
  factor: number | null;
  timed: boolean;
  viaEnzyme: boolean;
  skipped: string;
  pooled: string;
}

export interface PooledEnzyme {
  enzyme: string;
  multiplied: number;
  combined: number;
}

export interface RegimenSeries {
  drugId: string;
  name: string;
  unit: string;
  mec: number;
  mtc: number;
  doseMg: number;
  intervalH: number;
  concentrations: number[];
  alone: number[];
  exposureFactor: number;
  tdmStatus: string;
  applied: AppliedInteraction[];
  drivers: string[];
  pooled: PooledEnzyme[];
}

export interface Clash {
  effect: string;
  participants: string[];
  names: string[];
  overlapH: number;
  firstH: number | null;
}

export interface RegimenResult {
  timeGrid: number[];
  series: RegimenSeries[];
  enzymeActivity: Record<string, number[]>;
  enzymeDrivers: Record<string, string[]>;
  clashes: Clash[];
  pdUnassessed: string[];
  unmatched: string[];
}

async function getJson<T>(path: string): Promise<T | null> {
  try {
    const res = await fetchWithTimeout(`${API_BASE}${path}`);
    return res.ok ? ((await res.json()) as T) : null;
  } catch {
    return null;
  }
}

export async function presetsApi(): Promise<PresetSummary[] | null> {
  const d = await getJson<any[]>("/api/v1/pk/presets");
  return d
    ? d.map((p) => ({
        id: p.id, name: p.name, unit: p.concentration_unit, mec: p.mec, mtc: p.mtc,
        defaultDoseMg: p.default_dose_mg, defaultIntervalH: p.default_interval_h,
      }))
    : null;
}

export async function perpetratorsApi(): Promise<Perpetrator[] | null> {
  const d = await getJson<any[]>("/api/v1/pk/perpetrators");
  return d
    ? d.map((p) => ({
        name: p.name, nameCs: p.name_cs, timed: p.timed, enzyme: p.enzyme, induces: p.induces,
        victims: p.victims,
      }))
    : null;
}

export async function multiDrugApi(req: {
  drugs: RegimenDrugIn[];
  substances: SubstanceIn[];
  durationH: number;
  bodyWeightKg: number;
}): Promise<RegimenResult | null> {
  const d = await postJson<any>("/api/v1/pk/multi-drug", {
    drugs: [
      ...req.drugs.map((x) => ({ drug_id: x.drugId, dose_mg: x.doseMg, interval_h: x.intervalH })),
      ...req.substances.map((s) => ({ perpetrator: s.name, start_h: s.startH, stop_h: s.stopH })),
    ],
    duration_h: req.durationH,
    body_weight_kg: req.bodyWeightKg,
    output_points: 800,
  }, MULTI_DRUG_TIMEOUT_MS);
  if (!d) return null;
  return {
    timeGrid: d.time_grid,
    series: d.series.map((s: any) => ({
      drugId: s.drug_id, name: s.drug_name, unit: s.unit, mec: s.mec, mtc: s.mtc, doseMg: s.dose_mg,
      intervalH: s.interval_h, concentrations: s.concentrations, alone: s.concentrations_alone,
      exposureFactor: s.exposure_factor, tdmStatus: s.tdm.status,
      applied: s.applied.map((a: any) => ({
        perpetrator: a.perpetrator, site: a.site, factor: a.factor, timed: a.timed,
        viaEnzyme: a.via_enzyme, skipped: a.skipped, pooled: a.pooled ?? "",
      })),
      drivers: s.drivers,
      pooled: s.pooled ?? [],
    })),
    enzymeActivity: d.enzyme_activity,
    enzymeDrivers: d.enzyme_drivers,
    clashes: d.clashes.map((c: any) => ({
      effect: c.effect, participants: c.participants, names: c.names, overlapH: c.overlap_h, firstH: c.first_h,
    })),
    pdUnassessed: d.pd_unassessed,
    unmatched: d.unmatched_perpetrators,
  };
}
