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

async function fetchWithTimeout(url: string, options: RequestInit = {}): Promise<Response> {
  const controller = new AbortController();
  const id = setTimeout(() => controller.abort(), REQUEST_TIMEOUT_MS);
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
