"""Client for querying the Czech State Institute for Drug Control (SÚKL) open API.

Endpoints:
- https://prehledy.sukl.cz/prehledy/v1/dlprc (registered medications by query)
- https://prehledy.sukl.cz/prehledy/v1/dlp (general DLP catalog)
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.request
from dataclasses import dataclass
from typing import Any

TIMEOUT_SECONDS = 5.0
USER_AGENT = "DrugLevelModelingApp/2.0 (HealthTech PK Simulator)"

# Mapping of common SÚKL ATC active substances to internal drug preset IDs
ATC_TO_PRESET_MAP: dict[str, str] = {
    "LISDEXAMFETAMIN": "livizux",
    "LISDEXAMFETAMINE": "livizux",
    "LEVETIRACETAM": "levetiracetam",
    "PARACETAMOL": "paracetamol",
    "IBUPROFEN": "ibuprofen",
    "AMOXICILIN": "amoxicillin",
    "AMOXICILLIN": "amoxicillin",
    "LITHIUM": "lithium",
    "DIAZEPAM": "diazepam",
    "PREGABALIN": "pregabalin",
    "PROPRANOLOL": "propranolol",
    "KLONAZEPAM": "rivotril",
    "CLONAZEPAM": "rivotril",
    "MIRTAZAPIN": "mirtazapine",
    "MIRTAZAPINE": "mirtazapine",
    "RISPERIDON": "risperidone",
    "RISPERIDONE": "risperidone",
}

# Default official SÚKL DLP codes for clinical presets
PRESET_DEFAULT_SUKL_CODES: dict[str, str] = {
    "livizux": "0265888",        # LIVIZUX 70MG CPS DUR 100
    "rivotril": "0014957",       # RIVOTRIL 0,5MG TBL NOB 50
    "pregabalin": "0211488",     # PRAGIOLA 75MG CPS DUR 56
    "propranolol": "0284404",    # PROPRANOLOL 80MG CPS PRO 28
    "mirtazapine": "0127760",    # MIRZATEN 30MG TBL FLM 30
    "risperidone": "0102963",    # MEDORISPER 1MG TBL FLM 60
    "lithium": "0207946",        # LITHIUM CARBONICUM SLOVAKOFARMA 300MG
    "levetiracetam": "0267073",  # LEVELANZ 250MG TBL FLM 56 II
    "ibuprofen": "0254292",      # IBALGIN 200MG TBL FLM 10
    "paracetamol": "0254292",    # PARALEN 500 500MG TBL NOB 24
    "amoxicillin": "0203097",    # AMOKSIKLAV 1G TBL FLM 21
    "diazepam": "0002479",       # DIAZEPAM SLOVAKOFARMA 5MG TBL NOB 20
}


@dataclass(frozen=True)
class SuklMedication:
    kod_sukl: str
    nazev: str
    doplnek: str
    sila: str
    lekova_forma: str
    cesta_podani: str
    atc_kod: str
    atc_nazev: str
    parsed_dose_mg: float | None
    matched_preset_id: str | None
    sukl_url: str

    @property
    def display_name(self) -> str:
        dop = f" ({self.doplnek})" if self.doplnek else ""
        return f"{self.nazev}{dop} [SÚKL: {self.kod_sukl}]"


@dataclass(frozen=True)
class SuklDrugDetail:
    kod_sukl: str
    nazev: str
    doplnek: str
    registracni_cislo: str
    drzitel: str
    zpusob_vydeje_kod: str      # "rx", "otc", "blue_stripe"
    zpusob_vydeje_text: str     # Human readable description
    je_dodavka: bool            # Currently marketed and supplied
    dostupnost_status: str      # "active", "interrupted", "ended", "unknown"
    dostupnost_text: str        # Human readable supply description
    cena_lek: float | None      # End / indicative price
    max_cena_lek: float | None  # Maximum pharmacy price
    uhrada_zp: float | None     # Health insurance reimbursement
    max_doplatek: float | None  # Maximum patient copay
    zapocitatelny_doplatek: float | None  # Copay counting towards annual limit
    uhrada_rezim: str           # "full", "partial", "none"
    indikacni_omezeni: str | None
    preskripcni_omezeni: str | None
    spc_url: str | None
    pil_url: str | None
    sukl_url: str


def parse_dose_mg(text: str) -> float | None:
    """Extract numeric dosage in milligrams from strength or supplement text.

    Examples:
    '70MG' -> 70.0
    '500 MG' -> 500.0
    '2,5MG/ML' -> 2.5
    '0.5MG' -> 0.5
    '1 G' -> 1000.0
    """
    if not text:
        return None

    # Check for grams first (e.g. 1G)
    g_match = re.search(r"(\d+(?:[.,]\d+)?)\s*G(?![A-Za-z])", text, re.I)
    if g_match:
        try:
            return float(g_match.group(1).replace(",", ".")) * 1000.0
        except ValueError:
            pass

    # Check for milligrams (e.g. 70MG, 100 MG)
    mg_match = re.search(r"(\d+(?:[.,]\d+)?)\s*MG", text, re.I)
    if mg_match:
        try:
            return float(mg_match.group(1).replace(",", "."))
        except ValueError:
            pass

    # Check for micrograms (e.g. 500MCG, 500UG)
    ug_match = re.search(r"(\d+(?:[.,]\d+)?)\s*(?:MCG|UG)", text, re.I)
    if ug_match:
        try:
            return float(ug_match.group(1).replace(",", ".")) / 1000.0
        except ValueError:
            pass

    return None


def match_preset_id(atc_nazev: str, nazev_lp: str) -> str | None:
    """Match SÚKL item to an internal drug preset ID."""
    clean_atc = atc_nazev.upper().strip()
    for key, preset_id in ATC_TO_PRESET_MAP.items():
        if key in clean_atc:
            return preset_id

    clean_nazev = nazev_lp.upper().strip()
    for key, preset_id in ATC_TO_PRESET_MAP.items():
        if key in clean_nazev:
            return preset_id

    return None


def _clean_str(val: Any) -> str:
    if val is None:
        return ""
    if isinstance(val, dict):
        return str(val.get("cs") or val.get("en") or "")
    return str(val).strip()


def query_sukl_api(query: str, max_results: int = 10) -> list[SuklMedication]:
    """Query SÚKL official public API for registered medications matching query string."""
    query = query.strip()
    if len(query) < 2:
        return []

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": USER_AGENT,
    }

    # Try /dlprc first (matches commercial brands and active substances), then fallback to /dlp
    endpoints = ["dlprc", "dlp"]
    raw_items: list[dict] = []

    for ep in endpoints:
        payload = json.dumps({"filtr": query, "pocet": max_results, "stranka": 1}).encode("utf-8")
        req = urllib.request.Request(
            f"https://prehledy.sukl.cz/prehledy/v1/{ep}",
            data=payload,
            headers=headers,
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT_SECONDS) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                items = data.get("data", [])
                if items:
                    raw_items = items
                    break
        except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError):
            continue

    results: list[SuklMedication] = []
    seen_codes: set[str] = set()

    for item in raw_items:
        kod = _clean_str(item.get("kodSUKL"))
        if not kod or kod in seen_codes:
            continue
        seen_codes.add(kod)

        nazev = _clean_str(item.get("nazevLP"))
        doplnek = _clean_str(item.get("doplnekNazvu"))
        sila = _clean_str(item.get("sila"))
        lekova_forma = _clean_str(item.get("lekovaForma"))
        cesta = _clean_str(item.get("cestaPodani"))

        atc_obj = item.get("atc") or item.get("ATCskupina") or {}
        if isinstance(atc_obj, dict):
            atc_kod = _clean_str(atc_obj.get("kod"))
            atc_nazev = _clean_str(atc_obj.get("nazev"))
        else:
            atc_kod = str(atc_obj)
            atc_nazev = _clean_str(item.get("ATCnazev"))

        dose = parse_dose_mg(sila) or parse_dose_mg(doplnek)
        preset = match_preset_id(atc_nazev, nazev)
        url = f"https://prehledy.sukl.cz/prehled_leciv.html#/detail/dlp/{kod}"

        results.append(
            SuklMedication(
                kod_sukl=kod,
                nazev=nazev,
                doplnek=doplnek,
                sila=sila,
                lekova_forma=lekova_forma,
                cesta_podani=cesta,
                atc_kod=atc_kod,
                atc_nazev=atc_nazev,
                parsed_dose_mg=dose,
                matched_preset_id=preset,
                sukl_url=url,
            )
        )

    return results


def fetch_sukl_detail(kod_sukl: str) -> SuklDrugDetail | None:
    """Fetch complete regulatory, pricing, reimbursement, and shortage data for a SÚKL code.

    Endpoint: https://prehledy.sukl.cz/prehledy/v1/dlp/{kodSUKL}
    """
    clean_kod = str(kod_sukl).strip().zfill(7)
    if len(clean_kod) < 4 or not clean_kod.isdigit():
        return None

    url = f"https://prehledy.sukl.cz/prehledy/v1/dlp/{clean_kod}"
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_SECONDS) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError):
        return None

    if not isinstance(data, dict) or not data.get("kodSUKL"):
        return None

    nazev = _clean_str(data.get("nazevLP"))
    doplnek = _clean_str(data.get("doplnekNazvu"))
    reg_cislo = _clean_str(data.get("registracniCisloDisplay") or data.get("registracniCislo"))
    drzitel_obj = data.get("registrace", {}).get("drzitelRozhodnuti") if isinstance(data.get("registrace"), dict) else {}
    drzitel = _clean_str(drzitel_obj.get("nazev") if isinstance(drzitel_obj, dict) else "")

    # Vydej (Rx / OTC / Modry pruh)
    vydej_obj = data.get("zpusobVydeje") or {}
    vydej_nazev = _clean_str(vydej_obj.get("nazevDlw") or vydej_obj.get("nazev"))
    vydej_kod_dlw = _clean_str(vydej_obj.get("kodDlw") or vydej_obj.get("kod")).upper()
    if "modr" in vydej_nazev.lower():
        zpusob_vydeje_kod = "blue_stripe"
    elif "bez" in vydej_nazev.lower() or vydej_kod_dlw in ["F", "OTC"]:
        zpusob_vydeje_kod = "otc"
    else:
        zpusob_vydeje_kod = "rx"

    # Dostupnost & dodavky
    je_dodavka = bool(data.get("jeDodavka"))
    dostupnost_raw = str(data.get("dostupnost") or "").upper()
    if dostupnost_raw in ["D", "A"] or je_dodavka:
        dostupnost_status = "active"
        dostupnost_text = "Běžně dodáváno na trh ČR"
    elif dostupnost_raw == "P":
        dostupnost_status = "interrupted"
        dostupnost_text = "Přerušení dodávek (výpadek)"
    elif dostupnost_raw == "U":
        dostupnost_status = "ended"
        dostupnost_text = "Uvádění na trh ukončeno"
    else:
        dostupnost_status = "active" if je_dodavka else "unknown"
        dostupnost_text = "Běžně dodáváno" if je_dodavka else "Bez hlášení o dodávkách"

    # Ceny a uhrady
    uhrady_obj = data.get("uhrady") if isinstance(data.get("uhrady"), dict) else {}
    zaklad_uhrada = uhrady_obj.get("zakladniUhrada") if isinstance(uhrady_obj.get("zakladniUhrada"), dict) else {}

    def _safe_float(v: Any) -> float | None:
        if v is None:
            return None
        try:
            return float(v)
        except (ValueError, TypeError):
            return None

    cena_lek = _safe_float(data.get("cena"))
    max_cena_lek = _safe_float(uhrady_obj.get("maxCenaLek")) or cena_lek
    uhrada_zp = _safe_float(zaklad_uhrada.get("maximalniUhrada"))
    max_doplatek = _safe_float(zaklad_uhrada.get("maximalniDoplatek"))
    zapocitatelny_doplatek = _safe_float(zaklad_uhrada.get("zapocitatelnyDoplatek"))

    # Rezim uhrady
    if uhrada_zp is None or uhrada_zp <= 0:
        uhrada_rezim = "none"
        if max_cena_lek and max_doplatek is None:
            max_doplatek = max_cena_lek
    elif max_doplatek is not None and max_doplatek > 0:
        uhrada_rezim = "partial"
    else:
        uhrada_rezim = "full"

    indikacni_omezeni = _clean_str(zaklad_uhrada.get("indikacniOmezeni")) or None
    preskripcni_omezeni = _clean_str(zaklad_uhrada.get("specializacePredepisujicihoLekare")) or None

    # Dokumenty SPC a PIL
    dokumenty = data.get("dokumenty") if isinstance(data.get("dokumenty"), list) else []
    spc_url, pil_url = None, None
    for doc in dokumenty:
        typ = str(doc.get("typDokumentu") or "").upper()
        nazev_doc = _clean_str(doc.get("nazev"))
        if typ == "SPC" and nazev_doc and not spc_url:
            spc_url = f"https://www.sukl.cz/modules/medication/download.php?file={nazev_doc}&type=spc"
        elif typ == "PIL" and nazev_doc and not pil_url:
            pil_url = f"https://www.sukl.cz/modules/medication/download.php?file={nazev_doc}&type=pil"

    sukl_url = f"https://prehledy.sukl.cz/prehled_leciv.html#/detail/dlp/{clean_kod}"

    return SuklDrugDetail(
        kod_sukl=clean_kod,
        nazev=nazev,
        doplnek=doplnek,
        registracni_cislo=reg_cislo,
        drzitel=drzitel,
        zpusob_vydeje_kod=zpusob_vydeje_kod,
        zpusob_vydeje_text=vydej_nazev or ("Na lékařský předpis (Rx)" if zpusob_vydeje_kod == "rx" else "Volně prodejný (OTC)"),
        je_dodavka=je_dodavka,
        dostupnost_status=dostupnost_status,
        dostupnost_text=dostupnost_text,
        cena_lek=cena_lek,
        max_cena_lek=max_cena_lek,
        uhrada_zp=uhrada_zp,
        max_doplatek=max_doplatek,
        zapocitatelny_doplatek=zapocitatelny_doplatek,
        uhrada_rezim=uhrada_rezim,
        indikacni_omezeni=indikacni_omezeni,
        preskripcni_omezeni=preskripcni_omezeni,
        spc_url=spc_url,
        pil_url=pil_url,
        sukl_url=sukl_url,
    )

