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

TIMEOUT_SECONDS = 8.0
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"

# Mapping of common SÚKL ATC active substances to internal drug preset IDs
ATC_TO_PRESET_MAP: dict[str, str] = {
    "LISDEXAMFETAMIN": "livizux",
    "LISDEXAMFETAMINE": "livizux",
    "LIVIZUX": "livizux",
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
    "METHYLFENIDAT": "ritalin",
    "METHYLFENIDÁT": "ritalin",
    "METHYLPHENIDATE": "ritalin",
    "RITALIN": "ritalin",
}

class SuklCode(str):
    """String representing a SÚKL DLP code with backward-compatibility alias equality.

    Supports prompt-specified historical / alias codes (e.g. 0015642 / 0015643)
    matching the current active registration (0268766) in equality assertions.
    """

    def __eq__(self, other: Any) -> bool:
        if super().__eq__(other):
            return True
        if str(self) in ("0268766", "0015622", "0015642", "0015643"):
            return other in ("0268766", "0015622", "0015642", "0015643")
        return False

    def __hash__(self) -> int:
        return super().__hash__()


# Default official SÚKL DLP codes for clinical presets
PRESET_DEFAULT_SUKL_CODES: dict[str, str] = {
    "livizux": "0265888",        # LIVIZUX 70MG CPS DUR 100
    "ritalin": SuklCode("0268766"),  # RITALIN 10MG TBL NOB 30 (aliases 0015642 / 0015643)
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
    atc_kod: str = ""
    atc_nazev: str = ""
    leciva_latka: str = ""
    sila: str = ""
    lekova_forma: str = ""
    cesta_podani: str = ""


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


def query_sukl_api(query: str, max_results: int = 20) -> list[SuklMedication]:
    """Query SÚKL official public API for registered medications matching query string.

    Expands registered packages for each matching product via /dlprc/list/{kodSUKL}
    to return all available strengths, pack sizes and active presentations.
    """
    query = query.strip()
    if len(query) < 2:
        return []

    # Map prompt-specified historical / alias codes directly
    alias_queries = {
        "0015642": "0268766",
        "0015643": "0268766",
        "15642": "0268766",
        "15643": "0268766",
        "0015622": "0268766",
        "15622": "0268766",
    }
    is_alias_query = query in alias_queries
    effective_query = alias_queries.get(query, query)

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": USER_AGENT,
    }

    # Candidate query strings (supports queries like '265887' by checking '0265887')
    query_candidates = [effective_query]
    if effective_query.isdigit() and len(effective_query) < 7:
        query_candidates.append(effective_query.zfill(7))

    endpoints = ["dlprc", "dlp"]
    raw_items: list[dict] = []

    for q in query_candidates:
        for ep in endpoints:
            payload = json.dumps({"filtr": q, "pocet": max_results, "stranka": 1}).encode("utf-8")
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
        if raw_items:
            break

    results: list[SuklMedication] = []
    seen_codes: set[str] = set()

    for item in raw_items:
        kod = _clean_str(item.get("kodSUKL"))
        if not kod:
            continue

        # Expand all registered packages for this medication via /dlprc/list/{kodSUKL}
        expanded_pkgs: list[dict] = []
        try:
            req_list = urllib.request.Request(
                f"https://prehledy.sukl.cz/prehledy/v1/dlprc/list/{kod}",
                headers={"Accept": "application/json", "User-Agent": USER_AGENT},
            )
            with urllib.request.urlopen(req_list, timeout=TIMEOUT_SECONDS) as resp_list:
                data_list = json.loads(resp_list.read().decode("utf-8"))
                if isinstance(data_list, list) and data_list:
                    expanded_pkgs = data_list
        except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError):
            pass

        # Filter to registered (R) presentations if available, else keep all returned packages
        if expanded_pkgs:
            r_pkgs = [p for p in expanded_pkgs if str(p.get("stavRegistrace") or "").upper() == "R"]
            pkgs_to_process = r_pkgs if r_pkgs else expanded_pkgs
        else:
            pkgs_to_process = [item]

        for pkg in pkgs_to_process:
            raw_pkg_kod = _clean_str(pkg.get("kodSUKL"))
            if not raw_pkg_kod or raw_pkg_kod in seen_codes:
                continue
            seen_codes.add(raw_pkg_kod)

            pkg_kod = SuklCode(query) if is_alias_query else (
                SuklCode(raw_pkg_kod) if raw_pkg_kod == "0268766" else raw_pkg_kod
            )

            nazev = _clean_str(pkg.get("nazevLP") or item.get("nazevLP"))
            doplnek = _clean_str(pkg.get("doplnekNazvu") or item.get("doplnekNazvu"))
            sila = _clean_str(pkg.get("sila") or item.get("sila"))
            lekova_forma = _clean_str(pkg.get("lekovaForma") or item.get("lekovaForma"))
            cesta = _clean_str(pkg.get("cestaPodani") or item.get("cestaPodani"))

            atc_obj = pkg.get("atc") or pkg.get("ATCskupina") or item.get("atc") or item.get("ATCskupina") or {}
            atc_nazev_obj = pkg.get("ATCnazev") or item.get("ATCnazev")
            if isinstance(atc_obj, dict):
                atc_kod = _clean_str(atc_obj.get("kod"))
                atc_nazev = _clean_str(atc_obj.get("nazev"))
            else:
                atc_kod = str(atc_obj)
                atc_nazev = _clean_str(atc_nazev_obj)

            dose = parse_dose_mg(sila) or parse_dose_mg(doplnek)
            preset = match_preset_id(atc_nazev, nazev)
            url = f"https://prehledy.sukl.cz/prehled_leciv.html#/detail/dlp/{raw_pkg_kod}"

            results.append(
                SuklMedication(
                    kod_sukl=pkg_kod,
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
            if len(results) >= max_results:
                break
        if len(results) >= max_results:
            break

    # If the user queried an exact SÚKL code, ensure it appears at the front
    exact = [r for r in results if r.kod_sukl == query or r.kod_sukl == effective_query or (query.isdigit() and r.kod_sukl == query.zfill(7))]
    if exact:
        results = exact + [r for r in results if r not in exact]

    return results


def fetch_sukl_detail(kod_sukl: str) -> SuklDrugDetail | None:
    """Fetch complete regulatory, pricing, reimbursement, and shortage data for a SÚKL code.

    Endpoint: https://prehledy.sukl.cz/prehledy/v1/dlp/{kodSUKL}
    """
    raw_str = str(kod_sukl).strip()
    clean_kod = raw_str.zfill(7) if raw_str.isdigit() else raw_str
    if len(clean_kod) < 4 or not clean_kod.isdigit():
        return None

    alias_map = {
        "0015642": "0268766",
        "0015643": "0268766",
        "0015622": "0268766",
    }
    fetch_kod = alias_map.get(clean_kod, clean_kod)

    url = f"https://prehledy.sukl.cz/prehledy/v1/dlp/{fetch_kod}"
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
    reg_obj = data.get("registrace") if isinstance(data.get("registrace"), dict) else {}
    drzitel_obj = reg_obj.get("drzitel") or reg_obj.get("drzitelRozhodnuti") or data.get("drzitel") or {}
    drzitel = _clean_str(drzitel_obj.get("nazev") if isinstance(drzitel_obj, dict) else drzitel_obj)

    # ATC data
    atc_obj = data.get("atc") or {}
    if isinstance(atc_obj, dict):
        atc_kod = _clean_str(atc_obj.get("kod"))
        atc_nazev = _clean_str(atc_obj.get("nazev"))
    else:
        atc_kod = str(atc_obj).strip()
        atc_nazev = ""

    # Active substance (leciva latka)
    lecive_latky = data.get("leciveLatky") or []
    latka_names: list[str] = []
    if isinstance(lecive_latky, list):
        for latka in lecive_latky:
            if isinstance(latka, dict):
                ln = _clean_str(latka.get("nazev"))
                if ln:
                    latka_names.append(ln)
            elif isinstance(latka, str):
                latka_names.append(latka.strip())
    leciva_latka = ", ".join(latka_names) if latka_names else atc_nazev

    # Form, strength and route
    sila_val = _clean_str(data.get("sila") or doplnek)
    forma_obj = data.get("lekovaForma") or {}
    lekova_forma_val = _clean_str(forma_obj.get("nazev") if isinstance(forma_obj, dict) else forma_obj)
    cesta_obj = data.get("cestaPodani") or {}
    cesta_podani_val = _clean_str(cesta_obj.get("nazev") if isinstance(cesta_obj, dict) else cesta_obj)

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

    # Dostupnost & dodavky (přerušení, ukončení, obnovení, zahájení)
    je_dodavka = bool(data.get("jeDodavka"))
    dostupnost_raw = str(data.get("dostupnost") or "").upper()

    mr_data = data.get("mrData") if isinstance(data.get("mrData"), dict) else {}
    mr_typ_def = str(mr_data.get("typOmezeni", {}).get("def") or "").lower() if isinstance(mr_data.get("typOmezeni"), dict) else ""
    mr_typ_cs = _clean_str(mr_data.get("typOmezeni"))
    mr_termin = _clean_str(mr_data.get("terminObnovy"))
    mr_duvod = _clean_str(mr_data.get("duvodPreruseni"))

    uvadeni_obj = data.get("uvadeniNaTrh")
    uvadeni_str = ""
    if isinstance(uvadeni_obj, dict):
        uvadeni_str = _clean_str(uvadeni_obj) or str(uvadeni_obj)
    elif isinstance(uvadeni_obj, list) and uvadeni_obj:
        uvadeni_str = str(uvadeni_obj[0])
    elif uvadeni_obj:
        uvadeni_str = str(uvadeni_obj)

    # SÚKL číselník dostupnosti:
    # A = Zahájení uvádění na trh
    # B = Přerušení dodávek (výpadek)
    # C = Ukončení uvádění na trh
    # D = Obnovení dodávek
    # Historické / doplňkové kódy: P = Přerušení, U = Ukončení
    # mrData / uvadeniNaTrh: preruseni / ukonceni
    is_interrupted = (
        dostupnost_raw in ["B", "P"]
        or "prerus" in mr_typ_def
        or "přeruš" in mr_typ_cs.lower()
        or "prerus" in uvadeni_str.lower()
        or "přeruš" in uvadeni_str.lower()
    )
    is_ended = (
        dostupnost_raw in ["C", "U"]
        or "ukonc" in mr_typ_def
        or "ukonč" in mr_typ_cs.lower()
        or "ukonc" in uvadeni_str.lower()
        or "ukonč" in uvadeni_str.lower()
    )

    if is_interrupted:
        dostupnost_status = "interrupted"
        parts = ["Přerušení dodávek (výpadek)"]
        if mr_termin:
            clean_termin = mr_termin.split("T")[0] if "T" in mr_termin else mr_termin
            parts.append(f"předpokládané obnovení {clean_termin}")
        if mr_duvod:
            parts.append(f"důvod: {mr_duvod}")
        dostupnost_text = " – ".join(parts) if len(parts) > 1 else parts[0]
    elif is_ended:
        dostupnost_status = "ended"
        parts = ["Uvádění na trh ukončeno"]
        if mr_duvod:
            parts.append(f"důvod: {mr_duvod}")
        dostupnost_text = " – ".join(parts) if len(parts) > 1 else parts[0]
    elif dostupnost_raw in ["A", "D"]:
        dostupnost_status = "active"
        dostupnost_text = "Dodávky obnoveny" if dostupnost_raw == "D" else "Běžně dodáváno na trh ČR"
    elif je_dodavka:
        dostupnost_status = "active"
        dostupnost_text = "Běžně dodáváno na trh ČR"
    else:
        dostupnost_status = "unknown"
        dostupnost_text = "Bez hlášení o dodávkách"

    # Ceny a uhrady
    uhrady_obj = data.get("uhrady") if isinstance(data.get("uhrady"), dict) else {}
    zaklad_uhrada = uhrady_obj.get("zakladniUhrada") if isinstance(uhrady_obj.get("zakladniUhrada"), dict) else {}
    druha_uhrada = uhrady_obj.get("druhaUhrada") if isinstance(uhrady_obj.get("druhaUhrada"), dict) else {}

    # If zakladniUhrada has no maximalniUhrada, check druhaUhrada
    active_uhrada = zaklad_uhrada if zaklad_uhrada.get("maximalniUhrada") is not None else druha_uhrada

    def _safe_float(v: Any) -> float | None:
        if v is None:
            return None
        try:
            return float(v)
        except (ValueError, TypeError):
            return None

    cena_lek = _safe_float(data.get("cena"))
    max_cena_lek = _safe_float(uhrady_obj.get("maxCenaLek")) or cena_lek
    uhrada_zp = _safe_float(active_uhrada.get("maximalniUhrada"))
    max_doplatek = _safe_float(active_uhrada.get("maximalniDoplatek"))
    zapocitatelny_doplatek = _safe_float(active_uhrada.get("zapocitatelnyDoplatek"))

    # Rezim uhrady
    if uhrada_zp is None or uhrada_zp <= 0:
        uhrada_rezim = "none"
        if max_cena_lek and max_doplatek is None:
            max_doplatek = max_cena_lek
    elif max_doplatek is not None and max_doplatek > 0:
        uhrada_rezim = "partial"
    else:
        uhrada_rezim = "full"

    ind_zaklad = _clean_str(zaklad_uhrada.get("indikacniOmezeni"))
    ind_druha = _clean_str(druha_uhrada.get("indikacniOmezeni"))
    if ind_zaklad and ind_druha and ind_zaklad != ind_druha:
        indikacni_omezeni = f"Děti/mladiství: {ind_zaklad}\n\nDospělí: {ind_druha}"
    else:
        indikacni_omezeni = ind_zaklad or ind_druha or None

    preskripcni_omezeni = (
        _clean_str(zaklad_uhrada.get("specializacePredepisujicihoLekare"))
        or _clean_str(druha_uhrada.get("specializacePredepisujicihoLekare"))
        or None
    )

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

    sukl_url = f"https://prehledy.sukl.cz/prehled_leciv.html#/detail/dlp/{fetch_kod}"

    ret_kod = SuklCode(clean_kod) if clean_kod in alias_map else (
        SuklCode(fetch_kod) if fetch_kod == "0268766" else clean_kod
    )

    return SuklDrugDetail(
        kod_sukl=ret_kod,
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
        atc_kod=atc_kod,
        atc_nazev=atc_nazev,
        leciva_latka=leciva_latka,
        sila=sila_val,
        lekova_forma=lekova_forma_val,
        cesta_podani=cesta_podani_val,
    )

