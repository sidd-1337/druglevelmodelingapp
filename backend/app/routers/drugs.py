from __future__ import annotations

import dataclasses
import unicodedata
from fastapi import APIRouter, HTTPException, Query

from backend.app.schemas.drugs import (
    AlternativeDrugItem,
    DrugDetailResponse,
    SuklSearchResultItem,
)
from core.drug_database import load_drugs, get_drug
from core.sukl_client import (
    PRESET_DEFAULT_SUKL_CODES,
    fetch_sukl_detail,
    query_sukl_api,
    match_preset_id,
    SuklMedication,
    SuklDrugDetail,
)

router = APIRouter()


def _strip_accents(text: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFD", text.lower()) if unicodedata.category(c) != "Mn")


@router.get("/search", response_model=list[SuklSearchResultItem])
def search_drugs(
    q: str = Query(..., min_length=1, description="Search query for SÚKL registry or clinical presets"),
    limit: int = Query(20, ge=1, le=100, description="Maximum number of results to return"),
) -> list[SuklSearchResultItem]:
    """Search registered medications in the Czech SÚKL registry and clinical presets."""
    query_str = q.strip()
    results: list[SuklSearchResultItem] = []

    # 1. Query official SÚKL registry
    try:
        sukl_results = query_sukl_api(query_str, max_results=limit)
        for med in sukl_results:
            results.append(
                SuklSearchResultItem(
                    kod_sukl=str(med.kod_sukl),
                    nazev=med.nazev,
                    doplnek=med.doplnek,
                    sila=med.sila,
                    lekova_forma=med.lekova_forma,
                    cesta_podani=med.cesta_podani,
                    atc_kod=med.atc_kod,
                    atc_nazev=med.atc_nazev,
                    parsed_dose_mg=med.parsed_dose_mg,
                    matched_preset_id=med.matched_preset_id,
                    sukl_url=med.sukl_url,
                )
            )
    except Exception:
        pass

    # 2. If results are few or empty, match against internal clinical presets
    if len(results) < limit:
        q_norm = _strip_accents(query_str)
        for drug in load_drugs():
            d_norm = _strip_accents(f"{drug.id} {drug.name} {drug.category}")
            if q_norm in d_norm:
                # Avoid duplicate if already found in SÚKL
                if any(r.matched_preset_id == drug.id for r in results):
                    continue
                results.append(
                    SuklSearchResultItem(
                        kod_sukl=f"PRESET_{drug.id.upper()}",
                        nazev=drug.name,
                        doplnek=f"{drug.category} (Clinical preset)",
                        sila=f"{drug.default_dose_mg:g} {drug.dose_basis}",
                        lekova_forma="POR TBL",
                        cesta_podani="Perorální podání",
                        atc_kod="",
                        atc_nazev=drug.name,
                        parsed_dose_mg=drug.default_dose_mg,
                        matched_preset_id=drug.id,
                        sukl_url="https://prehledy.sukl.cz/",
                    )
                )
            if len(results) >= limit:
                break

    return results


@router.get("/{kod_sukl}", response_model=DrugDetailResponse)
def get_drug_by_code(kod_sukl: str) -> DrugDetailResponse:
    """Retrieve full regulatory, pricing, reimbursement, and shortage status for a SÚKL code."""
    clean_code = kod_sukl.strip()

    # 1. Check if SÚKL detail exists
    detail: SuklDrugDetail | None = None
    try:
        detail = fetch_sukl_detail(clean_code)
    except Exception:
        detail = None

    if detail is not None:
        matched_preset_dict = None
        matched_preset_id = None
        # Try matching active substance or ATC to clinical preset
        found_pid = match_preset_id(detail.atc_nazev, detail.nazev) or match_preset_id(detail.leciva_latka, detail.nazev)
        if not found_pid:
            combined_low = f"{detail.nazev} {detail.atc_nazev} {detail.leciva_latka} {detail.doplnek}".lower()
            for drug in load_drugs():
                if drug.id in combined_low or drug.name.lower() in combined_low or drug.id in clean_code.lower():
                    found_pid = drug.id
                    break
        if found_pid:
            try:
                drug_obj = get_drug(found_pid)
                matched_preset_dict = dataclasses.asdict(drug_obj)
                matched_preset_id = drug_obj.id
            except KeyError:
                pass

        return DrugDetailResponse(
            kod_sukl=str(detail.kod_sukl),
            nazev=detail.nazev,
            doplnek=detail.doplnek,
            registracni_cislo=detail.registracni_cislo,
            drzitel=detail.drzitel,
            zpusob_vydeje_kod=detail.zpusob_vydeje_kod,
            zpusob_vydeje_text=detail.zpusob_vydeje_text,
            je_dodavka=detail.je_dodavka,
            dostupnost_status=detail.dostupnost_status,
            dostupnost_text=detail.dostupnost_text,
            cena_lek=detail.cena_lek,
            max_cena_lek=detail.max_cena_lek,
            uhrada_zp=detail.uhrada_zp,
            max_doplatek=detail.max_doplatek,
            zapocitatelny_doplatek=detail.zapocitatelny_doplatek,
            uhrada_rezim=detail.uhrada_rezim,
            indikacni_omezeni=detail.indikacni_omezeni,
            preskripcni_omezeni=detail.preskripcni_omezeni,
            spc_url=detail.spc_url,
            pil_url=detail.pil_url,
            sukl_url=detail.sukl_url,
            atc_kod=detail.atc_kod,
            atc_nazev=detail.atc_nazev,
            leciva_latka=detail.leciva_latka,
            sila=detail.sila,
            lekova_forma=detail.lekova_forma,
            cesta_podani=detail.cesta_podani,
            matched_preset_id=matched_preset_id,
            matched_preset=matched_preset_dict,
        )

    # 2. Check if clean_code matches a preset ID directly (e.g. "ritalin", "PRESET_RITALIN", "ibuprofen") or default SÚKL code
    preset_id = clean_code.lower().removeprefix("preset_")
    for pid, code in PRESET_DEFAULT_SUKL_CODES.items():
        if clean_code == str(code):
            preset_id = pid
            break
    try:
        drug = get_drug(preset_id)
        return DrugDetailResponse(
            kod_sukl=clean_code,
            nazev=drug.name,
            doplnek=f"{drug.category} (Clinical preset)",
            registracni_cislo="CLINICAL/PRESET",
            drzitel="National Formulary / AGNP",
            zpusob_vydeje_kod="rx",
            zpusob_vydeje_text="Na lékařský předpis (Rx)",
            je_dodavka=True,
            dostupnost_status="active",
            dostupnost_text="Dostupný klinický model",
            cena_lek=None,
            max_cena_lek=None,
            uhrada_zp=None,
            max_doplatek=None,
            zapocitatelny_doplatek=None,
            uhrada_rezim="none",
            indikacni_omezeni=None,
            preskripcni_omezeni=None,
            spc_url=None,
            pil_url=None,
            sukl_url=f"https://prehledy.sukl.cz/",
            atc_kod="",
            atc_nazev=drug.name,
            leciva_latka=drug.name,
            sila=f"{drug.default_dose_mg:g} {drug.dose_basis}",
            lekova_forma="POR TBL",
            cesta_podani="Perorální podání",
            matched_preset_id=drug.id,
            matched_preset=dataclasses.asdict(drug),
        )
    except KeyError:
        pass

    raise HTTPException(status_code=404, detail=f"Medication with code {kod_sukl!r} not found in SÚKL or presets.")


@router.get("/{kod_sukl}/alternatives", response_model=list[AlternativeDrugItem])
def get_drug_alternatives(
    kod_sukl: str,
    limit: int = Query(20, ge=1, le=100, description="Max alternative medications to return"),
) -> list[AlternativeDrugItem]:
    """Retrieve therapeutic and generic alternatives (same active substance / ATC group) for a medication."""
    clean_code = kod_sukl.strip()
    alternatives: list[AlternativeDrugItem] = []

    # 1. Retrieve the base drug detail to obtain name or ATC code
    base_name = ""
    base_preset_id = ""
    try:
        detail = fetch_sukl_detail(clean_code)
        if detail:
            # First word of the name is usually the active substance or brand root
            base_name = detail.nazev.split()[0]
    except Exception:
        pass

    if not base_name:
        preset_id = clean_code.lower().removeprefix("preset_")
        for pid, code in PRESET_DEFAULT_SUKL_CODES.items():
            if clean_code == str(code):
                preset_id = pid
                break
        try:
            drug = get_drug(preset_id)
            base_name = drug.name.split()[0]
            base_preset_id = drug.id
        except KeyError:
            pass

    if not base_name:
        base_name = clean_code

    # 2. Query SÚKL for the same substance
    try:
        matching = query_sukl_api(base_name, max_results=limit + 5)
        for med in matching:
            med_code = str(med.kod_sukl)
            if med_code == clean_code:
                continue
            alternatives.append(
                AlternativeDrugItem(
                    kod_sukl=med_code,
                    nazev=med.nazev,
                    doplnek=med.doplnek,
                    sila=med.sila,
                    lekova_forma=med.lekova_forma,
                    cesta_podani=med.cesta_podani,
                    atc_kod=med.atc_kod,
                    atc_nazev=med.atc_nazev,
                    parsed_dose_mg=med.parsed_dose_mg,
                    matched_preset_id=med.matched_preset_id,
                    sukl_url=med.sukl_url,
                )
            )
            if len(alternatives) >= limit:
                break
    except Exception:
        pass

    # 3. If no external results, check presets in the same category
    if not alternatives and base_preset_id:
        try:
            base_drug = get_drug(base_preset_id)
            for d in load_drugs():
                if d.id != base_drug.id and (d.category == base_drug.category or base_name.lower() in d.name.lower()):
                    alternatives.append(
                        AlternativeDrugItem(
                            kod_sukl=f"PRESET_{d.id.upper()}",
                            nazev=d.name,
                            doplnek=d.category,
                            sila=f"{d.default_dose_mg:g} mg",
                            lekova_forma="POR TBL",
                            cesta_podani="Perorální podání",
                            atc_kod="",
                            atc_nazev=d.name,
                            parsed_dose_mg=d.default_dose_mg,
                            matched_preset_id=d.id,
                            sukl_url="https://prehledy.sukl.cz/",
                        )
                    )
        except Exception:
            pass

    return alternatives
