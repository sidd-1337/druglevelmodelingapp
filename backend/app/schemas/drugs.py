from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class SuklSearchResultItem(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    kod_sukl: str
    nazev: str
    doplnek: str
    sila: str
    lekova_forma: str
    cesta_podani: str
    atc_kod: str
    atc_nazev: str
    parsed_dose_mg: float | None = None
    matched_preset_id: str | None = None
    sukl_url: str


class DrugDetailResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    kod_sukl: str
    nazev: str
    doplnek: str
    registracni_cislo: str
    drzitel: str
    zpusob_vydeje_kod: str
    zpusob_vydeje_text: str
    je_dodavka: bool
    dostupnost_status: str
    dostupnost_text: str
    cena_lek: float | None = None
    max_cena_lek: float | None = None
    uhrada_zp: float | None = None
    max_doplatek: float | None = None
    zapocitatelny_doplatek: float | None = None
    uhrada_rezim: str
    indikacni_omezeni: str | None = None
    preskripcni_omezeni: str | None = None
    spc_url: str | None = None
    pil_url: str | None = None
    sukl_url: str
    atc_kod: str = ""
    atc_nazev: str = ""
    leciva_latka: str = ""
    sila: str = ""
    lekova_forma: str = ""
    cesta_podani: str = ""
    matched_preset_id: str | None = None
    matched_preset: dict | None = None


class AlternativeDrugItem(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    kod_sukl: str
    nazev: str
    doplnek: str
    sila: str
    lekova_forma: str
    cesta_podani: str
    atc_kod: str
    atc_nazev: str
    parsed_dose_mg: float | None = None
    matched_preset_id: str | None = None
    sukl_url: str
