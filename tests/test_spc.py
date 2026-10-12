from core.spc import extract_section_52, pk_statements

CZ_SPC = """
5.1 Farmakodynamické vlastnosti
Pregabalin je analog GABA.
5.2 Farmakokinetické vlastnosti
Absorpce
Pregabalin se rychle vstřebává. Maximálních plazmatických koncentrací je dosaženo během 1 hodiny.
Biologická dostupnost pregabalinu je ≥ 90 % a nezávisí na dávce. Zdánlivý distribuční objem je
přibližně 0,56 l/kg. Průměrný eliminační poločas je 6,3 hodiny. Clearance pregabalinu je přímo
úměrná clearance kreatininu.
5.3 Předklinické údaje vztahující se k bezpečnosti
Toxicita 10 mg/kg.
"""

EN_SPC = """5.2 Pharmacokinetic properties
Peak plasma concentrations (Cmax 3.9 ug/mL) occur within 1.5 hours (Tmax). The AUC was 27.9 ug.h/mL.
5.3 Preclinical safety data
"""


def test_extracts_czech_section_and_classifies_sentences():
    section = extract_section_52(CZ_SPC)
    assert section is not None and "Toxicita" not in section and "GABA" not in section
    st = pk_statements(section)
    assert any("0,56 l/kg" in s for s in st["volume"])
    assert any("6,3 hodiny" in s for s in st["half_life"])
    assert any("90 %" in s for s in st["bioavailability"])
    assert any("1 hodiny" in s for s in st["tmax"])


def test_extracts_english_section():
    st = pk_statements(extract_section_52(EN_SPC))
    assert st["cmax"] and st["auc"] and st["tmax"]


def test_missing_section_returns_none():
    assert extract_section_52("4.8 Nežádoucí účinky\nnic") is None
