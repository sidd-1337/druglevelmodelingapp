"""Pull pharmacokinetic statements out of SmPC/SPC text (section 5.2) for manual checking."""
from __future__ import annotations

import re

_SECTION_52 = re.compile(
    r"5\.2\.?\s+(?:Farmakokinetick[éeá]\s+vlastnosti|Pharmacokinetic\s+properties)(?P<body>.*?)"
    r"(?=\n\s*5\.3\.?\s+(?:P[řr]edklinick|Preclinical)|\Z)",
    re.IGNORECASE | re.DOTALL,
)

PK_PATTERNS: dict[str, re.Pattern] = {
    "cmax": re.compile(r"\bc\s?max\b|maxim[áa]ln[íi]\w*\s+(?:plazmatick\w*\s+)?koncentrac|peak\s+(?:plasma\s+)?conc", re.I),
    "tmax": re.compile(r"\bt\s?max\b|time\s+to\s+peak|maxim\w+\s+(?:plazmatick\w+\s+)?koncentrac\w+\s+(?:je\s+)?dosa[žz]",
                      re.I),
    "auc": re.compile(r"\bAUC|plocha\s+pod\s+k[řr]ivkou|area\s+under", re.I),
    "half_life": re.compile(r"polo[čc]as|half-?life|t\s?1/2|t½", re.I),
    "volume": re.compile(r"distribu[čc]n[íi]\w*\s+objem|volume\s+of\s+distribution|\bV\s?d\b|\bVss\b", re.I),
    "bioavailability": re.compile(r"biologick\w+\s+dostupnost|bioavailab", re.I),
    "clearance": re.compile(r"clearance|\bCL\b", re.I),
}

# Sentence ends at ". " followed by an upper-case letter; decimals like 0.5 or 0,5 stay intact.
_SENTENCE_SPLIT = re.compile(r"(?<=[.;])\s+(?=[A-ZÁČĎÉĚÍŇÓŘŠŤÚŮÝŽ])")


def extract_section_52(text: str) -> str | None:
    """Body of section 5.2 (pharmacokinetic properties), or None if not found."""
    m = _SECTION_52.search(text)
    return m.group("body").strip() if m else None


def pk_statements(section: str) -> dict[str, list[str]]:
    """Sentences mentioning each PK quantity and containing a number."""
    flat = re.sub(r"\s+", " ", section)
    sentences = [s.strip() for s in _SENTENCE_SPLIT.split(flat) if re.search(r"\d", s)]
    return {key: [s for s in sentences if pat.search(s)] for key, pat in PK_PATTERNS.items()}
