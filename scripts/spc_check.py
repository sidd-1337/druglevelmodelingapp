"""Download each preset's SÚKL SPC and list its section 5.2 PK statements next to our reference values.

Usage:  pip install pypdf && python scripts/spc_check.py [report.md]

Needs network access to prehledy.sukl.cz and www.sukl.cz. The report is a checklist
for a person: confirm or correct each value in data/reference_pk.json and set
"verified": true once it matches the SPC.
"""
from __future__ import annotations

import io
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.spc import extract_section_52, pk_statements  # noqa: E402
from core.sukl_client import PRESET_DEFAULT_SUKL_CODES, USER_AGENT, fetch_sukl_detail  # noqa: E402
from core.validation import load_references  # noqa: E402


def _pdf_text(url: str) -> str:
    from pypdf import PdfReader  # optional dependency, only needed for this script

    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=30) as resp:
        reader = PdfReader(io.BytesIO(resp.read()))
    return "\n".join(page.extract_text() or "" for page in reader.pages)


def main(out_path: str = "spc_check_report.md") -> None:
    refs = defaultdict(list)
    for r in load_references():
        refs[r["drug_id"]].append(r)
    lines = ["# SPC section 5.2 check", ""]
    for drug_id, code in PRESET_DEFAULT_SUKL_CODES.items():
        lines.append(f"## {drug_id} (SÚKL {code})")
        detail = fetch_sukl_detail(code)
        if detail is None or not detail.spc_url:
            lines += ["SPC not found (network blocked or no SPC listed).", ""]
            continue
        lines.append(f"Product: {detail.nazev} {detail.doplnek} - [SPC]({detail.spc_url})")
        try:
            section = extract_section_52(_pdf_text(detail.spc_url))
        except Exception as exc:  # network or PDF parsing failure: report and continue
            lines += [f"Download/parse failed: {exc}", ""]
            continue
        if not section:
            lines += ["Section 5.2 not found in the PDF text.", ""]
            continue
        lines += ["", "Our reference values:"]
        lines += [f"- {r['metric']} {r['dose_mg']} mg: {r['value']} ({r['source']})" for r in refs.get(drug_id, [])]
        for key, sentences in pk_statements(section).items():
            if sentences:
                lines += ["", f"**{key}**"] + [f"- {s}" for s in sentences]
        lines.append("")
    Path(out_path).write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main(*sys.argv[1:])
