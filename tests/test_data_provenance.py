"""A value marked verified must say where it was verified - and that place must be primary.

"verified": true means the value was checked against the primary document: an SPC or
regulator label, the original study's full text (or its abstract, where the note says
so), or a registry's posted results. Drug-information aggregators and label mirrors
restate those documents and were the source of several wrong values here, so a
verified entry may not cite one.
"""
import json
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent.parent / "data"
SECONDARY = ("drugs.com", "medlibrary.org", "rxreasoner", "drugbank.com", "ebmconsult", "patientcareonline",
             "sigmaaldrich", "uptodate.com", "citedrive", "erowid", "kobv.de", "drugcentral")


def _flagged(node, path=""):
    """Every dict carrying a 'verified' key, with a readable path."""
    if isinstance(node, dict):
        if "verified" in node:
            yield path, node
        for key, value in node.items():
            yield from _flagged(value, f"{path}.{key}" if path else key)
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _flagged(value, f"{path}[{i}]")


ENTRIES = [(f"{name}:{path}", entry)
           for name in ("reference_pk.json", "drugs.json", "enzymes.json", "interactions.json")
           for path, entry in _flagged(json.loads((DATA / name).read_text(encoding="utf-8")))]


def test_flags_are_found():
    assert len(ENTRIES) >= 100


@pytest.mark.parametrize("where, entry", ENTRIES, ids=[w for w, _ in ENTRIES])
def test_verified_entries_cite_a_primary_document(where, entry):
    assert isinstance(entry["verified"], bool)
    if entry["verified"]:
        assert entry.get("source", "").startswith("http"), where
        assert not any(domain in entry["source"] for domain in SECONDARY), f"{where} cites {entry['source']}"
        assert entry.get("note"), f"{where} is verified but its note does not say against what"
