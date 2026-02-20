#!/usr/bin/env python3
"""
add_ktype_notes.py — Add 'K type:KLxxx' notes to G1/G4 source features.

Kaptive reads the best-match *type* (the 'Best match type' column) from
the /note="K type:XXX" qualifier in each GenBank source feature, using the
regex r'(?<=type:)\\w+|(?<=type: ).*'.  Without this note, Kaptive writes
'unknown (KLxxx)' for every locus it finds.

The Gladstone Group 2 & 3 records already include 'K type:Kx' notes, so
they report clean type names.  Our G1/G4 records (KL300–KL423, K24, K96)
currently lack this note, hence the 'unknown' prefix.

This script:
  1. Reads the v0.5 G1/G4 GenBank.
  2. For each record that does not already have a 'K type:' note, inserts
       /note="K type:{locus_name}"
     immediately after the existing 'K locus:' note.
  3. Rewrites the G1/G4 GenBank in-place.
  4. Rewrites the combined all-groups GenBank in-place (only G1/G4 records
     are modified; G2/G3 records already have type annotations).
"""
import re
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
REPO_DIR = Path(__file__).resolve().parent.parent
DB_DIR   = REPO_DIR / "DB"
G1G4_GBK = DB_DIR / "EC-K-typing_group1and4_v0.5.gbk"
ALL_GBK  = DB_DIR / "EC-K-typing_all_groups_v0.5.gbk"
# ---------------------------------------------------------------------------


_LOCUS_RE = re.compile(r"K locus:\s*(\S+)")
_TYPE_RE  = re.compile(r"K type:", re.IGNORECASE)


def get_locus_name(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = _LOCUS_RE.search(note)
                if m:
                    return m.group(1)
    # Fallback: first underscore-delimited token of the record name
    return rec.name.split("_")[0]


def has_type_note(rec) -> bool:
    for f in rec.features:
        if f.type == "source":
            return any(_TYPE_RE.search(n) for n in f.qualifiers.get("note", []))
    return False


def insert_type_note(rec, locus_name: str) -> bool:
    """
    Insert 'K type:{locus_name}' immediately after the 'K locus:' note.
    Returns True if the record was modified.
    """
    for f in rec.features:
        if f.type == "source":
            notes     = f.qualifiers.get("note", [])
            new_notes = []
            inserted  = False
            for note in notes:
                new_notes.append(note)
                if _LOCUS_RE.search(note) and not inserted:
                    new_notes.append(f"K type:{locus_name}")
                    inserted = True
            if not inserted:
                new_notes.append(f"K type:{locus_name}")
            f.qualifiers["note"] = new_notes
            return True
    return False


def annotate_records(records: list, label: str) -> int:
    n = 0
    for rec in records:
        if not has_type_note(rec):
            locus = get_locus_name(rec)
            insert_type_note(rec, locus)
            print(f"  [{label}] Added K type:{locus}  ({rec.name})")
            n += 1
    return n


# ---------------------------------------------------------------------------
# G1/G4 database
# ---------------------------------------------------------------------------
print(f"\n[1] Processing G1/G4 database: {G1G4_GBK.name}")
recs_g1g4 = list(SeqIO.parse(str(G1G4_GBK), "genbank"))
n_mod = annotate_records(recs_g1g4, "G1G4")
print(f"    Modified {n_mod}/{len(recs_g1g4)} records.")
SeqIO.write(recs_g1g4, str(G1G4_GBK), "genbank")
print(f"    Written: {G1G4_GBK.name}")

# ---------------------------------------------------------------------------
# Combined all-groups database
# ---------------------------------------------------------------------------
print(f"\n[2] Processing combined database: {ALL_GBK.name}")
recs_all = list(SeqIO.parse(str(ALL_GBK), "genbank"))
n_mod_all = annotate_records(recs_all, "ALL")
print(f"    Modified {n_mod_all}/{len(recs_all)} records "
      f"({len(recs_all)-n_mod_all} G2/G3 records already annotated).")
SeqIO.write(recs_all, str(ALL_GBK), "genbank")
print(f"    Written: {ALL_GBK.name}")

print("\nDone. Run kaptive to verify 'Best match type' now shows KL names directly.")
