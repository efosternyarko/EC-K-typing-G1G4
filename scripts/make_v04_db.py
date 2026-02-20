#!/usr/bin/env python3
"""
make_v04_db.py — Build v0.4 database by stripping conserved CDS from G1/G4 loci.

The 5 remaining self-typing failures in v0.3.1 (KL300, KL301, KL303, KL306, KL307)
are all KX01 loci that share identical conserved flanking/export genes with KL302.
These conserved genes (wza, wzb, wzc, galF, gnd, ugd) contribute equal bitscore to
every KX01 reference, adding shared background noise that partially masks the variable
region signal.

Fix: remove conserved CDS features from all 93 G1/G4 GenBank records, leaving only
the variable biosynthetic region. Kaptive then scores exclusively on discriminating genes.

Usage
-----
    python scripts/make_v04_db.py

Outputs (written to DB/)
------------------------
    EC-K-typing_group1and4_v0.4.gbk     — G1/G4 only, conserved genes stripped
    EC-K-typing_all_groups_v0.4.gbk     — combined (G2/G3 unchanged + stripped G1/G4)
"""

import re
import sys
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR = Path(__file__).resolve().parent.parent
DB_DIR   = REPO_DIR / "DB"

INPUT_G1G4  = DB_DIR / "EC-K-typing_group1and4_v0.3.1.gbk"
INPUT_G2G3  = DB_DIR / "EC-K-typing_group2and3_v3.0.0.gbk"
OUTPUT_G1G4 = DB_DIR / "EC-K-typing_group1and4_v0.4.gbk"
OUTPUT_ALL  = DB_DIR / "EC-K-typing_all_groups_v0.4.gbk"

# ---------------------------------------------------------------------------
# Conserved genes to strip
# ---------------------------------------------------------------------------
CONSERVED_GENE_NAMES = {'galF', 'galF_2', 'gnd', 'ugd', 'wza', 'wzb', 'wzc'}


def get_locus_name(rec) -> str:
    """Extract K locus name from GenBank source feature note."""
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def strip_conserved_cds(rec):
    """Remove CDS features whose gene name is in CONSERVED_GENE_NAMES."""
    new_features = [
        f for f in rec.features
        if not (f.type == "CDS" and
                f.qualifiers.get("gene", [""])[0] in CONSERVED_GENE_NAMES)
    ]
    rec.features = new_features
    return rec


def main():
    # Verify inputs exist
    for path in (INPUT_G1G4, INPUT_G2G3):
        if not path.exists():
            print(f"ERROR: input not found: {path}", file=sys.stderr)
            sys.exit(1)

    print("=" * 65)
    print("make_v04_db.py — Strip conserved CDS from G1/G4 loci")
    print("=" * 65)
    print(f"\nConserved genes to remove: {sorted(CONSERVED_GENE_NAMES)}")
    print(f"\nInput  G1/G4: {INPUT_G1G4.name}")
    print(f"Input  G2/G3: {INPUT_G2G3.name}")

    # ------------------------------------------------------------------
    # Process G1/G4 records
    # ------------------------------------------------------------------
    print("\n" + "-" * 65)
    print(f"{'Locus':<12}  {'Before':>6}  {'Removed':>7}  {'After':>5}")
    print("-" * 65)

    stripped_records = []
    total_before = 0
    total_removed = 0

    for rec in SeqIO.parse(str(INPUT_G1G4), "genbank"):
        lname = get_locus_name(rec)
        before = sum(1 for f in rec.features if f.type == "CDS")
        rec    = strip_conserved_cds(rec)
        after  = sum(1 for f in rec.features if f.type == "CDS")
        removed = before - after
        total_before  += before
        total_removed += removed
        stripped_records.append(rec)
        print(f"{lname:<12}  {before:>6}  {removed:>7}  {after:>5}")

    print("-" * 65)
    print(f"{'TOTAL':<12}  {total_before:>6}  {total_removed:>7}  {total_before-total_removed:>5}")
    print(f"\nExpected ~528 removed (16.0% of 3,304 CDS across 93 loci)")

    # ------------------------------------------------------------------
    # Write G1/G4-only output
    # ------------------------------------------------------------------
    SeqIO.write(stripped_records, str(OUTPUT_G1G4), "genbank")
    print(f"\nWritten: {OUTPUT_G1G4.name}  ({len(stripped_records)} loci)")

    # ------------------------------------------------------------------
    # Load G2/G3 records (unchanged)
    # ------------------------------------------------------------------
    g2g3_records = list(SeqIO.parse(str(INPUT_G2G3), "genbank"))
    print(f"Loaded:  {INPUT_G2G3.name}  ({len(g2g3_records)} loci, unchanged)")

    # ------------------------------------------------------------------
    # Write combined all-groups output
    # ------------------------------------------------------------------
    all_records = g2g3_records + stripped_records
    SeqIO.write(all_records, str(OUTPUT_ALL), "genbank")
    print(f"Written: {OUTPUT_ALL.name}  ({len(all_records)} loci total)")

    print("\nDone. Use DB/EC-K-typing_all_groups_v0.4.gbk with Kaptive or")
    print("type_normalized.py --db DB/EC-K-typing_all_groups_v0.4.gbk --suffix v0.4norm")


if __name__ == "__main__":
    main()
