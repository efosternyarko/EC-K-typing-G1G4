#!/usr/bin/env python3
"""
06_build_v1.3.py

Build the v1.3 database from v1.2 by:
  1. Removing IS-variant and deletion-variant loci (keeping only primary references)
  2. Stripping IS CDS annotations from all retained loci
  3. Writing the merged GenBank + updated mapping TSV

Usage:
  python3 scripts/06_build_v1.3.py [--dry-run]

With --dry-run: print what would be removed without writing anything.
"""

import csv
import re
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

REPO_DIR = Path(__file__).resolve().parent.parent
IS_DIR = REPO_DIR / "DB" / "is_analysis"
DB_DIR = REPO_DIR / "DB"

GBK_IN      = DB_DIR / "EC-K-typing_group1and4_v1.2.gbk"
GBK_OUT     = DB_DIR / "EC-K-typing_group1and4_v1.3.gbk"
MAPPING_IN  = DB_DIR / "KL_G1G4_mapping.tsv"
MAPPING_OUT = DB_DIR / "KL_G1G4_mapping_v1.3.tsv"

import os
_merge_env = os.environ.get("PANAROO_MERGE")
MERGE_TSV = Path(_merge_env) if _merge_env else IS_DIR / "panaroo_loci_to_merge.tsv"
if not MERGE_TSV.exists():
    MERGE_TSV = IS_DIR / "loci_to_merge.tsv"   # BLASTp fallback
IS_IDS_TXT  = IS_DIR / "is_protein_ids.txt"
FAM_TSV     = IS_DIR / "protein_families.tsv"

IS_PATTERN = re.compile(
    r"(IS\d|transpos|tnp[ABCDRSab]|integrase|resolvase|insertion seq)", re.IGNORECASE
)

DRY_RUN = "--dry-run" in sys.argv

# ── 1. Load merge decisions ───────────────────────────────────────────────────
print("Loading variant merge decisions...")
loci_to_remove = {}  # locus_name -> (primary, reason)

if MERGE_TSV.exists():
    with open(MERGE_TSV) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            loci_to_remove[row["locus_to_remove"]] = (row["primary_locus"], row["reason"])
    print(f"  {len(loci_to_remove)} loci flagged for removal")
else:
    print(f"  WARNING: {MERGE_TSV} not found — no loci will be removed")

# ── 2. Load IS protein IDs and IS-tagged protein families ────────────────────
print("Loading IS-tagged protein families...")
is_families = set()

if FAM_TSV.exists():
    with open(FAM_TSV) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["is_family"] == "True":
                is_families.add(int(row["family_id"]))
    print(f"  {len(is_families)} IS-tagged protein families")

# Also load per-protein IS flags for direct CDS matching
is_protein_ids = set()
if IS_IDS_TXT.exists():
    with open(IS_IDS_TXT) as fh:
        is_protein_ids = {l.strip() for l in fh if l.strip()}
    print(f"  {len(is_protein_ids)} individual IS-flagged proteins")

# Helper: is a CDS feature IS-related?
def is_cds_is_related(feat, locus_name, cds_idx):
    gene = feat.qualifiers.get("gene", ["?"])[0]
    product = feat.qualifiers.get("product", [""])[0]
    # Pattern match on gene/product
    if IS_PATTERN.search(gene) or IS_PATTERN.search(product):
        return True
    # Check against hmmsearch-flagged protein IDs
    pid = f"{locus_name}|{cds_idx}|{gene}"
    if pid in is_protein_ids:
        return True
    return False

# ── 3. Process each GenBank record ───────────────────────────────────────────
print("\nProcessing GenBank records...")
records_in = list(SeqIO.parse(GBK_IN, "genbank"))
print(f"  Loaded {len(records_in)} records from v1.2")

records_out = []
n_removed_loci = 0
n_stripped_is_cds = 0
n_retained = 0

removal_log = []
stripping_log = []

for rec in records_in:
    locus = rec.name

    # Check if this locus is flagged for removal
    if locus in loci_to_remove:
        primary, reason = loci_to_remove[locus]
        removal_log.append(f"{locus}\t{primary}\t{reason}")
        n_removed_loci += 1
        continue

    # Strip IS CDS from retained loci
    new_features = []
    cds_idx = 0
    n_is_in_locus = 0
    for feat in rec.features:
        if feat.type != "CDS":
            new_features.append(feat)
            continue
        if is_cds_is_related(feat, locus, cds_idx):
            n_is_in_locus += 1
            n_stripped_is_cds += 1
        else:
            new_features.append(feat)
        cds_idx += 1

    if n_is_in_locus > 0:
        stripping_log.append(f"{locus}\t{n_is_in_locus}")

    rec.features = new_features
    records_out.append(rec)
    n_retained += 1

print(f"\n  Loci removed (IS/deletion variants): {n_removed_loci}")
print(f"  Loci retained:                       {n_retained}")
print(f"  IS CDS annotations stripped:         {n_stripped_is_cds}")

# ── 4. Write outputs ──────────────────────────────────────────────────────────
if DRY_RUN:
    print("\n[DRY RUN] Would write:")
    print(f"  {GBK_OUT}  ({n_retained} records)")
    print(f"\nRemoval log preview (first 20):")
    for line in removal_log[:20]:
        print(f"  {line}")
    print(f"\nIS-stripping log preview (first 20):")
    for line in stripping_log[:20]:
        print(f"  {line}")
    sys.exit(0)

print(f"\nWriting {GBK_OUT}...")
with open(GBK_OUT, "w") as fh:
    SeqIO.write(records_out, fh, "genbank")
print(f"  Wrote {n_retained} records")

# Write removal log
log_path = IS_DIR / "v1.3_removal_log.tsv"
with open(log_path, "w") as fh:
    fh.write("removed_locus\tprimary_locus\treason\n")
    fh.write("\n".join(removal_log) + "\n")
print(f"  Removal log: {log_path}")

# Write IS stripping log
strip_log_path = IS_DIR / "v1.3_is_strip_log.tsv"
with open(strip_log_path, "w") as fh:
    fh.write("locus\tn_is_cds_stripped\n")
    fh.write("\n".join(stripping_log) + "\n")
print(f"  IS stripping log: {strip_log_path}")

# ── 5. Update mapping file ────────────────────────────────────────────────────
if MAPPING_IN.exists():
    print(f"\nUpdating mapping file...")
    kept_rows = []
    with open(MAPPING_IN) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            kl = row.get("KL", list(row.values())[0])
            # Check if any record in records_out matches this KL
            if any(kl in r.name for r in records_out):
                kept_rows.append(row)

    with open(MAPPING_OUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept_rows)
    print(f"  {len(kept_rows)} mapping entries retained → {MAPPING_OUT}")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("v1.3 BUILD COMPLETE")
print(f"{'='*60}")
print(f"  Input loci (v1.2):  {len(records_in)}")
print(f"  Output loci (v1.3): {n_retained}")
print(f"  Loci removed:       {n_removed_loci}")
print(f"  IS CDS stripped:    {n_stripped_is_cds}")
print(f"  Output GenBank:     {GBK_OUT}")
