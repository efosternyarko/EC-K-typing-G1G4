#!/usr/bin/env python3
"""
08_clinker_is_variants.py

Prepare GenBank files and run clinker for IS-variant groups (9 true IS-variant
groups where primary vs variant differ in IS gene family count).

For each group, writes per-locus GenBK slices from v1.2 (so IS annotations
are still present) and runs clinker to produce an HTML comparison.

Usage:
  python3 scripts/08_clinker_is_variants.py [--only-prep]

--only-prep  Write GBK files without running clinker (for debugging)

Requires:
  clinker from /Users/lshef4/lshef4_to_copy/anaconda3/bin/
"""

import csv
import sys
import subprocess
from pathlib import Path
from Bio import SeqIO

REPO = Path(__file__).resolve().parent.parent
IS_DIR = REPO / "DB" / "is_analysis"
GBK_V12 = REPO / "DB" / "EC-K-typing_group1and4_v1.2.gbk"
CLINKER_OUT = REPO / "DB" / "clinker_is_variants"
CLINKER_PY = Path("/Users/lshef4/lshef4_to_copy/anaconda3/bin/python3")

ONLY_PREP = "--only-prep" in sys.argv

IS_VAR_TSV = IS_DIR / "panaroo_is_variant_groups.tsv"

print("Loading IS-variant groups...")
true_is_groups = []
with open(IS_VAR_TSV) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        primary_is = int(row["primary_is_count"])
        var_is = [int(x) for x in row["variant_is_counts"].split(";")]
        # True IS-variant: some member has IS (count > 0)
        if primary_is > 0 or any(v > 0 for v in var_is):
            true_is_groups.append(row)

print(f"  {len(true_is_groups)} true IS-variant groups to visualize")

# Load all v1.2 records indexed by name
print("Loading v1.2 GenBank...")
v12_records = {rec.name: rec for rec in SeqIO.parse(GBK_V12, "genbank")}
print(f"  {len(v12_records)} records")

CLINKER_OUT.mkdir(parents=True, exist_ok=True)
gbk_dir = CLINKER_OUT / "gbk"
gbk_dir.mkdir(exist_ok=True)

# For each group, write GBK files and run clinker
for group in true_is_groups:
    primary = group["primary"]
    variants = group["variants"].split(";")
    members = [primary] + variants

    # Short names for the group
    primary_kl = primary.split("_")[0]
    variant_kls = [v.split("_")[0] for v in variants]
    group_name = f"{primary_kl}_vs_{'_'.join(variant_kls)}"

    print(f"\nGroup: {group_name}")
    print(f"  primary IS count: {group['primary_is_count']}")
    print(f"  variant IS counts: {group['variant_is_counts']}")

    # Write per-member GBK
    gbk_paths = []
    all_found = True
    for member in members:
        if member not in v12_records:
            print(f"  WARNING: {member} not in v1.2 GBK, skipping group")
            all_found = False
            break
        rec = v12_records[member]
        out_gbk = gbk_dir / f"{member}.gbk"
        with open(out_gbk, "w") as fh:
            SeqIO.write(rec, fh, "genbank")
        gbk_paths.append(str(out_gbk))
        kl = member.split("_")[0]
        is_count = group["primary_is_count"] if member == primary else \
            group["variant_is_counts"].split(";")[variants.index(member)]
        n_cds = sum(1 for f in rec.features if f.type == "CDS")
        print(f"  {kl}: {n_cds} CDS, {is_count} IS families")

    if not all_found:
        continue

    if ONLY_PREP:
        print(f"  [prep only] GBK files written")
        continue

    # Run clinker
    html_out = CLINKER_OUT / f"{group_name}.html"
    matrix_out = CLINKER_OUT / f"{group_name}_matrix.tsv"

    print(f"  Running clinker → {html_out.name}")
    cmd = [
        str(CLINKER_PY), "-c",
        f"""
import sys
sys.argv = ['clinker'] + {gbk_paths!r} + [
    '-p', '{html_out}',
    '-mo', '{matrix_out}',
    '-i', '0.3',
    '-j', '2',
    '-f',
]
from clinker.main import main
main()
"""
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-300:]}")
    else:
        print(f"  Done: {html_out}")

print(f"\n{'='*60}")
print(f"CLINKER IS-VARIANT VISUALIZATION")
print(f"{'='*60}")
print(f"  Groups processed: {len(true_is_groups)}")
print(f"  Output directory: {CLINKER_OUT}/")
print(f"  Open HTML files in a browser to view comparisons")
