#!/usr/bin/env python3
"""
04b_panaroo_find_variants.py

Parse Panaroo's gene_presence_absence.csv to find IS-variant and
deletion-variant locus groups.

This replaces the BLASTp-based approach in 04/05 when Panaroo output is
available, providing more reliable gene family assignments (graph-based,
accounts for paralogues and fragmented genes).

IS gene tagging:
  - hmmsearch hits from 02_detect_is_proteins.sh
  - Pattern matching on Panaroo gene/annotation columns

Outputs:
  - panaroo_is_variant_groups.tsv
  - panaroo_deletion_variant_groups.tsv
  - panaroo_loci_to_merge.tsv
  - panaroo_is_families.txt       gene families flagged as IS
  - panaroo_presence_absence_nois.tsv  matrix with IS families removed
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent.parent
IS_DIR = REPO_DIR / "DB" / "is_analysis"
PANAROO_DIR = REPO_DIR / "DB" / "panaroo_out"

PA_CSV      = PANAROO_DIR / "gene_presence_absence.csv"
HMMER_TBL   = IS_DIR / "hmmsearch_is_hits.tbl"

IS_VAR_OUT  = IS_DIR / "panaroo_is_variant_groups.tsv"
DEL_VAR_OUT = IS_DIR / "panaroo_deletion_variant_groups.tsv"
MERGE_OUT   = IS_DIR / "panaroo_loci_to_merge.tsv"
IS_FAM_OUT  = IS_DIR / "panaroo_is_families.txt"
PA_NOIS_OUT     = IS_DIR / "panaroo_presence_absence_nois.tsv"
SMALL_LOCI_OUT  = IS_DIR / "panaroo_small_locus_subsets.tsv"

IS_PATTERN = re.compile(
    r"(IS\d|transpos|tnp[ABCDRSab]|integrase|resolvase|insertion seq)",
    re.IGNORECASE
)

CORE_ASSEMBLY_GENES = {
    "wza", "wzb", "wzc", "wzy", "wzx", "kpsM", "kpsMT",
    "wza_2", "wzb_2", "wzc_2",
}

# Deletion-variant thresholds (tier 2: expanded to named biosynthetic genes)
MAX_DELETION_GENES = 5    # max genes that can be missing for a deletion-variant call
SMALL_LOCUS_THRESHOLD = 15  # ≥ this many genes missing → "small locus subset" class (not auto-merged)


# ── 1. Load Panaroo presence/absence ─────────────────────────────────────────
print("Loading Panaroo gene_presence_absence.csv...")

if not PA_CSV.exists():
    raise FileNotFoundError(
        f"{PA_CSV} not found.\n"
        "Run:\n"
        "  python3 scripts/02b_gbk_to_gff3.py\n"
        "  bash scripts/02c_run_panaroo.sh"
    )

with open(PA_CSV) as fh:
    reader = csv.reader(fh)
    header = next(reader)

# Identify locus columns: everything after fixed metadata columns
META_COLS = {"Gene", "Non-unique Gene name", "Annotation",
             "No. isolates", "No. sequences", "Avg sequences per isolate",
             "Genome Fragment", "Order within Fragment",
             "Accessory Fragment", "Accessory Order with Fragment",
             "QC", "Min group size nuc", "Max group size nuc",
             "Avg group size nuc"}

locus_cols = [c for c in header if c not in META_COLS]
meta_col_idx = {c: i for i, c in enumerate(header) if c in META_COLS}
locus_col_idx = {c: i for i, c in enumerate(header) if c in set(locus_cols)}

print(f"  {len(locus_cols)} loci  ×  N gene families")

rows = []
with open(PA_CSV) as fh:
    reader = csv.reader(fh)
    next(reader)
    for row in reader:
        rows.append(row)

gene_col = meta_col_idx.get("Gene", 0)
annot_col = meta_col_idx.get("Annotation", 2)

gene_families = []       # list of family names (row order)
family_genes = []        # gene name strings
family_annotations = []  # annotation strings

for row in rows:
    gene_families.append(row[gene_col])
    family_genes.append(row[gene_col] if gene_col < len(row) else "")
    family_annotations.append(row[annot_col] if annot_col < len(row) else "")

n_families = len(gene_families)
print(f"  {n_families} gene families in Panaroo output")


# ── 2. Load hmmsearch IS hits to tag IS families ──────────────────────────────
print("Loading hmmsearch IS hits...")
is_hit_proteins = set()

if HMMER_TBL.exists():
    with open(HMMER_TBL) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 5 and float(parts[4]) <= 1e-5:
                is_hit_proteins.add(parts[0])
    print(f"  {len(is_hit_proteins)} IS-flagged proteins from hmmsearch")

# Build set of IS-flagged locus tags from hmmsearch
# protein IDs are: locus|cds_index|gene  →  extract locus + gene for matching
is_locus_tags = set()
for pid in is_hit_proteins:
    parts = pid.split("|")
    if len(parts) >= 3:
        locus, idx, gene = parts[0], parts[1], parts[2]
        kl_prefix = locus.split("_")[0]   # KL387_ESC_... -> KL387
        is_locus_tags.add(f"{kl_prefix}_{int(idx)+1:05d}")

# Tag gene families as IS
is_family_names = set()
with open(PA_CSV) as fh:
    reader = csv.reader(fh)
    next(reader)
    for row in rows:
        fam_name = row[gene_col] if gene_col < len(row) else ""
        annot = row[annot_col] if annot_col < len(row) else ""

        # Pattern-match gene/annotation
        is_by_name = IS_PATTERN.search(fam_name) or IS_PATTERN.search(annot)

        # Check if any member locus_tag matches hmmsearch IS hits
        is_by_hmmer = False
        for locus in locus_cols:
            cell = row[locus_col_idx[locus]] if locus_col_idx[locus] < len(row) else ""
            if cell and cell != "":
                # Panaroo uses locus_tags in the cells
                for tag in cell.split(";"):
                    tag = tag.strip()
                    if tag in is_locus_tags:
                        is_by_hmmer = True
                        break
            if is_by_hmmer:
                break

        if is_by_name or is_by_hmmer:
            is_family_names.add(fam_name)

print(f"  {len(is_family_names)} IS-tagged gene families")

with open(IS_FAM_OUT, "w") as fh:
    for name in sorted(is_family_names):
        fh.write(name + "\n")


# ── 3. Build per-locus gene content sets ─────────────────────────────────────
print("Building per-locus gene content sets...")

locus_biosyn = defaultdict(set)   # locus -> frozenset of biosynthetic family names
locus_is = defaultdict(set)       # locus -> set of IS family names

with open(PA_CSV) as fh:
    reader = csv.reader(fh)
    next(reader)
    for row in rows:
        fam_name = row[gene_col] if gene_col < len(row) else ""
        is_fam = fam_name in is_family_names

        for locus in locus_cols:
            cell = row[locus_col_idx[locus]] if locus_col_idx[locus] < len(row) else ""
            if cell and cell.strip() not in ("", "nan"):
                if is_fam:
                    locus_is[locus].add(fam_name)
                else:
                    locus_biosyn[locus].add(fam_name)

# Freeze biosynthetic profiles for comparison
biosyn_profiles = {l: frozenset(locus_biosyn[l]) for l in locus_cols}


# ── 4. IS counts per locus ───────────────────────────────────────────────────
is_counts = {l: len(locus_is[l]) for l in locus_cols}
n_with_is = sum(1 for c in is_counts.values() if c > 0)
print(f"  {n_with_is} loci carry ≥1 IS gene family")


# ── 5. IS-variant groups ─────────────────────────────────────────────────────
print("Finding IS-variant groups...")

profile_to_loci = defaultdict(list)
for locus, profile in biosyn_profiles.items():
    profile_to_loci[profile].append(locus)

is_variant_groups = []
for profile, group in profile_to_loci.items():
    if len(group) < 2:
        continue
    group_sorted = sorted(group, key=lambda l: is_counts[l])
    primary = group_sorted[0]
    variants = group_sorted[1:]
    is_variant_groups.append({
        "primary": primary,
        "variants": variants,
        "n_biosyn_families": len(profile),
        "primary_is_count": is_counts[primary],
        "variant_is_counts": [is_counts[v] for v in variants],
    })

print(f"  {len(is_variant_groups)} IS-variant groups")
total_is_variants = sum(len(g["variants"]) for g in is_variant_groups)
print(f"  {total_is_variants} loci are IS-variants")

with open(IS_VAR_OUT, "w") as fh:
    fh.write("primary\tvariants\tn_members\tn_biosyn_families\t"
             "primary_is_count\tvariant_is_counts\n")
    for g in sorted(is_variant_groups, key=lambda x: x["primary"]):
        fh.write(
            f"{g['primary']}\t{';'.join(g['variants'])}\t"
            f"{1+len(g['variants'])}\t{g['n_biosyn_families']}\t"
            f"{g['primary_is_count']}\t{';'.join(str(c) for c in g['variant_is_counts'])}\n"
        )


# ── 6. Deletion-variant groups ───────────────────────────────────────────────
# Two-tier detection:
#   Tier 1 (core-assembly): B ⊂ A, all missing genes ⊆ CORE_ASSEMBLY_GENES
#   Tier 2 (biosynthetic):  B ⊂ A, |A-B| ≤ MAX_DELETION_GENES, no IS in diff,
#                           B has ZERO unique biosynthetic genes (pb - pa = ∅)
# Both tiers add to the merge table.  Loci with ≥ SMALL_LOCUS_THRESHOLD missing
# genes (the ~15-17 kb "small-locus" class) are flagged in a separate report
# but NOT auto-merged — they need curator review.
print("Finding deletion-variant groups...")

del_variant_groups = []
small_locus_subsets = []     # flagged but not merged
assigned = set()

loci_by_size = sorted(locus_cols, key=lambda l: -len(biosyn_profiles[l]))

for primary in loci_by_size:
    if primary in assigned:
        continue
    profile_a = biosyn_profiles[primary]
    deletion_vars = []

    for locus_b in loci_by_size:
        if locus_b == primary or locus_b in assigned:
            continue
        profile_b = biosyn_profiles[locus_b]
        if not (profile_b < profile_a):   # B must be strict subset of A
            continue
        diff = profile_a - profile_b       # genes in A but not B
        n_diff = len(diff)
        unique_b = profile_b - profile_a   # genes unique to B (should be ∅ for strict subset)

        # Large gene deficit → "small locus subset" class; flag but don't merge
        if n_diff >= SMALL_LOCUS_THRESHOLD:
            small_locus_subsets.append({
                "locus": locus_b,
                "parent": primary,
                "n_missing": n_diff,
                "n_unique": len(unique_b),
                "missing_named": sorted(
                    g for g in diff if not g.startswith("group_")
                )[:6],
            })
            continue

        # Check for IS genes in the missing set
        is_in_diff = any(IS_PATTERN.search(g) for g in diff if g)
        if is_in_diff:
            continue

        # Tier 1: all missing genes are core assembly genes
        named_diff = {g for g in diff if g and not IS_PATTERN.search(g)}
        if named_diff.issubset(CORE_ASSEMBLY_GENES) or not named_diff:
            deletion_vars.append((locus_b, "core_assembly"))
            assigned.add(locus_b)
            continue

        # Tier 2: small deletion (≤ MAX_DELETION_GENES) of any named/unnamed
        # biosynthetic genes, with B having zero unique genes
        if n_diff <= MAX_DELETION_GENES and len(unique_b) == 0:
            deletion_vars.append((locus_b, "biosynthetic_deletion"))
            assigned.add(locus_b)

    if deletion_vars:
        del_variant_groups.append({
            "primary": primary,
            "deletion_variants": [(v, t) for v, t in deletion_vars],
            "primary_gene_count": len(profile_a),
        })
        assigned.add(primary)

n_core = sum(
    sum(1 for _, t in g["deletion_variants"] if t == "core_assembly")
    for g in del_variant_groups
)
n_biosyn = sum(
    sum(1 for _, t in g["deletion_variants"] if t == "biosynthetic_deletion")
    for g in del_variant_groups
)
print(f"  {len(del_variant_groups)} deletion-variant groups")
print(f"  {n_core} core-assembly deletion-variants")
print(f"  {n_biosyn} biosynthetic deletion-variants")
total_del = n_core + n_biosyn
print(f"  {total_del} loci are deletion-variants (to be removed)")
print(f"  {len(small_locus_subsets)} small-locus subsets (flagged, not auto-merged)")

with open(DEL_VAR_OUT, "w") as fh:
    fh.write("primary\tdeletion_variants\tn_members\tprimary_gene_count\tdeletion_types\n")
    for g in sorted(del_variant_groups, key=lambda x: x["primary"]):
        variants = [v for v, _ in g["deletion_variants"]]
        types = [t for _, t in g["deletion_variants"]]
        fh.write(
            f"{g['primary']}\t{';'.join(variants)}\t"
            f"{1+len(variants)}\t{g['primary_gene_count']}\t"
            f"{';'.join(types)}\n"
        )

with open(SMALL_LOCI_OUT, "w") as fh:
    fh.write("locus\tparent\tn_missing\tn_unique\tmissing_named_sample\n")
    seen = set()
    for s in sorted(small_locus_subsets, key=lambda x: x["locus"]):
        key = (s["locus"], s["parent"])
        if key in seen:
            continue
        seen.add(key)
        fh.write(
            f"{s['locus']}\t{s['parent']}\t{s['n_missing']}\t"
            f"{s['n_unique']}\t{','.join(s['missing_named'])}\n"
        )
print(f"  Small-locus report: {SMALL_LOCI_OUT.name}")


# ── 7. Combined merge table ───────────────────────────────────────────────────
merge_rows = []
for g in is_variant_groups:
    for v in g["variants"]:
        merge_rows.append((v, g["primary"], "IS_variant",
                           is_counts[v], is_counts[g["primary"]]))
for g in del_variant_groups:
    for v, dtype in g["deletion_variants"]:
        merge_rows.append((v, g["primary"], dtype,
                           len(biosyn_profiles[v]), len(biosyn_profiles[g["primary"]])))
merge_rows.sort()

with open(MERGE_OUT, "w") as fh:
    fh.write("locus_to_remove\tprimary_locus\treason\t"
             "variant_metric\tprimary_metric\n")
    for row in merge_rows:
        fh.write("\t".join(str(x) for x in row) + "\n")


# ── 8. Write non-IS presence/absence matrix ──────────────────────────────────
all_biosyn_families = sorted(set(
    fam for profile in biosyn_profiles.values() for fam in profile
))

with open(PA_NOIS_OUT, "w") as fh:
    fh.write("locus\t" + "\t".join(all_biosyn_families) + "\n")
    for locus in locus_cols:
        row = ["1" if f in biosyn_profiles[locus] else "0"
               for f in all_biosyn_families]
        fh.write(locus + "\t" + "\t".join(row) + "\n")


# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("PANAROO VARIANT DETECTION SUMMARY")
print(f"{'='*60}")
print(f"  Total loci analysed:          {len(locus_cols)}")
print(f"  Total gene families:          {n_families}")
print(f"  IS-tagged families:           {len(is_family_names)}")
print(f"  Loci with IS:                 {n_with_is}")
print(f"  IS-variant groups:            {len(is_variant_groups)}")
print(f"  IS-variants to resolve:       {total_is_variants}")
print(f"  Deletion-variant groups:      {len(del_variant_groups)}")
print(f"  Deletion-variants to resolve: {total_del} ({n_core} core-assembly, {n_biosyn} biosynthetic)")
print(f"  Small-locus subsets (flagged):{len(set(s['locus'] for s in small_locus_subsets))}")
print(f"  Total loci to merge/remove:   {len(merge_rows)}")
print(f"\nOutputs in {IS_DIR}/")
for p in [IS_VAR_OUT, DEL_VAR_OUT, MERGE_OUT, IS_FAM_OUT, PA_NOIS_OUT, SMALL_LOCI_OUT]:
    print(f"  {p.name}")
