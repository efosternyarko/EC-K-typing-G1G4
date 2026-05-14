#!/usr/bin/env python3
"""
05_find_variants.py

Identify IS-variant and deletion-variant locus groups from the presence/absence matrix.

IS-variant definition (Kat Holt):
  Two loci are IS-variants if their non-IS biosynthetic gene content is identical.
  Primary reference = the member with fewest IS families (preferring 0).

Deletion-variant definition:
  Locus B is a deletion-variant of locus A if:
    - B's non-IS gene set is a proper subset of A's
    - The genes in A but not B are ALL core assembly genes (wza/wzb/wzc/wzy/wzx)

Outputs:
  - is_variant_groups.tsv       one row per group (primary + IS-variant members)
  - deletion_variant_groups.tsv one row per group (full-length primary + deletion variants)
  - loci_to_merge.tsv           all loci flagged for merging with their primary + reason
  - is_counts_per_locus.tsv     number of IS families per locus
"""

import csv
from collections import defaultdict
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent.parent
IS_DIR = REPO_DIR / "DB" / "is_analysis"

PA_NOIS    = IS_DIR / "presence_absence_nois.tsv"
FAM_TSV    = IS_DIR / "protein_families.tsv"
PA_ALL     = IS_DIR / "presence_absence.tsv"

IS_VARIANT_OUT  = IS_DIR / "is_variant_groups.tsv"
DEL_VARIANT_OUT = IS_DIR / "deletion_variant_groups.tsv"
MERGE_OUT       = IS_DIR / "loci_to_merge.tsv"
IS_COUNT_OUT    = IS_DIR / "is_counts_per_locus.tsv"

CORE_ASSEMBLY_GENES = {"wza", "wzb", "wzc", "wzy", "wzx", "kpsM", "kpsMT"}

# ── 1. Load non-IS presence/absence matrix ────────────────────────────────────
print("Loading non-IS presence/absence matrix...")
loci = []
biosyn_profiles = {}  # locus -> frozenset of family IDs

with open(PA_NOIS) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    fam_ids = [int(h[1:]) for h in header[1:]]  # strip leading 'f'
    for row in reader:
        locus = row[0]
        loci.append(locus)
        present = frozenset(fam_ids[i] for i, v in enumerate(row[1:]) if v == "1")
        biosyn_profiles[locus] = present

print(f"  {len(loci)} loci loaded")

# ── 2. Load IS family counts per locus ───────────────────────────────────────
print("Loading IS family counts per locus...")
is_counts = defaultdict(int)

with open(PA_ALL) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    is_col_idx = [i for i, h in enumerate(header[1:], 1) if h.endswith("*")]
    for row in reader:
        locus = row[0]
        is_counts[locus] = sum(1 for i in is_col_idx if row[i] == "1")

with open(IS_COUNT_OUT, "w") as fh:
    fh.write("locus\tn_is_families\n")
    for locus in loci:
        fh.write(f"{locus}\t{is_counts[locus]}\n")

n_with_is = sum(1 for l in loci if is_counts[l] > 0)
print(f"  {n_with_is} loci carry at least one IS family")

# ── 3. Identify IS-variant groups ────────────────────────────────────────────
print("Finding IS-variant groups...")

# Group loci by their biosynthetic gene profile
profile_to_loci = defaultdict(list)
for locus, profile in biosyn_profiles.items():
    profile_to_loci[profile].append(locus)

is_variant_groups = []
for profile, group in profile_to_loci.items():
    if len(group) < 2:
        continue
    # Sort: IS-free first (fewest IS families = primary)
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

print(f"  {len(is_variant_groups)} IS-variant groups found")
total_redundant_is = sum(len(g["variants"]) for g in is_variant_groups)
print(f"  {total_redundant_is} loci are IS-variants (to be merged into their primary)")

with open(IS_VARIANT_OUT, "w") as fh:
    fh.write("primary\tvariants\tn_members\tn_biosyn_families\t"
             "primary_is_count\tvariant_is_counts\n")
    for g in sorted(is_variant_groups, key=lambda x: x["primary"]):
        variants_str = ";".join(g["variants"])
        var_counts_str = ";".join(str(c) for c in g["variant_is_counts"])
        fh.write(
            f"{g['primary']}\t{variants_str}\t{1+len(g['variants'])}\t"
            f"{g['n_biosyn_families']}\t{g['primary_is_count']}\t{var_counts_str}\n"
        )

# ── 4. Identify deletion-variant groups ──────────────────────────────────────
print("\nFinding deletion-variant groups...")

# Build gene-name lookup: family_id -> set of gene names across all proteins in family
print("  Loading protein family gene annotations...")
fam_genes = defaultdict(set)
with open(FAM_TSV) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        fam = int(row["family_id"])
        gene = row["gene"]
        if gene and gene != "?":
            fam_genes[fam].add(gene)

# For each locus pair (A, B) where B ⊂ A (proper subset):
# check if A - B is entirely core assembly genes
del_variant_groups = []
assigned_del = set()

for i, locus_a in enumerate(loci):
    if locus_a in assigned_del:
        continue
    profile_a = biosyn_profiles[locus_a]
    group = [locus_a]

    for locus_b in loci[i+1:]:
        if locus_b in assigned_del:
            continue
        profile_b = biosyn_profiles[locus_b]

        # B is a proper subset of A?
        if profile_b < profile_a:
            diff_fams = profile_a - profile_b
            # Are all diff families core assembly gene families?
            diff_genes = set()
            for fam in diff_fams:
                diff_genes.update(fam_genes.get(fam, set()))
            # Allow unnamed families in diff only if all named ones are core
            named_diff = diff_genes - {"?", ""}
            if named_diff and named_diff.issubset(CORE_ASSEMBLY_GENES):
                group.append(locus_b)
                assigned_del.add(locus_b)

        # A is a proper subset of B?
        elif profile_a < profile_b:
            diff_fams = profile_b - profile_a
            diff_genes = set()
            for fam in diff_fams:
                diff_genes.update(fam_genes.get(fam, set()))
            named_diff = diff_genes - {"?", ""}
            if named_diff and named_diff.issubset(CORE_ASSEMBLY_GENES):
                # B is the full-length one; make it primary
                # We'll handle re-ordering when writing output
                pass  # handled in second pass below

    if len(group) > 1:
        # Sort: most genes first (full-length = primary)
        group.sort(key=lambda l: -len(biosyn_profiles[l]))
        primary = group[0]
        del_variant_groups.append({
            "primary": primary,
            "deletion_variants": group[1:],
            "primary_gene_count": len(biosyn_profiles[primary]),
        })
        assigned_del.update(group)

# Second pass: catch cases where A ⊂ B (B is full-length, A is deletion)
for i, locus_a in enumerate(loci):
    if locus_a in assigned_del:
        continue
    profile_a = biosyn_profiles[locus_a]
    group = [locus_a]

    for locus_b in loci[i+1:]:
        if locus_b in assigned_del:
            continue
        profile_b = biosyn_profiles[locus_b]
        if profile_a < profile_b:
            diff_fams = profile_b - profile_a
            diff_genes = set()
            for fam in diff_fams:
                diff_genes.update(fam_genes.get(fam, set()))
            named_diff = diff_genes - {"?", ""}
            if named_diff and named_diff.issubset(CORE_ASSEMBLY_GENES):
                group.append(locus_b)
                assigned_del.add(locus_b)

    if len(group) > 1:
        group.sort(key=lambda l: -len(biosyn_profiles[l]))
        primary = group[0]
        del_variant_groups.append({
            "primary": primary,
            "deletion_variants": group[1:],
            "primary_gene_count": len(biosyn_profiles[primary]),
        })
        assigned_del.update(group)

print(f"  {len(del_variant_groups)} deletion-variant groups found")
total_del = sum(len(g["deletion_variants"]) for g in del_variant_groups)
print(f"  {total_del} loci are deletion-variants")

with open(DEL_VARIANT_OUT, "w") as fh:
    fh.write("primary\tdeletion_variants\tn_members\tprimary_gene_count\n")
    for g in sorted(del_variant_groups, key=lambda x: x["primary"]):
        vars_str = ";".join(g["deletion_variants"])
        fh.write(
            f"{g['primary']}\t{vars_str}\t{1+len(g['deletion_variants'])}\t"
            f"{g['primary_gene_count']}\n"
        )

# ── 5. Combined merge table ───────────────────────────────────────────────────
print("\nWriting combined loci_to_merge.tsv...")
merge_rows = []
for g in is_variant_groups:
    for v in g["variants"]:
        merge_rows.append((v, g["primary"], "IS_variant", is_counts[v], is_counts[g["primary"]]))

for g in del_variant_groups:
    for v in g["deletion_variants"]:
        merge_rows.append((v, g["primary"], "deletion_variant",
                           len(biosyn_profiles[v]), len(biosyn_profiles[g["primary"]])))

merge_rows.sort(key=lambda r: r[0])

with open(MERGE_OUT, "w") as fh:
    fh.write("locus_to_remove\tprimary_locus\treason\t"
             "variant_metric\tprimary_metric\n")
    for row in merge_rows:
        fh.write("\t".join(str(x) for x in row) + "\n")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("VARIANT DETECTION SUMMARY")
print(f"{'='*60}")
print(f"  Total loci analysed:          {len(loci)}")
print(f"  Loci carrying IS elements:    {n_with_is}")
print(f"  IS-variant groups:            {len(is_variant_groups)}")
print(f"  IS-variants to merge/remove:  {total_redundant_is}")
print(f"  Deletion-variant groups:      {len(del_variant_groups)}")
print(f"  Deletion-variants to resolve: {total_del}")
print(f"  Total loci flagged:           {len(merge_rows)}")
print(f"\nOutputs:")
for p in [IS_VARIANT_OUT, DEL_VARIANT_OUT, MERGE_OUT, IS_COUNT_OUT]:
    print(f"  {p.name}")
