#!/usr/bin/env python3
"""
Apply the gene-content deduplication merging plan to produce database v1.0.

Reads merging_plan.tsv produced by investigate_duplicates.py and:
  1. Removes merged (redundant) loci from the GenBank database
  2. Updates novel_kl_summary.tsv and novel_rep_kl_map.tsv
  3. Writes a changelog

Input DB:  EC-K-typing_group1and4_v0.9.gbk  (651 loci)
Output DB: EC-K-typing_group1and4_v1.0.gbk  (628 loci)

Usage:
    python3 apply_merging_plan.py [--plan PATH] [--db_dir PATH] [--out_dir PATH]
"""

import argparse
import csv
from datetime import date
from pathlib import Path

from Bio import SeqIO


# ── Defaults ──────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parents[1]
DB_DIR      = REPO_ROOT / "DB"
DEFAULT_PLAN = Path("/tmp/duplicate_investigation/merging_plan.tsv")
IN_GBK      = DB_DIR / "EC-K-typing_group1and4_v0.9.gbk"
OUT_GBK     = DB_DIR / "EC-K-typing_group1and4_v1.0.gbk"
NOVEL_SUMMARY = DB_DIR / "novel_kl_summary.tsv"
NOVEL_MAP     = DB_DIR / "novel_rep_kl_map.tsv"


# ── Parse merging plan ────────────────────────────────────────────────────────
def load_merging_plan(plan_path: Path) -> tuple[set[str], dict[str, str]]:
    """
    Returns:
        loci_to_remove: set of KL types to drop from the database
        merge_map: dict mapping each removed locus → its retained representative
    """
    loci_to_remove = set()
    merge_map = {}
    with open(plan_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            keep = row["keep"]
            for merged in row["merge_into"].split(","):
                merged = merged.strip()
                if merged:
                    loci_to_remove.add(merged)
                    merge_map[merged] = keep
    return loci_to_remove, merge_map


# ── Extract KL type from a GenBank record ─────────────────────────────────────
def get_kl_type(record) -> str:
    for feat in record.features:
        if feat.type == "source":
            for note in feat.qualifiers.get("note", []):
                if note.startswith("K type:"):
                    return note.split("K type:")[-1].strip()
    return record.name.split("_")[0]


# ── Step 1: Filter GenBank ────────────────────────────────────────────────────
def filter_genbank(
    in_gbk: Path,
    out_gbk: Path,
    loci_to_remove: set[str],
) -> tuple[int, int]:
    kept = []
    removed_found = set()

    for record in SeqIO.parse(in_gbk, "genbank"):
        kl_type = get_kl_type(record)
        if kl_type in loci_to_remove:
            removed_found.add(kl_type)
        else:
            kept.append(record)

    # Warn if any planned removals weren't found
    not_found = loci_to_remove - removed_found
    if not_found:
        print(f"  WARNING: {len(not_found)} loci in plan not found in GenBank: "
              f"{sorted(not_found)}")

    with open(out_gbk, "w") as fh:
        SeqIO.write(kept, fh, "genbank")

    return len(kept), len(removed_found)


# ── Step 2: Update TSV mapping files ─────────────────────────────────────────
def filter_tsv(
    in_path: Path,
    out_path: Path,
    loci_to_remove: set[str],
    kl_col: str,
) -> tuple[int, int]:
    """Remove rows where kl_col value is in loci_to_remove."""
    kept_rows = []
    removed = 0
    with open(in_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            if row[kl_col] in loci_to_remove:
                removed += 1
            else:
                kept_rows.append(row)

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept_rows)

    return len(kept_rows), removed


# ── Step 3: Write changelog ───────────────────────────────────────────────────
def write_changelog(
    plan_path: Path,
    merge_map: dict[str, str],
    loci_to_remove: set[str],
    out_path: Path,
    n_input_loci: int = 0,
) -> None:
    # Reconstruct groups for the changelog
    groups: dict[str, list[str]] = {}
    for removed, kept in merge_map.items():
        groups.setdefault(kept, []).append(removed)

    lines = [
        f"EC-K-typing G1/G4 Database — v0.9 → v1.0 Changelog",
        f"Date: {date.today().isoformat()}",
        f"",
        f"CHANGE: Gene-content deduplication implementing Kaptive locus-definition principle",
        f"        'A unique locus is a unique set of genes, defined at 95% amino acid",
        f"        identity / 80% bidirectional coverage (CD-HIT)'",
        f"",
        f"METHOD: All CDS protein sequences were clustered with CD-HIT (95% aa / 80% cov).",
        f"        Each locus was represented as a frozenset of gene-cluster IDs.",
        f"        Loci with identical gene-content sets were merged; the lowest KL number",
        f"        (corresponding to the largest ATB cluster) was retained.",
        f"",
        f"RESULT: {len(loci_to_remove)} redundant loci removed across {len(groups)} merge groups.",
        f"        Database reduced from {n_input_loci} → {n_input_loci - len(loci_to_remove)} loci.",
        f"",
        f"MERGE GROUPS:",
    ]
    for kept, removed_list in sorted(groups.items()):
        lines.append(f"  KEEP {kept:<10} REMOVE {', '.join(sorted(removed_list))}")

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n".join(lines))


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plan",    default=DEFAULT_PLAN, type=Path)
    parser.add_argument("--db_dir",  default=DB_DIR,       type=Path)
    parser.add_argument("--out_dir", default=DB_DIR,       type=Path)
    args = parser.parse_args()

    in_gbk        = args.db_dir / "EC-K-typing_group1and4_v0.9.gbk"
    out_gbk       = args.out_dir / "EC-K-typing_group1and4_v1.0.gbk"
    novel_summary = args.db_dir / "novel_kl_summary.tsv"
    novel_map     = args.db_dir / "novel_rep_kl_map.tsv"
    changelog     = args.out_dir / "CHANGELOG_v0.9_to_v1.0.txt"

    print(f"Loading merging plan: {args.plan}")
    loci_to_remove, merge_map = load_merging_plan(args.plan)
    print(f"  {len(loci_to_remove)} loci to remove across "
          f"{len(set(merge_map.values()))} groups")

    print(f"\nStep 1: Filtering GenBank {in_gbk.name} → {out_gbk.name}")
    n_kept, n_removed = filter_genbank(in_gbk, out_gbk, loci_to_remove)
    print(f"  Kept: {n_kept}  |  Removed: {n_removed}")

    print(f"\nStep 2: Updating novel_kl_summary.tsv")
    kept_rows, removed_rows = filter_tsv(
        novel_summary, novel_summary.with_name("novel_kl_summary_v1.0.tsv"),
        loci_to_remove, kl_col="KL_type"
    )
    print(f"  Kept: {kept_rows}  |  Removed: {removed_rows}")

    print(f"\nStep 3: Updating novel_rep_kl_map.tsv")
    kept_rows2, removed_rows2 = filter_tsv(
        novel_map, novel_map.with_name("novel_rep_kl_map_v1.0.tsv"),
        loci_to_remove, kl_col="KL"
    )
    print(f"  Kept: {kept_rows2}  |  Removed: {removed_rows2}")

    print(f"\nStep 4: Writing changelog → {changelog}")
    n_input_loci = n_kept + n_removed
    write_changelog(args.plan, merge_map, loci_to_remove, changelog, n_input_loci)

    print(f"\n{'='*60}")
    print(f"v1.0 database: {out_gbk}")
    print(f"  Loci: {n_input_loci} → {n_kept}")
    print(f"  Duplicates removed: {n_removed}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
