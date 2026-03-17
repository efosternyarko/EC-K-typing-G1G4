#!/usr/bin/env python3
"""
remove_oversized_loci.py — Remove aberrantly large novel loci from the v0.9 GenBank.

Novel ATB loci with >60 CDS captured extensive flanking chromosomal sequence during
extraction.  These universally outscore legitimate loci in Kaptive because their
conserved chromosomal genes are present in every assembly, giving a normalized score
close to the maximum regardless of K-locus identity.

The original 93 BSI loci have a maximum of 53 CDS (KL308).  A cutoff of 60 CDS
cleanly separates legitimate K loci from extraction artefacts in the novel set.

Usage:
    python scripts/remove_oversized_loci.py [--cds-threshold N] [--dry-run]
"""

import argparse
import re
from pathlib import Path

from Bio import SeqIO

REPO_DIR       = Path(__file__).resolve().parent.parent
DB_DIR         = REPO_DIR / "DB"
INPUT_GBK      = DB_DIR / "EC-K-typing_group1and4_v0.9.gbk"
OUTPUT_GBK     = DB_DIR / "EC-K-typing_group1and4_v0.9.gbk"   # overwrite
NOVEL_SUMMARY  = DB_DIR / "novel_kl_summary.tsv"
NOVEL_MAP      = DB_DIR / "novel_rep_kl_map.tsv"
NOVEL_FASTAS   = DB_DIR / "novel_rep_fastas"


def get_locus_name(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cds-threshold", type=int, default=60,
                        help="Remove loci with more than this many CDS (default: 60)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be removed without modifying files")
    args = parser.parse_args()

    # 1. Identify oversized loci
    to_remove = set()
    kept = []
    for rec in SeqIO.parse(str(INPUT_GBK), "genbank"):
        lname = get_locus_name(rec)
        n_cds = len([f for f in rec.features if f.type == "CDS"])
        if n_cds > args.cds_threshold:
            to_remove.add(lname)
            print(f"  REMOVE {lname:10s}  {n_cds:4d} CDS  {len(rec):,} bp")
        else:
            kept.append(rec)

    print(f"\nRemoving {len(to_remove)} loci; keeping {len(kept)}")

    if args.dry_run:
        print("(dry run — no files modified)")
        return

    # 2. Write filtered GenBank (overwrite in-place, using tmp file)
    tmp_gbk = INPUT_GBK.with_suffix(".tmp.gbk")
    SeqIO.write(kept, str(tmp_gbk), "genbank")
    tmp_gbk.replace(INPUT_GBK)
    print(f"Written: {INPUT_GBK}  ({len(kept)} records)")

    # 3. Remove from novel_kl_summary.tsv
    if NOVEL_SUMMARY.exists():
        lines = NOVEL_SUMMARY.read_text().splitlines()
        header = lines[0]
        kept_lines = [l for l in lines[1:] if l.split("\t")[0] not in to_remove]
        NOVEL_SUMMARY.write_text("\n".join([header] + kept_lines) + "\n")
        print(f"Updated: {NOVEL_SUMMARY}  ({len(kept_lines)} entries)")

    # 4. Remove from novel_rep_kl_map.tsv
    if NOVEL_MAP.exists():
        lines = NOVEL_MAP.read_text().splitlines()
        header = lines[0]
        kept_lines = [l for l in lines[1:] if l.split("\t")[0] not in to_remove]
        NOVEL_MAP.write_text("\n".join([header] + kept_lines) + "\n")
        print(f"Updated: {NOVEL_MAP}  ({len(kept_lines)} entries)")

    print("Done.")


if __name__ == "__main__":
    main()
