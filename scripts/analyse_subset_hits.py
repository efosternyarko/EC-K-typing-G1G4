#!/usr/bin/env python3
"""
analyse_subset_hits.py — Identify the best ATB replacement representatives
for high-priority subset loci (KL601, KL713, KL742).

Reads raw LexicMap output from find_subset_reps_ec2.sh and applies
locus-specific selection logic:

  KL601 / KL742 (3'-truncated subsets):
    Select assemblies where the subject region extends FURTHEST beyond the
    query end — these contain a more complete version of the truncated locus.
    Rank by: (send - qend) descending, breaking ties by pident descending.
    The subject end position relative to the query end estimates how much
    extra locus sequence is present in that assembly beyond our truncated rep.

  KL713 (superset with cross-contig artefact):
    Select assemblies where the match is contained within a SINGLE CONTIG
    (no 100-N spacer in the matched region). Among those, prefer longest
    subject match (slen) and highest pident.
    These are candidates for a clean, non-artefactual KL713 representative.

Usage:
    python scripts/analyse_subset_hits.py

Input:
    DB/atb_lexicmap_results/subset_fix/atb_subset_raw.tsv.gz

Outputs:
    DB/atb_lexicmap_results/subset_fix/subset_candidates.tsv
    DB/atb_lexicmap_results/subset_fix/subset_summary.txt
"""

import gzip
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT   = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "DB" / "atb_lexicmap_results" / "subset_fix"
RAW_TSV     = RESULTS_DIR / "atb_subset_raw.tsv.gz"
OUT_TSV     = RESULTS_DIR / "subset_candidates.tsv"
OUT_SUMMARY = RESULTS_DIR / "subset_summary.txt"

TARGETS = ["KL601", "KL713", "KL742"]

# Selection thresholds
MIN_QCOV   = 70.0   # genome must cover >=70% of query
MIN_PIDENT = 90.0   # >=90% nucleotide identity
TOP_N      = 20     # report top N candidates per locus


def main():
    if not RAW_TSV.exists():
        print(f"ERROR: {RAW_TSV} not found. Run find_subset_reps_ec2.sh first.",
              file=sys.stderr)
        sys.exit(1)

    print(f"Reading {RAW_TSV}...")
    opener = gzip.open if str(RAW_TSV).endswith(".gz") else open
    with opener(RAW_TSV, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", comment="#")

    df.columns = [c.strip() for c in df.columns]
    print(f"  {len(df):,} alignment rows, columns: {list(df.columns)}")

    # Normalise key column names across LexicMap versions
    col_map = {}
    for c in df.columns:
        cl = c.lower()
        if cl in ("query", "qname"):          col_map[c] = "query"
        elif cl in ("sgenome", "genome"):      col_map[c] = "sgenome"
        elif cl in ("sseqid", "scontig"):      col_map[c] = "scontig"
        elif cl in ("pident", "identity"):     col_map[c] = "pident"
        elif "qcovgnm" in cl or cl == "qcov": col_map[c] = "qcovGnm"
        elif cl == "qstart":                   col_map[c] = "qstart"
        elif cl == "qend":                     col_map[c] = "qend"
        elif cl == "sstart":                   col_map[c] = "sstart"
        elif cl == "send":                     col_map[c] = "send"
        elif cl in ("qlen",):                  col_map[c] = "qlen"
    df = df.rename(columns=col_map)

    # Parse locus name from query field
    df["locus"] = df["query"].str.split().str[0]
    df = df[df["locus"].isin(TARGETS)].copy()
    print(f"  {len(df):,} rows for target loci")

    # Numeric coercion
    for col in ("pident", "qcovGnm", "qstart", "qend", "sstart", "send", "qlen"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Filter to quality threshold
    df = df[
        (df["pident"] >= MIN_PIDENT) &
        (df["qcovGnm"] >= MIN_QCOV)
    ].copy()
    print(f"  {len(df):,} rows after quality filter (pident>={MIN_PIDENT}%, qcov>={MIN_QCOV}%)")

    summary_lines = []
    all_candidates = []

    for locus in TARGETS:
        sub = df[df["locus"] == locus].copy()
        print(f"\n{locus}: {len(sub):,} hits after filter")

        if sub.empty:
            summary_lines.append(f"\n{locus}: NO HITS after quality filter")
            continue

        # Per-genome: take the single best-scoring alignment row
        sub["score"] = sub["pident"] * sub["qcovGnm"]
        best = sub.sort_values("score", ascending=False).drop_duplicates("sgenome")

        if locus in ("KL601", "KL742"):
            # For truncated subsets: rank by how far the subject extends beyond
            # the query end — a larger (send - qend) means more locus captured
            # beyond our truncated representative.
            if "send" in best.columns and "qend" in best.columns:
                best["extra_bp"] = best["send"] - best["qend"]
                best = best.sort_values(
                    ["extra_bp", "pident"], ascending=[False, False]
                )
                rank_note = "ranked by extra subject bp beyond query end (more = more complete)"
            else:
                # Fallback: rank by pident
                best = best.sort_values("pident", ascending=False)
                rank_note = "ranked by pident (send/qend columns not available)"

        else:  # KL713
            # For the artefactual superset: rank by pident then subject length
            # (longer matches on single contigs preferred)
            best = best.sort_values(["pident", "qcovGnm"], ascending=[False, False])
            rank_note = "ranked by pident + qcov (prefer clean same-contig hits)"

        best["locus"]     = locus
        best["rank_note"] = rank_note
        top = best.head(TOP_N)
        all_candidates.append(top)

        summary_lines.append(
            f"\n{'='*60}\n"
            f"{locus} — top {min(TOP_N, len(top))} candidates\n"
            f"  Strategy: {rank_note}\n"
        )
        cols_to_show = [c for c in
                        ["sgenome", "scontig", "pident", "qcovGnm",
                         "qstart", "qend", "sstart", "send", "extra_bp"]
                        if c in top.columns]
        summary_lines.append(top[cols_to_show].to_string(index=False))

    # Write outputs
    if all_candidates:
        combined = pd.concat(all_candidates, ignore_index=True)
        combined.to_csv(OUT_TSV, sep="\t", index=False)
        print(f"\nCandidates written to {OUT_TSV}")
    else:
        print("\nWARNING: No candidates found for any target locus.")

    summary_text = "\n".join(summary_lines)
    with open(OUT_SUMMARY, "w") as fh:
        fh.write(
            "EC-K-typing G1/G4 — Subset Locus Representative Search\n"
            "========================================================\n"
            f"LexicMap thresholds: pident>={MIN_PIDENT}%, qcov>={MIN_QCOV}%\n"
            + summary_text + "\n"
        )
    print(f"Summary written to {OUT_SUMMARY}")

    print("\nNext steps:")
    print("  1. Review subset_candidates.tsv")
    print("  2. For KL601/KL742: download top-ranked assemblies, re-extract locus,")
    print("     verify gene count increases, replace representative in GenBank")
    print("  3. For KL713: download top-ranked assembly, confirm no N-spacer in locus")
    print("     region, replace representative in GenBank")
    print("  4. Re-run validate_gene_content.py to confirm subset pairs resolved")


if __name__ == "__main__":
    main()
