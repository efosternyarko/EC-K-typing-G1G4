#!/usr/bin/env python3
"""
analyse_atb_results.py — Rank ATB matches to identify better KL300/303/306/307 reps.

Reads the raw LexicMap output from the ATB search and ranks candidate assemblies
by two criteria:

  1. High coverage of the query locus (>= 80%) — confirms it really has this locus type
  2. Low identity to KL302 variable region — means it would discriminate better

Since we can't directly measure criterion 2 from LexicMap output alone, we use a
proxy: among assemblies that match the query locus well, prefer those with the
LOWEST overall identity to KL302. These are candidates for being "more diverged"
representatives that might self-type correctly.

Usage
-----
    python scripts/analyse_atb_results.py

Outputs
-------
    DB/atb_lexicmap_results/atb_candidates.tsv  — ranked candidate assemblies
    DB/atb_lexicmap_results/atb_summary.txt     — per-locus summary
"""

import gzip
import sys
from pathlib import Path

import pandas as pd

REPO_DIR    = Path(__file__).resolve().parent.parent
DB_DIR      = REPO_DIR / "DB"
RESULTS_DIR = DB_DIR / "atb_lexicmap_results"
RAW_TSV     = RESULTS_DIR / "atb_lexicmap_raw.tsv.gz"
OUT_TSV     = RESULTS_DIR / "atb_candidates.tsv"
OUT_SUMMARY = RESULTS_DIR / "atb_summary.txt"

TARGETS = ["KL300", "KL303", "KL306", "KL307"]

# Thresholds
MIN_QCOV   = 70.0   # % of query covered by this genome
MIN_PIDENT = 90.0   # % identity


def main():
    if not RAW_TSV.exists():
        print(f"ERROR: {RAW_TSV} not found. Run run_atb_lexicmap_ec2.sh first.",
              file=sys.stderr)
        sys.exit(1)

    print(f"Reading {RAW_TSV}...")
    opener = gzip.open if str(RAW_TSV).endswith(".gz") else open
    with opener(RAW_TSV, "rt") as fh:
        # LexicMap output columns (v0.8):
        # query  qlen  hits  sgenome  sseqid  qcovGnm  hsp  qcovHSP
        # alenHSP  alenSeg  pident  gaps  qstart  qend  sstart  send  sstr  cigar
        df = pd.read_csv(fh, sep="\t", comment="#")

    print(f"  {len(df):,} alignment rows")
    print(f"  Columns: {list(df.columns)}")

    # Normalise column names (LexicMap versions may vary slightly)
    df.columns = [c.strip() for c in df.columns]

    # Identify relevant columns
    query_col  = [c for c in df.columns if c.lower() in ("query", "qname")][0]
    genome_col = [c for c in df.columns if c.lower() in ("sgenome", "genome", "subject")][0]
    pident_col = [c for c in df.columns if c.lower() in ("pident", "identity")][0]
    qcov_col   = [c for c in df.columns if "qcovgnm" in c.lower() or "qcov" in c.lower()][0]

    # Parse query locus name (strip length/cds annotations)
    df["locus"] = df[query_col].str.split().str[0]

    # Filter to target loci
    df = df[df["locus"].isin(TARGETS)].copy()
    print(f"  {len(df):,} rows for target loci {TARGETS}")

    # Per genome: take best hit (highest qcovGnm * pident)
    df[pident_col] = pd.to_numeric(df[pident_col], errors="coerce")
    df[qcov_col]   = pd.to_numeric(df[qcov_col],   errors="coerce")
    df["score"]    = df[pident_col] * df[qcov_col]

    best = (df.sort_values("score", ascending=False)
              .drop_duplicates(subset=["locus", genome_col])
              .copy())

    # Filter by thresholds
    candidates = best[
        (best[pident_col] >= MIN_PIDENT) &
        (best[qcov_col]   >= MIN_QCOV)
    ].copy()

    print(f"\n  Candidates (pident>={MIN_PIDENT}%, qcov>={MIN_QCOV}%): "
          f"{len(candidates):,}")

    # Rank within each locus: highest coverage first, then highest identity
    candidates = candidates.sort_values(
        [qcov_col, pident_col], ascending=[False, False]
    )

    # Tidy output
    out_cols = ["locus", genome_col, qcov_col, pident_col]
    out_cols += [c for c in candidates.columns
                 if c not in out_cols and c not in ("score",)]
    candidates[out_cols].to_csv(OUT_TSV, sep="\t", index=False)
    print(f"  Written: {OUT_TSV}")

    # Summary
    lines = []
    lines.append("ATB LexicMap search — candidate assemblies per locus")
    lines.append("=" * 60)
    for locus in TARGETS:
        sub = candidates[candidates["locus"] == locus]
        lines.append(f"\n{locus}:  {len(sub)} candidates")
        if not sub.empty:
            lines.append(f"  {'Assembly':<30}  {'qcov%':>6}  {'pident%':>7}")
            lines.append(f"  {'-'*30}  {'------':>6}  {'-------':>7}")
            for _, row in sub.head(20).iterrows():
                lines.append(
                    f"  {str(row[genome_col]):<30}  "
                    f"{row[qcov_col]:>6.1f}  "
                    f"{row[pident_col]:>7.2f}"
                )

    summary_text = "\n".join(lines)
    print("\n" + summary_text)
    with open(OUT_SUMMARY, "w") as fh:
        fh.write(summary_text + "\n")
    print(f"\nWritten: {OUT_SUMMARY}")
    print("\nNext steps:")
    print("  1. Review atb_candidates.tsv — pick 5–10 candidates per locus")
    print("  2. Download assemblies from ATB / ENA")
    print("  3. Run through kaptive + type_normalized.py to check self-typing")


if __name__ == "__main__":
    main()
