#!/usr/bin/env python3
"""
validate_db_v2.py — Validate the v2.0 EC-K-typing G1/G4 database using Kaptive.

Runs Kaptive on all 565 assemblies from which a locus was successfully extracted,
then evaluates:

  1. Self-typing:    the 125 representative assemblies should each return their own KL locus
  2. Typeability:    fraction of all 565 extracted assemblies that are confidently typed
  3. Utilisation:    how many of the 93 filtered loci are correctly self-typed

Outputs
-------
  DB/kaptive_validation_results_v2.tsv  — raw Kaptive output
  DB/kaptive_validation_summary_v2.tsv  — per-assembly table with expected vs observed KL

Usage
-----
  python scripts/validate_db_v2.py [--threads N]
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR          = Path(__file__).resolve().parent.parent
DB                = REPO_DIR / "DB" / "EC-K-typing_all_groups_v2.0.gbk"
MAPPING_FILE      = REPO_DIR / "DB" / "KL_G1G4_mapping.tsv"
FILTERED_MAPPING  = REPO_DIR / "DB" / "KL_G1G4_mapping_filtered.tsv"
EXTRACTION_SUMMARY = Path(
    "/Users/LSHEF4/Dropbox/work_2025/e_coli_db/db_build_v2/extraction_summary.tsv"
)
GENOMES_DIR       = Path(
    "/Users/LSHEF4/Dropbox/work_2025/e_coli_db/db_build_v2/nohits_genomes"
)
RESULTS_TSV       = REPO_DIR / "DB" / "kaptive_validation_results_v2.tsv"
SUMMARY_TSV       = REPO_DIR / "DB" / "kaptive_validation_summary_v2.tsv"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_kaptive(db: Path, genome_files: list[Path], output: Path, threads: int) -> None:
    cmd = (
        ["kaptive", "assembly", str(db)]
        + [str(g) for g in genome_files]
        + ["-o", str(output), "-t", str(threads)]
    )
    print(f"\n[kaptive] Running on {len(genome_files)} assemblies → {output.name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("ERROR running kaptive:", file=sys.stderr)
        print(result.stderr[:1000], file=sys.stderr)
        sys.exit(1)
    if result.stderr:
        print(result.stderr[:500])


def parse_results(tsv: Path, expected_kl: dict, filtered_kls: set) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t")

    # Kaptive uses the filename stem as the assembly ID
    asm_col = df.columns[0]
    df["assembly_stem"] = df[asm_col].apply(lambda x: Path(x).stem)

    df["expected_KL"]      = df["assembly_stem"].map(expected_kl)
    df["is_representative"] = df["expected_KL"].notna()
    df["in_filtered_set"]  = df["expected_KL"].apply(
        lambda x: x in filtered_kls if pd.notna(x) else False
    )

    best_col = [c for c in df.columns if "best match locus" in c.lower()][0]
    df["best_match_locus"] = df[best_col]
    df["correct_type"]     = df.apply(
        lambda r: (r["best_match_locus"] == r["expected_KL"])
        if r["is_representative"] else pd.NA,
        axis=1,
    )

    conf_col = next(
        (c for c in df.columns if "confidence" in c.lower()), None
    )
    if conf_col:
        df["typeable"] = (
            df[conf_col].str.contains("Typeable", na=False)
            & ~df[conf_col].str.contains("Untypeable", na=False)
        )
    else:
        df["typeable"] = pd.NA

    return df


def print_summary(df: pd.DataFrame) -> None:
    reps          = df[df["is_representative"]]
    reps_filtered = df[df["in_filtered_set"]]

    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY — v2.0 database (Kaptive v3.1.0)")
    print("=" * 60)

    # 1. Self-typing
    n_all       = len(reps)
    n_all_ok    = int(reps["correct_type"].sum())
    n_filt      = len(reps_filtered)
    n_filt_ok   = int(reps_filtered["correct_type"].sum())

    print(f"\n1. Self-typing (representative assemblies)")
    print(f"   All 125 loci:          {n_all_ok}/{n_all}  "
          f"({100*n_all_ok/n_all:.1f}%)")
    print(f"   Filtered 93 loci:      {n_filt_ok}/{n_filt}  "
          f"({100*n_filt_ok/n_filt:.1f}%)")

    # Wrong types
    wrong = reps[reps["correct_type"] == False][
        ["assembly_stem", "expected_KL", "best_match_locus"]
    ]
    if not wrong.empty:
        print(f"\n   Incorrect self-types ({len(wrong)}):")
        for _, row in wrong.iterrows():
            print(f"     {row['assembly_stem']:30s}  expected {row['expected_KL']:10s}  "
                  f"got {row['best_match_locus']}")

    # 2. Typeability
    if "typeable" in df.columns and df["typeable"].notna().any():
        n_total     = len(df)
        n_typeable  = int(df["typeable"].sum())
        print(f"\n2. Typeability (all {n_total} extracted assemblies)")
        print(f"   Typeable:              {n_typeable}/{n_total}  "
              f"({100*n_typeable/n_total:.1f}%)")

        # Confidence breakdown
        conf_col = next((c for c in df.columns if "confidence" in c.lower()), None)
        if conf_col:
            print(f"\n   Confidence breakdown:")
            for val, cnt in df[conf_col].value_counts().items():
                print(f"     {val:30s}  {cnt}")

    # 3. Locus utilisation
    utilised = set(reps_filtered[reps_filtered["correct_type"] == True]["expected_KL"])
    not_util = set(reps_filtered["expected_KL"]) - utilised
    print(f"\n3. Locus utilisation (filtered 93-locus set)")
    print(f"   Correctly self-typed:  {len(utilised)}/93")
    if not_util:
        print(f"   Not self-typed ({len(not_util)}):  {', '.join(sorted(not_util))}")

    print(f"\nFull Kaptive output:  {RESULTS_TSV}")
    print(f"Per-assembly summary: {SUMMARY_TSV}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-t", "--threads", type=int, default=0,
        help="Kaptive alignment threads (0 = all available, default: 0)"
    )
    parser.add_argument(
        "--skip-kaptive", action="store_true",
        help="Skip running kaptive and parse existing results file"
    )
    args = parser.parse_args()

    # Load mapping tables
    mapping          = pd.read_csv(MAPPING_FILE, sep="\t")
    filtered_mapping = pd.read_csv(FILTERED_MAPPING, sep="\t")
    extraction       = pd.read_csv(EXTRACTION_SUMMARY, sep="\t")

    # expected_kl: assembly_stem -> KL (e.g. "ESC_GC2276AA_AS" -> "KL300")
    expected_kl  = dict(zip(mapping["source_assembly"], mapping["KL"]))
    filtered_kls = set(filtered_mapping["KL"])

    # Assemblies to run: those where a locus was extracted
    extracted_asms = extraction[extraction["status"] == "extracted"]["assembly"].tolist()
    print(f"Extracted assemblies in summary: {len(extracted_asms)}")

    # Resolve to genome file paths
    genome_files, missing = [], []
    for asm in extracted_asms:
        fa = GENOMES_DIR / asm
        if fa.exists():
            genome_files.append(fa)
        else:
            missing.append(asm)

    if missing:
        print(f"WARNING: {len(missing)} genome file(s) not found on disk — skipping",
              file=sys.stderr)
        for m in missing[:5]:
            print(f"  {m}", file=sys.stderr)
        if len(missing) > 5:
            print(f"  ... and {len(missing)-5} more", file=sys.stderr)

    print(f"Genome files found:              {len(genome_files)}")

    # Run Kaptive (or skip if --skip-kaptive)
    if not args.skip_kaptive:
        if not DB.exists():
            print(f"ERROR: database not found: {DB}", file=sys.stderr)
            sys.exit(1)
        run_kaptive(DB, genome_files, RESULTS_TSV, args.threads)
    else:
        if not RESULTS_TSV.exists():
            print(f"ERROR: --skip-kaptive set but results file not found: {RESULTS_TSV}",
                  file=sys.stderr)
            sys.exit(1)
        print(f"\n[kaptive] Skipping run, reading existing: {RESULTS_TSV}")

    # Parse and summarise
    df = parse_results(RESULTS_TSV, expected_kl, filtered_kls)
    df.to_csv(SUMMARY_TSV, sep="\t", index=False)

    print_summary(df)


if __name__ == "__main__":
    main()
