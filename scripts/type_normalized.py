#!/usr/bin/env python3
"""
type_normalized.py — Normalised K-locus typing using Kaptive scores.

Kaptive's default scoring accumulates raw alignment score (AS) across all
genes found in a reference locus.  Large references (e.g. KL302: 40 CDS /
37.8 kb; KL304: 44 CDS / 46.5 kb) systematically dominate smaller same-KX-type
loci even at lower per-base identity, causing 23/93 filtered loci to fail
self-typing in v3.0.

Fix
---
Run Kaptive with `--scores` to obtain the full per-locus score matrix for every
assembly, then re-rank each assembly's loci by:

    normalised_score  =  AS  /  total_expected_gene_bp

where `total_expected_gene_bp` is the sum of CDS sequence lengths in the
reference GenBank record.  This converts raw bitscore into "alignment score per
expected reference base" — a per-base identity metric that is independent of
locus size.

Loci that match all their genes at 100 % identity receive the maximum possible
normalised score regardless of how many genes they have.  Only loci that either
miss genes (lower AS numerator) or have lower per-base identity score worse.

Usage
-----
    python scripts/type_normalized.py [--threads N] [--suffix v3norm]
                                      [--db PATH] [--skip-kaptive]

    --threads N       Kaptive alignment threads (0 = all; default 0)
    --suffix SUFFIX   Output file suffix (default: v3norm)
    --db PATH         Kaptive database (default: EC-K-typing_all_groups_v3.0.gbk)
    --skip-kaptive    Re-use an existing scores TSV instead of re-running Kaptive

Outputs (written to DB/)
------------------------
    kaptive_scores_{suffix}.tsv          — raw Kaptive scores matrix (183 loci × N assemblies)
    kaptive_validation_results_{suffix}.tsv — normalised typing results (Kaptive-format)
    kaptive_validation_summary_{suffix}.tsv — per-assembly comparison table
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR           = Path(__file__).resolve().parent.parent
DB_DIR             = REPO_DIR / "DB"
_DEFAULT_DB        = DB_DIR / "EC-K-typing_all_groups_v3.0.gbk"
MAPPING_FILE       = DB_DIR / "KL_G1G4_mapping.tsv"
FILTERED_MAPPING   = DB_DIR / "KL_G1G4_mapping_filtered.tsv"
EXTRACTION_SUMMARY = Path(
    "/Users/LSHEF4/Dropbox/work_2025/e_coli_db/db_build_v2/extraction_summary.tsv"
)
GENOMES_DIR        = Path(
    "/Users/LSHEF4/Dropbox/work_2025/e_coli_db/db_build_v2/nohits_genomes"
)

# Typeable threshold: fraction of expected genes that must be found
MIN_GENE_COVERAGE = 0.50


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_locus_name(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def load_locus_stats(db: Path) -> dict:
    """Return {locus: total_expected_gene_bp} and {locus: n_genes} from GenBank."""
    locus_total_bp = {}
    locus_n_genes  = {}
    for rec in SeqIO.parse(str(db), "genbank"):
        lname = get_locus_name(rec)
        cds   = [f for f in rec.features if f.type == "CDS"]
        locus_n_genes[lname]  = len(cds)
        locus_total_bp[lname] = sum(len(f) for f in cds)
    return locus_total_bp, locus_n_genes


def run_kaptive_scores(db: Path, genome_files: list, scores_tsv: Path, threads: int) -> None:
    cmd = (
        ["kaptive", "assembly", str(db)]
        + [str(g) for g in genome_files]
        + ["--scores", str(scores_tsv), "-t", str(threads)]
    )
    print(f"\n[kaptive] Running --scores on {len(genome_files)} assemblies → {scores_tsv.name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("ERROR running kaptive:", file=sys.stderr)
        print(result.stderr[:1000], file=sys.stderr)
        sys.exit(1)
    if result.stderr:
        print(result.stderr[:300])


def normalise_and_rank(scores_tsv: Path, locus_total_bp: dict) -> pd.DataFrame:
    """
    Read the Kaptive scores matrix TSV and re-rank each assembly by normalised score.

    Returns a DataFrame with one row per assembly:
      Assembly, Best match locus, Best match confidence,
      Genes found, Genes expected, Gene coverage %,
      Raw AS, Norm AS, Norm rank (1 = best)
    """
    df = pd.read_csv(scores_tsv, sep="\t")

    # Column names from Kaptive v3.1 --scores output
    # Assembly  Locus  AS  mlen  blen  q_len  genes_found  genes_expected
    asm_col = df.columns[0]   # "Assembly"

    df["total_expected_bp"] = df["Locus"].map(locus_total_bp)
    # Avoid div-by-zero for loci not in our GenBank (e.g. G2/G3)
    # For G2/G3 loci, use q_len as proxy for total_expected_bp
    mask_missing = df["total_expected_bp"].isna() | (df["total_expected_bp"] == 0)
    df.loc[mask_missing, "total_expected_bp"] = df.loc[mask_missing, "q_len"].clip(lower=1)

    df["AS_norm"] = df["AS"] / df["total_expected_bp"]

    rows = []
    for asm, grp in df.groupby(asm_col):
        active = grp[grp["AS"] > 0]
        if active.empty:
            rows.append({
                "Assembly":             asm,
                "Best match locus":     "none",
                "Best match confidence": "Untypeable",
                "Genes found":          0,
                "Genes expected":       0,
                "Gene coverage":        "0.0%",
                "Raw AS":               0,
                "Norm AS":              0.0,
            })
            continue

        # Primary rank: AS_norm descending.
        # Near-tie tiebreaker: among loci scoring within 0.5 % of the best
        # normalised score, prefer the one with more genes found (more absolute
        # evidence of a correct match), then fall back to AS_norm.
        top_norm = active["AS_norm"].max()
        candidates = active[active["AS_norm"] >= top_norm * 0.995]
        best = candidates.sort_values(
            ["genes_found", "AS_norm"], ascending=[False, False]
        ).iloc[0]
        genes_found    = int(best["genes_found"])
        genes_expected = int(best["genes_expected"])
        coverage       = genes_found / genes_expected if genes_expected > 0 else 0
        conf           = "Typeable" if coverage >= MIN_GENE_COVERAGE else "Untypeable"

        rows.append({
            "Assembly":             asm,
            "Best match locus":     best["Locus"],
            "Best match confidence": conf,
            "Genes found":          genes_found,
            "Genes expected":       genes_expected,
            "Gene coverage":        f"{100*coverage:.1f}%",
            "Raw AS":               int(best["AS"]),
            "Norm AS":              round(best["AS_norm"], 4),
        })

    return pd.DataFrame(rows)


def parse_results(results_df: pd.DataFrame, expected_kl: dict, filtered_kls: set) -> pd.DataFrame:
    df = results_df.copy()
    df["assembly_stem"] = df["Assembly"].apply(lambda x: Path(x).stem)
    df["expected_KL"]      = df["assembly_stem"].map(expected_kl)
    df["is_representative"] = df["expected_KL"].notna()
    df["in_filtered_set"]  = df["expected_KL"].apply(
        lambda x: x in filtered_kls if pd.notna(x) else False
    )
    df["best_match_locus"] = df["Best match locus"]
    df["correct_type"] = df.apply(
        lambda r: (r["best_match_locus"] == r["expected_KL"])
        if r["is_representative"] else pd.NA,
        axis=1,
    )
    df["typeable"] = df["Best match confidence"].str.contains("Typeable", na=False) & \
                     ~df["Best match confidence"].str.contains("Untypeable", na=False)
    return df


def print_summary(df: pd.DataFrame) -> None:
    reps          = df[df["is_representative"]]
    reps_filtered = df[df["in_filtered_set"]]

    print("\n" + "=" * 60)
    print("NORMALISED TYPING SUMMARY  (AS / total_expected_gene_bp)")
    print("=" * 60)

    n_all     = len(reps)
    n_all_ok  = int(reps["correct_type"].sum())
    n_filt    = len(reps_filtered)
    n_filt_ok = int(reps_filtered["correct_type"].sum())

    print(f"\n1. Self-typing (representative assemblies)")
    print(f"   All 125 loci:          {n_all_ok}/{n_all}  "
          f"({100*n_all_ok/n_all:.1f}%)")
    print(f"   Filtered 93 loci:      {n_filt_ok}/{n_filt}  "
          f"({100*n_filt_ok/n_filt:.1f}%)")

    wrong = reps[reps["correct_type"] == False][
        ["assembly_stem", "expected_KL", "best_match_locus"]
    ]
    if not wrong.empty:
        print(f"\n   Incorrect self-types ({len(wrong)}):")
        for _, row in wrong.iterrows():
            print(f"     {row['assembly_stem']:30s}  expected {row['expected_KL']:10s}  "
                  f"got {row['best_match_locus']}")

    n_total    = len(df)
    n_typeable = int(df["typeable"].sum())
    print(f"\n2. Typeability (all {n_total} extracted assemblies)")
    print(f"   Typeable:              {n_typeable}/{n_total}  "
          f"({100*n_typeable/n_total:.1f}%)")

    utilised = set(reps_filtered[reps_filtered["correct_type"] == True]["expected_KL"])
    not_util = set(reps_filtered["expected_KL"]) - utilised
    n_filt   = len(set(reps_filtered["expected_KL"]))
    print(f"\n3. Locus utilisation (filtered {n_filt}-locus set)")
    print(f"   Correctly self-typed:  {len(utilised)}/{n_filt}")
    if not_util:
        print(f"   Not self-typed ({len(not_util)}):  {', '.join(sorted(not_util))}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-t", "--threads", type=int, default=0)
    parser.add_argument("--suffix",      default="v3norm",
                        help="Output file suffix (default: v3norm)")
    parser.add_argument("--db",          type=Path, default=None,
                        help="Kaptive database (default: v3.0 all-groups)")
    parser.add_argument("--skip-kaptive", action="store_true",
                        help="Skip Kaptive; re-use existing scores TSV")
    args = parser.parse_args()

    db         = args.db if args.db else _DEFAULT_DB
    suffix     = args.suffix
    scores_tsv = DB_DIR / f"kaptive_scores_{suffix}.tsv"
    results_tsv= DB_DIR / f"kaptive_validation_results_{suffix}.tsv"
    summary_tsv= DB_DIR / f"kaptive_validation_summary_{suffix}.tsv"

    # Load reference locus statistics
    print(f"[1] Loading locus stats from {db.name}...")
    locus_total_bp, locus_n_genes = load_locus_stats(db)
    print(f"    {len(locus_total_bp)} loci; "
          f"total expected bp range: "
          f"{min(locus_total_bp.values())}–{max(locus_total_bp.values())} bp")

    # Load mapping tables
    mapping          = pd.read_csv(MAPPING_FILE, sep="\t")
    filtered_mapping = pd.read_csv(FILTERED_MAPPING, sep="\t")
    extraction       = pd.read_csv(EXTRACTION_SUMMARY, sep="\t")

    expected_kl  = dict(zip(mapping["source_assembly"], mapping["KL"]))
    filtered_kls = set(filtered_mapping["KL"])

    extracted_asms = extraction[extraction["status"] == "extracted"]["assembly"].tolist()
    genome_files, missing = [], []
    for asm in extracted_asms:
        fa = GENOMES_DIR / asm
        if fa.exists():
            genome_files.append(fa)
        else:
            missing.append(asm)
    if missing:
        print(f"WARNING: {len(missing)} genome file(s) not found", file=sys.stderr)
    print(f"[2] Genome files: {len(genome_files)}")

    # Run Kaptive --scores (or skip)
    if not args.skip_kaptive:
        if not db.exists():
            print(f"ERROR: database not found: {db}", file=sys.stderr)
            sys.exit(1)
        run_kaptive_scores(db, genome_files, scores_tsv, args.threads)
    else:
        if not scores_tsv.exists():
            print(f"ERROR: --skip-kaptive but {scores_tsv} not found", file=sys.stderr)
            sys.exit(1)
        print(f"\n[kaptive] Skipping run, reading: {scores_tsv}")

    # Normalise and report
    print(f"\n[3] Normalising scores (AS / total_expected_gene_bp)...")
    results_df = normalise_and_rank(scores_tsv, locus_total_bp)
    results_df.to_csv(results_tsv, sep="\t", index=False)
    print(f"    Written: {results_tsv.name}")

    df = parse_results(results_df, expected_kl, filtered_kls)
    df.to_csv(summary_tsv, sep="\t", index=False)

    print_summary(df)

    print(f"\n   Raw scores:   {scores_tsv}")
    print(f"   Results:      {results_tsv}")
    print(f"   Summary:      {summary_tsv}")


if __name__ == "__main__":
    main()
