#!/usr/bin/env python3
"""
test_blast_candidates.py — Download NCBI BLAST candidates and test with Kaptive.

For each candidate assembly found by blast_ncbi_candidates.py:
  1. Download the chromosome FASTA via NCBI Entrez
  2. Run Kaptive (normalised scoring) against the v0.4 database
  3. Report which candidates type correctly as their expected locus (not KL302)

The best candidates — those that score higher as their own type than as KL302 —
are replacement representative sequences for the failing loci.

Usage
-----
    python scripts/test_blast_candidates.py [--top-n 10] [--threads 8]

Outputs
-------
    DB/blast_ncbi_results/candidate_genomes/   — downloaded FASTAs
    DB/blast_ncbi_results/candidate_kaptive_scores.tsv
    DB/blast_ncbi_results/candidate_kaptive_results.tsv
    DB/blast_ncbi_results/candidate_typing_summary.tsv  ← key output
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "lshef4@lshtm.ac.uk"

REPO_DIR     = Path(__file__).resolve().parent.parent
DB_DIR       = REPO_DIR / "DB"
RESULTS_DIR  = DB_DIR / "blast_ncbi_results"
GENOMES_DIR  = RESULTS_DIR / "candidate_genomes"
CANDIDATES   = RESULTS_DIR / "candidates_summary.tsv"   # rebuilt from per-locus files
DATABASE     = DB_DIR / "EC-K-typing_all_groups_v0.4.gbk"

TARGETS = ["KL300", "KL303", "KL306", "KL307"]


def pick_candidates(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """
    Select top_n candidates per locus.
    Strategy: highest qcov first, then pident.
    For KL306 (only 6 hits), take all of them.
    """
    rows = []
    for locus in TARGETS:
        sub = df[df["locus"] == locus].sort_values(
            ["qcov", "pident"], ascending=[False, False]
        )
        n = len(sub) if locus == "KL306" else min(top_n, len(sub))
        rows.append(sub.head(n))
        print(f"  {locus}: {len(sub)} hits → selecting top {n}")
    return pd.concat(rows, ignore_index=True)


def download_fasta(accession: str, out_fasta: Path) -> bool:
    """Download a nucleotide record from NCBI as FASTA."""
    if out_fasta.exists() and out_fasta.stat().st_size > 1000:
        return True
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="fasta", retmode="text"
        )
        seq = handle.read()
        handle.close()
        if len(seq) < 100:
            return False
        out_fasta.write_text(seq)
        return True
    except Exception as e:
        print(f"    ERROR downloading {accession}: {e}", file=sys.stderr)
        return False


def run_kaptive_scores(db: Path, genome_files: list, scores_tsv: Path, threads: int) -> bool:
    """Run kaptive --scores on a list of genome files."""
    cmd = (
        ["kaptive", "assembly", str(db)]
        + [str(g) for g in genome_files]
        + ["--scores", str(scores_tsv), "-t", str(threads)]
    )
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("ERROR running kaptive:", result.stderr[:500], file=sys.stderr)
        return False
    return True


def normalise_scores(scores_tsv: Path, db: Path) -> pd.DataFrame:
    """Re-rank scores by AS / total_expected_gene_bp."""
    import re
    # Load locus CDS lengths from GenBank
    locus_bp = {}
    for rec in SeqIO.parse(str(db), "genbank"):
        name = None
        for f in rec.features:
            if f.type == "source":
                for note in f.qualifiers.get("note", []):
                    m = re.search(r"K locus:\s*(\S+)", note)
                    if m:
                        name = m.group(1)
        if name is None:
            name = rec.name.split("_")[0]
        locus_bp[name] = sum(len(f) for f in rec.features if f.type == "CDS")

    df = pd.read_csv(scores_tsv, sep="\t")
    asm_col = df.columns[0]
    df["total_bp"] = df["Locus"].map(locus_bp)
    mask = df["total_bp"].isna() | (df["total_bp"] == 0)
    df.loc[mask, "total_bp"] = df.loc[mask, "q_len"].clip(lower=1)
    df["AS_norm"] = df["AS"] / df["total_bp"]

    rows = []
    for asm, grp in df.groupby(asm_col):
        active = grp[grp["AS"] > 0]
        if active.empty:
            rows.append({"Assembly": asm, "Best_locus": "none", "AS_norm": 0})
            continue
        best = active.loc[active["AS_norm"].idxmax()]
        rows.append({
            "Assembly":   asm,
            "Best_locus": best["Locus"],
            "AS_norm":    round(best["AS_norm"], 4),
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-n",   type=int, default=10,
                        help="Candidates to test per locus (default 10; KL306 uses all)")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-kaptive",  action="store_true")
    args = parser.parse_args()

    GENOMES_DIR.mkdir(parents=True, exist_ok=True)

    # Load candidates — rebuild from per-locus files (the combined summary
    # gets overwritten when loci are run in parallel as separate processes)
    per_locus_files = [RESULTS_DIR / f"{l}_blast_hits.tsv" for l in TARGETS]
    missing = [f for f in per_locus_files if not f.exists()]
    if missing:
        print(f"ERROR: missing hit files: {missing}\nRun blast_ncbi_candidates.py first.",
              file=sys.stderr)
        sys.exit(1)

    all_candidates = pd.concat(
        [pd.read_csv(f, sep="\t") for f in per_locus_files],
        ignore_index=True
    )
    # Overwrite combined summary with the correct merged version
    all_candidates.to_csv(CANDIDATES, sep="\t", index=False)
    print(f"Loaded {len(all_candidates)} candidates across {TARGETS}")
    # Deduplicate: same accession may appear in multiple loci (shouldn't, but be safe)
    selected = pick_candidates(all_candidates, args.top_n)
    print(f"\nTotal candidates to test: {len(selected)}")

    # Download FASTAs
    print(f"\n[1] Downloading {len(selected)} genome FASTAs from NCBI...")
    genome_files, acc_to_locus = [], {}
    for _, row in selected.iterrows():
        acc   = row["accession"]
        locus = row["locus"]
        out   = GENOMES_DIR / f"{locus}_{acc}.fasta"

        if not args.skip_download:
            ok = download_fasta(acc, out)
            time.sleep(0.4)   # NCBI rate limit: ~3 requests/sec
        else:
            ok = out.exists()

        if ok and out.exists():
            genome_files.append(out)
            acc_to_locus[out.stem] = locus
            print(f"  {locus}  {acc}  ✓  ({out.stat().st_size // 1024} kb)")
        else:
            print(f"  {locus}  {acc}  FAILED")

    if not genome_files:
        print("No genome files to process.", file=sys.stderr)
        sys.exit(1)

    # Run Kaptive --scores
    scores_tsv  = RESULTS_DIR / "candidate_kaptive_scores.tsv"
    print(f"\n[2] Running Kaptive --scores on {len(genome_files)} assemblies...")
    if not args.skip_kaptive:
        ok = run_kaptive_scores(DATABASE, genome_files, scores_tsv, args.threads)
        if not ok:
            sys.exit(1)
    else:
        if not scores_tsv.exists():
            print("ERROR: --skip-kaptive but scores TSV not found", file=sys.stderr)
            sys.exit(1)

    # Normalise and assess
    print(f"\n[3] Normalising scores...")
    results = normalise_scores(scores_tsv, DATABASE)

    # Add expected locus from filename
    results["stem"] = results["Assembly"].apply(lambda x: Path(x).stem)
    results["expected_locus"] = results["stem"].map(acc_to_locus)
    results["correct"] = results["Best_locus"] == results["expected_locus"]

    out_tsv = RESULTS_DIR / "candidate_typing_summary.tsv"
    results.to_csv(out_tsv, sep="\t", index=False)

    # Report
    print("\n" + "=" * 65)
    print("CANDIDATE TYPING RESULTS  (v0.4 database, normalised scoring)")
    print("=" * 65)

    for locus in TARGETS:
        sub = results[results["expected_locus"] == locus]
        n_correct = int(sub["correct"].sum())
        print(f"\n{locus}:  {n_correct}/{len(sub)} type correctly")
        print(f"  {'Assembly':<35}  {'Called':>8}  {'Norm AS':>8}  {'OK?'}")
        print(f"  {'-'*35}  {'------':>8}  {'-------':>8}  {'---'}")
        for _, row in sub.sort_values("AS_norm", ascending=False).iterrows():
            ok = "✓ CORRECT" if row["correct"] else f"✗ → {row['Best_locus']}"
            print(f"  {row['stem']:<35}  {row['Best_locus']:>8}  "
                  f"{row['AS_norm']:>8.4f}  {ok}")

    total_correct = int(results["correct"].sum())
    print(f"\n{'='*65}")
    print(f"Summary: {total_correct}/{len(results)} candidates type correctly")
    print(f"\nAssemblies that type correctly = candidate replacement representatives")
    print(f"Written: {out_tsv}")


if __name__ == "__main__":
    main()
