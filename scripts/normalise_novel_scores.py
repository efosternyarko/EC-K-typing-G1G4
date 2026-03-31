#!/usr/bin/env python3
"""
normalise_novel_scores.py

Normalised self-typing validation for novel ATB loci (KL392-KL968).

Reads the Kaptive --scores matrix and re-ranks each assembly by:
    normalised_score = AS / total_expected_gene_bp

Then checks whether each representative assembly correctly types to
its expected KL type.

Usage:
    python3 scripts/normalise_novel_scores.py \
        --scores  DB/kaptive_scores_novel_v0.8norm.tsv \
        --db      DB/EC-K-typing_group1and4_v0.8.gbk \
        --mapping DB/novel_rep_kl_map.tsv \
        --out     DB/novel_validation_summary_v0.8norm.tsv
"""

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

MIN_GENE_COVERAGE = 0.50


def get_locus_name(rec):
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def load_locus_stats(db: Path):
    locus_total_bp = {}
    for rec in SeqIO.parse(str(db), "genbank"):
        ln = get_locus_name(rec)
        cds = [f for f in rec.features if f.type == "CDS"]
        locus_total_bp[ln] = sum(len(f) for f in cds)
    return locus_total_bp


def normalise_and_rank(scores_tsv: Path, locus_total_bp: dict) -> pd.DataFrame:
    df = pd.read_csv(scores_tsv, sep="\t")
    asm_col = df.columns[0]

    df["total_expected_bp"] = df["Locus"].map(locus_total_bp)
    mask = df["total_expected_bp"].isna() | (df["total_expected_bp"] == 0)
    df.loc[mask, "total_expected_bp"] = df.loc[mask, "q_len"].clip(lower=1)
    df["AS_norm"] = df["AS"] / df["total_expected_bp"]

    rows = []
    for asm, grp in df.groupby(asm_col):
        active = grp[grp["AS"] > 0]
        if active.empty:
            rows.append({"Assembly": asm, "Best_locus": "none",
                         "Confidence": "Untypeable", "Norm_AS": 0.0,
                         "Genes_found": 0, "Genes_expected": 0})
            continue
        top_norm = active["AS_norm"].max()
        candidates = active[active["AS_norm"] >= top_norm * 0.995]
        best = candidates.sort_values(
            ["genes_found", "AS_norm"], ascending=[False, False]
        ).iloc[0]
        genes_found    = int(best["genes_found"])
        genes_expected = int(best["genes_expected"])
        coverage = genes_found / genes_expected if genes_expected > 0 else 0
        conf = "Typeable" if coverage >= MIN_GENE_COVERAGE else "Untypeable"
        rows.append({
            "Assembly":       asm,
            "Best_locus":     best["Locus"],
            "Confidence":     conf,
            "Norm_AS":        round(best["AS_norm"], 4),
            "Genes_found":    genes_found,
            "Genes_expected": genes_expected,
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores",  required=True, type=Path)
    parser.add_argument("--db",      required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--out",     required=True, type=Path)
    args = parser.parse_args()

    print(f"Loading locus stats from {args.db.name}...")
    locus_total_bp = load_locus_stats(args.db)
    print(f"  {len(locus_total_bp)} loci loaded")

    # Load expected KL per assembly
    mapping = pd.read_csv(args.mapping, sep="\t")
    # Assembly column in scores TSV uses the full path stem
    mapping["stem"] = mapping["fasta_path"].apply(lambda x: Path(x).stem)
    expected_kl = dict(zip(mapping["stem"], mapping["KL_type"]))

    print(f"Normalising scores from {args.scores.name}...")
    results = normalise_and_rank(args.scores, locus_total_bp)

    # Match assembly stem to expected KL
    results["stem"] = results["Assembly"].apply(lambda x: Path(x).stem)
    results["Expected_KL"] = results["stem"].map(expected_kl)
    results["Correct"] = results.apply(
        lambda r: r["Best_locus"] == r["Expected_KL"]
        if pd.notna(r["Expected_KL"]) else None,
        axis=1
    )
    results["Typeable"] = results["Confidence"] == "Typeable"

    results.to_csv(args.out, sep="\t", index=False)
    print(f"Written: {args.out}")

    # Summary
    reps = results[results["Expected_KL"].notna()]
    n_total   = len(reps)
    n_correct = int(reps["Correct"].sum())
    n_typeable = int(results["Typeable"].sum())

    print(f"\n{'='*55}")
    print(f"NOVEL LOCI SELF-TYPING SUMMARY (v0.8 normalised scoring)")
    print(f"{'='*55}")
    print(f"  Novel loci validated:    {n_total}")
    print(f"  Correctly self-typed:    {n_correct}/{n_total}  ({100*n_correct/n_total:.1f}%)")
    print(f"  Typeable:                {n_typeable}/{len(results)}  ({100*n_typeable/len(results):.1f}%)")

    wrong = reps[reps["Correct"] == False][
        ["stem", "Expected_KL", "Best_locus", "Norm_AS"]
    ].sort_values("Expected_KL")
    if not wrong.empty:
        print(f"\n  Incorrect ({len(wrong)}):")
        for _, row in wrong.iterrows():
            print(f"    {row['stem']:<30}  expected {row['Expected_KL']:<10}  "
                  f"got {row['Best_locus']}")
    else:
        print("\n  All novel loci self-typed correctly.")


if __name__ == "__main__":
    main()
