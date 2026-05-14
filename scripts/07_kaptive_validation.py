#!/Users/lshef4/lshef4_to_copy/anaconda3/envs/kleborate_test/bin/python3
"""
07_kaptive_validation.py

Self-typing validation: type each v1.3 locus FASTA back against the v1.3
database and check it gets the correct call.

Usage:
  python3 scripts/07_kaptive_validation.py

Requires:
  - kaptive (python3 from kleborate_test env)
  - minimap2 in PATH (from kleborate_test env)

Outputs:
  - DB/kaptive_validation/kaptive_self_typing.tsv   full kaptive table
  - DB/kaptive_validation/self_typing_summary.tsv   pass/fail per locus
"""

import sys
import os
import csv
from pathlib import Path

KAPTIVE_ENV = Path("/Users/lshef4/lshef4_to_copy/anaconda3/envs/kleborate_test")
os.environ["PATH"] = str(KAPTIVE_ENV / "bin") + ":" + os.environ.get("PATH", "")

from kaptive.database import load_database
from kaptive.assembly import typing_pipeline, write_headers

REPO = Path(__file__).resolve().parent.parent
DB_PATH = REPO / "DB" / "EC-K-typing_group1and4_v1.3.gbk"
FASTA_DIR = REPO / "DB" / "kaptive_validation" / "loci_fasta"
VAL_DIR = REPO / "DB" / "kaptive_validation"
OUT_TSV = VAL_DIR / "kaptive_self_typing.tsv"
SUMMARY_TSV = VAL_DIR / "self_typing_summary.tsv"

VAL_DIR.mkdir(parents=True, exist_ok=True)

print(f"Loading v1.3 database...")
db = load_database(str(DB_PATH), load_locus_seqs=True, verbose=False)
print(f"  {len(db.loci)} loci loaded")

fasta_files = sorted(FASTA_DIR.glob("*.fna"))
print(f"  {len(fasta_files)} FASTA files to type")

results = []
n_pass = 0
n_fail = 0
n_no_result = 0

with open(OUT_TSV, "w") as tsv_out:
    write_headers(tsv_out, no_header=False, scores=None)

    for i, fasta in enumerate(fasta_files, 1):
        if i % 50 == 0 or i == 1:
            print(f"  [{i}/{len(fasta_files)}] {fasta.stem}...")

        result = typing_pipeline(str(fasta), db, threads=4, verbose=False)

        if result is None:
            n_no_result += 1
            results.append({
                "locus": fasta.stem,
                "best_match": "NO_RESULT",
                "expected_kl": fasta.stem.split("_")[0],
                "pass": "FAIL",
                "confidence": "",
                "identity": "",
                "coverage": "",
            })
            continue

        result.write(tsv_out, None, None, None, None, None, None)

        # Parse result — best_locus is the KL name in the database
        expected_kl = fasta.stem.split("_")[0]   # e.g. KL300 from KL300_ESC_...
        best_match = str(result.best_match) if result.best_match else "NONE"

        passed = best_match == expected_kl
        if passed:
            n_pass += 1
        else:
            n_fail += 1

        results.append({
            "locus": fasta.stem,
            "best_match": best_match,
            "expected_kl": expected_kl,
            "pass": "PASS" if passed else "FAIL",
            "confidence": str(result.confidence),
            "identity": f"{result.percent_identity:.2f}%",
            "coverage": f"{result.percent_coverage:.2f}%",
        })

print(f"\nWriting results...")
with open(SUMMARY_TSV, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["locus","expected_kl","best_match","pass",
                                             "confidence","identity","coverage"],
                            delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

print(f"\n{'='*60}")
print("KAPTIVE SELF-TYPING VALIDATION SUMMARY")
print(f"{'='*60}")
print(f"  Loci typed:       {len(fasta_files)}")
print(f"  PASS (correct):   {n_pass}  ({100*n_pass/len(fasta_files):.1f}%)")
print(f"  FAIL (wrong KL):  {n_fail}  ({100*n_fail/len(fasta_files):.1f}%)")
print(f"  No result:        {n_no_result}")
print(f"\nOutputs:")
print(f"  {OUT_TSV}")
print(f"  {SUMMARY_TSV}")

if n_fail > 0:
    print(f"\nFailed loci:")
    for r in results:
        if r["pass"] == "FAIL":
            print(f"  {r['locus']}  →  {r['best_match']}  (expected {r['expected_kl']})")
