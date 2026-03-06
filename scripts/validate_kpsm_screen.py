#!/usr/bin/env python3
"""
validate_kpsm_screen.py — Batch kpsM screen for validation.

Runs a fast kpsM BLAST screen across a directory of genome assemblies
(plain FASTA or .gz) using blastn -subject (no makeblastdb needed per
genome). Parallelised with multiprocessing.

Expected results for BSI bacteraemia genomes (pre-screened as G2/G3
no-hits):
  - All should be kpsM-negative (Group 1/4 candidates)
  - Any kpsM-positive hits are unexpected and flagged for inspection.

Usage:
    python validate_kpsm_screen.py \\
        --genome-dir /path/to/genomes/ \\
        --kpsm-ref DB/kpsM_reference.fasta \\
        --output results/kpsm_validation.tsv \\
        --workers 8
"""

import argparse
import gzip
import multiprocessing
import os
import subprocess
import sys
import tempfile
from pathlib import Path

KPSM_PIDENT_THRESH = 90.0
KPSM_QCOV_THRESH   = 80.0
BLAST_EVALUE       = 1e-10
BLAST_BIN          = os.path.expanduser("~/.local/bin/blastn")


def decompress_to_tmp(gz_path: Path, tmpdir: str) -> str:
    """Decompress a .gz FASTA to a temp file; return path. Pass-through if not .gz."""
    if str(gz_path).endswith('.gz'):
        out = os.path.join(tmpdir, gz_path.stem)
        with gzip.open(gz_path, 'rb') as fin, open(out, 'wb') as fout:
            fout.write(fin.read())
        return out
    return str(gz_path)


def screen_one(args):
    """Worker: screen one genome. Returns result dict."""
    genome_path, kpsm_ref, pident_thresh, qcov_thresh = args
    genome_path = Path(genome_path)
    name = genome_path.name.replace('.fa.gz', '').replace('.fasta.gz', '').replace('.fna.gz', '')

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            fa = decompress_to_tmp(genome_path, tmpdir)
        except Exception as e:
            return {"assembly": name, "kpsM_hit": False,
                    "gene": "ERROR", "pident": "NA", "qcov": "NA",
                    "group": f"ERROR:{e}"}

        result = subprocess.run(
            [
                BLAST_BIN,
                "-query", str(kpsm_ref),
                "-subject", fa,
                "-evalue", str(BLAST_EVALUE),
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send qlen",
                "-max_target_seqs", "5",
                "-perc_identity", "70",   # low pre-filter; apply thresholds in Python
            ],
            capture_output=True, text=True, timeout=60
        )

    best_pident, best_qcov, best_gene = None, None, None
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        gene    = parts[0]
        pident  = float(parts[2])
        aln_len = int(parts[3])
        qlen    = int(parts[8])
        qcov    = 100.0 * aln_len / qlen

        if best_pident is None or pident > best_pident:
            best_pident = pident
            best_qcov   = qcov
            best_gene   = gene

    if best_pident is None:
        return {"assembly": name, "kpsM_hit": False,
                "gene": "NA", "pident": "NA", "qcov": "NA", "group": "G1_G4_candidate"}

    hit = (best_pident >= pident_thresh) and (best_qcov >= qcov_thresh)
    return {
        "assembly": name,
        "kpsM_hit": hit,
        "gene":     best_gene,
        "pident":   f"{best_pident:.1f}",
        "qcov":     f"{best_qcov:.1f}",
        "group":    "Group2_3" if hit else "G1_G4_candidate",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Batch kpsM screen for validation against large assembly collections."
    )
    parser.add_argument("--genome-dir", type=Path, required=True,
                        help="Directory of genome FASTAs (.fa, .fa.gz, .fasta, .fna).")
    parser.add_argument("--kpsm-ref", type=Path,
                        default=Path(__file__).parent.parent / "DB" / "kpsM_reference.fasta",
                        help="kpsM/kpsM_3 reference FASTA.")
    parser.add_argument("--output", type=Path, default=Path("kpsm_validation.tsv"),
                        help="Output TSV path.")
    parser.add_argument("--workers", type=int, default=8,
                        help="Parallel workers (default: 8).")
    parser.add_argument("--pident", type=float, default=KPSM_PIDENT_THRESH)
    parser.add_argument("--qcov",   type=float, default=KPSM_QCOV_THRESH)
    args = parser.parse_args()

    pident_thresh = args.pident
    qcov_thresh   = args.qcov

    # Collect genomes
    exts = ('.fa', '.fa.gz', '.fasta', '.fasta.gz', '.fna', '.fna.gz')
    genomes = [p for p in sorted(args.genome_dir.iterdir())
               if p.name.endswith(exts) and not p.name.startswith('.')]

    if not genomes:
        sys.exit(f"No genome files found in {args.genome_dir}")

    print(f"kpsM batch validation")
    print(f"  Genomes      : {len(genomes)}")
    print(f"  kpsM ref     : {args.kpsm_ref}")
    print(f"  Thresholds   : pident >= {KPSM_PIDENT_THRESH}%, qcov >= {KPSM_QCOV_THRESH}%")
    print(f"  Workers      : {args.workers}")
    print(f"  Output       : {args.output}\n")

    tasks = [(str(g), str(args.kpsm_ref), pident_thresh, qcov_thresh) for g in genomes]

    results = []
    done = 0
    with multiprocessing.Pool(processes=args.workers) as pool:
        for res in pool.imap_unordered(screen_one, tasks, chunksize=20):
            results.append(res)
            done += 1
            if done % 500 == 0 or done == len(genomes):
                n_pos = sum(1 for r in results if r["kpsM_hit"])
                print(f"  [{done}/{len(genomes)}] kpsM+ so far: {n_pos}", flush=True)

    # Write TSV
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cols = ["assembly", "kpsM_hit", "gene", "pident", "qcov", "group"]
    with open(args.output, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in sorted(results, key=lambda x: x["assembly"]):
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # Summary
    n_pos  = sum(1 for r in results if r["kpsM_hit"])
    n_neg  = len(results) - n_pos
    print(f"\n{'='*50}")
    print(f"  Total screened        : {len(results)}")
    print(f"  kpsM+ (Group 2/3)     : {n_pos}  {'<-- UNEXPECTED, investigate' if n_pos > 0 else '(expected: 0)'}")
    print(f"  kpsM- (G1/G4 cand.)  : {n_neg}")
    print(f"  Results written to    : {args.output}")

    if n_pos > 0:
        print(f"\n  kpsM+ assemblies:")
        for r in results:
            if r["kpsM_hit"]:
                print(f"    {r['assembly']}  gene={r['gene']}  pident={r['pident']}  qcov={r['qcov']}")


if __name__ == "__main__":
    main()
