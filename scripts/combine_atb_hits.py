#!/usr/bin/env python3
"""
combine_atb_hits.py — Combine per-batch screening results and identify G1/G4 candidates.

Reads all batch_*_hits.tsv files from the results directory, filters for
genomes with both galF AND gnd hits, and writes:
  - atb_g1g4_candidates.tsv  : genomes with both flanking markers (locus extraction targets)
  - atb_all_hits.tsv         : all hits across all batches (for QC)
  - atb_screen_summary.txt   : summary statistics

Usage:
    python combine_atb_hits.py
"""

import csv
import os
from collections import defaultdict
from pathlib import Path

RESULTS_DIR  = "/home/ebenezef/js66_scratch/ebenn/atb_screen/results"
OUT_DIR      = "/home/ebenezef/js66_scratch/ebenn/atb_screen"
ARCHIVES_DIR = "/scratch2/js66/atb_archives"
N_BATCHES    = 888

# Require both of these to call a genome a G1/G4 candidate
REQUIRED_GENES = {"galF", "gnd"}


def main():
    results_dir = Path(RESULTS_DIR)
    hit_files   = sorted(results_dir.glob("batch_*_hits.tsv"))

    print(f"Found {len(hit_files)} completed batch result files (of {N_BATCHES} total).")

    if len(hit_files) < N_BATCHES:
        missing = N_BATCHES - len(hit_files)
        print(f"  WARNING: {missing} batches not yet processed.")

    # genome_id -> {gene -> {pident, qcov, contig, batch}}
    all_hits = defaultdict(dict)

    for hit_file in hit_files:
        batch_num = hit_file.stem.replace("batch_", "").replace("_hits", "")
        with open(hit_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                genome = row["genome"]
                gene   = row["gene"]
                if gene not in all_hits[genome] or float(row["pident"]) > all_hits[genome][gene]["pident"]:
                    all_hits[genome][gene] = {
                        "pident": float(row["pident"]),
                        "qcov":   float(row["qcov"]),
                        "contig": row["contig"],
                        "batch":  batch_num,
                    }

    # Filter: genomes with both galF AND gnd
    candidates = {
        gid: genes for gid, genes in all_hits.items()
        if REQUIRED_GENES.issubset(genes.keys())
    }

    # Write candidates TSV
    cand_out = os.path.join(OUT_DIR, "atb_g1g4_candidates.tsv")
    with open(cand_out, "w") as f:
        f.write("genome\tbatch\tgalF_pident\tgalF_qcov\tgalF_contig\tgnd_pident\tgnd_qcov\tgnd_contig\tsame_contig\n")
        for genome_id in sorted(candidates):
            genes = candidates[genome_id]
            galf  = genes.get("galF", {})
            gnd   = genes.get("gnd", {})
            same  = galf.get("contig") == gnd.get("contig")
            f.write(
                f"{genome_id}\t{galf.get('batch','')}\t"
                f"{galf.get('pident','')}\t{galf.get('qcov','')}\t{galf.get('contig','')}\t"
                f"{gnd.get('pident','')}\t{gnd.get('qcov','')}\t{gnd.get('contig','')}\t"
                f"{same}\n"
            )

    # Write all-hits TSV
    all_out = os.path.join(OUT_DIR, "atb_all_hits.tsv")
    with open(all_out, "w") as f:
        f.write("genome\tgene\tcontig\tpident\tqcov\tbatch\n")
        for genome_id in sorted(all_hits):
            for gene, h in sorted(all_hits[genome_id].items()):
                f.write(f"{genome_id}\t{gene}\t{h['contig']}\t{h['pident']}\t{h['qcov']}\t{h['batch']}\n")

    # Summary
    n_any      = len(all_hits)
    n_cands    = len(candidates)
    n_same_ctg = sum(
        1 for g in candidates
        if candidates[g].get("galF", {}).get("contig") == candidates[g].get("gnd", {}).get("contig")
    )

    summary = [
        f"AllTheBacteria G1/G4 screen summary",
        f"  Batches processed     : {len(hit_files)} / {N_BATCHES}",
        f"  Genomes with any hit  : {n_any}",
        f"  G1/G4 candidates      : {n_cands}  (galF + gnd both present)",
        f"  Both on same contig   : {n_same_ctg}  (highest confidence)",
        f"  Candidates file       : {cand_out}",
        f"  All hits file         : {all_out}",
    ]
    print("\n".join(summary))
    with open(os.path.join(OUT_DIR, "atb_screen_summary.txt"), "w") as f:
        f.write("\n".join(summary) + "\n")

    print(f"\nNext step: run extract_atb_loci.py on {cand_out}")


if __name__ == "__main__":
    main()
