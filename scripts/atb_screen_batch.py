#!/usr/bin/env python3
"""
atb_screen_batch.py — Screen one AllTheBacteria batch for G1/G4 flanking genes.

For each genome in the batch, detects presence of galF and gnd (the conserved
flanking genes that bracket the E. coli Group 1/4 cps locus). Genomes with
both markers are G1/G4 candidates for locus extraction.

Usage:
    python atb_screen_batch.py <batch_number>

Output:
    <OUT_DIR>/batch_<N>_hits.tsv
    Columns: genome, gene, contig, pident, qcov

Post-processing (combine_atb_hits.py) filters for genomes with both
galF AND gnd hits, then runs locus extraction on candidates.
"""

import os
import subprocess
import sys
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
ARCHIVES_DIR = "/scratch2/js66/atb_archives"
FLANKING_FA  = "/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/flanking_genes/flanking_genes.fasta"
OUT_DIR      = "/home/ebenezef/js66_scratch/ebenn/atb_screen/results"

PIDENT_THRESH = 80.0   # % nucleotide identity
QCOV_THRESH   = 70.0   # % query coverage
BLAST_EVALUE  = 1e-10
# ─────────────────────────────────────────────────────────────────────────────


def build_combined_fasta(batch_dir: Path, out_fa: str) -> dict:
    """
    Concatenate all genomes in the batch into a single FASTA, prefixing each
    contig ID with its genome ID so we can map hits back after BLAST.
    Returns {prefixed_contig_id: genome_id}.
    """
    genome_map = {}
    with open(out_fa, "w") as out:
        for fa_path in sorted(batch_dir.glob("*.fa")):
            genome_id = fa_path.stem
            with open(fa_path) as f:
                for line in f:
                    if line.startswith(">"):
                        contig = line.strip()[1:].split()[0]
                        prefixed = f"{genome_id}___{contig}"
                        genome_map[prefixed] = (genome_id, contig)
                        out.write(f">{prefixed}\n")
                    else:
                        out.write(line)
    return genome_map


def run_blast(query: str, db: str) -> str:
    result = subprocess.run(
        [
            "blastn",
            "-query", query,
            "-db", db,
            "-evalue", str(BLAST_EVALUE),
            "-outfmt", "6 qseqid sseqid pident length qlen",
            "-perc_identity", "75",        # pre-filter; thresholds applied below
            "-max_target_seqs", "100000",  # high cap — we want all hits
        ],
        capture_output=True, text=True
    )
    return result.stdout


def parse_hits(blast_output: str, genome_map: dict) -> dict:
    """
    Parse BLAST tabular output. Returns:
    {genome_id: {gene: {"pident": float, "qcov": float, "contig": str}}}
    Keeps the best (highest pident) hit per gene per genome.
    """
    hits = defaultdict(dict)
    for line in blast_output.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        gene    = parts[0]
        seq_id  = parts[1]
        pident  = float(parts[2])
        aln_len = int(parts[3])
        qlen    = int(parts[4])
        qcov    = 100.0 * aln_len / qlen

        if pident < PIDENT_THRESH or qcov < QCOV_THRESH:
            continue

        genome_id, contig = genome_map.get(seq_id, (seq_id, seq_id))

        if gene not in hits[genome_id] or pident > hits[genome_id][gene]["pident"]:
            hits[genome_id][gene] = {"pident": pident, "qcov": qcov, "contig": contig}

    return hits


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: atb_screen_batch.py <batch_number>")

    batch_num = sys.argv[1]
    archive   = f"{ARCHIVES_DIR}/atb.assembly.r0.2.batch.{batch_num}.tar.xz"

    if not os.path.exists(archive):
        sys.exit(f"Archive not found: {archive}")

    os.makedirs(OUT_DIR, exist_ok=True)
    out_file = os.path.join(OUT_DIR, f"batch_{batch_num}_hits.tsv")

    # Skip if already done (allows safe re-submission)
    if os.path.exists(out_file):
        print(f"Batch {batch_num}: already processed, skipping.")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Batch {batch_num}: extracting {archive} ...", flush=True)
        with tarfile.open(archive, "r:xz") as tar:
            tar.extractall(tmpdir)

        batch_dir = next(Path(tmpdir).iterdir())
        n_genomes = len(list(batch_dir.glob("*.fa")))
        print(f"Batch {batch_num}: {n_genomes} genomes extracted.", flush=True)

        # Build combined FASTA and BLAST DB
        combined_fa = os.path.join(tmpdir, "combined.fa")
        genome_map  = build_combined_fasta(batch_dir, combined_fa)

        db = os.path.join(tmpdir, "batch_db")
        subprocess.run(
            ["makeblastdb", "-in", combined_fa, "-dbtype", "nucl", "-out", db],
            check=True, capture_output=True
        )

        # Run BLAST (one query against combined DB — much faster than per-genome)
        blast_output = run_blast(FLANKING_FA, db)

        hits = parse_hits(blast_output, genome_map)

    # Write results
    with open(out_file, "w") as f:
        f.write("genome\tgene\tcontig\tpident\tqcov\n")
        for genome_id in sorted(hits):
            for gene, h in sorted(hits[genome_id].items()):
                f.write(f"{genome_id}\t{gene}\t{h['contig']}\t{h['pident']:.1f}\t{h['qcov']:.1f}\n")

    n_hits = len(hits)
    n_both = sum(1 for g in hits if "galF" in hits[g] and "gnd" in hits[g])
    print(f"Batch {batch_num}: {n_hits} genomes with any hit, {n_both} with both galF+gnd (G1/G4 candidates).")


if __name__ == "__main__":
    main()
