#!/usr/bin/env python3
"""
extract_novel_loci.py — Extract G1/G4 cps loci from novel ATB candidates.

For each same-contig G1/G4 candidate in a given batch:
  1. Screen contig against existing 93 G1/G4 reference loci (blastn)
     → if ≥95% identity AND ≥80% query coverage → known type, skip
  2. For novel candidates: BLAST flanking genes vs contig to get coordinates
  3. Extract locus sequence between galF_end and gnd_start (±500 bp flanking)
  4. Write to output FASTA

Output per batch:
  <OUT_DIR>/batch_<N>_novel_loci.fasta   — extracted novel loci
  <OUT_DIR>/batch_<N>_extract_summary.tsv — per-genome result

Usage:
    python extract_novel_loci.py <batch_number>
"""

import csv
import os
import subprocess
import sys
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
ARCHIVES_DIR  = "/scratch2/js66/atb_archives"
CANDIDATES    = "/home/ebenezef/js66_scratch/ebenn/atb_screen/atb_g1g4_candidates.tsv"
FLANKING_FA   = "/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/flanking_genes/flanking_genes.fasta"
NOVELTY_DB    = "/home/ebenezef/js66_scratch/ebenn/atb_screen/novelty_db/g1g4_refs"
OUT_DIR       = "/home/ebenezef/js66_scratch/ebenn/atb_screen/novel_loci"

# Novelty thresholds — match at BOTH means "known type, skip"
KNOWN_PIDENT  = 95.0   # % nucleotide identity
KNOWN_QCOV    = 80.0   # % query coverage

# Flanking gene BLAST thresholds for coordinate detection
FLANK_PIDENT  = 80.0
FLANK_QCOV    = 70.0

# Locus extraction padding
FLANK_PAD     = 500    # bp added each side of galF..gnd region
MIN_LOCUS_LEN = 15000  # discard extractions shorter than this (fragments)
# ─────────────────────────────────────────────────────────────────────────────


def load_candidates(batch_num: str) -> list:
    """Load same-contig candidates for this batch from the combined TSV."""
    candidates = []
    with open(CANDIDATES) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["batch"] == batch_num and row["same_contig"] == "True":
                candidates.append(row)
    return candidates


def is_known_type(contig_fa: str, genome_id: str) -> tuple:
    """
    BLAST contig vs existing G1/G4 references.
    Returns (is_known: bool, best_match: str, pident: float, qcov: float).
    """
    result = subprocess.run(
        [
            "blastn",
            "-query", contig_fa,
            "-db", NOVELTY_DB,
            "-evalue", "1e-10",
            "-outfmt", "6 sseqid pident length qlen",
            "-max_target_seqs", "1",
            "-perc_identity", "85",
        ],
        capture_output=True, text=True
    )
    best_match, best_pident, best_qcov = "none", 0.0, 0.0
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        ref_id  = parts[0]
        pident  = float(parts[1])
        aln_len = int(parts[2])
        qlen    = int(parts[3])
        qcov    = 100.0 * aln_len / qlen
        if pident > best_pident:
            best_pident = pident
            best_qcov   = qcov
            best_match  = ref_id
    is_known = (best_pident >= KNOWN_PIDENT) and (best_qcov >= KNOWN_QCOV)
    return is_known, best_match, best_pident, best_qcov


def get_flanking_coords(contig_fa: str) -> dict:
    """
    BLAST flanking genes vs contig to get galF and gnd coordinates.
    Returns {gene: {"start": int, "end": int, "strand": str}}.
    Coordinates are 1-based, on the contig.
    """
    result = subprocess.run(
        [
            "blastn",
            "-query", FLANKING_FA,
            "-subject", contig_fa,
            "-evalue", "1e-10",
            "-outfmt", "6 qseqid pident length qlen sstart send",
            "-perc_identity", "75",
            "-max_target_seqs", "3",
        ],
        capture_output=True, text=True
    )
    coords = {}
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        gene    = parts[0]
        pident  = float(parts[1])
        aln_len = int(parts[2])
        qlen    = int(parts[3])
        sstart  = int(parts[4])
        send    = int(parts[5])
        qcov    = 100.0 * aln_len / qlen
        if pident < FLANK_PIDENT or qcov < FLANK_QCOV:
            continue
        # Keep best hit per gene
        if gene not in coords or pident > coords[gene]["pident"]:
            start  = min(sstart, send)
            end    = max(sstart, send)
            strand = "+" if send >= sstart else "-"
            coords[gene] = {"start": start, "end": end,
                            "strand": strand, "pident": pident}
    return coords


def extract_locus(contig_fa: str, coords: dict, genome_id: str) -> str:
    """
    Extract locus sequence from galF_end to gnd_start with FLANK_PAD padding.
    Returns the extracted sequence as a string, or "" if coordinates are invalid.
    """
    if "galF" not in coords or "gnd" not in coords:
        return ""

    # Determine locus boundaries regardless of strand
    galF_end   = coords["galF"]["end"]
    gnd_start  = coords["gnd"]["start"]

    # Ensure galF is upstream of gnd
    locus_start = min(galF_end, gnd_start)
    locus_end   = max(galF_end, gnd_start)

    if locus_end - locus_start < 5000:
        return ""  # suspiciously short — likely same gene hit twice

    # Read contig sequence
    seq = ""
    with open(contig_fa) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()

    contig_len = len(seq)
    extract_start = max(0, locus_start - FLANK_PAD - 1)  # to 0-based
    extract_end   = min(contig_len, locus_end + FLANK_PAD)

    locus_seq = seq[extract_start:extract_end]

    if len(locus_seq) < MIN_LOCUS_LEN:
        return ""

    return locus_seq


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: extract_novel_loci.py <batch_number>")

    batch_num = sys.argv[1]
    archive   = f"{ARCHIVES_DIR}/atb.assembly.r0.2.batch.{batch_num}.tar.xz"

    if not os.path.exists(archive):
        print(f"Archive not found: {archive}")
        return

    os.makedirs(OUT_DIR, exist_ok=True)

    out_fasta   = os.path.join(OUT_DIR, f"batch_{batch_num}_novel_loci.fasta")
    out_summary = os.path.join(OUT_DIR, f"batch_{batch_num}_extract_summary.tsv")

    if os.path.exists(out_fasta):
        print(f"Batch {batch_num}: already processed, skipping.")
        return

    candidates = load_candidates(batch_num)
    if not candidates:
        print(f"Batch {batch_num}: no same-contig candidates, skipping.")
        # Write empty outputs so batch is marked done
        open(out_fasta, "w").close()
        with open(out_summary, "w") as f:
            f.write("genome\tcontig\tstatus\tbest_ref\tpident\tqcov\tlocus_len\n")
        return

    print(f"Batch {batch_num}: {len(candidates)} same-contig candidates to process.", flush=True)

    # Build set of genomes+contigs we need
    needed = {c["genome"]: c["galF_contig"] for c in candidates}

    n_novel = 0
    n_known = 0
    n_failed = 0
    summary_rows = []

    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract only the needed genome FASTAs from the archive
        print(f"Batch {batch_num}: extracting {len(needed)} genomes from archive...", flush=True)
        with tarfile.open(archive, "r:xz") as tar:
            members = tar.getmembers()
            to_extract = [m for m in members
                          if any(m.name.endswith(f"{gid}.fa") for gid in needed)]
            tar.extractall(tmpdir, members=to_extract)

        batch_dir = next(Path(tmpdir).iterdir())

        with open(out_fasta, "w") as fa_out, open(out_summary, "w") as sum_out:
            sum_out.write("genome\tcontig\tstatus\tbest_ref\tpident\tqcov\tlocus_len\n")

            for candidate in candidates:
                genome_id    = candidate["genome"]
                target_contig = candidate["galF_contig"]  # same as gnd_contig

                genome_fa = batch_dir / f"{genome_id}.fa"
                if not genome_fa.exists():
                    sum_out.write(f"{genome_id}\t{target_contig}\tarchive_missing\t\t\t\t\n")
                    n_failed += 1
                    continue

                # Write target contig to temp file
                contig_fa = os.path.join(tmpdir, f"{genome_id}_contig.fa")
                found_contig = False
                with open(genome_fa) as fin, open(contig_fa, "w") as fout:
                    write = False
                    for line in fin:
                        if line.startswith(">"):
                            contig_name = line.strip()[1:].split()[0]
                            write = (contig_name == target_contig)
                            if write:
                                found_contig = True
                        if write:
                            fout.write(line)

                if not found_contig:
                    sum_out.write(f"{genome_id}\t{target_contig}\tcontig_not_found\t\t\t\t\n")
                    n_failed += 1
                    continue

                # Novelty check
                is_known, best_ref, pident, qcov = is_known_type(contig_fa, genome_id)

                if is_known:
                    sum_out.write(f"{genome_id}\t{target_contig}\tknown\t{best_ref}\t{pident:.1f}\t{qcov:.1f}\t\n")
                    n_known += 1
                    continue

                # Get flanking coordinates for novel candidates
                coords = get_flanking_coords(contig_fa)
                locus_seq = extract_locus(contig_fa, coords, genome_id)

                if not locus_seq:
                    sum_out.write(f"{genome_id}\t{target_contig}\textract_failed\t{best_ref}\t{pident:.1f}\t{qcov:.1f}\t\n")
                    n_failed += 1
                    continue

                # Write to output FASTA
                galF_pi = coords.get("galF", {}).get("pident", 0)
                gnd_pi  = coords.get("gnd",  {}).get("pident", 0)
                header  = (f">{genome_id} novel_locus {len(locus_seq)}bp "
                           f"[galF({galF_pi:.1f}%)--gnd({gnd_pi:.1f}%)] "
                           f"best_ref={best_ref}({pident:.1f}%id,{qcov:.1f}%cov)")
                fa_out.write(f"{header}\n{locus_seq}\n")
                sum_out.write(f"{genome_id}\t{target_contig}\tnovel\t{best_ref}\t{pident:.1f}\t{qcov:.1f}\t{len(locus_seq)}\n")
                n_novel += 1

    print(f"Batch {batch_num}: {n_novel} novel loci extracted, "
          f"{n_known} known types skipped, {n_failed} failed.")


if __name__ == "__main__":
    main()
