#!/usr/bin/env python3
"""
replace_subset_reps.py — Replace GenBank records for KL601, KL713, KL742
with better representatives identified from AllTheBacteria.

Strategy:
  KL601: BLAST KL414 (superset, defines full locus extent) against SAMN12943639
         to find the full boundary, then extract and verify gene count > 17 CDS.
  KL713: BLAST current KL713 query against SAMN07152258; extract the single-contig
         region (no N-spacer) and verify.
  KL742: BLAST KL394 (superset, defines full locus extent) against SAMN41228912.

For each locus:
  1. BLAST reference sequence against candidate assembly to define boundaries
  2. Extract precise locus sequence (no flanking buffer)
  3. Annotate CDS with pyrodigal
  4. Build a new GenBank record matching the existing record format
  5. Write updated GenBank: DB/EC-K-typing_group1and4_v1.1.gbk

Usage:
    ~/miniforge3/bin/python3 scripts/replace_subset_reps.py

Inputs:
    DB/EC-K-typing_group1and4_v1.0.gbk
    DB/subset_candidates/KL601/SAMN12943639_genomic.fna.gz
    DB/subset_candidates/KL713/SAMN07152258_genomic.fna.gz
    DB/subset_candidates/KL742/SAMN41228912_genomic.fna.gz

Output:
    DB/EC-K-typing_group1and4_v1.1.gbk
    DB/subset_candidates/<KL>/<KL>_new_rep.fasta   (extracted locus FASTAs)
"""

import csv
import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pyrodigal
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

REPO_ROOT   = Path(__file__).resolve().parents[1]
GBK_IN      = REPO_ROOT / "DB" / "EC-K-typing_group1and4_v1.0.gbk"
GBK_OUT     = REPO_ROOT / "DB" / "EC-K-typing_group1and4_v1.1.gbk"
CAND_DIR    = REPO_ROOT / "DB" / "subset_candidates"
QUERY_DIR   = REPO_ROOT / "DB" / "subset_queries"

REPLACEMENTS = [
    dict(
        locus       = "KL601",
        accession   = "SAMN12943639",
        fasta_gz    = CAND_DIR / "KL601" / "SAMN12943639_genomic.fna.gz",
        ref_locus   = "KL414",       # superset — defines full locus boundary
        min_cds     = 18,            # must have at least this many (current rep = 17)
    ),
    dict(
        locus       = "KL713",
        accession   = "SAMN07152258",
        fasta_gz    = CAND_DIR / "KL713" / "SAMN07152258_genomic.fna.gz",
        ref_locus   = "KL713",       # use self — just need clean single-contig version
        min_cds     = 14,            # current rep has 15 but 1 is artefact at N-spacer
        boundary_buf = 150,          # extend merged span by 150 bp each side
    ),
    dict(
        locus       = "KL742",
        accession   = "SAMN41228912",
        fasta_gz    = CAND_DIR / "KL742" / "SAMN41228912_genomic.fna.gz",
        ref_locus   = "KL394",       # superset — defines full locus boundary
        min_cds     = 16,            # must have at least this many (current rep = 15)
    ),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_gbk_records(gbk: Path) -> dict:
    """Return {kl_name: SeqRecord} for all records in the GenBank file."""
    records = {}
    for rec in SeqIO.parse(gbk, "genbank"):
        kl = None
        for f in rec.features:
            if f.type == "source":
                for n in f.qualifiers.get("note", []):
                    if n.startswith("K type:") or n.startswith("K locus:"):
                        kl = n.split(":")[-1].strip()
        if kl is None:
            kl = rec.name.split("_")[0]
        records[kl] = rec
    return records


def write_tmp_fasta(seq: str | Seq, name: str, tmp_dir: Path) -> Path:
    path = tmp_dir / f"{name}.fasta"
    with open(path, "w") as fh:
        fh.write(f">{name}\n{seq}\n")
    return path


def decompress_fasta(fasta_gz: Path, out_dir: Path) -> Path:
    out = out_dir / fasta_gz.stem
    if not out.exists():
        with gzip.open(fasta_gz, "rt") as fin, open(out, "w") as fout:
            for line in fin:
                fout.write(line)
    return out


def make_blastdb(fasta: Path, out_dir: Path) -> Path:
    db = out_dir / "blastdb"
    if not (out_dir / "blastdb.nsi").exists() and not (out_dir / "blastdb.nsq").exists():
        subprocess.run(
            ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
             "-out", str(db), "-parse_seqids"],
            capture_output=True, check=True,
        )
    return db


def blast_locus(query_fasta: Path, db: Path, out_dir: Path,
                tag: str) -> list:
    """
    Run blastn and return list of hits sorted by qcovs desc.
    Each hit: (sseqid_bare, sstart, send, strand, pident, qcovs, qstart, qend)
    """
    out = out_dir / f"blast_{tag}.tsv"
    subprocess.run(
        ["blastn", "-query", str(query_fasta), "-db", str(db),
         "-out", str(out),
         "-outfmt", "6 sseqid sstart send sstrand pident qcovs qstart qend",
         "-perc_identity", "85", "-max_target_seqs", "5",
         "-max_hsps", "20", "-num_threads", "4"],
        capture_output=True, check=True,
    )
    hits = []
    with open(out) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row:
                continue
            sseqid, sstart, send, strand, pident, qcovs, qstart, qend = row
            bare = sseqid.split("|")[1] if "|" in sseqid else sseqid
            hits.append((bare, int(sstart), int(send), strand,
                         float(pident), float(qcovs), int(qstart), int(qend)))
    hits.sort(key=lambda x: x[5], reverse=True)
    return hits


def merged_span(hits: list) -> tuple:
    """
    For queries with N-spacers (like KL713), BLAST returns multiple HSPs on the
    same contig. Merge all HSPs on the top contig into a single span covering
    [min(sstart,send) .. max(sstart,send)] across all HSPs.
    Returns (contig, lo, hi, strand, mean_pident, total_span_bp).
    """
    if not hits:
        return None
    top_contig = hits[0][0]
    contig_hits = [h for h in hits if h[0] == top_contig]
    positions = []
    for h in contig_hits:
        positions += [h[1], h[2]]
    lo = min(positions) - 1   # 0-based
    hi = max(positions)
    strand = hits[0][3]
    mean_pident = sum(h[4] for h in contig_hits) / len(contig_hits)
    return top_contig, lo, hi, strand, mean_pident, hi - lo


def extract_region(records: dict, contig: str, lo: int, hi: int) -> Seq:
    """Extract seq[lo:hi] (0-based, half-open) from records dict."""
    rec = records.get(contig)
    if rec is None:
        raise ValueError(f"Contig {contig} not found (available: {list(records)[:3]})")
    return rec.seq[lo:hi]


def annotate_pyrodigal(seq: Seq | str, locus: str) -> list:
    """
    Return list of CDS SeqFeature objects annotated by pyrodigal metagenome mode.
    """
    orf_finder = pyrodigal.GeneFinder(meta=True)
    genes = orf_finder.find_genes(bytes(str(seq).upper(), "ascii"))

    features = []
    for i, gene in enumerate(genes, 1):
        start = min(gene.begin, gene.end) - 1   # 0-based
        end   = max(gene.begin, gene.end)
        strand = 1 if gene.strand == 1 else -1
        loc = SeqFeature.FeatureLocation(start, end, strand=strand)
        feat = SeqFeature.SeqFeature(
            loc, type="CDS",
            qualifiers={
                "locus_tag": [f"{locus}_{i:05d}"],
                "product":   ["hypothetical protein"],
            },
        )
        features.append(feat)
    return features


def build_gbk_record(locus: str, accession: str, seq: Seq,
                     cds_features: list) -> SeqRecord:
    """Construct a GenBank SeqRecord in the same style as existing DB records."""
    rec_id   = f"{locus}_{accession}"
    rec_name = rec_id[:16]     # GenBank NAME field max 16 chars

    source_feat = SeqFeature.SeqFeature(
        SeqFeature.FeatureLocation(0, len(seq), strand=1),
        type="source",
        qualifiers={
            "organism":  ["Escherichia coli"],
            "mol_type":  ["genomic DNA"],
            "note":      [f"K locus: {locus}"],
        },
    )

    rec = SeqRecord(
        seq,
        id=rec_id,
        name=rec_name,
        description=f"{locus} K-locus representative, {accession}",
        features=[source_feat] + cds_features,
        annotations={
            "molecule_type": "DNA",
            "topology":      "linear",
            "data_file_division": "BCT",
        },
    )
    return rec


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"Loading {GBK_IN.name}...")
    all_records = load_gbk_records(GBK_IN)
    print(f"  {len(all_records)} loci loaded")

    new_reps = {}   # locus → new SeqRecord

    for t in REPLACEMENTS:
        locus     = t["locus"]
        accession = t["accession"]
        fasta_gz  = t["fasta_gz"]
        ref_locus = t["ref_locus"]
        min_cds   = t["min_cds"]

        print(f"\n{'='*60}")
        print(f"{locus}  →  {accession}  (ref: {ref_locus})")
        print(f"{'='*60}")

        work_dir = CAND_DIR / locus
        work_dir.mkdir(exist_ok=True)

        # 1. Decompress and index assembly
        asm_fasta = decompress_fasta(fasta_gz, work_dir)
        asm_records = {r.id: r for r in SeqIO.parse(asm_fasta, "fasta")}
        db = make_blastdb(asm_fasta, work_dir)

        # 2. Get reference locus sequence from GenBank
        ref_rec = all_records.get(ref_locus)
        if ref_rec is None:
            print(f"  ERROR: {ref_locus} not found in GenBank — skipping")
            continue

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            ref_fasta = write_tmp_fasta(ref_rec.seq, ref_locus, tmp_dir)

            # 3. BLAST reference locus against assembly
            print(f"  BLASTing {ref_locus} ({len(ref_rec.seq)} bp) against {accession}...")
            hits = blast_locus(ref_fasta, db, work_dir, f"{ref_locus}_vs_{accession}")

        if not hits:
            print(f"  ERROR: no BLAST hit for {ref_locus} — skipping")
            continue

        best = hits[0]
        contig, sstart, send, strand, pident, qcovs, qstart, qend = best
        print(f"  Best HSP: {contig} {sstart}-{send} ({strand}) "
              f"pident={pident:.1f}% qcovs={qcovs:.1f}% "
              f"qstart={qstart} qend={qend}")

        # 4. Define locus boundaries — merge all HSPs on top contig
        # (needed for queries with N-spacers, e.g. KL713, which produce
        # multiple HSPs that together define the full locus span)
        span = merged_span(hits)
        contig, lo, hi, strand, mean_pident, span_bp = span
        print(f"  Merged span ({len([h for h in hits if h[0]==contig])} HSPs): "
              f"{lo+1}–{hi} ({span_bp} bp, mean pident={mean_pident:.1f}%)")

        # Apply per-target boundary buffer if specified
        buf = t.get("boundary_buf", 0)
        if buf:
            lo = max(0, lo - buf)
            hi = hi + buf
            print(f"  + boundary buffer ±{buf} bp → {lo+1}–{hi} ({hi-lo} bp)")

        print(f"  Locus boundaries: {lo+1}–{hi} ({hi-lo} bp) on {contig}")

        # 5. Extract sequence
        try:
            locus_seq = extract_region(asm_records, contig, lo, hi)
        except ValueError as e:
            print(f"  ERROR: {e} — skipping")
            continue

        if strand == "minus":
            locus_seq = locus_seq.reverse_complement()

        # Check for N-spacer
        n_spacer = bool(re.search(r"N{20,}", str(locus_seq)))
        print(f"  N-spacer (≥20 N): {n_spacer}")
        if n_spacer and locus == "KL713":
            print(f"  WARNING: N-spacer present in KL713 candidate — check other hits")

        # Save extracted FASTA
        out_fasta = work_dir / f"{locus}_new_rep.fasta"
        with open(out_fasta, "w") as fh:
            fh.write(f">{locus}_{accession} len={len(locus_seq)} "
                     f"blast_hit={contig}:{lo+1}-{hi}({strand}) "
                     f"pident={pident:.1f} qcovs={qcovs:.1f}\n"
                     f"{locus_seq}\n")
        print(f"  Saved: {out_fasta.name}")

        # 6. Annotate
        print(f"  Annotating with pyrodigal...")
        cds_features = annotate_pyrodigal(locus_seq, locus)
        print(f"  CDS found: {len(cds_features)}  (current rep: "
              f"{sum(1 for f in all_records[locus].features if f.type == 'CDS')}  "
              f"| minimum required: {min_cds})")

        if len(cds_features) < min_cds:
            print(f"  WARNING: CDS count {len(cds_features)} < {min_cds} — "
                  f"check boundaries")

        # 7. Build new GenBank record
        new_rec = build_gbk_record(locus, accession, locus_seq, cds_features)
        new_reps[locus] = new_rec
        print(f"  New record: {new_rec.id}, {len(locus_seq)} bp, "
              f"{len(cds_features)} CDS")

    # 8. Write updated GenBank
    if not new_reps:
        print("\nERROR: no replacements built — aborting")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Writing {GBK_OUT.name}...")
    written = 0
    with open(GBK_OUT, "w") as fh:
        for kl, rec in all_records.items():
            if kl in new_reps:
                SeqIO.write(new_reps[kl], fh, "genbank")
                print(f"  Replaced: {kl}")
            else:
                SeqIO.write(rec, fh, "genbank")
            written += 1
    print(f"  {written} records written ({len(new_reps)} replaced)")
    print(f"\nDone → {GBK_OUT}")
    print("\nNext steps:")
    print("  1. python scripts/validate_gene_content.py   "
          "# re-run validation on v1.1")
    print("  2. Confirm KL601/KL713/KL742 no longer appear in subset pairs")
    print("  3. Run type_normalized.py self-typing test on v1.1")


if __name__ == "__main__":
    main()
