#!/usr/bin/env python3
"""
replace_kl300_kl303.py

Given validated assembly FASTAs for KL300 and KL303, this script:
  1. Aligns the existing KL300/KL303 locus sequences (from v0.8 GenBank) to
     each assembly using minimap2 to locate the locus region
  2. Extracts the locus region (plus 500 bp flanks, trimmed to contig bounds)
  3. Annotates with pyrodigal (metagenomic) + BLASTp vs Klebsiella reference
  4. Replaces the KL300 and KL303 records in v0.8 GenBank → v0.9

Usage:
    python3 scripts/replace_kl300_kl303.py \
        --kl300_fasta DB/kl300_303_fastas/SAMD00053151.fa \
        --kl303_fasta DB/kl300_303_fastas/SAMEA6656333.fa

    # Or let the script auto-select from DB/kl300_303_fastas/ based on
    # the validation summary (picks first correctly self-typed assembly):
    python3 scripts/replace_kl300_kl303.py --auto
"""

import argparse
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR   = Path(__file__).resolve().parent.parent
DB_DIR     = REPO_DIR / "DB"
INPUT_GBK  = DB_DIR / "EC-K-typing_group1and4_v0.8.gbk"
OUTPUT_GBK = DB_DIR / "EC-K-typing_group1and4_v0.9.gbk"
FASTA_DIR  = DB_DIR / "kl300_303_fastas"
VALID_TSV  = DB_DIR / "kl300_303_validation_summary.tsv"

KLEB_REF_GBK = Path(
    "/Users/lshef4/Library/Python/3.9/lib/python/site-packages"
    "/kaptive/data/Klebsiella_k_locus_primary_reference.gbk"
)
BLASTP      = "/Users/lshef4/lshef4_to_copy/anaconda3/envs/my_python.env/bin/blastp"
MAKEBLASTDB = "/Users/lshef4/lshef4_to_copy/anaconda3/envs/my_python.env/bin/makeblastdb"
MINIMAP2 = "/Users/lshef4/lshef4_to_copy/anaconda3/envs/kleborate_test/bin/minimap2"
FLANK    = 500   # bp flanking each side of the locus region

MIN_GENE_ID  = 30.0
MIN_GENE_COV = 50.0


# ---------------------------------------------------------------------------
# Helpers shared with annotate_and_merge_atb_loci.py
# ---------------------------------------------------------------------------

def get_locus_name(rec):
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def extract_klebsiella_proteins(kleb_gbk: Path, out_faa: Path):
    with open(out_faa, "w") as fh:
        for rec in SeqIO.parse(str(kleb_gbk), "genbank"):
            ln = get_locus_name(rec)
            for f in rec.features:
                if f.type != "CDS":
                    continue
                gene = f.qualifiers.get("gene", ["unknown"])[0]
                aa = f.qualifiers.get("translation", [None])[0]
                if aa:
                    fh.write(f">{ln}__{gene}\n{aa}\n")


def make_blastdb(faa: Path, db_path: Path):
    subprocess.run(
        [MAKEBLASTDB, "-in", str(faa), "-dbtype", "prot", "-out", str(db_path)],
        check=True, capture_output=True
    )


def blastp_annotate(query_faa: Path, db_path: Path):
    """Returns dict: seq_id -> gene_name (best hit >= thresholds)."""
    result = subprocess.run(
        [BLASTP, "-query", str(query_faa), "-db", str(db_path),
         "-outfmt", "6 qseqid sseqid pident qcovs",
         "-max_target_seqs", "1", "-num_threads", "4"],
        capture_output=True, text=True, check=True
    )
    hits = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        qid, sid, pident, qcovs = parts[0], parts[1], float(parts[2]), float(parts[3])
        if pident >= MIN_GENE_ID and qcovs >= MIN_GENE_COV:
            gene = sid.split("__")[-1]
            if qid not in hits:
                hits[qid] = gene
    return hits


def predict_and_annotate(locus_seq: str, locus_name: str, blast_db: Path):
    """Predict CDS, BLASTp, return list of SeqFeature."""
    orf_finder = pyrodigal.GeneFinder(meta=True)
    genes = orf_finder.find_genes(locus_seq.encode())

    with tempfile.NamedTemporaryFile(suffix=".faa", mode="w", delete=False) as fh:
        faa_path = Path(fh.name)
        for i, gene in enumerate(genes):
            fh.write(f">gene_{i+1}\n{gene.translate()}\n")

    gene_names = blastp_annotate(faa_path, blast_db)
    faa_path.unlink()

    features = []
    for i, gene in enumerate(genes):
        gid = f"gene_{i+1}"
        gname = gene_names.get(gid, "")
        loc = FeatureLocation(gene.begin - 1, gene.end,
                              strand=gene.strand)
        quals = {
            "product": [gname if gname else "hypothetical protein"],
            "translation": [gene.translate()],
        }
        if gname:
            quals["gene"] = [gname]
        features.append(SeqFeature(loc, type="CDS", qualifiers=quals))

    return features


def build_gbk_record(locus_name: str, seq: str, accession: str,
                     features: list) -> SeqRecord:
    source = SeqFeature(
        FeatureLocation(0, len(seq), strand=1),
        type="source",
        qualifiers={"note": [f"K locus: {locus_name}"],
                    "organism": ["Escherichia coli"],
                    "mol_type": ["genomic DNA"]}
    )
    rec = SeqRecord(
        Seq(seq),
        id=f"{locus_name}_ATB",
        name=f"{locus_name}_ATB",
        description=f"E. coli {locus_name} K-locus, extracted from ATB assembly {accession}",
        features=[source] + features,
        annotations={"molecule_type": "DNA"}
    )
    return rec


# ---------------------------------------------------------------------------
# Locus extraction via minimap2
# ---------------------------------------------------------------------------

def extract_locus_from_assembly(locus_seq: str, assembly_fasta: Path,
                                locus_name: str) -> Optional[str]:
    """
    Align locus_seq to assembly using minimap2 (asm5 preset).
    Find the best-covered contig region, extract it with FLANK bp on each side.
    Returns the extracted sequence string, or None if alignment fails.
    """
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as fh:
        fh.write(f">{locus_name}\n{locus_seq}\n")
        query_fa = Path(fh.name)

    try:
        result = subprocess.run(
            [MINIMAP2, "-c", "--cs", "-x", "asm5",
             str(assembly_fasta), str(query_fa)],
            capture_output=True, text=True, check=True
        )
    finally:
        query_fa.unlink()

    if not result.stdout.strip():
        print(f"    WARNING: no minimap2 alignments for {locus_name}")
        return None

    # Parse PAF: find the alignment with the highest number of matching bases
    best = None
    for line in result.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        t_name  = cols[5]
        t_len   = int(cols[6])
        t_start = int(cols[7])
        t_end   = int(cols[8])
        n_match = int(cols[9])
        if best is None or n_match > best[4]:
            best = (t_name, t_len, t_start, t_end, n_match)

    if best is None:
        return None

    t_name, t_len, t_start, t_end, _ = best
    extract_start = max(0, t_start - FLANK)
    extract_end   = min(t_len, t_end + FLANK)

    # Load the assembly and extract the region
    assembly = SeqIO.to_dict(SeqIO.parse(str(assembly_fasta), "fasta"))
    if t_name not in assembly:
        print(f"    WARNING: contig {t_name} not in assembly dict")
        return None

    region = str(assembly[t_name].seq[extract_start:extract_end])
    print(f"    Extracted {len(region):,} bp from {t_name}:{extract_start}-{extract_end}")
    return region


# ---------------------------------------------------------------------------
# Auto-select: pick first correctly self-typed assembly per locus
# ---------------------------------------------------------------------------

def auto_select_fastas():
    import csv
    selected = {}
    if not VALID_TSV.exists():
        raise FileNotFoundError(f"Validation summary not found: {VALID_TSV}\n"
                                "Run validate_kl300_kl303_m3.sh first and scp the results.")
    with open(VALID_TSV) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            kl = row.get("Expected_KL", "")
            correct = row.get("Correct", "")
            asm = Path(row.get("Assembly", "")).stem
            if correct == "True" and kl in ("KL300", "KL303") and kl not in selected:
                fa = FASTA_DIR / f"{asm}.fa"
                if fa.exists():
                    selected[kl] = (fa, asm)
                    print(f"    Auto-selected for {kl}: {asm}")
    return selected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kl300_fasta", type=Path)
    parser.add_argument("--kl303_fasta", type=Path)
    parser.add_argument("--auto", action="store_true",
                        help="Auto-select from validation summary")
    args = parser.parse_args()

    if args.auto:
        selected = auto_select_fastas()
        kl300_fa, kl300_acc = selected.get("KL300", (None, None))
        kl303_fa, kl303_acc = selected.get("KL303", (None, None))
    else:
        kl300_fa  = args.kl300_fasta
        kl303_fa  = args.kl303_fasta
        kl300_acc = kl300_fa.stem if kl300_fa else None
        kl303_acc = kl303_fa.stem if kl303_fa else None

    to_replace = {}
    if kl300_fa:
        to_replace["KL300"] = (kl300_fa, kl300_acc)
    if kl303_fa:
        to_replace["KL303"] = (kl303_fa, kl303_acc)

    if not to_replace:
        raise ValueError("Need at least one FASTA (--kl300_fasta, --kl303_fasta, or --auto).")

    for kl, (fa, acc) in to_replace.items():
        print(f"\n{kl} assembly: {acc}  ({fa})")

    # ── Build BLASTp DB from Klebsiella reference ─────────────────────────────
    print("\n[1] Building BLASTp reference DB...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        kleb_faa  = tmpdir / "kleb_ref.faa"
        blast_db  = tmpdir / "kleb_ref"
        extract_klebsiella_proteins(KLEB_REF_GBK, kleb_faa)
        make_blastdb(kleb_faa, blast_db)
        print(f"    BLASTp DB built.")

        # ── Load existing KL300 and KL303 locus sequences from v0.8 ──────────
        print("\n[2] Loading existing KL300/KL303 sequences from v0.8 GenBank...")
        existing_loci = {}
        for rec in SeqIO.parse(str(INPUT_GBK), "genbank"):
            ln = get_locus_name(rec)
            if ln in ("KL300", "KL303"):
                existing_loci[ln] = str(rec.seq)
                print(f"    {ln}: {len(rec.seq):,} bp")

        # ── Extract locus regions from validated assemblies ───────────────────
        print("\n[3] Extracting locus regions from assemblies via minimap2...")
        new_seqs = {}
        for kl, (fa, acc) in to_replace.items():
            seq = extract_locus_from_assembly(existing_loci[kl], fa, kl)
            if seq is None:
                raise RuntimeError(f"Failed to extract {kl} locus from {acc}")
            new_seqs[kl] = (seq, acc)

        # ── Annotate extracted loci ───────────────────────────────────────────
        print("\n[4] Annotating extracted loci...")
        new_records = {}
        for kl, (seq, acc) in new_seqs.items():
            print(f"    {kl} ({len(seq):,} bp)...")
            features = predict_and_annotate(seq, kl, blast_db)
            n_ann = sum(1 for f in features if f.qualifiers.get("gene"))
            print(f"      {len(features)} CDS predicted, {n_ann} annotated")
            new_records[kl] = build_gbk_record(kl, seq, acc, features)

    # ── Build v0.9 GenBank: replace KL300/KL303, keep everything else ─────────
    print(f"\n[5] Writing v0.9 GenBank: {OUTPUT_GBK}")
    written = 0
    replaced = set()
    with open(OUTPUT_GBK, "w") as out:
        for rec in SeqIO.parse(str(INPUT_GBK), "genbank"):
            ln = get_locus_name(rec)
            if ln in new_records:
                SeqIO.write(new_records[ln], out, "genbank")
                replaced.add(ln)
                print(f"    Replaced {ln}")
            else:
                SeqIO.write(rec, out, "genbank")
            written += 1

    print(f"\n{'='*55}")
    print(f"v0.9 GenBank written: {OUTPUT_GBK}")
    print(f"  Total records: {written}")
    print(f"  Replaced:      {', '.join(sorted(replaced))}")
    print(f"\nNext: run self-typing validation on v0.9")
    print(f"  python3 scripts/type_normalized.py")


if __name__ == "__main__":
    main()
