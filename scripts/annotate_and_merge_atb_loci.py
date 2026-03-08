#!/usr/bin/env python3
"""
annotate_and_merge_atb_loci.py

Annotate novel G1/G4 K-loci from AllTheBacteria clustering and merge
with the existing v0.6 reference database to produce v0.8.

Steps:
  1. Screen novel loci against existing DB: discard any novel locus where
     >= 80% of the novel sequence aligns at >= 95% identity to an existing locus
     (bidirectional check prevents partial captures of existing loci slipping through)
  2. Extract Klebsiella K-locus reference proteins (for BLASTp annotation)
  3. For each novel locus in atb_clustering/novel_kl_types.fasta:
       - Predict CDS with pyrodigal (metagenomic mode)
       - BLASTp predicted proteins vs Klebsiella reference
       - Transfer gene names (>=30% identity, >=50% coverage)
       - Build Kaptive-compatible GenBank record
  4. Merge annotated novel loci with existing v0.6 GenBank -> v0.8
  5. Write annotation summary TSV

Usage:
    python3 scripts/annotate_and_merge_atb_loci.py

Requirements:
    pip install pyrodigal biopython
    BLAST+ (blastp, makeblastdb)
"""

import os
import re
import subprocess
import tempfile
from collections import OrderedDict
from datetime import date
from pathlib import Path

import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR    = Path(__file__).resolve().parent.parent
DB_DIR      = REPO_DIR / "DB"
ATB_DIR     = REPO_DIR.parent / "atb_clustering"

NOVEL_FASTA = ATB_DIR / "novel_kl_types.fasta"
EXISTING_GBK = DB_DIR / "EC-K-typing_group1and4_v0.6.gbk"
OUTPUT_GBK  = DB_DIR / "EC-K-typing_group1and4_v0.8.gbk"
SUMMARY_TSV = DB_DIR / "annotation_summary_atb_v0.8.tsv"

KLEB_REF_GBK = Path(
    "/Users/lshef4/Library/Python/3.9/lib/python/site-packages"
    "/kaptive/data/Klebsiella_k_locus_primary_reference.gbk"
)

BLASTP      = "/Users/lshef4/lshef4_to_copy/anaconda3/envs/my_python.env/bin/blastp"
MAKEBLASTDB = "/Users/lshef4/lshef4_to_copy/anaconda3/envs/my_python.env/bin/makeblastdb"

# BLASTp thresholds
MIN_PIDENT = 30.0
MIN_QCOV   = 50.0

# ---------------------------------------------------------------------------


def build_kleb_blast_db():
    """
    Extract CDS protein sequences from Klebsiella K-locus reference GenBank,
    write to a temp FASTA, and build a BLASTp database.

    Returns:
        (db_path, uid_to_gene) where uid_to_gene maps 'ref_N' -> (gene_name, product)
    """
    print(f"[1] Building Klebsiella reference BLAST database from {KLEB_REF_GBK.name}...")
    gene_info = {}    # gene_name -> product
    proteins  = []   # (gene_name, prot_seq, uid)

    idx = 0
    for rec in SeqIO.parse(str(KLEB_REF_GBK), "genbank"):
        for f in rec.features:
            if f.type != "CDS":
                continue
            if "gene" not in f.qualifiers:
                continue
            gene    = f.qualifiers["gene"][0]
            product = f.qualifiers.get("product", ["hypothetical protein"])[0]
            nuc_seq = f.location.extract(rec.seq)
            if len(nuc_seq) % 3 != 0:
                continue
            try:
                prot = str(nuc_seq.translate(to_stop=True))
            except Exception:
                continue
            if len(prot) < 20:
                continue
            if gene not in gene_info:
                gene_info[gene] = product
            proteins.append((gene, prot, f"ref_{idx}"))
            idx += 1

    # Write protein FASTA
    ref_fa = tempfile.NamedTemporaryFile(
        mode="w", suffix=".faa", delete=False, prefix="kleb_ref_"
    )
    for gene, prot, uid in proteins:
        ref_fa.write(f">{uid} {gene}\n{prot}\n")
    ref_fa.close()

    # Build BLAST db
    subprocess.run(
        [MAKEBLASTDB, "-in", ref_fa.name, "-dbtype", "prot"],
        capture_output=True, check=True
    )

    uid_to_gene = {uid: (gene, gene_info[gene]) for gene, _, uid in proteins}
    print(f"    {len(proteins)} reference proteins indexed ({len(gene_info)} unique gene names)")
    return ref_fa.name, uid_to_gene


def predict_and_annotate(sequence, ref_db, uid_to_gene):
    """
    Predict CDS with pyrodigal (metagenomic mode) and annotate by BLASTp
    against the Klebsiella reference.

    Returns list of (start0, end, strand, gene_name_or_None, product)
    """
    orf_finder = pyrodigal.GeneFinder(meta=True)
    genes = orf_finder.find_genes(bytes(sequence))

    # Write predicted proteins
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as pf:
        prot_file = pf.name
        for i, gene in enumerate(genes):
            prot = gene.translate()
            if prot and prot[-1] == "*":
                prot = prot[:-1]
            pf.write(f">cds_{i}\n{prot}\n")

    annotations = {}  # gene_idx -> (gene_name, product, pident, bitscore)

    try:
        result = subprocess.run(
            [BLASTP, "-query", prot_file, "-db", ref_db,
             "-outfmt", "6 qseqid sseqid pident qcovs evalue bitscore",
             "-evalue", "1e-5", "-max_target_seqs", "5",
             "-num_threads", "4"],
            capture_output=True, text=True
        )
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            qid     = parts[0]
            sid     = parts[1]
            pident  = float(parts[2])
            qcovs   = float(parts[3])
            bitscore = float(parts[5])
            gene_idx = int(qid.split("_")[1])
            ref_gene, ref_product = uid_to_gene.get(sid, (None, None))
            if ref_gene is None:
                continue
            if pident >= MIN_PIDENT and qcovs >= MIN_QCOV:
                if gene_idx not in annotations or bitscore > annotations[gene_idx][3]:
                    annotations[gene_idx] = (ref_gene, ref_product, pident, bitscore)
    finally:
        os.unlink(prot_file)

    result_genes = []
    for i, gene in enumerate(genes):
        start = gene.begin - 1   # 0-based
        end   = gene.end
        strand = gene.strand

        # Trim CDS to multiple of 3
        cds_len = end - start
        if cds_len % 3 != 0:
            excess = cds_len % 3
            if strand == 1:
                end -= excess
            else:
                start += excess

        if end - start < 60:
            continue

        if i in annotations:
            gene_name, product, _, _ = annotations[i]
        else:
            gene_name = None
            product   = "hypothetical protein"

        result_genes.append((start, end, strand, gene_name, product))

    return result_genes


def make_genbank_record(locus_name, source_assembly, sequence, gene_list):
    """Create a Kaptive-compatible GenBank record."""
    record = SeqRecord(
        Seq(str(sequence)),
        id=f"{locus_name}_{source_assembly}",
        name=f"{locus_name}_{source_assembly}"[:16],
        description=(
            f"Escherichia coli capsular polysaccharide synthesis "
            f"gene cluster, {locus_name}"
        ),
    )
    record.annotations.update({
        "molecule_type": "DNA",
        "topology": "linear",
        "data_file_division": "BCT",
        "date": date.today().strftime("%d-%b-%Y").upper(),
        "organism": "Escherichia coli",
        "taxonomy": [
            "Bacteria", "Pseudomonadota", "Gammaproteobacteria",
            "Enterobacterales", "Enterobacteriaceae", "Escherichia"
        ],
        "source": "Escherichia coli",
    })

    # Source feature
    src_quals = OrderedDict([
        ("organism", ["Escherichia coli"]),
        ("mol_type", ["genomic DNA"]),
        ("note",     [f"K locus: {locus_name}"]),
    ])
    record.features.append(
        SeqFeature(FeatureLocation(0, len(sequence)), type="source",
                   qualifiers=src_quals)
    )

    # CDS features
    gene_name_counts = {}
    for i, (start, end, strand, gene_name, product) in enumerate(gene_list):
        locus_tag    = f"{locus_name}_{str(i+1).zfill(5)}"
        display_gene = gene_name
        if gene_name:
            if gene_name in gene_name_counts:
                gene_name_counts[gene_name] += 1
                display_gene = f"{gene_name}_{gene_name_counts[gene_name]}"
            else:
                gene_name_counts[gene_name] = 1

        quals = OrderedDict([("locus_tag", [locus_tag])])
        if display_gene:
            quals["gene"]    = [display_gene]
        quals["product"] = [product]

        record.features.append(
            SeqFeature(FeatureLocation(start, end, strand=strand),
                       type="CDS", qualifiers=quals)
        )

    return record


def screen_against_existing(novel_records, existing_gbk):
    """
    Remove novel loci that are partial captures of existing loci.

    A novel locus is rejected if >= 80% of its length aligns to any existing
    locus at >= 95% nucleotide identity (one-directional: novel as query).
    This catches truncated extractions that slipped through bidirectional
    clustering (e.g. KL486/KL562/KL812 in v0.8 build).

    Returns: (kept_records, rejected_names)
    """
    print("\n[Screening] Checking novel loci against existing DB for partial captures...")

    # Write existing loci to temp FASTA
    existing_fa = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, prefix="existing_loci_"
    )
    for rec in SeqIO.parse(str(existing_gbk), "genbank"):
        existing_fa.write(f">{rec.id}\n{str(rec.seq)}\n")
    existing_fa.close()

    # Write novel loci to temp FASTA
    novel_fa = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, prefix="novel_loci_"
    )
    for rec in novel_records:
        novel_fa.write(f">{rec.id}\n{str(rec.seq)}\n")
    novel_fa.close()

    # Build BLAST db from existing loci
    db_path = existing_fa.name + "_db"
    subprocess.run(
        [MAKEBLASTDB, "-in", existing_fa.name, "-dbtype", "nucl", "-out", db_path],
        capture_output=True, check=True
    )

    # BLAST novel (query) vs existing (subject): qcovs = coverage of the novel sequence
    result = subprocess.run(
        ["/Users/lshef4/.local/bin/blastn",
         "-query", novel_fa.name, "-db", db_path,
         "-outfmt", "6 qseqid sseqid pident qcovs",
         "-perc_identity", "95", "-qcov_hsp_perc", "80",
         "-max_target_seqs", "1"],
        capture_output=True, text=True
    )

    # Collect rejected loci (any hit = partial capture)
    rejected = set()
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        novel_id, existing_id, pident, qcov = parts[0], parts[1], float(parts[2]), float(parts[3])
        rejected.add(novel_id)
        print(f"    REJECT {novel_id} → {existing_id} ({pident:.1f}% id, {qcov:.0f}% qcov)")

    # Cleanup
    os.unlink(existing_fa.name)
    os.unlink(novel_fa.name)
    for ext in [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto"]:
        try:
            os.unlink(db_path + ext)
        except FileNotFoundError:
            pass

    kept = [r for r in novel_records if r.id not in rejected]
    print(f"    Kept {len(kept)}/{len(novel_records)} novel loci ({len(rejected)} rejected as partial captures)")
    return kept, rejected


def cleanup_blast_db(db_path):
    os.unlink(db_path)
    for ext in [".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto"]:
        try:
            os.unlink(db_path + ext)
        except FileNotFoundError:
            pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    # Step 1: Screen novel loci against existing DB
    novel_records_raw = list(SeqIO.parse(str(NOVEL_FASTA), "fasta"))
    novel_records_input, _ = screen_against_existing(novel_records_raw, EXISTING_GBK)

    # Step 2: Build Klebsiella reference BLAST DB
    ref_db, uid_to_gene = build_kleb_blast_db()

    # Step 3: Annotate novel ATB loci
    print(f"\n[3] Annotating novel loci from {NOVEL_FASTA.name}...")

    print(f"    {len(novel_records_input)} novel loci to annotate (after screening)")

    novel_gbk_records = []
    summary_rows = []

    for i, rec in enumerate(novel_records_input, 1):
        # Parse header: ">KL392 SAMN30284532 19219bp cluster_size=86028"
        parts = rec.description.split()
        locus_name      = parts[0]
        source_assembly = parts[1] if len(parts) > 1 else "unknown"

        gene_list = predict_and_annotate(rec.seq, ref_db, uid_to_gene)

        n_cds = len(gene_list)
        n_ann = sum(1 for g in gene_list if g[3] is not None)
        pct   = (100 * n_ann / n_cds) if n_cds else 0

        if i % 50 == 0 or i <= 5:
            print(f"    [{i:4d}/{len(novel_records_input)}] {locus_name} "
                  f"({len(rec.seq):,} bp): {n_cds} CDS, {n_ann} annotated ({pct:.0f}%)")

        gbk_rec = make_genbank_record(locus_name, source_assembly, rec.seq, gene_list)
        novel_gbk_records.append(gbk_rec)

        summary_rows.append({
            "KL_type":        locus_name,
            "source_assembly": source_assembly,
            "seq_len":        len(rec.seq),
            "n_cds":          n_cds,
            "n_annotated":    n_ann,
            "pct_annotated":  f"{pct:.1f}",
        })

    cleanup_blast_db(ref_db)
    print(f"\n    Annotated {len(novel_gbk_records)} novel loci")

    # Step 4: Load existing v0.6 and merge
    print(f"\n[4] Loading existing database {EXISTING_GBK.name}...")
    existing_records = list(SeqIO.parse(str(EXISTING_GBK), "genbank"))
    print(f"    {len(existing_records)} existing loci (v0.6)")

    all_records = existing_records + novel_gbk_records
    print(f"    Total after merge: {len(all_records)} loci")

    # Write v0.8 GenBank
    print(f"\n[5] Writing {OUTPUT_GBK.name}...")
    SeqIO.write(all_records, str(OUTPUT_GBK), "genbank")
    print(f"    Done.")

    # Write annotation summary
    with open(str(SUMMARY_TSV), "w") as f:
        f.write("KL_type\tsource_assembly\tseq_len\tn_cds\tn_annotated\tpct_annotated\n")
        for row in summary_rows:
            f.write(
                f"{row['KL_type']}\t{row['source_assembly']}\t{row['seq_len']}\t"
                f"{row['n_cds']}\t{row['n_annotated']}\t{row['pct_annotated']}\n"
            )
    print(f"    Annotation summary: {SUMMARY_TSV.name}")

    # Step 5: Basic validation
    print(f"\n[6] Validation...")
    parsed = list(SeqIO.parse(str(OUTPUT_GBK), "genbank"))
    locus_re = re.compile(r"K locus:\s*(\S+)")
    fails = []

    for rec in parsed:
        src = [f for f in rec.features if f.type == "source"]
        has_locus = any(
            locus_re.search(note)
            for f in src
            for note in f.qualifiers.get("note", [])
        )
        n_cds = sum(1 for f in rec.features if f.type == "CDS")
        if not has_locus or n_cds == 0:
            fails.append(rec.id)

    total_cds = sum(
        sum(1 for f in r.features if f.type == "CDS") for r in parsed
    )
    total_ann = sum(
        sum(1 for f in r.features if f.type == "CDS" and "gene" in f.qualifiers)
        for r in parsed
    )

    print(f"    Records parsed:   {len(parsed)}")
    print(f"    Total CDS:        {total_cds}")
    print(f"    Annotated CDS:    {total_ann} ({100*total_ann/total_cds:.1f}%)")
    if fails:
        print(f"    FAILURES ({len(fails)}): {fails[:10]}")
    else:
        print(f"    All records passed validation.")

    print(f"\nDone. v0.8 database: {OUTPUT_GBK}")
    print(f"Next: run self-typing validation with Kaptive normalised scoring.")


if __name__ == "__main__":
    main()
