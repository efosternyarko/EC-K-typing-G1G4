#!/usr/bin/env python3
"""
swap_reps_v04.py — Swap KL302 and KL305 to size-matched representatives (v0.4).

Root cause of 23 self-typing failures in v0.3:
  KL302 (40 CDS, 41 kb) and KL305 (50 CDS, 51 kb) are oversized representatives
  that accumulate more Kaptive bitscore than smaller same-KX-type loci even at
  100% identity. Fix: replace with size-matched cluster members.

Swaps
-----
  KL302: ESC_GB6443AA_AS (41,151 bp, 40 CDS) → ESC_CC1376AA_AS (KX01, 31,433 bp)
  KL305: ESC_BA8240AA_AS (50,958 bp, 50 CDS) → ESC_TA7291AA_AS (KX34, 32,468 bp)

The new sequences are re-annotated with pyrodigal + Klebsiella K-locus BLASTp,
then spliced into the v2.0 GenBank to create EC-K-typing_group1and4_v4.0_base.gbk.
Positional naming is applied by name_loci_positional.py --suffix v4 in the next step.

Usage
-----
    python scripts/swap_reps_v04.py
"""

import os
import re
import subprocess
import sys
import tempfile
from collections import OrderedDict
from datetime import date
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR       = Path(__file__).resolve().parent.parent
DB_DIR         = REPO_DIR / "DB"
INPUT_GBK      = DB_DIR / "EC-K-typing_group1and4_v2.0.gbk"
OUTPUT_GBK     = DB_DIR / "EC-K-typing_group1and4_v4.0_base.gbk"
MAPPING_FILE   = DB_DIR / "KL_G1G4_mapping.tsv"
FILT_MAPPING   = DB_DIR / "KL_G1G4_mapping_filtered.tsv"
EXTRACTED_DIR  = Path("/Users/LSHEF4/Dropbox/work_2025/e_coli_db/db_build_v2/extracted_loci")
KLEB_REF_GBK   = Path("/tmp/kleb_k_ref.gbk")

# ---------------------------------------------------------------------------
# Swaps: KL name → new assembly info
# ---------------------------------------------------------------------------
SWAPS = {
    "KL302": {
        "new_asm": "ESC_CC1376AA_AS",
        "new_kx":  "KX01",
        "fasta":   EXTRACTED_DIR / "KX01__ESC_CC1376AA_AS.fa",
    },
    "KL305": {
        "new_asm": "ESC_TA7291AA_AS",
        "new_kx":  "KX34",
        "fasta":   EXTRACTED_DIR / "KX34__ESC_TA7291AA_AS.fa",
    },
}


# ---------------------------------------------------------------------------
# Annotation helpers  (same logic as annotate_loci.py)
# ---------------------------------------------------------------------------

def build_kleb_blastdb():
    """Extract proteins from KLEB_REF_GBK and build a BLAST protein database.

    Returns (db_prefix, uid_to_gene_product).
    """
    import pyrodigal  # noqa: F401  (just to check it is available)

    gene_info = {}
    proteins  = []
    idx       = 0
    for rec in SeqIO.parse(str(KLEB_REF_GBK), "genbank"):
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers:
                continue
            gene    = f.qualifiers["gene"][0]
            product = f.qualifiers.get("product", ["hypothetical protein"])[0]
            nuc     = f.location.extract(rec.seq)
            if len(nuc) % 3 != 0:
                continue
            try:
                prot = str(nuc.translate(to_stop=True))
            except Exception:
                continue
            if len(prot) < 20:
                continue
            if gene not in gene_info:
                gene_info[gene] = product
            proteins.append((gene, prot, f"ref_{idx}"))
            idx += 1

    faa = tempfile.NamedTemporaryFile(
        mode="w", suffix=".faa", delete=False, prefix="kleb_ref_"
    )
    for gene, prot, uid in proteins:
        faa.write(f">{uid} {gene}\n{prot}\n")
    faa.close()

    subprocess.run(
        ["makeblastdb", "-in", faa.name, "-dbtype", "prot"],
        capture_output=True, check=True,
    )

    uid_to_gene = {uid: (gene, gene_info[gene]) for gene, _, uid in proteins}
    return faa.name, uid_to_gene


def predict_and_annotate(sequence, ref_db, uid_to_gene):
    """Run pyrodigal + BLASTp; return list of (start, end, strand, gene_name, product)."""
    import pyrodigal

    orf_finder = pyrodigal.GeneFinder(meta=True)
    genes      = orf_finder.find_genes(bytes(sequence))

    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as pf:
        prot_file = pf.name
        for i, gene in enumerate(genes):
            prot = gene.translate()
            if prot and prot[-1] == "*":
                prot = prot[:-1]
            pf.write(f">cds_{i}\n{prot}\n")

    annotations = {}
    try:
        result = subprocess.run(
            [
                "blastp", "-query", prot_file, "-db", ref_db,
                "-outfmt", "6 qseqid sseqid pident qcovs evalue bitscore",
                "-evalue", "1e-5", "-max_target_seqs", "5",
            ],
            capture_output=True, text=True,
        )
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            parts     = line.split("\t")
            qid       = parts[0]
            sid       = parts[1]
            pident    = float(parts[2])
            qcovs     = float(parts[3])
            bitscore  = float(parts[5])
            gene_idx  = int(qid.split("_")[1])
            ref_gene, ref_product = uid_to_gene.get(sid, (None, None))
            if ref_gene is None:
                continue
            if pident >= 30 and qcovs >= 50:
                if gene_idx not in annotations or bitscore > annotations[gene_idx][3]:
                    annotations[gene_idx] = (ref_gene, ref_product, pident, bitscore)
    finally:
        os.unlink(prot_file)

    result_genes = []
    for i, gene in enumerate(genes):
        start   = gene.begin - 1
        end     = gene.end
        strand  = gene.strand
        cds_len = end - start
        if cds_len % 3 != 0:
            excess = cds_len % 3
            if strand == 1:
                end   -= excess
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


def create_genbank_record(locus_name, sequence, source_assembly, gene_list):
    """Build a Kaptive-compatible GenBank record (same format as annotate_loci.py)."""
    record = SeqRecord(
        Seq(str(sequence)),
        id          = f"{locus_name}_{source_assembly}",
        name        = f"{locus_name}_{source_assembly}",
        description = (
            f"Escherichia coli capsular polysaccharide synthesis "
            f"gene cluster, {locus_name}"
        ),
    )
    record.annotations["molecule_type"]       = "DNA"
    record.annotations["topology"]            = "linear"
    record.annotations["data_file_division"]  = "BCT"
    record.annotations["date"]                = date.today().strftime("%d-%b-%Y").upper()
    record.annotations["organism"]            = "Escherichia coli"
    record.annotations["taxonomy"]            = [
        "Bacteria", "Pseudomonadota", "Gammaproteobacteria",
        "Enterobacterales", "Enterobacteriaceae", "Escherichia",
    ]
    record.annotations["source"] = "Escherichia coli"

    src_q = OrderedDict()
    src_q["organism"] = ["Escherichia coli"]
    src_q["mol_type"] = ["genomic DNA"]
    src_q["note"]     = [f"K locus: {locus_name}"]
    record.features.append(
        SeqFeature(FeatureLocation(0, len(sequence)), type="source", qualifiers=src_q)
    )

    gene_name_counts = {}
    for i, (start, end, strand, gene_name, product) in enumerate(gene_list):
        gene_num  = str(i + 1).zfill(5)
        locus_tag = f"{locus_name}_{gene_num}"

        display_gene = gene_name
        if gene_name:
            if gene_name in gene_name_counts:
                gene_name_counts[gene_name] += 1
                display_gene = f"{gene_name}_{gene_name_counts[gene_name]}"
            else:
                gene_name_counts[gene_name] = 1

        q = OrderedDict()
        q["locus_tag"] = [locus_tag]
        q["gene"]      = [display_gene if display_gene else locus_tag]
        q["product"]   = [product]

        record.features.append(
            SeqFeature(FeatureLocation(start, end, strand=strand), type="CDS", qualifiers=q)
        )

    return record


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Swapping KL302 and KL305 representatives  —  v0.4")
    print("=" * 70)

    if not KLEB_REF_GBK.exists():
        print(f"ERROR: Klebsiella reference not found: {KLEB_REF_GBK}", file=sys.stderr)
        sys.exit(1)

    # Build Klebsiella BLAST database
    print("\n[1] Building Klebsiella reference BLAST database...")
    ref_db, uid_to_gene = build_kleb_blastdb()
    print(f"    {len(uid_to_gene)} reference proteins")

    # Annotate new representatives
    new_records = {}
    for kl_name, swap in SWAPS.items():
        asm   = swap["new_asm"]
        fasta = swap["fasta"]
        print(f"\n[2] Annotating {kl_name} → {asm}  ({fasta.name})")
        if not fasta.exists():
            print(f"    ERROR: {fasta} not found", file=sys.stderr)
            sys.exit(1)
        seq_rec   = next(SeqIO.parse(str(fasta), "fasta"))
        gene_list = predict_and_annotate(seq_rec.seq, ref_db, uid_to_gene)
        n_cds     = len(gene_list)
        n_ann     = sum(1 for g in gene_list if g[3] is not None)
        print(f"    {n_cds} CDS, {n_ann} annotated ({100*n_ann/n_cds:.0f}%)")
        gbk_rec = create_genbank_record(kl_name, seq_rec.seq, asm, gene_list)
        new_records[kl_name] = gbk_rec

    # Clean up BLAST database
    os.unlink(ref_db)
    for ext in [".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto"]:
        try:
            os.unlink(ref_db + ext)
        except FileNotFoundError:
            pass

    # Splice into v2.0 GenBank
    print(f"\n[3] Loading {INPUT_GBK.name} and splicing in new records...")
    records_out = []
    n_replaced  = 0
    for rec in SeqIO.parse(str(INPUT_GBK), "genbank"):
        kl = rec.name.split("_")[0]
        if kl in new_records:
            records_out.append(new_records[kl])
            old_cds = sum(1 for f in rec.features if f.type == "CDS")
            new_cds = sum(1 for f in new_records[kl].features if f.type == "CDS")
            print(f"    Replaced {kl}: {len(rec.seq)} bp/{old_cds} CDS  →  "
                  f"{len(new_records[kl].seq)} bp/{new_cds} CDS")
            n_replaced += 1
        else:
            records_out.append(rec)

    print(f"    {n_replaced} loci replaced, {len(records_out)} total records")

    # Write output
    print(f"\n[4] Writing {OUTPUT_GBK.name}...")
    with open(OUTPUT_GBK, "w") as fh:
        SeqIO.write(records_out, fh, "genbank")
    print(f"    Written: {OUTPUT_GBK}")

    # Update mapping files
    print("\n[5] Updating mapping files...")
    for mfile in [MAPPING_FILE, FILT_MAPPING]:
        df = pd.read_csv(mfile, sep="\t")
        for kl_name, swap in SWAPS.items():
            mask = df["KL"] == kl_name
            if mask.any():
                old_asm = df.loc[mask, "source_assembly"].iloc[0]
                old_len = df.loc[mask, "length_bp"].iloc[0]
                new_len = len(next(SeqIO.parse(str(swap["fasta"]), "fasta")).seq)
                df.loc[mask, "source_assembly"] = swap["new_asm"]
                df.loc[mask, "length_bp"]       = new_len
                print(f"    {mfile.name}  {kl_name}: {old_asm} ({old_len} bp) "
                      f"→ {swap['new_asm']} ({new_len} bp)")
        df.to_csv(mfile, sep="\t", index=False)

    print(f"\n{'=' * 70}")
    print("SWAP COMPLETE")
    print(f"{'=' * 70}")
    for rec in records_out:
        kl = rec.name.split("_")[0]
        if kl in SWAPS:
            n_cds = sum(1 for f in rec.features if f.type == "CDS")
            print(f"  {kl}: {len(rec.seq)} bp, {n_cds} CDS  [{rec.name}]")
    print(f"\nNext step:")
    print(f"  python scripts/name_loci_positional.py --suffix v4")


if __name__ == "__main__":
    main()
