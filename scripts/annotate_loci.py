#!/usr/bin/env python3
"""
Annotate Group 1/4 K-locus reference sequences using Klebsiella K-locus gene models.

Uses pyrodigal for gene prediction, then BLASTs predicted CDS against the
Klebsiella K-locus reference database CDS to transfer gene names.
This provides much richer annotation than flanking-genes-only, enabling
Kaptive to discriminate between loci effectively.
"""

import os
import sys
import subprocess
import tempfile
from datetime import date
from collections import OrderedDict

import pyrodigal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Paths
DB_DIR = "/Users/LSHEF4/Dropbox/work_2025/e_coli_db/db_build/database"
FASTA_FILE = os.path.join(DB_DIR, "EC-K-typing_group1and4_v1.0.fasta")
MAPPING_FILE = os.path.join(DB_DIR, "KL_G1G4_mapping.tsv")
OUTPUT_GBK = os.path.join(DB_DIR, "EC-K-typing_group1and4_v1.0.gbk")

KLEB_REF_GBK = "/tmp/kleb_k_ref.gbk"


def load_mapping():
    """Load KL mapping to get source assembly info."""
    mapping = {}
    with open(MAPPING_FILE) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                mapping[parts[0]] = {
                    'kx_origin': parts[1],
                    'source_assembly': parts[2],
                    'length': int(parts[3])
                }
    return mapping


def extract_kleb_reference_proteins():
    """Extract all CDS protein sequences from Klebsiella K-locus reference database.

    Returns dict: gene_name -> (product, [protein_sequences])
    Also writes a combined protein FASTA for BLAST.
    """
    gene_info = {}  # gene_name -> product
    proteins = []   # list of (gene_name, protein_seq, unique_id)

    idx = 0
    for rec in SeqIO.parse(KLEB_REF_GBK, 'genbank'):
        for f in rec.features:
            if f.type == 'CDS' and 'gene' in f.qualifiers:
                gene = f.qualifiers['gene'][0]
                product = f.qualifiers.get('product', ['hypothetical protein'])[0]

                # Extract and translate
                nuc_seq = f.location.extract(rec.seq)
                if len(nuc_seq) % 3 != 0:
                    continue
                try:
                    prot_seq = str(nuc_seq.translate(to_stop=True))
                except Exception:
                    continue

                if len(prot_seq) < 20:
                    continue

                if gene not in gene_info:
                    gene_info[gene] = product

                proteins.append((gene, prot_seq, f"ref_{idx}"))
                idx += 1

    # Write protein FASTA
    ref_prot_file = tempfile.NamedTemporaryFile(
        mode='w', suffix='.faa', delete=False, prefix='kleb_ref_'
    )
    for gene, prot, uid in proteins:
        ref_prot_file.write(f">{uid} {gene}\n{prot}\n")
    ref_prot_file.close()

    # Build BLAST database
    subprocess.run(
        ['makeblastdb', '-in', ref_prot_file.name, '-dbtype', 'prot'],
        capture_output=True
    )

    # Create lookup: uid -> (gene_name, product)
    uid_to_gene = {uid: (gene, gene_info[gene]) for gene, _, uid in proteins}

    return ref_prot_file.name, uid_to_gene


def predict_and_annotate(sequence, ref_db, uid_to_gene):
    """Predict genes with pyrodigal and annotate via BLAST against Klebsiella reference.

    Returns list of (start, end, strand, gene_name, product).
    """
    orf_finder = pyrodigal.GeneFinder(meta=True)
    genes = orf_finder.find_genes(bytes(sequence))

    # Write predicted proteins to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as pf:
        prot_file = pf.name
        for i, gene in enumerate(genes):
            prot = gene.translate()
            if prot and prot[-1] == '*':
                prot = prot[:-1]
            pf.write(f">cds_{i}\n{prot}\n")

    # BLAST against Klebsiella reference
    annotations = {}  # gene_index -> (gene_name, product, identity, coverage)

    try:
        result = subprocess.run(
            ['blastp', '-query', prot_file, '-db', ref_db,
             '-outfmt', '6 qseqid sseqid pident qcovs evalue bitscore',
             '-evalue', '1e-5', '-max_target_seqs', '5'],
            capture_output=True, text=True
        )

        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            qid = parts[0]  # cds_N
            sid = parts[1]  # ref_N
            pident = float(parts[2])
            qcovs = float(parts[3])
            bitscore = float(parts[5])

            gene_idx = int(qid.split('_')[1])
            ref_gene, ref_product = uid_to_gene.get(sid, (None, None))

            if ref_gene is None:
                continue

            # Accept hits with >=30% identity and >=50% coverage
            if pident >= 30 and qcovs >= 50:
                # Keep best hit by bitscore
                if gene_idx not in annotations or bitscore > annotations[gene_idx][3]:
                    annotations[gene_idx] = (ref_gene, ref_product, pident, bitscore)

    finally:
        os.unlink(prot_file)

    # Build result list
    result_genes = []
    for i, gene in enumerate(genes):
        start = gene.begin - 1  # 0-based
        end = gene.end
        strand = gene.strand

        # Fix CDS length
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
            gene_name, product, pident, _ = annotations[i]
        else:
            gene_name = None
            product = 'hypothetical protein'

        result_genes.append((start, end, strand, gene_name, product))

    return result_genes


def create_genbank_record(locus_name, sequence, source_assembly, gene_list):
    """Create a Kaptive-compatible GenBank record."""
    record = SeqRecord(
        Seq(str(sequence)),
        id=f"{locus_name}_{source_assembly}",
        name=f"{locus_name}_{source_assembly}",
        description=f"Escherichia coli capsular polysaccharide synthesis gene cluster, {locus_name}",
    )
    record.annotations['molecule_type'] = 'DNA'
    record.annotations['topology'] = 'linear'
    record.annotations['data_file_division'] = 'BCT'
    record.annotations['date'] = date.today().strftime('%d-%b-%Y').upper()
    record.annotations['organism'] = 'Escherichia coli'
    record.annotations['taxonomy'] = [
        'Bacteria', 'Pseudomonadota', 'Gammaproteobacteria',
        'Enterobacterales', 'Enterobacteriaceae', 'Escherichia'
    ]
    record.annotations['source'] = 'Escherichia coli'

    # Source feature
    source_qualifiers = OrderedDict()
    source_qualifiers['organism'] = ['Escherichia coli']
    source_qualifiers['mol_type'] = ['genomic DNA']
    source_qualifiers['note'] = [f'K locus: {locus_name}']

    source_feature = SeqFeature(
        FeatureLocation(0, len(sequence)),
        type='source',
        qualifiers=source_qualifiers
    )
    record.features.append(source_feature)

    # CDS features
    gene_name_counts = {}

    for i, (start, end, strand, gene_name, product) in enumerate(gene_list):
        gene_num = str(i + 1).zfill(5)
        locus_tag = f"{locus_name}_{gene_num}"

        # Ensure unique gene names within this locus
        display_gene = gene_name
        if gene_name:
            if gene_name in gene_name_counts:
                gene_name_counts[gene_name] += 1
                display_gene = f"{gene_name}_{gene_name_counts[gene_name]}"
            else:
                gene_name_counts[gene_name] = 1

        qualifiers = OrderedDict()
        qualifiers['locus_tag'] = [locus_tag]
        if display_gene:
            qualifiers['gene'] = [display_gene]
        qualifiers['product'] = [product]

        cds_feature = SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type='CDS',
            qualifiers=qualifiers
        )
        record.features.append(cds_feature)

    return record


def main():
    print("Loading Klebsiella K-locus reference proteins...")
    ref_db, uid_to_gene = extract_kleb_reference_proteins()
    print(f"  Built BLAST database with {len(uid_to_gene)} reference proteins")

    mapping = load_mapping()

    print(f"\nReading reference loci from {FASTA_FILE}...")
    records = list(SeqIO.parse(FASTA_FILE, 'fasta'))
    print(f"  Found {len(records)} loci")

    genbank_records = []

    for rec in records:
        locus_name = rec.id

        if locus_name in mapping:
            source_assembly = mapping[locus_name]['source_assembly']
        else:
            source_assembly = 'unknown'

        print(f"  Annotating {locus_name} ({len(rec.seq)} bp)...", end='')

        gene_list = predict_and_annotate(rec.seq, ref_db, uid_to_gene)

        n_cds = len(gene_list)
        n_annotated = sum(1 for g in gene_list if g[3] is not None)
        print(f" {n_cds} CDS, {n_annotated} annotated ({100*n_annotated/n_cds:.0f}%)")

        gbk_record = create_genbank_record(locus_name, rec.seq, source_assembly, gene_list)
        genbank_records.append(gbk_record)

    # Clean up temp files
    os.unlink(ref_db)
    for ext in ['.phr', '.pin', '.psq', '.pdb', '.pot', '.ptf', '.pto']:
        try:
            os.unlink(ref_db + ext)
        except FileNotFoundError:
            pass

    print(f"\nWriting GenBank file to {OUTPUT_GBK}...")
    with open(OUTPUT_GBK, 'w') as out:
        SeqIO.write(genbank_records, out, 'genbank')

    # Validate
    print("\nValidation:")
    import re
    _LOCUS_REGEX = re.compile(r'(?<=locus:)\w+|(?<=locus: ).*')

    parsed = list(SeqIO.parse(OUTPUT_GBK, 'genbank'))
    all_ok = True
    total_cds = 0
    total_annotated = 0

    for rec in parsed:
        n_cds = sum(1 for f in rec.features if f.type == 'CDS')
        n_ann = sum(1 for f in rec.features if f.type == 'CDS' and 'gene' in f.qualifiers)
        total_cds += n_cds
        total_annotated += n_ann

        source = [f for f in rec.features if f.type == 'source']
        has_locus = False
        if source and 'note' in source[0].qualifiers:
            for note in source[0].qualifiers['note']:
                if _LOCUS_REGEX.search(note):
                    has_locus = True

        if not has_locus or n_cds == 0:
            print(f"  [FAIL] {rec.id}")
            all_ok = False

    print(f"  {len(parsed)} records parsed back")
    print(f"  {total_cds} total CDS, {total_annotated} with gene names ({100*total_annotated/total_cds:.1f}%)")
    print(f"  All checks passed: {all_ok}")


if __name__ == '__main__':
    main()
