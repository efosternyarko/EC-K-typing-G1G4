#!/usr/bin/env python3
"""
Annotate Group 1/4 K-locus reference sequences and generate Kaptive-compatible GenBank file.

Uses pyrodigal for gene prediction and BLAST for known gene annotation.
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
FLANKING_GENES = "/Users/LSHEF4/Dropbox/work_2025/e_coli_db/db_build/flanking_genes/flanking_genes.fasta"
OUTPUT_GBK = os.path.join(DB_DIR, "EC-K-typing_group1and4_v1.0.gbk")
MAPPING_FILE = os.path.join(DB_DIR, "KL_G1G4_mapping.tsv")


def load_mapping():
    """Load KL mapping to get source assembly info."""
    mapping = {}
    with open(MAPPING_FILE) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                kl, kx_origin, source_asm, length = parts[0], parts[1], parts[2], parts[3]
                mapping[kl] = {
                    'kx_origin': kx_origin,
                    'source_assembly': source_asm,
                    'length': int(length)
                }
    return mapping


def predict_genes(sequence):
    """Predict genes using pyrodigal in metagenomic mode."""
    orf_finder = pyrodigal.GeneFinder(meta=True)
    genes = orf_finder.find_genes(bytes(sequence))
    return genes


def blast_annotate(predicted_proteins, flanking_fasta):
    """BLAST predicted proteins against known capsule genes to annotate them.

    Uses tblastn: protein query vs nucleotide subject (flanking genes).
    Returns dict mapping protein_index -> (gene_name, product, identity).
    """
    annotations = {}

    # Write predicted proteins to temp FASTA
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as prot_f:
        prot_file = prot_f.name
        for i, (prot_seq, start, end, strand) in enumerate(predicted_proteins):
            prot_f.write(f">gene_{i}\n{prot_seq}\n")

    # Known gene products
    gene_products = {
        'galF': 'UTP--glucose-1-phosphate uridylyltransferase',
        'gnd': 'UDP-glucose 6-dehydrogenase Gnd',
        'ugd': 'UDP-glucose 6-dehydrogenase Ugd',
        'wza': 'Polysaccharide export outer membrane protein Wza',
        'wzc': 'Tyrosine-protein kinase Wzc',
    }

    try:
        # Run tblastn: protein query vs nucleotide subject
        result = subprocess.run(
            ['tblastn', '-query', prot_file, '-subject', flanking_fasta,
             '-outfmt', '6 qseqid sseqid pident qcovs evalue',
             '-evalue', '1e-10'],
            capture_output=True, text=True
        )

        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            qid, sid, pident, qcovs, evalue = parts
            gene_idx = int(qid.split('_')[1])
            pident = float(pident)
            qcovs = float(qcovs)

            # Accept hits with >=50% identity and >=70% coverage
            if pident >= 50 and qcovs >= 70:
                gene_name = sid  # flanking gene names are the FASTA headers
                product = gene_products.get(gene_name, gene_name)

                # Keep best hit per gene
                if gene_idx not in annotations or pident > annotations[gene_idx][2]:
                    annotations[gene_idx] = (gene_name, product, pident)

    finally:
        os.unlink(prot_file)

    return annotations


def create_genbank_record(locus_name, sequence, source_assembly, kx_origin):
    """Create a Kaptive-compatible GenBank record for one locus."""

    # Predict genes
    genes = predict_genes(sequence)

    # Collect predicted proteins for BLAST annotation
    predicted_proteins = []
    for gene in genes:
        prot = gene.translate()
        # Remove trailing stop if present
        if prot and prot[-1] == '*':
            prot = prot[:-1]
        predicted_proteins.append((prot, gene.begin, gene.end, gene.strand))

    # BLAST annotate against known genes
    annotations = blast_annotate(predicted_proteins, FLANKING_GENES)

    # Create SeqRecord
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

    # Source feature (mandatory for Kaptive)
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
    gene_name_counts = {}  # Track gene name usage for uniqueness

    for i, gene in enumerate(genes):
        # pyrodigal uses 1-based coordinates; BioPython uses 0-based
        start = gene.begin - 1  # Convert to 0-based
        end = gene.end
        strand = gene.strand  # 1 or -1

        # Check CDS length is multiple of 3
        cds_len = end - start
        if cds_len % 3 != 0:
            # Trim to nearest multiple of 3
            excess = cds_len % 3
            if strand == 1:
                end -= excess
            else:
                start += excess

        # Skip if too short after trimming
        if end - start < 60:
            continue

        # Gene number (1-based, zero-padded)
        gene_num = str(i + 1).zfill(5)
        locus_tag = f"{locus_name}_{gene_num}"

        # Get annotation from BLAST
        if i in annotations:
            gene_name, product, _ = annotations[i]
            # Ensure unique gene names
            if gene_name in gene_name_counts:
                gene_name_counts[gene_name] += 1
                gene_name = f"{gene_name}_{gene_name_counts[gene_name]}"
            else:
                gene_name_counts[gene_name] = 1
        else:
            gene_name = None
            product = 'hypothetical protein'

        qualifiers = OrderedDict()
        qualifiers['locus_tag'] = [locus_tag]
        if gene_name:
            qualifiers['gene'] = [gene_name]
        qualifiers['product'] = [product]

        cds_feature = SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type='CDS',
            qualifiers=qualifiers
        )
        record.features.append(cds_feature)

    return record


def main():
    print("Loading K-locus mapping...")
    mapping = load_mapping()

    print(f"Reading reference loci from {FASTA_FILE}...")
    records = list(SeqIO.parse(FASTA_FILE, 'fasta'))
    print(f"  Found {len(records)} loci")

    genbank_records = []

    for rec in records:
        # Parse locus name from FASTA header: "K96 (K96/KX17) 46066bp from ERR4035267"
        locus_name = rec.id  # e.g., "K96"

        # Get metadata from mapping
        if locus_name in mapping:
            source_assembly = mapping[locus_name]['source_assembly']
            kx_origin = mapping[locus_name]['kx_origin']
        else:
            source_assembly = 'unknown'
            kx_origin = 'UNK'

        print(f"  Annotating {locus_name} ({len(rec.seq)} bp)...")

        gbk_record = create_genbank_record(locus_name, rec.seq, source_assembly, kx_origin)
        genbank_records.append(gbk_record)

        # Count annotated vs hypothetical
        n_cds = sum(1 for f in gbk_record.features if f.type == 'CDS')
        n_annotated = sum(1 for f in gbk_record.features
                         if f.type == 'CDS' and 'gene' in f.qualifiers)
        print(f"    -> {n_cds} CDS predicted, {n_annotated} with known gene names")

    # Write GenBank file
    print(f"\nWriting GenBank file to {OUTPUT_GBK}...")
    with open(OUTPUT_GBK, 'w') as out:
        SeqIO.write(genbank_records, out, 'genbank')

    print(f"Done! {len(genbank_records)} loci annotated.")

    # Validate: check the output can be parsed back
    print("\nValidation:")
    parsed = list(SeqIO.parse(OUTPUT_GBK, 'genbank'))
    print(f"  Parsed back {len(parsed)} records from GenBank file")

    for rec in parsed:
        n_cds = sum(1 for f in rec.features if f.type == 'CDS')
        has_source = any(f.type == 'source' for f in rec.features)
        source_notes = []
        for f in rec.features:
            if f.type == 'source' and 'note' in f.qualifiers:
                source_notes = f.qualifiers['note']
        locus_match = any('locus:' in n or 'locus: ' in n for n in source_notes)
        status = 'OK' if (has_source and locus_match and n_cds > 0) else 'ISSUE'
        print(f"  [{status}] {rec.id}: {n_cds} CDS, source={'yes' if has_source else 'NO'}, "
              f"locus_note={'yes' if locus_match else 'NO'}")


if __name__ == '__main__':
    main()
