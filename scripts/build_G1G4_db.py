#!/usr/bin/env python3
"""
Build E. coli Group 1 & 4 K-locus reference database.

Pipeline:
  1. Select best representative assemblies per KX type from FastKaptive output
  2. Extract those genomes from the tar.gz archive
  3. BLAST for capsule locus flanking genes (galF/gnd) to locate the cps region
  4. Extract the inter-flanking K-locus region (including wza-wzc upstream)
  5. Cluster extracted loci by sequence identity
  6. Output representative loci as multi-FASTA database
"""

import csv
import os
import subprocess
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# ── Configuration ──
BASE_DIR = Path("/Users/LSHEF4/Dropbox/work_2025/e_coli_db")
SUMMARY_CSV = BASE_DIR / "G1_and_G4_fastkaptive_summary.csv"
TARBALL = BASE_DIR / "bacteraemia_genomes.dir.tar.gz"
WORK_DIR = BASE_DIR / "db_build"
REPRESENTATIVES_DIR = WORK_DIR / "representative_genomes"
LOCI_DIR = WORK_DIR / "extracted_loci"
DB_DIR = WORK_DIR / "database"
FLANKING_FA = WORK_DIR / "flanking_genes" / "flanking_genes.fasta"

REPS_PER_KX = 10  # request many to compensate for truncated archive
MIN_LOCUS_LEN = 5000
MAX_LOCUS_LEN = 50000
BLAST_EVALUE = 1e-10


def select_representatives(summary_csv, n_per_type=10):
    """Select best representative assemblies per KX type."""
    df = pd.read_csv(summary_csv)
    df = df[df['best.match'].notna() & (df['best.match'] != 'NA')]
    df['best.match.cov'] = pd.to_numeric(df['best.match.cov'], errors='coerce')
    df['no.ref.genes.missing'] = pd.to_numeric(df['no.ref.genes.missing'], errors='coerce')

    representatives = {}
    for kx_type, group in df.groupby('best.match'):
        ranked = group.sort_values(
            by=['best.match.cov', 'no.ref.genes.missing'],
            ascending=[False, True]
        )
        reps = ranked.head(n_per_type)
        representatives[kx_type] = list(reps['assembly'])
    return representatives


def extract_genomes(tarball, assembly_names, output_dir):
    """Extract specific assemblies from the tar.gz archive."""
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = set(f"bacteraemia_genomes.dir/{name}" for name in assembly_names)

    extracted = {}
    print(f"  Extracting {len(targets)} genomes from {tarball.name}...")
    try:
        with tarfile.open(tarball, 'r:gz') as tar:
            for member in tar:
                if member.name in targets:
                    basename = os.path.basename(member.name)
                    member.name = basename
                    tar.extract(member, path=output_dir)
                    extracted[basename] = output_dir / basename
                    if len(extracted) == len(targets):
                        break
    except EOFError:
        print(f"  NOTE: tar.gz truncated; extracted {len(extracted)} before EOF")

    missing = set(assembly_names) - set(extracted.keys())
    if missing:
        print(f"  {len(missing)} assemblies not found (archive may be truncated)")
    print(f"  Extracted {len(extracted)} genomes")
    return extracted


def blast_flanking(genome_fa, flanking_fa):
    """BLAST flanking gene references directly against a genome assembly."""
    result = subprocess.run(
        ['blastn', '-query', str(flanking_fa), '-subject', str(genome_fa),
         '-evalue', str(BLAST_EVALUE),
         '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen'],
        capture_output=True, text=True
    )
    hits = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        hits.append({
            'gene': parts[0],
            'contig': parts[1],
            'pident': float(parts[2]),
            'length': int(parts[3]),
            'qstart': int(parts[4]),
            'qend': int(parts[5]),
            'sstart': int(parts[6]),
            'send': int(parts[7]),
            'evalue': float(parts[8]),
            'bitscore': float(parts[9]),
            'slen': int(parts[10]),
        })
    return hits


def find_cps_locus(hits):
    """Find the cps locus region between wza and gnd (the full G1/G4 locus).

    The G1/G4 capsule locus spans: wza-wzb-wzc-...-galF-<biosynthesis>-gnd-ugd
    We extract from wza to ugd/gnd to capture the full locus.
    """
    by_contig = defaultdict(list)
    for h in hits:
        if h['pident'] >= 80 and h['length'] >= 200:
            by_contig[h['contig']].append(h)

    best_region = None
    best_score = 0

    for contig, contig_hits in by_contig.items():
        # We need at least galF and gnd on the same contig
        galf_hits = [h for h in contig_hits if h['gene'] == 'galF']
        gnd_hits = [h for h in contig_hits if h['gene'] == 'gnd']
        wza_hits = [h for h in contig_hits if h['gene'] == 'wza']

        if not galf_hits or not gnd_hits:
            continue

        for gf in galf_hits:
            for gd in gnd_hits:
                # Determine the full locus extent
                positions = [gf['sstart'], gf['send'], gd['sstart'], gd['send']]

                # Also include wza if on same contig (upstream of the locus)
                wza_str = ""
                for wz in wza_hits:
                    positions.extend([wz['sstart'], wz['send']])
                    wza_str = f", wza({wz['pident']:.1f}%)"

                start = min(positions)
                end = max(positions)
                locus_len = end - start

                if MIN_LOCUS_LEN <= locus_len <= MAX_LOCUS_LEN:
                    score = gf['bitscore'] + gd['bitscore']
                    if score > best_score:
                        best_score = score
                        best_region = {
                            'contig': contig,
                            'start': start,
                            'end': end,
                            'length': locus_len,
                            'type': 'cps',
                            'flanking': (
                                f"galF({gf['pident']:.1f}%)--"
                                f"gnd({gd['pident']:.1f}%)"
                                f"{wza_str}"
                            ),
                            'contig_len': gf['slen'],
                        }

    return best_region


def extract_locus_sequence(genome_fa, region):
    """Extract the locus region sequence from the genome."""
    for record in SeqIO.parse(genome_fa, 'fasta'):
        if record.id == region['contig'] or record.name == region['contig']:
            # Add 500bp flanking on each side for context
            start = max(0, region['start'] - 500)
            end = min(len(record.seq), region['end'] + 500)
            seq = record.seq[start:end]
            return str(seq)
    return None


def cluster_loci(loci_fasta, identity_threshold=0.95):
    """Cluster locus sequences by pairwise BLAST identity."""
    if not os.path.exists(loci_fasta) or os.path.getsize(loci_fasta) == 0:
        return {}

    db_prefix = str(loci_fasta) + '_db'
    subprocess.run(
        ['makeblastdb', '-in', str(loci_fasta), '-dbtype', 'nucl', '-out', db_prefix],
        capture_output=True, check=True
    )

    result = subprocess.run(
        ['blastn', '-query', str(loci_fasta), '-db', db_prefix,
         '-evalue', '1e-20', '-outfmt', '6 qseqid sseqid pident qcovs',
         '-max_target_seqs', '100'],
        capture_output=True, text=True
    )

    pairs = defaultdict(dict)
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        q, s = parts[0], parts[1]
        pident = float(parts[2])
        qcovs = float(parts[3])
        if q != s and pident >= identity_threshold * 100 and qcovs >= 80:
            pairs[q][s] = pident
            pairs[s][q] = pident

    # Parse sequences for length info
    all_seqs = list(SeqIO.parse(loci_fasta, 'fasta'))
    seq_lens = {r.id: len(r.seq) for r in all_seqs}
    seq_ids_sorted = sorted(seq_lens, key=lambda x: -seq_lens[x])

    # Greedy clustering
    assigned = set()
    clusters = {}
    cluster_id = 0

    for sid in seq_ids_sorted:
        if sid in assigned:
            continue
        cluster_id += 1
        cluster = [sid]
        assigned.add(sid)

        for other in seq_ids_sorted:
            if other in assigned:
                continue
            if other in pairs.get(sid, {}):
                cluster.append(other)
                assigned.add(other)

        clusters[cluster_id] = cluster

    # Clean up BLAST db files
    for ext in ['.ndb', '.nhr', '.nin', '.njs', '.not', '.nsq', '.ntf', '.nto']:
        p = Path(db_prefix + ext)
        if p.exists():
            p.unlink()

    return clusters


def main():
    print("=" * 70)
    print("E. coli Group 1 & 4 K-locus Database Builder")
    print("=" * 70)

    # Ensure flanking gene FASTA exists
    if not FLANKING_FA.exists():
        print(f"ERROR: {FLANKING_FA} not found. Run NCBI download step first.")
        sys.exit(1)

    # ── Step 1: Select representatives ──
    print("\n[Step 1] Selecting representative assemblies per KX type...")
    representatives = select_representatives(SUMMARY_CSV, n_per_type=REPS_PER_KX)

    all_assemblies = []
    kx_to_assemblies = {}
    for kx, asm_list in sorted(representatives.items()):
        print(f"  {kx}: {len(asm_list)} candidates")
        all_assemblies.extend(asm_list)
        kx_to_assemblies[kx] = asm_list

    print(f"  Total: {len(all_assemblies)} candidates")

    # ── Step 2: Extract genomes ──
    print(f"\n[Step 2] Extracting genomes from {TARBALL.name}...")
    genome_paths = extract_genomes(TARBALL, all_assemblies, REPRESENTATIVES_DIR)

    if not genome_paths:
        print("  ERROR: No genomes could be extracted. Check the tar.gz archive.")
        sys.exit(1)

    # Map assemblies to KX types
    assembly_to_kx = {}
    for kx, asm_list in kx_to_assemblies.items():
        for asm in asm_list:
            assembly_to_kx[asm] = kx

    # ── Step 3: Extract K-locus regions ──
    print("\n[Step 3] Extracting K-locus regions via BLAST...")
    LOCI_DIR.mkdir(parents=True, exist_ok=True)

    extraction_results = []
    all_loci_records = []
    kx_extracted = defaultdict(list)  # track which KX types have extractions

    for asm_name, genome_path in sorted(genome_paths.items()):
        kx = assembly_to_kx.get(asm_name, 'unknown')

        hits = blast_flanking(genome_path, FLANKING_FA)
        region = find_cps_locus(hits)

        if region is None:
            extraction_results.append({
                'assembly': asm_name, 'kx_type': kx,
                'status': 'no_locus', 'locus_type': '-',
                'length': 0, 'flanking': '-',
            })
            continue

        seq = extract_locus_sequence(genome_path, region)
        if seq is None or len(seq) < MIN_LOCUS_LEN:
            extraction_results.append({
                'assembly': asm_name, 'kx_type': kx,
                'status': 'too_short', 'locus_type': region['type'],
                'length': len(seq) if seq else 0, 'flanking': region['flanking'],
            })
            continue

        asm_base = asm_name.replace('.fa', '').replace('_AS', '')
        locus_id = f"{kx}__{asm_base}"
        record = SeqRecord(
            Seq(seq),
            id=locus_id,
            name=locus_id,
            description=(
                f"{region['type']} locus {region['length']}bp "
                f"from {asm_name} [{region['flanking']}]"
            )
        )
        all_loci_records.append(record)
        kx_extracted[kx].append(locus_id)

        locus_fa = LOCI_DIR / f"{locus_id}.fa"
        SeqIO.write([record], locus_fa, 'fasta')

        print(f"  {asm_name} ({kx}): {region['length']}bp {region['type']} locus "
              f"[{region['flanking']}]")

        extraction_results.append({
            'assembly': asm_name, 'kx_type': kx,
            'status': 'extracted', 'locus_type': region['type'],
            'length': region['length'], 'flanking': region['flanking'],
        })

    # Summary of extraction
    results_df = pd.DataFrame(extraction_results)
    results_df.to_csv(WORK_DIR / 'extraction_summary.tsv', sep='\t', index=False)

    n_extracted = len(all_loci_records)
    print(f"\n  Extracted {n_extracted} loci from {len(genome_paths)} genomes")
    print(f"  KX types with extractions: {len(kx_extracted)}")
    for kx in sorted(kx_extracted):
        print(f"    {kx}: {len(kx_extracted[kx])} loci")

    kx_missing = set(kx_to_assemblies.keys()) - set(kx_extracted.keys())
    if kx_missing:
        print(f"  KX types without extractions: {sorted(kx_missing)}")

    if n_extracted == 0:
        print("\n  ERROR: No loci extracted!")
        sys.exit(1)

    # Save combined FASTA
    combined_fasta = WORK_DIR / 'all_extracted_loci.fasta'
    SeqIO.write(all_loci_records, combined_fasta, 'fasta')

    # ── Step 4: Cluster loci ──
    print("\n[Step 4] Clustering extracted loci by sequence identity (95%)...")
    clusters = cluster_loci(combined_fasta, identity_threshold=0.95)

    print(f"  Found {len(clusters)} clusters:")
    cluster_reps = []
    cluster_info = []
    for cid, members in sorted(clusters.items()):
        rep = members[0]
        kx_types = set(m.split('__')[0] for m in members)
        kx_label = '/'.join(sorted(kx_types))
        print(f"    Cluster {cid} ({kx_label}): {len(members)} members, rep={rep}")
        cluster_reps.append(rep)
        cluster_info.append({
            'cluster': cid, 'representative': rep,
            'kx_types': kx_label, 'n_members': len(members),
            'members': ';'.join(members),
        })

    # ── Step 5: Build reference database ──
    print("\n[Step 5] Building reference database...")
    DB_DIR.mkdir(parents=True, exist_ok=True)

    loci_dict = {r.id: r for r in all_loci_records}
    db_records = []

    # Assign KL nomenclature (G1/G4 types start from KL300)
    kl_counter = 300
    nomenclature = []

    for rep_id in cluster_reps:
        if rep_id not in loci_dict:
            continue
        record = loci_dict[rep_id]
        kx = rep_id.split('__')[0]

        if kx.startswith('K') and not kx.startswith('KX'):
            kl_name = kx  # direct K-type call
        else:
            kl_name = f"KL{kl_counter}"
            kl_counter += 1

        asm_source = rep_id.split('__')[1] if '__' in rep_id else rep_id
        locus_len = len(record.seq)

        # Create clean record for database
        db_record = SeqRecord(
            record.seq,
            id=kl_name,
            name=kl_name,
            description=f"{kl_name} ({kx}) {locus_len}bp from {asm_source}"
        )
        db_records.append(db_record)

        nomenclature.append({
            'KL': kl_name,
            'KX_origin': kx,
            'source_assembly': asm_source,
            'length_bp': locus_len,
            'group': 'G1/G4',
        })

    # Write database FASTA
    db_fasta = DB_DIR / "EC-K-typing_group1and4_v1.0.fasta"
    SeqIO.write(db_records, db_fasta, 'fasta')

    # Write nomenclature mapping
    mapping_df = pd.DataFrame(nomenclature)
    mapping_file = DB_DIR / "KL_G1G4_mapping.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)

    # Write cluster info
    cluster_df = pd.DataFrame(cluster_info)
    cluster_df.to_csv(DB_DIR / "cluster_info.tsv", sep='\t', index=False)

    # ── Summary ──
    print(f"\n  Database: {db_fasta}")
    print(f"  Mapping:  {mapping_file}")
    print(f"  Total reference loci: {len(db_records)}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  KX types in input:       {len(kx_to_assemblies)}")
    print(f"  Genomes extracted:       {len(genome_paths)}")
    print(f"  Loci successfully found: {n_extracted}")
    print(f"  Sequence clusters:       {len(clusters)}")
    print(f"  Final DB references:     {len(db_records)}")

    print(f"\n  Output directory: {WORK_DIR}")

    print("\n  Reference loci in database:")
    for entry in nomenclature:
        print(f"    {entry['KL']:8s}  ({entry['KX_origin']:6s})  "
              f"{entry['length_bp']:6d}bp  from {entry['source_assembly']}")

    if kx_missing:
        print(f"\n  NOTE: {len(kx_missing)} KX types could not be extracted")
        print(f"  (archive truncated). Missing: {sorted(kx_missing)}")
        print(f"  To complete the database, re-run with the full tar.gz archive")
        print(f"  or provide the missing genomes in {REPRESENTATIVES_DIR}/")

    print("\n" + "-" * 70)
    print("NEXT STEPS:")
    print(f"  1. Annotate with Prokka:")
    print(f"     for f in {LOCI_DIR}/*.fa; do")
    print(f"       prokka --kingdom Bacteria --genus Escherichia \\")
    print(f"         --outdir ${{f%.fa}}_prokka --prefix $(basename ${{f%.fa}}) \"$f\"")
    print(f"     done")
    print(f"  2. Convert Prokka .gbk to Kaptive-compatible GenBank format")
    print(f"  3. Merge with Gladstone G2/G3 database:")
    print(f"     cat EC-K-typing_group2and3_v3.0.0.gbk EC-K-typing_G1G4_v1.0.gbk \\")
    print(f"       > EC-K-typing_all_groups_v1.0.gbk")
    print(f"  4. Validate: re-type the full collection with the merged database")
    print("-" * 70)


if __name__ == '__main__':
    main()
