#!/usr/bin/env python3
"""
Build E. coli Group 1 & 4 K-locus reference database.

Uses ALL 1,112 no-hits genomes from EnteroBase (bloodstream infection isolates
with no hit to the Gladstone Group 2/3 database).

Pipeline:
  1. BLAST flanking genes (galF, gnd, ugd, wza, wzc) against all genomes
  2. Extract the cps locus region between flanking genes
  3. Handle fragmented assemblies (galF/gnd on different contigs)
  4. Cluster extracted loci at 95% identity / 80% coverage
  5. Select representative (longest) per cluster
  6. Output FASTA database with KL nomenclature

Paths are hard-coded for the original build environment; update BASE_DIR
and WORK_DIR for local use.
"""

import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# ── Configuration ──
BASE_DIR = Path("/Users/LSHEF4/Dropbox/work_2025/e_coli_db")
WORK_DIR = BASE_DIR / "db_build_v2"
GENOMES_DIR = WORK_DIR / "nohits_genomes"
LOCI_DIR = WORK_DIR / "extracted_loci"
DB_DIR = WORK_DIR / "database"
FLANKING_FA = BASE_DIR / "db_build" / "flanking_genes" / "flanking_genes.fasta"
SUMMARY_CSV = BASE_DIR / "G1_and_G4_fastkaptive_summary.csv"

MIN_LOCUS_LEN = 5000
MAX_LOCUS_LEN = 60000  # increased from 50k to capture larger loci
BLAST_EVALUE = 1e-10
CLUSTER_IDENTITY = 0.95
CLUSTER_COVERAGE = 80


def load_fastkaptive_types(summary_csv):
    """Load FastKaptive KX type assignments for all assemblies."""
    df = pd.read_csv(summary_csv)
    # Map assembly name to KX type
    type_map = {}
    for _, row in df.iterrows():
        asm = str(row['assembly']).strip()
        kx = str(row['best.match']).strip() if pd.notna(row['best.match']) else 'unknown'
        type_map[asm] = kx
    return type_map


def blast_flanking(genome_fa, flanking_fa):
    """BLAST flanking gene references against a genome assembly."""
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
    """Find the cps locus region.

    Strategy:
    1. Try same-contig: galF + gnd on same contig (ideal case)
    2. Fallback: galF and gnd on different contigs (fragmented assembly)

    Returns dict with region info or None.
    """
    by_contig = defaultdict(list)
    for h in hits:
        if h['pident'] >= 80 and h['length'] >= 200:
            by_contig[h['contig']].append(h)

    # ── Try same-contig extraction ──
    best_region = None
    best_score = 0

    for contig, contig_hits in by_contig.items():
        galf_hits = [h for h in contig_hits if h['gene'] == 'galF']
        gnd_hits = [h for h in contig_hits if h['gene'] == 'gnd']
        wza_hits = [h for h in contig_hits if h['gene'] == 'wza']

        if not galf_hits or not gnd_hits:
            continue

        for gf in galf_hits:
            for gd in gnd_hits:
                positions = [gf['sstart'], gf['send'], gd['sstart'], gd['send']]

                for wz in wza_hits:
                    positions.extend([wz['sstart'], wz['send']])

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
                            'type': 'same_contig',
                            'galF_pident': gf['pident'],
                            'gnd_pident': gd['pident'],
                            'contig_len': gf['slen'],
                        }

    if best_region:
        return best_region

    # ── Fallback: cross-contig extraction ──
    all_galf = [h for h in hits if h['gene'] == 'galF' and h['pident'] >= 80 and h['length'] >= 200]
    all_gnd = [h for h in hits if h['gene'] == 'gnd' and h['pident'] >= 80 and h['length'] >= 200]
    all_wza = [h for h in hits if h['gene'] == 'wza' and h['pident'] >= 80 and h['length'] >= 200]

    if not all_galf or not all_gnd:
        return None

    # Pick best galF and gnd hits on different contigs
    best_gf = max(all_galf, key=lambda h: h['bitscore'])
    best_gd = max(all_gnd, key=lambda h: h['bitscore'])

    if best_gf['contig'] == best_gd['contig']:
        return None  # same contig but didn't pass length filter above

    return {
        'contig': best_gf['contig'],
        'contig2': best_gd['contig'],
        'start': min(best_gf['sstart'], best_gf['send']),
        'end': max(best_gf['sstart'], best_gf['send']),
        'start2': min(best_gd['sstart'], best_gd['send']),
        'end2': max(best_gd['sstart'], best_gd['send']),
        'type': 'cross_contig',
        'galF_pident': best_gf['pident'],
        'gnd_pident': best_gd['pident'],
        'contig_len': best_gf['slen'],
        'contig2_len': best_gd['slen'],
        'length': 0,  # unknown for cross-contig
    }


def extract_locus_sequence(genome_fa, region):
    """Extract the locus region sequence from the genome."""
    records = {r.id: r for r in SeqIO.parse(genome_fa, 'fasta')}

    if region['type'] == 'same_contig':
        record = records.get(region['contig'])
        if not record:
            return None
        start = max(0, region['start'] - 500)
        end = min(len(record.seq), region['end'] + 500)
        return str(record.seq[start:end])

    elif region['type'] == 'cross_contig':
        rec1 = records.get(region['contig'])
        rec2 = records.get(region['contig2'])
        if not rec1 or not rec2:
            return None

        # Extract from galF contig: galF position to end of contig
        gf_pos = min(region['start'], region['end'])
        start1 = max(0, gf_pos - 500)
        seq1 = str(rec1.seq[start1:])

        # Extract from gnd contig: start of contig to gnd position
        gd_pos = max(region['start2'], region['end2'])
        end2 = min(len(rec2.seq), gd_pos + 500)
        seq2 = str(rec2.seq[:end2])

        # Concatenate with N-spacer
        combined = seq1 + 'N' * 100 + seq2
        if len(combined) < MIN_LOCUS_LEN or len(combined) > MAX_LOCUS_LEN:
            return None
        return combined

    return None


def cluster_loci(loci_fasta, identity_threshold=0.95, coverage_threshold=80):
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
         '-evalue', '1e-20',
         '-outfmt', '6 qseqid sseqid pident qcovs',
         '-max_target_seqs', '5000',
         '-num_threads', '4'],
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
        if q != s and pident >= identity_threshold * 100 and qcovs >= coverage_threshold:
            pairs[q][s] = pident
            pairs[s][q] = pident

    # Parse sequences for length info
    all_seqs = list(SeqIO.parse(loci_fasta, 'fasta'))
    seq_lens = {r.id: len(r.seq) for r in all_seqs}
    seq_ids_sorted = sorted(seq_lens, key=lambda x: -seq_lens[x])

    # Greedy clustering (longest first)
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
    print("Using ALL 1,112 no-hits genomes from EnteroBase")
    print("=" * 70)

    # Check dependencies
    if not FLANKING_FA.exists():
        print(f"ERROR: {FLANKING_FA} not found")
        sys.exit(1)

    # Load FastKaptive type assignments
    print("\n[Step 0] Loading FastKaptive type assignments...")
    type_map = load_fastkaptive_types(SUMMARY_CSV)
    print(f"  Loaded types for {len(type_map)} assemblies")

    # Get list of genomes
    genome_files = sorted(GENOMES_DIR.glob("*.fa"))
    print(f"\n[Step 1] Found {len(genome_files)} genomes in {GENOMES_DIR.name}/")

    # ── Step 2: BLAST and extract loci ──
    print(f"\n[Step 2] Screening all genomes for G1/G4 loci...")
    LOCI_DIR.mkdir(parents=True, exist_ok=True)

    extraction_results = []
    all_loci_records = []
    kx_extracted = defaultdict(list)
    n_same_contig = 0
    n_cross_contig = 0
    n_no_locus = 0

    for i, genome_fa in enumerate(genome_files):
        asm_name = genome_fa.name  # e.g. ESC_AA8074AA_AS.fa
        asm_base = asm_name.replace('.fa', '')

        # Look up KX type
        kx = type_map.get(asm_name, type_map.get(asm_base, 'unknown'))

        hits = blast_flanking(genome_fa, FLANKING_FA)
        region = find_cps_locus(hits)

        if region is None:
            n_no_locus += 1
            extraction_results.append({
                'assembly': asm_name, 'kx_type': kx,
                'status': 'no_locus', 'extract_type': '-',
                'length': 0, 'galF_pident': 0, 'gnd_pident': 0,
            })
            continue

        seq = extract_locus_sequence(genome_fa, region)
        if seq is None or len(seq) < MIN_LOCUS_LEN:
            extraction_results.append({
                'assembly': asm_name, 'kx_type': kx,
                'status': 'too_short', 'extract_type': region['type'],
                'length': len(seq) if seq else 0,
                'galF_pident': region['galF_pident'],
                'gnd_pident': region['gnd_pident'],
            })
            continue

        if region['type'] == 'same_contig':
            n_same_contig += 1
        else:
            n_cross_contig += 1

        locus_id = f"{kx}__{asm_base}"
        record = SeqRecord(
            Seq(seq),
            id=locus_id,
            name=locus_id,
            description=f"{region['type']} locus {len(seq)}bp from {asm_name}"
        )
        all_loci_records.append(record)
        kx_extracted[kx].append(locus_id)

        # Save individual locus file
        locus_fa = LOCI_DIR / f"{locus_id}.fa"
        SeqIO.write([record], locus_fa, 'fasta')

        extraction_results.append({
            'assembly': asm_name, 'kx_type': kx,
            'status': 'extracted', 'extract_type': region['type'],
            'length': len(seq),
            'galF_pident': region['galF_pident'],
            'gnd_pident': region['gnd_pident'],
        })

        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(genome_files)} genomes... "
                  f"({len(all_loci_records)} loci extracted)")

    # Save extraction summary
    results_df = pd.DataFrame(extraction_results)
    results_df.to_csv(WORK_DIR / 'extraction_summary.tsv', sep='\t', index=False)

    n_extracted = len(all_loci_records)
    print(f"\n  Extraction complete:")
    print(f"    Total genomes:      {len(genome_files)}")
    print(f"    Loci extracted:     {n_extracted}")
    print(f"      Same-contig:      {n_same_contig}")
    print(f"      Cross-contig:     {n_cross_contig}")
    print(f"    No locus found:     {n_no_locus}")
    print(f"    KX types covered:   {len(kx_extracted)}")
    for kx in sorted(kx_extracted):
        print(f"      {kx}: {len(kx_extracted[kx])} loci")

    if n_extracted == 0:
        print("\n  ERROR: No loci extracted!")
        sys.exit(1)

    # Save combined FASTA
    combined_fasta = WORK_DIR / 'all_extracted_loci.fasta'
    SeqIO.write(all_loci_records, combined_fasta, 'fasta')
    print(f"\n  Combined FASTA: {combined_fasta}")

    # ── Step 3: Cluster loci ──
    print(f"\n[Step 3] Clustering {n_extracted} loci at {CLUSTER_IDENTITY*100:.0f}% identity, "
          f"{CLUSTER_COVERAGE}% coverage...")
    clusters = cluster_loci(combined_fasta, identity_threshold=CLUSTER_IDENTITY,
                           coverage_threshold=CLUSTER_COVERAGE)

    print(f"  Found {len(clusters)} clusters:")
    cluster_reps = []
    cluster_info = []
    for cid, members in sorted(clusters.items()):
        rep = members[0]  # longest sequence
        kx_types = set(m.split('__')[0] for m in members)
        kx_label = '/'.join(sorted(kx_types))
        print(f"    Cluster {cid}: {len(members)} members ({kx_label}), rep={rep}")
        cluster_reps.append(rep)
        cluster_info.append({
            'cluster': cid, 'representative': rep,
            'kx_types': kx_label, 'n_members': len(members),
            'members': ';'.join(members),
        })

    # ── Step 4: Build reference database ──
    print(f"\n[Step 4] Building reference database...")
    DB_DIR.mkdir(parents=True, exist_ok=True)

    loci_dict = {r.id: r for r in all_loci_records}
    db_records = []

    # Assign KL nomenclature
    kl_counter = 300
    nomenclature = []

    for rep_id in cluster_reps:
        if rep_id not in loci_dict:
            continue
        record = loci_dict[rep_id]
        kx = rep_id.split('__')[0]

        # Direct K-type calls keep their name
        if kx.startswith('K') and not kx.startswith('KX'):
            kl_name = kx
        else:
            kl_name = f"KL{kl_counter}"
            kl_counter += 1

        asm_source = rep_id.split('__')[1] if '__' in rep_id else rep_id
        locus_len = len(record.seq)

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
    db_fasta = DB_DIR / "EC-K-typing_group1and4_v2.0.fasta"
    SeqIO.write(db_records, db_fasta, 'fasta')

    # Write nomenclature mapping
    mapping_df = pd.DataFrame(nomenclature)
    mapping_file = DB_DIR / "KL_G1G4_mapping.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)

    # Write cluster info
    cluster_df = pd.DataFrame(cluster_info)
    cluster_df.to_csv(DB_DIR / "cluster_info.tsv", sep='\t', index=False)

    # ── Summary ──
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"  Input genomes:         {len(genome_files)}")
    print(f"  Loci extracted:        {n_extracted}")
    print(f"  Sequence clusters:     {len(clusters)}")
    print(f"  Final DB references:   {len(db_records)}")
    print(f"  Database FASTA:        {db_fasta}")
    print(f"  Mapping file:          {mapping_file}")

    print(f"\n  Reference loci:")
    for entry in nomenclature:
        print(f"    {entry['KL']:8s}  ({entry['KX_origin']:6s})  "
              f"{entry['length_bp']:6d}bp  from {entry['source_assembly']}")

    print(f"\n{'='*70}")


if __name__ == '__main__':
    main()
