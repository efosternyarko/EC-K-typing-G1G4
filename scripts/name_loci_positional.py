#!/usr/bin/env python3
"""
Systematic positional gene naming for G1/G4 K-locus references (v0.3).

Takes the annotated v2.0 GenBank, clusters all CDS proteins at 90% identity /
90% query+subject coverage, then assigns systematic gene names:

  - Conserved / functionally-annotated genes: keep the existing functional name
    (propagated to all cluster members that currently have only a locus_tag name)
  - Variable genes (unnamed clusters): KLxxx_N positional name, where xxx is the
    lowest-numbered locus that first encounters the family and N increments per locus

Proteins shared across loci receive the *same* gene name in every locus that
carries them, giving Kaptive the signal it needs to discriminate closely related
loci by their unique genes.

Usage
-----
    python scripts/name_loci_positional.py

Outputs
-------
    DB/EC-K-typing_group1and4_v3.0.gbk   -- 93 loci with positional names
    DB/EC-K-typing_all_groups_v3.0.gbk   -- 183-locus combined database
"""

import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR    = Path(__file__).resolve().parent.parent
DB_DIR      = REPO_DIR / "DB"
INPUT_GBK   = DB_DIR / "EC-K-typing_group1and4_v2.0.gbk"
OUTPUT_GBK  = DB_DIR / "EC-K-typing_group1and4_v3.0.gbk"
G23_GBK     = DB_DIR / "EC-K-typing_group2and3_v3.0.0.gbk"
MERGED_GBK  = DB_DIR / "EC-K-typing_all_groups_v3.0.gbk"

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
CLUSTER_IDENTITY = 90   # % amino acid identity
CLUSTER_COVERAGE = 90   # % coverage (both query and subject)

# Locus-tag gene name pattern: KL300_00001, K24_00005, KL423_00012, etc.
LOCUS_TAG_RE = re.compile(r'^(KL\d+|K\d+)_\d{5}$')


def is_locus_tag_name(name: str) -> bool:
    return bool(name) and LOCUS_TAG_RE.match(name) is not None


def get_locus_name(rec) -> str:
    """Extract K locus name from GenBank source feature /note."""
    for f in rec.features:
        if f.type == 'source' and 'note' in f.qualifiers:
            for note in f.qualifiers['note']:
                m = re.search(r'K locus:\s*(\S+)', note)
                if m:
                    return m.group(1)
    return rec.id.split('_')[0]


def kl_sort_key(name: str) -> int:
    """Numeric sort key for KL names (KL300 → 300, K24 → 24)."""
    m = re.search(r'(\d+)$', name)
    return int(m.group(1)) if m else 0


# ---------------------------------------------------------------------------
# Step 1 – load records and extract proteins
# ---------------------------------------------------------------------------

def load_and_extract(records):
    """
    Parse 93 loci.  For each CDS, translate and record:
      protein_id  = "{locus}__{cds_index:04d}"
      gene_name   = current /gene qualifier
      start       = genomic start (for ordering)

    Returns
    -------
    locus_cds   : {locus_name: [(start, cds_idx, pid, gene_name), ...]}  sorted by start
    proteins    : {pid: aa_seq}
    pid_to_gene : {pid: gene_name}
    pid_to_locus: {pid: locus_name}
    """
    locus_cds    = {}
    proteins     = {}
    pid_to_gene  = {}
    pid_to_locus = {}

    for rec in records:
        lname   = get_locus_name(rec)
        cds_idx = 0
        locus_cds[lname] = []

        for f in rec.features:
            if f.type != 'CDS':
                continue

            gene_name = f.qualifiers.get('gene', [''])[0]
            nuc       = f.location.extract(rec.seq)

            try:
                prot = str(nuc.translate(to_stop=True))
            except Exception:
                cds_idx += 1
                continue

            if len(prot) < 20:
                cds_idx += 1
                continue

            pid = f"{lname}__{cds_idx:04d}"
            proteins[pid]     = prot
            pid_to_gene[pid]  = gene_name
            pid_to_locus[pid] = lname

            start = int(f.location.start)
            locus_cds[lname].append((start, cds_idx, pid, gene_name))
            cds_idx += 1

        locus_cds[lname].sort()          # genomic order

    return locus_cds, proteins, pid_to_gene, pid_to_locus


# ---------------------------------------------------------------------------
# Step 2 – all-vs-all BLASTp
# ---------------------------------------------------------------------------

def allvsall_blastp(proteins: dict) -> dict:
    """
    Run all-vs-all BLASTp.

    Returns neighbors: {pid: set of pids} — pairs passing both pident ≥ 90
    and query+subject coverage ≥ 90 %.
    """
    # Write query FASTA
    fa = tempfile.NamedTemporaryFile(
        mode='w', suffix='.faa', delete=False, prefix='g14_all_'
    )
    for pid, seq in proteins.items():
        fa.write(f">{pid}\n{seq}\n")
    fa.close()

    db_pfx = fa.name + '_db'
    subprocess.run(
        ['makeblastdb', '-in', fa.name, '-dbtype', 'prot', '-out', db_pfx],
        capture_output=True, check=True
    )

    res = subprocess.run(
        ['blastp',
         '-query', fa.name, '-db', db_pfx,
         '-outfmt', '6 qseqid sseqid pident qcovs length qlen slen',
         '-evalue', '1e-5', '-num_threads', '4', '-max_target_seqs', '500'],
        capture_output=True, text=True
    )

    neighbors = defaultdict(set)
    for line in res.stdout.strip().split('\n'):
        if not line:
            continue
        parts    = line.split('\t')
        q, s     = parts[0], parts[1]
        if q == s:
            continue
        pident   = float(parts[2])
        qcovs    = float(parts[3])
        aln_len  = int(parts[4])
        slen     = int(parts[6])
        scovs    = aln_len / slen * 100 if slen else 0

        if pident >= CLUSTER_IDENTITY and qcovs >= CLUSTER_COVERAGE and scovs >= CLUSTER_COVERAGE:
            neighbors[q].add(s)
            neighbors[s].add(q)

    # Cleanup
    os.unlink(fa.name)
    for ext in ['.phr', '.pin', '.psq', '.pdb', '.pot', '.ptf', '.pto']:
        try:
            os.unlink(db_pfx + ext)
        except FileNotFoundError:
            pass

    return neighbors


# ---------------------------------------------------------------------------
# Step 3 – greedy clustering
# ---------------------------------------------------------------------------

def greedy_cluster(proteins: dict, neighbors: dict):
    """
    Greedy clustering (longest-sequence-first, CD-HIT style).

    Returns
    -------
    member_to_rep : {pid: representative_pid}
    rep_to_members: {rep_pid: [pid, ...]}
    """
    by_len = sorted(proteins, key=lambda p: -len(proteins[p]))

    member_to_rep  = {}
    rep_to_members = defaultdict(list)

    for pid in by_len:
        if pid in member_to_rep:
            continue
        rep_to_members[pid].append(pid)
        member_to_rep[pid] = pid

        for other in by_len:
            if other in member_to_rep:
                continue
            if other in neighbors.get(pid, set()):
                rep_to_members[pid].append(other)
                member_to_rep[other] = pid

    return member_to_rep, rep_to_members


# ---------------------------------------------------------------------------
# Steps 4+5 – name each protein family
# ---------------------------------------------------------------------------

def assign_family_names(
    rep_to_members: dict,
    member_to_rep: dict,
    locus_cds: dict,
    pid_to_gene: dict,
) -> dict:
    """
    Assign a gene name to every protein family (cluster).

    Pass 1  – functional names: if the cluster contains any protein whose
              gene_name is NOT a locus_tag, use the most-common such name.

    Pass 2  – positional names: iterate loci in ascending KL-number order;
              within each locus, CDS in genomic order.  The first locus to
              encounter an unnamed family becomes the "lead" and names it
              KL{N}_{counter}.

    Returns family_name: {representative_pid: gene_name}
    """
    family_name = {}

    # Pass 1 — propagate functional names
    for rep, members in rep_to_members.items():
        func = [
            pid_to_gene[m]
            for m in members
            if pid_to_gene.get(m) and not is_locus_tag_name(pid_to_gene[m])
        ]
        if func:
            family_name[rep] = Counter(func).most_common(1)[0][0]

    # Pass 2 — positional names for unnamed families
    all_loci_sorted = sorted(locus_cds.keys(), key=kl_sort_key)

    for lname in all_loci_sorted:
        kl_num     = re.search(r'(\d+)$', lname).group(1)
        pos_counter = 0

        for (_start, _idx, pid, _gene) in locus_cds[lname]:
            rep = member_to_rep.get(pid, pid)
            if rep not in family_name:
                pos_counter += 1
                family_name[rep] = f"KL{kl_num}_{pos_counter}"

    return family_name


# ---------------------------------------------------------------------------
# Step 6 – rewrite GenBank records
# ---------------------------------------------------------------------------

def rewrite_records(records, locus_cds, member_to_rep, family_name):
    """
    Return new GenBank records with updated /gene qualifiers.
    Within-locus duplicate family names get _2, _3 … suffixes.
    """
    new_records = []
    stats       = Counter()

    for rec in records:
        lname = get_locus_name(rec)

        new_rec              = SeqRecord(
            rec.seq,
            id          = rec.id,
            name        = rec.name,
            description = rec.description,
            annotations = rec.annotations.copy(),
            dbxrefs     = list(rec.dbxrefs),
        )

        cds_counter          = 0
        name_counts_in_locus = Counter()

        for f in rec.features:
            if f.type != 'CDS':
                new_rec.features.append(f)
                continue

            pid = f"{lname}__{cds_counter:04d}"
            cds_counter += 1

            if pid not in member_to_rep:
                # Short / untranslatable — keep original
                new_rec.features.append(f)
                stats['kept_short'] += 1
                continue

            rep       = member_to_rep[pid]
            base_name = family_name.get(rep, f.qualifiers.get('gene', [''])[0])

            # Deduplicate within locus
            name_counts_in_locus[base_name] += 1
            count      = name_counts_in_locus[base_name]
            final_name = base_name if count == 1 else f"{base_name}_{count}"

            old_name = f.qualifiers.get('gene', [''])[0]
            if is_locus_tag_name(old_name):
                stats['positional'] += 1
            else:
                stats['functional'] += 1

            new_q = OrderedDict()
            for key, val in f.qualifiers.items():
                new_q[key] = val
            new_q['gene'] = [final_name]

            new_rec.features.append(
                SeqFeature(location=f.location, type='CDS', qualifiers=new_q)
            )

        new_records.append(new_rec)

    return new_records, stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print('=' * 70)
    print('Systematic Positional Gene Naming  —  v0.3')
    print('=' * 70)

    # ── Step 1: load ────────────────────────────────────────────────────────
    print(f'\n[Step 1] Loading {INPUT_GBK.name}...')
    records = list(SeqIO.parse(INPUT_GBK, 'genbank'))
    locus_names = [get_locus_name(r) for r in records]
    print(f'  {len(records)} loci loaded: {locus_names[0]} … {locus_names[-1]}')

    locus_cds, proteins, pid_to_gene, pid_to_locus = load_and_extract(records)

    total_cds = sum(len(v) for v in locus_cds.values())
    unnamed   = sum(1 for g in pid_to_gene.values() if is_locus_tag_name(g))
    print(f'  {total_cds} CDS  |  {total_cds - unnamed} functional  |  {unnamed} locus_tag-named')

    # ── Step 2: all-vs-all BLASTp ───────────────────────────────────────────
    print(f'\n[Step 2] All-vs-all BLASTp ({len(proteins)} proteins)...')
    neighbors = allvsall_blastp(proteins)
    n_pairs   = sum(len(v) for v in neighbors.values()) // 2
    print(f'  {n_pairs} high-similarity pairs '
          f'(≥{CLUSTER_IDENTITY}% id, ≥{CLUSTER_COVERAGE}% cov)')

    # ── Step 3: cluster ─────────────────────────────────────────────────────
    print('\n[Step 3] Greedy clustering...')
    member_to_rep, rep_to_members = greedy_cluster(proteins, neighbors)
    n_fam     = len(rep_to_members)
    singleton = sum(1 for m in rep_to_members.values() if len(m) == 1)
    multi     = n_fam - singleton
    print(f'  {n_fam} protein families  |  {multi} multi-member  |  {singleton} singletons')

    # ── Steps 4+5: assign names ─────────────────────────────────────────────
    print('\n[Step 4] Assigning family names...')
    family_name = assign_family_names(rep_to_members, member_to_rep, locus_cds, pid_to_gene)

    POSITIONAL_RE = re.compile(r'^KL\d+_\d+$')
    n_pos  = sum(1 for n in family_name.values() if POSITIONAL_RE.match(n))
    n_func = len(family_name) - n_pos
    print(f'  {n_func} families with functional names')
    print(f'  {n_pos} families with positional names')

    # Multi-locus positional families
    rep_loci = defaultdict(set)
    for pid, rep in member_to_rep.items():
        rep_loci[rep].add(pid_to_locus[pid])

    shared = sorted(
        [(rep, loci) for rep, loci in rep_loci.items()
         if len(loci) > 1 and POSITIONAL_RE.match(family_name.get(rep, ''))],
        key=lambda x: -len(x[1])
    )
    print(f'  {len(shared)} multi-locus positional families  '
          f'(top 5 by locus count):')
    for rep, loci in shared[:5]:
        name = family_name[rep]
        loci_s = ', '.join(sorted(loci)[:6]) + ('…' if len(loci) > 6 else '')
        print(f'    {name}: {len(loci)} loci  ({loci_s})')

    # ── Step 6: rewrite records ──────────────────────────────────────────────
    print('\n[Step 5] Rewriting GenBank records...')
    new_records, stats = rewrite_records(
        records, locus_cds, member_to_rep, family_name
    )
    print(f'  {stats["functional"]} functional names kept')
    print(f'  {stats["positional"]} CDS given positional names')
    print(f'  {stats["kept_short"]} CDS kept unchanged (short/untranslatable)')

    # Quick sanity: no locus should have duplicate gene names at the same
    # position (only allowed via the _2/_3 suffix we apply above)
    for rec in new_records:
        names  = [f.qualifiers.get('gene', [''])[0]
                  for f in rec.features if f.type == 'CDS']
        dupes  = [n for n, c in Counter(names).items() if c > 1 and not n.endswith('_2')]
        if dupes:
            lname = get_locus_name(rec)
            print(f'  WARNING: {lname} has duplicated gene names: {dupes[:5]}')

    # ── Step 7: write G1/G4 GenBank ─────────────────────────────────────────
    print(f'\n[Step 6] Writing {OUTPUT_GBK.name}...')
    with open(OUTPUT_GBK, 'w') as fh:
        SeqIO.write(new_records, fh, 'genbank')
    print(f'  Written: {OUTPUT_GBK}')

    # ── Step 8: merge with G2/G3 ────────────────────────────────────────────
    print('\n[Step 7] Merging with Gladstone G2/G3 database...')
    if G23_GBK.exists():
        MERGED_GBK.write_text(G23_GBK.read_text() + '\n' + OUTPUT_GBK.read_text())
        merged = list(SeqIO.parse(MERGED_GBK, 'genbank'))
        print(f'  {len(merged)} records in merged database → {MERGED_GBK.name}')
    else:
        print(f'  WARNING: {G23_GBK} not found — skipping merge')

    # ── Step 9: summary statistics ──────────────────────────────────────────
    print(f'\n{"=" * 70}')
    print('GENE NAMING SUMMARY')
    print(f'{"=" * 70}')

    total_cds_out = func_total = pos_total = 0
    for rec in new_records:
        for f in rec.features:
            if f.type != 'CDS':
                continue
            total_cds_out += 1
            name = f.qualifiers.get('gene', [''])[0]
            if POSITIONAL_RE.match(name):
                pos_total  += 1
            else:
                func_total += 1

    print(f'  Total CDS      : {total_cds_out}')
    print(f'  Functional names: {func_total}  ({100*func_total/total_cds_out:.1f}%)')
    print(f'  Positional names: {pos_total}  ({100*pos_total/total_cds_out:.1f}%)')
    print(f'  Protein families: {n_fam}  (multi-locus positional: {len(shared)})')
    print(f'\n  Output files:')
    print(f'    G1/G4 v3.0    : {OUTPUT_GBK}')
    if G23_GBK.exists():
        print(f'    All-groups v3.0: {MERGED_GBK}')
    print(f'\nNext step:')
    print(f'  python scripts/validate_db_v2.py \\')
    print(f'    --db {MERGED_GBK} \\')
    print(f'    --suffix v3')


if __name__ == '__main__':
    main()
