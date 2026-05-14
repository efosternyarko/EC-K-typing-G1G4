#!/usr/bin/env python3
"""
04_build_presence_absence.py

Build a gene presence/absence matrix from:
  - BLASTp all-vs-all results (03_blastp_allvsall.sh)
  - hmmsearch IS protein hits (02_detect_is_proteins.sh)

Protein family clustering:
  - Two proteins are in the same family if pident >= 80% AND
    both query and subject coverage >= 80% (bidirectional).
  - Union-find greedy clustering on these pairs.

IS tagging:
  - Any protein with a hmmsearch hit (E-value <= 1e-5) is flagged as IS-related.
  - A protein family is tagged IS if ANY member is IS-flagged.

Outputs:
  - protein_families.tsv          protein_id -> family_id, is_is_family (bool)
  - presence_absence.tsv          locus x family_id (1/0), IS families marked
  - presence_absence_nois.tsv     same but IS families excluded
  - is_protein_ids.txt            all protein IDs flagged as IS
"""

import csv
from collections import defaultdict
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent.parent
IS_DIR = REPO_DIR / "DB" / "is_analysis"

BLAST_TSV  = IS_DIR / "blastp_allvsall.tsv"
HMMER_TBL  = IS_DIR / "hmmsearch_is_hits.tbl"
META_TSV   = IS_DIR / "protein_metadata.tsv"

FAM_OUT    = IS_DIR / "protein_families.tsv"
PA_OUT     = IS_DIR / "presence_absence.tsv"
PA_NOIS    = IS_DIR / "presence_absence_nois.tsv"
IS_IDS_OUT = IS_DIR / "is_protein_ids.txt"

PIDENT_THR = 80.0
COV_THR    = 80.0

# ── 1. Load protein metadata ──────────────────────────────────────────────────
print("Loading protein metadata...")
meta = {}  # protein_id -> (locus, gene, product)
all_proteins = set()
loci_order = []
seen_loci = set()

with open(META_TSV) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        locus = row["locus"]
        pid = f"{locus}|{row['cds_index']}|{row['gene']}"
        meta[pid] = (locus, row["gene"], row["product"])
        all_proteins.add(pid)
        if locus not in seen_loci:
            loci_order.append(locus)
            seen_loci.add(locus)

print(f"  {len(all_proteins)} proteins across {len(loci_order)} loci")

# ── 2. Parse hmmsearch hits to flag IS proteins ───────────────────────────────
print("Parsing hmmsearch IS hits...")
is_proteins = set()

if HMMER_TBL.exists():
    with open(HMMER_TBL) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            target_name = parts[0]  # protein_id is the target sequence
            evalue = float(parts[4])
            if evalue <= 1e-5 and target_name in all_proteins:
                is_proteins.add(target_name)
else:
    print(f"  WARNING: {HMMER_TBL} not found — no IS proteins will be flagged")

print(f"  {len(is_proteins)} IS-flagged proteins")

# ── 3. Parse BLASTp hits and cluster into protein families ───────────────────
print("Parsing BLASTp all-vs-all results...")
pairs = defaultdict(set)  # protein_id -> set of proteins in same family

n_hits = 0
with open(BLAST_TSV) as fh:
    for line in fh:
        parts = line.rstrip().split("\t")
        if len(parts) < 8:
            continue
        q, s = parts[0], parts[1]
        if q == s:
            continue
        pident = float(parts[2])
        qcovs = float(parts[3])
        scovs = float(parts[4])
        if pident >= PIDENT_THR and qcovs >= COV_THR and scovs >= COV_THR:
            pairs[q].add(s)
            pairs[s].add(q)
            n_hits += 1

print(f"  {n_hits} bidirectional hit pairs at ≥{PIDENT_THR}% id / ≥{COV_THR}% cov")

# Union-find clustering
parent = {p: p for p in all_proteins}

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra

for p, neighbours in pairs.items():
    for n in neighbours:
        union(p, n)

# Assign family IDs (root-based, sorted for determinism)
roots = sorted(set(find(p) for p in all_proteins))
root_to_fam = {r: i for i, r in enumerate(roots)}
protein_to_fam = {p: root_to_fam[find(p)] for p in all_proteins}

n_families = len(roots)
print(f"  {n_families} protein families")

# Tag families as IS if any member is IS-flagged
is_families = set()
for p in is_proteins:
    is_families.add(protein_to_fam[p])
print(f"  {len(is_families)} IS-tagged protein families")

# ── 4. Write protein families table ──────────────────────────────────────────
print("Writing protein_families.tsv...")
with open(FAM_OUT, "w") as fh:
    fh.write("protein_id\tlocus\tgene\tproduct\tfamily_id\tis_family\n")
    for p in sorted(all_proteins):
        locus, gene, product = meta[p]
        fam = protein_to_fam[p]
        is_fam = fam in is_families
        fh.write(f"{p}\t{locus}\t{gene}\t{product}\t{fam}\t{is_fam}\n")

# ── 5. Build presence/absence matrix ─────────────────────────────────────────
print("Building presence/absence matrix...")

# Families that appear in at least one locus
biosyn_families = sorted(f for f in range(n_families) if f not in is_families)
is_fam_list = sorted(is_families)

# Per-locus family presence
locus_biosyn = defaultdict(set)
locus_is = defaultdict(set)

for p, fam in protein_to_fam.items():
    locus = meta[p][0]
    if fam in is_families:
        locus_is[locus].add(fam)
    else:
        locus_biosyn[locus].add(fam)

# Write full matrix (all families)
all_families = sorted(set(range(n_families)))
with open(PA_OUT, "w") as fh:
    header_fams = [str(f) for f in all_families]
    # Tag IS families in header
    header = ["locus"] + [f"f{f}{'*' if f in is_families else ''}" for f in all_families]
    fh.write("\t".join(header) + "\n")
    for locus in loci_order:
        present = locus_biosyn[locus] | locus_is[locus]
        row = [locus] + ["1" if f in present else "0" for f in all_families]
        fh.write("\t".join(row) + "\n")

# Write non-IS matrix (biosynthetic families only)
with open(PA_NOIS, "w") as fh:
    header = ["locus"] + [f"f{f}" for f in biosyn_families]
    fh.write("\t".join(header) + "\n")
    for locus in loci_order:
        row = [locus] + ["1" if f in locus_biosyn[locus] else "0" for f in biosyn_families]
        fh.write("\t".join(row) + "\n")

# Write IS protein IDs
with open(IS_IDS_OUT, "w") as fh:
    for p in sorted(is_proteins):
        fh.write(p + "\n")

print(f"\nSummary:")
print(f"  Total protein families:      {n_families}")
print(f"  IS-tagged families:          {len(is_families)}")
print(f"  Biosynthetic families:       {len(biosyn_families)}")
print(f"  Loci with >=1 IS family:     {sum(1 for l in loci_order if locus_is[l])}")
print(f"\nOutputs in {IS_DIR}/")
