#!/usr/bin/env python3
"""Sanity check the v0.6 database."""
from Bio import SeqIO
import re

def get_locus_name(rec):
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]

gbk = "/Users/LSHEF4/Dropbox/work_2025/e_coli_db/EC-K-typing-G1G4/DB/EC-K-typing_all_groups_v0.6.gbk"
CONSERVED = {"galF", "galF_2", "gnd", "ugd", "wza", "wzb", "wzc"}

records = list(SeqIO.parse(gbk, "genbank"))
print(f"Total records: {len(records)}")

g1g4, g2g3 = [], []
for r in records:
    ln = get_locus_name(r)
    if re.match(r"KL[3-9]|^K24$|^K96$", ln):
        g1g4.append((ln, r))
    else:
        g2g3.append((ln, r))
print(f"G1/G4: {len(g1g4)}  G2/G3: {len(g2g3)}")

# ── 1. K type notes ──────────────────────────────────────────────────────────
missing_ktype = []
for ln, rec in g1g4 + g2g3:
    for f in rec.features:
        if f.type == "source":
            notes = f.qualifiers.get("note", [])
            has_ktype = any("K type:" in n for n in notes)
            if not has_ktype:
                missing_ktype.append(ln)
if missing_ktype:
    print(f"WARNING: Missing K type notes in {len(missing_ktype)} records: {missing_ktype[:5]}")
else:
    print(f"K type notes: all {len(g1g4) + len(g2g3)} records have K type notes  [PASS]")

# ── 2. Conserved gene counts ─────────────────────────────────────────────────
kl301_cons = None
zero_cons_others = []
for ln, rec in g1g4:
    n_cons = sum(
        1 for f in rec.features
        if f.type == "CDS" and f.qualifiers.get("gene", [""])[0] in CONSERVED
    )
    if ln == "KL301":
        kl301_cons = n_cons
    elif n_cons == 0:
        zero_cons_others.append(ln)

print(f"KL301 conserved CDS: {kl301_cons} (expected 0)  {'[PASS]' if kl301_cons == 0 else '[FAIL]'}")
if zero_cons_others:
    print(f"WARNING: Other G1/G4 loci with 0 conserved CDS: {zero_cons_others}")
else:
    print("All other 92 G1/G4 loci: have conserved CDS  [PASS]")

# ── 3. CDS lengths all multiples of 3 ────────────────────────────────────────
bad_cds = []
for ln, rec in g1g4:
    for f in rec.features:
        if f.type == "CDS":
            gene = f.qualifiers.get("gene", ["?"])[0]
            seq_len = len(f.extract(rec.seq))
            if seq_len % 3 != 0:
                bad_cds.append(f"{ln}/{gene} ({seq_len} bp)")
if bad_cds:
    print(f"WARNING: CDS not multiple of 3: {bad_cds}")
else:
    print("All G1/G4 CDS lengths are multiples of 3  [PASS]")

# ── 4. Total CDS counts ───────────────────────────────────────────────────────
total_cds = sum(
    sum(1 for f in r.features if f.type == "CDS")
    for _, r in g1g4
)
print(f"Total G1/G4 CDS: {total_cds}  (v0.3.1 had 3304, v0.5 had ~2776 stripped)")

# ── 5. KL306 / KL307 ugd check ───────────────────────────────────────────────
for locus_name in ("KL306", "KL307"):
    for ln, rec in g1g4:
        if ln == locus_name:
            for f in rec.features:
                if f.type == "CDS" and f.qualifiers.get("gene", [""])[0] == "ugd":
                    seq_len = len(f.extract(rec.seq))
                    print(f"{locus_name} ugd: {seq_len} bp, loc={f.location}  "
                          f"{'[PASS - partial OK]' if seq_len % 3 == 0 else '[FAIL - not mult of 3]'}")
