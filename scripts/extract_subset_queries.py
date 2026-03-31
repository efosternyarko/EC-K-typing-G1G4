#!/usr/bin/env python3
"""
Extract query FASTA sequences for the three high-priority subset loci,
ready for LexicMap search on ATB via find_subset_reps_ec2.sh.

Targets:
  KL713 — superset with 100-N cross-contig artefact
  KL601 — subset, 3'-end truncated (missing 3 genes vs KL414)
  KL742 — subset, 3'-end truncated (missing 8 genes vs KL394)

Output: DB/subset_queries/<KL>_query.fasta
"""
import re
from pathlib import Path
from Bio import SeqIO

REPO_ROOT = Path(__file__).resolve().parents[1]
GBK = REPO_ROOT / "DB" / "EC-K-typing_group1and4_v1.0.gbk"
OUT_DIR = REPO_ROOT / "DB" / "subset_queries"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGETS = {"KL713", "KL601", "KL742"}

for record in SeqIO.parse(GBK, "genbank"):
    kl = None
    for feat in record.features:
        if feat.type == "source":
            for note in feat.qualifiers.get("note", []):
                if note.startswith("K type:"):
                    kl = note.split("K type:")[-1].strip()
    if kl is None:
        kl = record.name.split("_")[0]
    if kl not in TARGETS:
        continue

    seq = str(record.seq)
    has_spacer = bool(re.search("N{20,}", seq))
    n_cds = sum(1 for f in record.features if f.type == "CDS")

    out_path = OUT_DIR / f"{kl}_query.fasta"
    with open(out_path, "w") as fh:
        fh.write(f">{kl} len={len(seq)} cds={n_cds} cross_contig={has_spacer}\n{seq}\n")
    print(f"  {kl}: {len(seq)} bp, {n_cds} CDS, cross-contig={has_spacer} → {out_path.name}")

print(f"\nQuery FASTAs written to {OUT_DIR}/")
