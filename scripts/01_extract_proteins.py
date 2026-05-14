#!/usr/bin/env python3
"""
01_extract_proteins.py

Extract all CDS protein sequences from the v1.2 GenBank database, translating
from coordinates since translations are not stored inline.

Outputs:
  - all_proteins.faa       one FASTA per CDS, header: >{locus}|{cds_index}|{gene}
  - protein_metadata.tsv   per-CDS metadata (locus, cds_index, gene, product, length)
"""

import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

REPO_DIR = Path(__file__).resolve().parent.parent
GBK = REPO_DIR / "DB" / "EC-K-typing_group1and4_v1.2.gbk"
OUT_DIR = REPO_DIR / "DB" / "is_analysis"
OUT_DIR.mkdir(exist_ok=True)

FAA_OUT = OUT_DIR / "all_proteins.faa"
META_OUT = OUT_DIR / "protein_metadata.tsv"


def main():
    records = list(SeqIO.parse(GBK, "genbank"))
    print(f"Loaded {len(records)} loci from {GBK.name}")

    protein_records = []
    meta_rows = []
    n_short = 0

    for rec in records:
        locus = rec.name
        cds_idx = 0
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            gene = feat.qualifiers.get("gene", ["?"])[0]
            product = feat.qualifiers.get("product", ["hypothetical protein"])[0]

            try:
                nuc = feat.extract(rec.seq)
                # Trim to multiple of 3
                trim = len(nuc) - (len(nuc) % 3)
                nuc = nuc[:trim]
                prot = nuc.translate(to_stop=True)
            except Exception:
                continue

            if len(prot) < 30:
                n_short += 1
                continue

            pid = f"{locus}|{cds_idx}|{gene}"
            protein_records.append(SeqRecord(prot, id=pid, description=""))
            meta_rows.append(f"{locus}\t{cds_idx}\t{gene}\t{product}\t{len(prot)}")
            cds_idx += 1

    with open(FAA_OUT, "w") as fh:
        SeqIO.write(protein_records, fh, "fasta")

    with open(META_OUT, "w") as fh:
        fh.write("locus\tcds_index\tgene\tproduct\tlength_aa\n")
        fh.write("\n".join(meta_rows) + "\n")

    print(f"Wrote {len(protein_records)} proteins to {FAA_OUT}")
    print(f"  ({n_short} CDS skipped: <30 aa)")


if __name__ == "__main__":
    main()
