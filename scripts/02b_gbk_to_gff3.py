#!/usr/bin/env python3
"""
02b_gbk_to_gff3.py

Convert v1.2 GenBank records to Prokka-style GFF3 files that Panaroo accepts.

Panaroo expects one GFF3 per "genome" with:
  - Standard GFF3 feature lines (seqid, source, type, start, end, ., strand, ., attributes)
  - A ##FASTA section at the end with the nucleotide sequence
  - Locus tags unique within each file

One file per KL locus → panaroo_inputs/<locus>.gff3
"""

import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

REPO_DIR = Path(__file__).resolve().parent.parent
GBK = REPO_DIR / "DB" / "EC-K-typing_group1and4_v1.2.gbk"
OUT_DIR = REPO_DIR / "DB" / "panaroo_inputs"
OUT_DIR.mkdir(exist_ok=True)


def feature_to_gff3_line(feat, seqid, locus_tag_prefix, cds_counter):
    """Convert a Biopython SeqFeature to a GFF3 line."""
    start = int(feat.location.start) + 1   # GFF3 is 1-based
    end = int(feat.location.end)
    strand = "+" if feat.location.strand >= 0 else "-"

    gene = feat.qualifiers.get("gene", [""])[0]
    product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
    locus_tag = f"{locus_tag_prefix}_{cds_counter:05d}"

    attrs = f"ID={locus_tag};locus_tag={locus_tag}"
    if gene and gene not in ("?", ""):
        attrs += f";gene={gene}"
    attrs += f";product={product}"

    return f"{seqid}\tpyrodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attrs}"


def gbk_record_to_gff3(rec, out_path):
    """Write one GenBank record as a Prokka-style GFF3 file."""
    seqid = rec.name
    locus_prefix = rec.name.split("_")[0]   # e.g. KL300

    lines = ["##gff-version 3"]

    cds_counter = 1
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        line = feature_to_gff3_line(feat, seqid, locus_prefix, cds_counter)
        lines.append(line)
        cds_counter += 1

    # Embed sequence (required by Panaroo)
    lines.append("##FASTA")
    lines.append(f">{seqid}")
    seq = str(rec.seq)
    # Wrap at 60 chars
    for i in range(0, len(seq), 60):
        lines.append(seq[i:i+60])

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    records = list(SeqIO.parse(GBK, "genbank"))
    print(f"Converting {len(records)} GenBank records to GFF3...")

    n_written = 0
    n_cds_total = 0
    for rec in records:
        out_path = OUT_DIR / f"{rec.name}.gff3"
        n_cds = sum(1 for f in rec.features if f.type == "CDS")
        if n_cds == 0:
            print(f"  WARNING: {rec.name} has 0 CDS features — skipping")
            continue
        gbk_record_to_gff3(rec, out_path)
        n_cds_total += n_cds
        n_written += 1

    print(f"Wrote {n_written} GFF3 files to {OUT_DIR}/")
    print(f"  Total CDS features: {n_cds_total}")
    print(f"  Average CDS per locus: {n_cds_total/n_written:.1f}")


if __name__ == "__main__":
    main()
