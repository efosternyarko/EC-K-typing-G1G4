#!/usr/bin/env python3
"""
Validate G1/G4 K-locus database using Kaptive gene-content principle.

Implements the core Kaptive locus-definition principle:
  "A unique locus is a unique set of genes, where genes are defined at
   threshold translated amino acid identity level."

Pipeline:
  1. Extract all CDS protein sequences from the GenBank database,
     labelled by locus (KL type)
  2. Cluster proteins with CD-HIT at 95% aa identity / 80% bidirectional
     coverage — each cluster = one "gene"
  3. Represent each locus as a frozenset of gene-cluster IDs
  4. Compare gene-content sets across all loci
  5. Report: duplicate loci (identical gene sets), unique loci, and
     any loci that are strict subsets of another

Usage:
    python3 validate_gene_content.py [--gbk PATH] [--cdhit PATH] [--out DIR]

Defaults:
    --gbk   DB/EC-K-typing_group1and4_v0.9.gbk
    --cdhit /tmp/cdhit/cd-hit
    --out   /tmp/gene_content_validation/
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


# ── Defaults ──────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GBK  = REPO_ROOT / "DB" / "EC-K-typing_group1and4_v0.9.gbk"
DEFAULT_CDHIT = Path("/tmp/cdhit/cd-hit")
DEFAULT_OUT  = Path("/tmp/gene_content_validation")

CDHIT_IDENTITY = 0.95   # 95% amino acid identity
CDHIT_COVERAGE = 0.80   # 80% bidirectional coverage
CDHIT_WORDSIZE = 5       # word size for 95% identity (CD-HIT recommendation)


# ── Step 1: Extract CDS translations ──────────────────────────────────────────
def extract_proteins(gbk_path: Path, out_fasta: Path) -> dict[str, list[str]]:
    """
    Parse GenBank file, extract CDS /translation qualifiers.

    Returns:
        locus_proteins: dict mapping KL-type → list of protein sequence IDs
    Also writes a combined FASTA of all proteins to out_fasta.
    """
    locus_proteins = defaultdict(list)
    records_written = 0

    with open(out_fasta, "w") as fh:
        for record in SeqIO.parse(gbk_path, "genbank"):
            # Extract KL type from the 'K type' note qualifier on the source feature
            kl_type = None
            for feat in record.features:
                if feat.type == "source":
                    for note in feat.qualifiers.get("note", []):
                        if note.startswith("K type:"):
                            kl_type = note.split("K type:")[-1].strip()
                            break
                if kl_type:
                    break

            if kl_type is None:
                # Fall back to record name (e.g. KL300_ESC_...)
                kl_type = record.name.split("_")[0]

            for feat in record.features:
                if feat.type != "CDS":
                    continue
                locus_tag = feat.qualifiers.get("locus_tag",
                                                [f"cds_{records_written}"])[0]
                protein_id = f"{kl_type}|{locus_tag}"

                # Prefer stored /translation; fall back to translating from nt
                if "translation" in feat.qualifiers:
                    aa_seq = feat.qualifiers["translation"][0]
                else:
                    try:
                        aa_seq = str(feat.extract(record.seq)
                                     .translate(to_stop=True))
                    except Exception:
                        continue
                if not aa_seq:
                    continue

                fh.write(f">{protein_id}\n{aa_seq}\n")
                locus_proteins[kl_type].append(protein_id)
                records_written += 1

    print(f"  Extracted {records_written} CDS from "
          f"{len(locus_proteins)} loci → {out_fasta}")
    return dict(locus_proteins)


# ── Step 2: Run CD-HIT ────────────────────────────────────────────────────────
def run_cdhit(cdhit_bin: Path, in_fasta: Path, out_prefix: Path) -> Path:
    """
    Cluster protein sequences with CD-HIT at 95% aa identity / 80% coverage.

    Returns path to the CD-HIT cluster file (.clstr).
    """
    cmd = [
        str(cdhit_bin),
        "-i", str(in_fasta),
        "-o", str(out_prefix),
        "-c", str(CDHIT_IDENTITY),
        "-aS", str(CDHIT_COVERAGE),   # coverage of shorter sequence
        "-aL", str(CDHIT_COVERAGE),   # coverage of longer sequence (bidirectional)
        "-n", str(CDHIT_WORDSIZE),
        "-d", "0",                     # unlimited description length in .clstr
        "-M", "2000",                  # memory in MB
        "-T", "4",                     # threads
    ]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("CD-HIT stderr:", result.stderr, file=sys.stderr)
        sys.exit(1)

    clstr_path = Path(str(out_prefix) + ".clstr")
    print(f"  CD-HIT complete → {clstr_path}")
    return clstr_path


# ── Step 3: Parse .clstr → protein → cluster ID mapping ──────────────────────
def parse_cdhit_clusters(clstr_path: Path) -> dict[str, int]:
    """
    Parse CD-HIT .clstr file.

    Returns:
        protein_to_cluster: dict mapping protein_id → integer cluster ID
    """
    protein_to_cluster = {}
    cluster_id = None

    with open(clstr_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1])
            elif line and cluster_id is not None:
                # Format: 0	241aa, >KL300|KL300_00001... *
                parts = line.split(">")
                if len(parts) < 2:
                    continue
                protein_id = parts[1].split("...")[0].strip()
                protein_to_cluster[protein_id] = cluster_id

    n_clusters = len(set(protein_to_cluster.values()))
    print(f"  Parsed {len(protein_to_cluster)} proteins into "
          f"{n_clusters} gene clusters")
    return protein_to_cluster


# ── Step 4: Build gene-content sets per locus ─────────────────────────────────
def build_gene_content_sets(
    locus_proteins: dict[str, list[str]],
    protein_to_cluster: dict[str, int],
) -> dict[str, frozenset]:
    """
    For each locus, build its gene-content set (frozenset of cluster IDs).
    """
    locus_gene_sets = {}
    missing = []

    for kl_type, proteins in locus_proteins.items():
        clusters = set()
        for pid in proteins:
            if pid in protein_to_cluster:
                clusters.add(protein_to_cluster[pid])
            else:
                missing.append(pid)
        locus_gene_sets[kl_type] = frozenset(clusters)

    if missing:
        print(f"  WARNING: {len(missing)} proteins not found in cluster file "
              f"(no /translation?): {missing[:5]}")
    return locus_gene_sets


# ── Step 5: Compare gene-content sets ─────────────────────────────────────────
def compare_gene_sets(
    locus_gene_sets: dict[str, frozenset],
    out_dir: Path,
) -> None:
    """
    Compare gene-content sets across all loci. Report:
      - Identical loci (same gene set → should not exist in a valid DB)
      - Subset loci (one locus gene set is a strict subset of another)
      - Summary counts
    """
    loci = list(locus_gene_sets.items())
    n = len(loci)

    # Group loci by gene-content set
    content_to_loci: dict[frozenset, list[str]] = defaultdict(list)
    for kl_type, gene_set in loci:
        content_to_loci[gene_set].append(kl_type)

    # Identical loci
    duplicates = {
        frozenset(v): v
        for v in content_to_loci.values()
        if len(v) > 1
    }

    # Subset relationships
    subset_pairs = []
    for i, (kl_a, set_a) in enumerate(loci):
        for kl_b, set_b in loci[i + 1:]:
            if set_a < set_b:
                subset_pairs.append((kl_a, kl_b, len(set_a), len(set_b)))
            elif set_b < set_a:
                subset_pairs.append((kl_b, kl_a, len(set_b), len(set_a)))

    # ── Write reports ──
    dup_path = out_dir / "duplicate_loci.tsv"
    with open(dup_path, "w") as fh:
        fh.write("loci_with_identical_gene_sets\tn_genes\n")
        for gene_set, loci_list in content_to_loci.items():
            if len(loci_list) > 1:
                fh.write(f"{','.join(sorted(loci_list))}\t{len(gene_set)}\n")

    sub_path = out_dir / "subset_loci.tsv"
    with open(sub_path, "w") as fh:
        fh.write("subset_locus\tsuperset_locus\tn_genes_subset\tn_genes_superset\n")
        for kl_a, kl_b, na, nb in sorted(subset_pairs):
            fh.write(f"{kl_a}\t{kl_b}\t{na}\t{nb}\n")

    summary_path = out_dir / "gene_content_summary.tsv"
    with open(summary_path, "w") as fh:
        fh.write("locus\tn_genes\n")
        for kl_type, gene_set in sorted(locus_gene_sets.items()):
            fh.write(f"{kl_type}\t{len(gene_set)}\n")

    # ── Print summary ──
    print("\n" + "=" * 60)
    print("GENE-CONTENT VALIDATION SUMMARY")
    print("=" * 60)
    print(f"  Total loci:              {n}")
    print(f"  Unique gene-content sets: {len(content_to_loci)}")
    print(f"  Duplicate loci (same set): "
          f"{sum(len(v) for v in duplicates.values())} "
          f"({len(duplicates)} groups)")
    print(f"  Subset relationships:     {len(subset_pairs)}")
    print()

    if duplicates:
        print("  ⚠ DUPLICATES (identical gene-content sets):")
        for loci_list in duplicates.values():
            print(f"    {', '.join(sorted(loci_list))}")
    else:
        print(f"  ✓ No duplicate loci — all {n} loci have unique gene-content sets")

    if subset_pairs:
        print(f"\n  ⚠ SUBSET RELATIONSHIPS ({len(subset_pairs)} pairs):")
        for kl_a, kl_b, na, nb in subset_pairs[:20]:
            print(f"    {kl_a} ({na} genes) ⊂ {kl_b} ({nb} genes)")
        if len(subset_pairs) > 20:
            print(f"    ... and {len(subset_pairs) - 20} more (see {sub_path})")
    else:
        print("  ✓ No strict subset relationships found")

    print()
    print(f"  Reports written to {out_dir}/")
    print(f"    {dup_path.name}")
    print(f"    {sub_path.name}")
    print(f"    {summary_path.name}")
    print("=" * 60)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gbk",   default=DEFAULT_GBK,   type=Path,
                        help="GenBank database file (default: %(default)s)")
    parser.add_argument("--cdhit", default=DEFAULT_CDHIT, type=Path,
                        help="Path to cd-hit binary (default: %(default)s)")
    parser.add_argument("--out",   default=DEFAULT_OUT,   type=Path,
                        help="Output directory (default: %(default)s)")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    print(f"\nValidating gene-content sets: {args.gbk.name}")
    print(f"CD-HIT:  {args.cdhit}")
    print(f"Output:  {args.out}\n")

    # Step 1
    print("Step 1: Extracting CDS protein sequences...")
    proteins_fasta = args.out / "all_proteins.faa"
    locus_proteins = extract_proteins(args.gbk, proteins_fasta)

    # Step 2
    print("\nStep 2: Clustering proteins with CD-HIT "
          f"({int(CDHIT_IDENTITY*100)}% aa identity, "
          f"{int(CDHIT_COVERAGE*100)}% bidirectional coverage)...")
    cdhit_prefix = args.out / "protein_clusters"
    clstr_path = run_cdhit(args.cdhit, proteins_fasta, cdhit_prefix)

    # Step 3
    print("\nStep 3: Parsing CD-HIT clusters...")
    protein_to_cluster = parse_cdhit_clusters(clstr_path)

    # Step 4
    print("\nStep 4: Building gene-content sets per locus...")
    locus_gene_sets = build_gene_content_sets(locus_proteins, protein_to_cluster)

    # Step 5
    print("\nStep 5: Comparing gene-content sets across all loci...")
    compare_gene_sets(locus_gene_sets, args.out)


if __name__ == "__main__":
    main()
