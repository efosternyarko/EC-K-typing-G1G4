#!/usr/bin/env python3
"""
Investigate duplicate and subset loci flagged by validate_gene_content.py.

For each duplicate group (identical gene-content sets at 95% aa / 80% cov):
  - Classify each locus as BSI (KL300-KL391) or novel ATB (KL392+)
  - Compute pairwise nucleotide identity between locus sequences
  - Propose which locus to retain (BSI > lowest novel KL number)
  - Report unique genes that differ between subset pairs

Outputs:
  merging_plan.tsv     — proposed merging actions for duplicate groups
  subset_analysis.tsv  — unique genes per subset pair
  investigation_report.txt — human-readable summary for collaborator meeting

Usage:
    python3 investigate_duplicates.py [--gbk PATH] [--validation_dir DIR] [--out DIR]
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


# ── Defaults ──────────────────────────────────────────────────────────────────
REPO_ROOT      = Path(__file__).resolve().parents[1]
DEFAULT_GBK    = REPO_ROOT / "DB" / "EC-K-typing_group1and4_v0.9.gbk"
DEFAULT_VAL    = Path("/tmp/gene_content_validation")
DEFAULT_OUT    = Path("/tmp/duplicate_investigation")

BSI_RANGE      = range(300, 392)   # KL300-KL391 (BSI reference loci)


# ── Helpers ───────────────────────────────────────────────────────────────────
def kl_num(kl_type: str) -> int:
    m = re.search(r"(\d+)$", kl_type)
    return int(m.group(1)) if m else 9999


def is_bsi(kl_type: str) -> bool:
    return kl_num(kl_type) in BSI_RANGE


def pairwise_nt_identity(seq_a: str, seq_b: str) -> float:
    """Global pairwise nucleotide identity (%) between two sequences."""
    # For speed, use a simple character comparison on aligned seqs
    # For long sequences use sampled comparison (every 10th base)
    step = max(1, min(len(seq_a), len(seq_b)) // 2000)
    a = seq_a[::step]
    b = seq_b[::step]
    min_len = min(len(a), len(b))
    if min_len == 0:
        return 0.0
    matches = sum(x == y for x, y in zip(a[:min_len], b[:min_len]))
    return round(matches / min_len * 100, 1)


# ── Load data ─────────────────────────────────────────────────────────────────
def load_locus_sequences(gbk_path: Path) -> dict[str, str]:
    """Return dict: KL_type -> full nucleotide sequence string."""
    locus_seqs = {}
    for record in SeqIO.parse(gbk_path, "genbank"):
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
            kl_type = record.name.split("_")[0]
        locus_seqs[kl_type] = str(record.seq)
    return locus_seqs


def load_locus_gene_sets(
    clstr_path: Path,
    gbk_path: Path,
) -> dict[str, frozenset]:
    """Re-derive gene-content sets (cluster IDs) per locus from .clstr file."""
    # Parse cluster file
    protein_to_cluster: dict[str, int] = {}
    cluster_id = None
    with open(clstr_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1])
            elif line and cluster_id is not None:
                parts = line.split(">")
                if len(parts) < 2:
                    continue
                pid = parts[1].split("...")[0].strip()
                protein_to_cluster[pid] = cluster_id

    # Build locus -> protein_ids from GenBank (same logic as step 1)
    locus_proteins: dict[str, list] = defaultdict(list)
    cds_counter = 0
    for record in SeqIO.parse(gbk_path, "genbank"):
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
            kl_type = record.name.split("_")[0]

        for feat in record.features:
            if feat.type != "CDS":
                continue
            lt = feat.qualifiers.get("locus_tag", [f"cds_{cds_counter}"])[0]
            locus_proteins[kl_type].append(f"{kl_type}|{lt}")
            cds_counter += 1

    # Build gene-content sets
    locus_gene_sets: dict[str, frozenset] = {}
    for kl_type, proteins in locus_proteins.items():
        clusters = frozenset(
            protein_to_cluster[pid]
            for pid in proteins
            if pid in protein_to_cluster
        )
        locus_gene_sets[kl_type] = clusters
    return locus_gene_sets


def parse_duplicate_groups(dup_tsv: Path) -> list[list[str]]:
    groups = []
    with open(dup_tsv) as fh:
        next(fh)  # header
        for line in fh:
            loci_str, _ = line.strip().split("\t")
            groups.append(sorted(loci_str.split(","), key=kl_num))
    return groups


def parse_subset_pairs(sub_tsv: Path) -> list[tuple]:
    pairs = []
    with open(sub_tsv) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split("\t")
            pairs.append((parts[0], parts[1], int(parts[2]), int(parts[3])))
    return pairs


# ── Analysis ──────────────────────────────────────────────────────────────────
def analyse_duplicates(
    groups: list[list[str]],
    locus_seqs: dict[str, str],
    locus_gene_sets: dict[str, frozenset],
) -> list[dict]:
    results = []
    for group in groups:
        bsi   = [k for k in group if is_bsi(k)]
        novel = [k for k in group if not is_bsi(k)]

        # Proposed keep: BSI locus if present, else smallest KL number (= largest cluster)
        keep = bsi[0] if bsi else novel[0]
        merge_into = [k for k in group if k != keep]

        # Pairwise nt identities vs representative
        nt_ids = {}
        ref_seq = locus_seqs.get(keep, "")
        for other in merge_into:
            other_seq = locus_seqs.get(other, "")
            nt_ids[other] = pairwise_nt_identity(ref_seq, other_seq)

        results.append({
            "group":      ",".join(group),
            "n_loci":     len(group),
            "n_genes":    len(locus_gene_sets.get(keep, set())),
            "bsi_loci":   ",".join(bsi) if bsi else "none",
            "novel_loci": ",".join(novel) if novel else "none",
            "keep":       keep,
            "merge":      ",".join(merge_into),
            "nt_identity_vs_keep": "; ".join(
                f"{k}:{v}%" for k, v in nt_ids.items()
            ),
            "justification": (
                "BSI reference locus retained" if bsi
                else "Lowest KL number retained (largest ATB cluster)"
            ),
        })
    return results


def analyse_subsets(
    pairs: list[tuple],
    locus_gene_sets: dict[str, frozenset],
    locus_seqs: dict[str, str],
) -> list[dict]:
    results = []
    for sub, sup, n_sub, n_sup in pairs:
        sub_genes = locus_gene_sets.get(sub, frozenset())
        sup_genes = locus_gene_sets.get(sup, frozenset())
        unique_in_sup = sup_genes - sub_genes  # genes in superset not in subset
        nt_id = pairwise_nt_identity(
            locus_seqs.get(sub, ""), locus_seqs.get(sup, "")
        )
        results.append({
            "subset_locus":     sub,
            "superset_locus":   sup,
            "n_genes_subset":   n_sub,
            "n_genes_superset": n_sup,
            "n_unique_in_superset": len(unique_in_sup),
            "nt_identity":      f"{nt_id}%",
            "sub_is_bsi":       is_bsi(sub),
            "sup_is_bsi":       is_bsi(sup),
        })
    return results


# ── Output ────────────────────────────────────────────────────────────────────
def write_merging_plan(results: list[dict], out_path: Path) -> None:
    header = [
        "group", "n_loci", "n_genes", "bsi_loci", "novel_loci",
        "keep", "merge_into", "nt_identity_vs_keep", "justification",
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in results:
            fh.write("\t".join(str(r[k]) for k in [
                "group", "n_loci", "n_genes", "bsi_loci", "novel_loci",
                "keep", "merge", "nt_identity_vs_keep", "justification",
            ]) + "\n")


def write_subset_analysis(results: list[dict], out_path: Path) -> None:
    header = [
        "subset_locus", "superset_locus", "n_genes_subset",
        "n_genes_superset", "n_unique_in_superset",
        "nt_identity", "sub_is_bsi", "sup_is_bsi",
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in results:
            fh.write("\t".join(str(r[k]) for k in header) + "\n")


def write_report(
    dup_results: list[dict],
    sub_results: list[dict],
    out_path: Path,
) -> None:
    lines = []
    lines.append("=" * 70)
    lines.append("EC-K-TYPING G1/G4 DATABASE — GENE-CONTENT VALIDATION REPORT")
    lines.append("Kaptive principle: unique locus = unique gene set")
    lines.append("Clustering: CD-HIT 95% aa identity / 80% bidirectional coverage")
    lines.append("=" * 70)

    lines.append(f"\nTotal loci in database: 651")
    lines.append(f"Unique gene-content sets: {651 - sum(r['n_loci'] - 1 for r in dup_results)}")
    lines.append(f"\nDUPLICATE GROUPS (identical gene-content sets): {len(dup_results)}")
    lines.append(f"  Loci to merge: {sum(len(r['merge'].split(',')) for r in dup_results)}")
    lines.append(f"\nSUBSET RELATIONSHIPS: {len(sub_results)} pairs")

    lines.append("\n" + "-" * 70)
    lines.append("PROPOSED MERGING PLAN")
    lines.append("-" * 70)
    for r in dup_results:
        lines.append(f"\n  Group: {r['group']}")
        lines.append(f"    Gene count:      {r['n_genes']}")
        lines.append(f"    BSI loci:        {r['bsi_loci']}")
        lines.append(f"    Novel loci:      {r['novel_loci']}")
        lines.append(f"    → KEEP:          {r['keep']}")
        lines.append(f"    → MERGE/REMOVE:  {r['merge']}")
        lines.append(f"    nt identity vs keep: {r['nt_identity_vs_keep']}")
        lines.append(f"    Justification:   {r['justification']}")

    lines.append("\n" + "-" * 70)
    lines.append("SUBSET RELATIONSHIPS SUMMARY")
    lines.append("(subset locus gene set ⊂ superset locus gene set)")
    lines.append("-" * 70)
    lines.append(
        f"{'Subset':<12} {'Superset':<12} {'Sub genes':>9} "
        f"{'Sup genes':>9} {'Unique in sup':>13} {'nt id':>7}"
    )
    for r in sub_results:
        bsi_flag = (
            "[BSI⊂BSI]" if r["sub_is_bsi"] and r["sup_is_bsi"]
            else "[BSI⊂nov]" if r["sub_is_bsi"]
            else "[nov⊂BSI]" if r["sup_is_bsi"]
            else ""
        )
        lines.append(
            f"  {r['subset_locus']:<10} {r['superset_locus']:<12} "
            f"{r['n_genes_subset']:>9} {r['n_genes_superset']:>9} "
            f"{r['n_unique_in_superset']:>13} {r['nt_identity']:>7}  {bsi_flag}"
        )

    lines.append("\n" + "=" * 70)
    lines.append("NOTE ON SUBSET RELATIONSHIPS")
    lines.append("=" * 70)
    lines.append(
        "Subset loci are not duplicates but may represent:\n"
        "  (a) Genuine biological variants (e.g. locus lacking an accessory gene)\n"
        "  (b) Incomplete locus capture during ATB screening (truncated sequence)\n"
        "  (c) Structural variants where flanking genes were missed\n"
        "Recommendation: investigate high-nt-identity subset pairs as potential\n"
        "incomplete captures; retain as distinct loci where biological basis is clear."
    )

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n".join(lines))


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gbk",            default=DEFAULT_GBK, type=Path)
    parser.add_argument("--validation_dir", default=DEFAULT_VAL, type=Path)
    parser.add_argument("--out",            default=DEFAULT_OUT, type=Path)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    clstr_path = args.validation_dir / "protein_clusters.clstr"
    dup_tsv    = args.validation_dir / "duplicate_loci.tsv"
    sub_tsv    = args.validation_dir / "subset_loci.tsv"

    print("Loading locus sequences...")
    locus_seqs = load_locus_sequences(args.gbk)
    print(f"  {len(locus_seqs)} loci loaded")

    print("Loading gene-content sets...")
    locus_gene_sets = load_locus_gene_sets(clstr_path, args.gbk)
    print(f"  {len(locus_gene_sets)} loci with gene sets")

    print("Parsing duplicate groups...")
    dup_groups = parse_duplicate_groups(dup_tsv)

    print("Parsing subset pairs...")
    sub_pairs = parse_subset_pairs(sub_tsv)

    print("Analysing duplicates...")
    dup_results = analyse_duplicates(dup_groups, locus_seqs, locus_gene_sets)

    print("Analysing subset pairs...")
    sub_results = analyse_subsets(sub_pairs, locus_gene_sets, locus_seqs)

    write_merging_plan(dup_results, args.out / "merging_plan.tsv")
    write_subset_analysis(sub_results, args.out / "subset_analysis.tsv")
    write_report(dup_results, sub_results, args.out / "investigation_report.txt")

    print(f"\nOutputs written to {args.out}/")


if __name__ == "__main__":
    main()
