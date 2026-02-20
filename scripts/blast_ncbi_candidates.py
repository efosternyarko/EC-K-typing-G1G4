#!/usr/bin/env python3
"""
blast_ncbi_candidates.py — Search NCBI for better KL300/303/306/307 representatives.

Submits each failing locus sequence to NCBI megaBLAST against the `nt` database
restricted to Escherichia coli (taxid:562). Collects assembly accessions from hits
with >=90% identity and >=70% query coverage, then ranks them as candidates for
replacement representatives.

Usage
-----
    python scripts/blast_ncbi_candidates.py [--loci KL300 KL303] [--max-hits 300]

Outputs
-------
    DB/blast_ncbi_results/<LOCUS>_blast_raw.xml   — raw NCBI XML
    DB/blast_ncbi_results/<LOCUS>_blast_hits.tsv  — parsed hits
    DB/blast_ncbi_results/candidates_summary.tsv  — ranked candidates across all loci
"""

import argparse
import re
import sys
import time
from pathlib import Path

import xml.etree.ElementTree as ET

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW

Entrez.email = "lshef4@lshtm.ac.uk"   # polite NCBI usage

REPO_DIR    = Path(__file__).resolve().parent.parent
DB_DIR      = REPO_DIR / "DB"
QUERY_FASTA = DB_DIR / "kl_failing_query.fasta"
OUT_DIR     = DB_DIR / "blast_ncbi_results"

# Thresholds for filtering hits
MIN_PIDENT    = 90.0
MIN_QCOV      = 70.0   # % of query length covered
MAX_HITS      = 300    # NCBI hits to request per query
E_VALUE       = 1e-50  # stringent — these are long sequences

TARGETS_DEFAULT = ["KL300", "KL303", "KL306", "KL307"]


def load_queries(fasta: Path, targets: list) -> dict:
    """Return {locus_name: SeqRecord} for each target."""
    queries = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        name = rec.id.split()[0]
        if name in targets:
            queries[name] = rec
    missing = [t for t in targets if t not in queries]
    if missing:
        print(f"WARNING: not found in {fasta.name}: {missing}", file=sys.stderr)
    return queries


def blast_locus(locus: str, seq: str, out_xml: Path, max_hits: int) -> None:
    """Submit one BLAST search and write raw XML."""
    print(f"  Submitting {locus} ({len(seq):,} bp) to NCBI megaBLAST...", flush=True)
    result = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=seq,
        hitlist_size=max_hits,
        entrez_query="Escherichia coli[Organism]",
        megablast=True,
        expect=E_VALUE,
        word_size=28,
    )
    xml_text = result.read()
    out_xml.write_bytes(xml_text if isinstance(xml_text, bytes) else xml_text.encode())
    print(f"    Saved raw XML: {out_xml.name}")


def parse_blast_xml(locus: str, xml_path: Path, qlen: int) -> list[dict]:
    """Parse BLAST XML directly with ElementTree — more robust than SearchIO."""
    rows = []
    try:
        tree = ET.parse(str(xml_path))
    except ET.ParseError as e:
        print(f"  ERROR parsing {xml_path.name}: {e}", file=sys.stderr)
        return rows

    root = tree.getroot()
    # Handle both bare BlastOutput and wrapped formats
    iterations = root.findall(".//Iteration")

    for iteration in iterations:
        for hit in iteration.findall(".//Hit"):
            hit_id   = hit.findtext("Hit_id",  "")
            hit_def  = hit.findtext("Hit_def", "")[:80]

            # Aggregate all HSPs for this hit
            covered_qranges = []   # list of (qstart, qend) to compute non-overlapping coverage
            best_pident     = 0.0

            for hsp in hit.findall(".//Hsp"):
                ident    = int(hsp.findtext("Hsp_identity",  "0"))
                aln_len  = int(hsp.findtext("Hsp_align-len", "1"))
                q_from   = int(hsp.findtext("Hsp_query-from","0"))
                q_to     = int(hsp.findtext("Hsp_query-to",  "0"))
                pident   = ident / aln_len * 100 if aln_len > 0 else 0
                best_pident = max(best_pident, pident)
                covered_qranges.append((min(q_from, q_to), max(q_from, q_to)))

            # Compute non-overlapping query coverage
            covered_qranges.sort()
            merged, total_covered = [], 0
            for start, end in covered_qranges:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append([start, end])
            total_covered = sum(e - s + 1 for s, e in merged)
            qcov = min(100.0, total_covered / qlen * 100)

            if best_pident < MIN_PIDENT or qcov < MIN_QCOV:
                continue

            # Extract clean accession (handles gi|xxx|gb|ACC| format)
            parts = hit_id.split("|")
            if len(parts) >= 4:
                acc = parts[3].split(".")[0]   # gb|ACC.1| → ACC
            elif len(parts) >= 2:
                acc = parts[1].split(".")[0]
            else:
                acc = hit_id.split(".")[0]

            rows.append({
                "locus":       locus,
                "accession":   acc,
                "hit_id":      hit_id,
                "description": hit_def,
                "pident":      round(best_pident, 2),
                "qcov":        round(qcov, 1),
                "n_hsps":      len(hit.findall(".//Hsp")),
            })

    return rows


def hits_to_tsv(rows: list[dict], out_tsv: Path) -> None:
    import csv
    if not rows:
        out_tsv.write_text("locus\taccession\thit_id\tdescription\tpident\tqcov\tn_hsps\n")
        return
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loci",     nargs="+", default=TARGETS_DEFAULT,
                        help="Loci to search (default: KL300 KL303 KL306 KL307)")
    parser.add_argument("--max-hits", type=int,  default=MAX_HITS)
    parser.add_argument("--skip-blast", action="store_true",
                        help="Re-use existing XML files (skip NCBI submission)")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("NCBI BLAST search for better locus representatives")
    print(f"Loci: {args.loci}")
    print(f"Database: nt  |  Organism: Escherichia coli  |  Megablast")
    print(f"Thresholds: pident>={MIN_PIDENT}%  qcov>={MIN_QCOV}%")
    print("=" * 60)

    queries = load_queries(QUERY_FASTA, args.loci)
    if not queries:
        print("No query sequences found. Check kl_failing_query.fasta", file=sys.stderr)
        sys.exit(1)

    all_hits = []

    for locus in args.loci:
        if locus not in queries:
            print(f"\n[{locus}] SKIP — not in query file")
            continue

        rec      = queries[locus]
        out_xml  = OUT_DIR / f"{locus}_blast_raw.xml"
        out_tsv  = OUT_DIR / f"{locus}_blast_hits.tsv"

        print(f"\n[{locus}] {len(rec.seq):,} bp")

        # BLAST (or re-use cached XML)
        if args.skip_blast and out_xml.exists():
            print(f"  Using cached XML: {out_xml.name}")
        else:
            blast_locus(locus, str(rec.seq), out_xml, args.max_hits)
            time.sleep(3)   # be polite to NCBI

        # Parse
        print(f"  Parsing results...")
        hits = parse_blast_xml(locus, out_xml, len(rec.seq))
        print(f"  Hits passing filters: {len(hits)}")

        # Sort: highest qcov then pident
        hits.sort(key=lambda r: (-r["qcov"], -r["pident"]))
        hits_to_tsv(hits, out_tsv)
        print(f"  Written: {out_tsv.name}")

        # Preview top 5
        if hits:
            print(f"  {'Accession':<15}  {'qcov':>6}  {'pident':>7}  Description")
            print(f"  {'-'*15}  {'------':>6}  {'-------':>7}  {'-'*30}")
            for h in hits[:5]:
                print(f"  {h['accession']:<15}  {h['qcov']:>6.1f}  "
                      f"{h['pident']:>7.2f}  {h['description'][:50]}")

        all_hits.extend(hits)

    # Write combined summary
    if all_hits:
        import csv
        summary = OUT_DIR / "candidates_summary.tsv"
        with open(summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=all_hits[0].keys(), delimiter="\t")
            w.writeheader()
            w.writerows(all_hits)
        print(f"\nCombined summary: {summary}  ({len(all_hits)} hits)")

    print("\n" + "=" * 60)
    print("Next steps:")
    print("  python scripts/download_candidates.py")
    print("  (downloads top candidate assemblies via NCBI Datasets)")
    print("=" * 60)


if __name__ == "__main__":
    main()
