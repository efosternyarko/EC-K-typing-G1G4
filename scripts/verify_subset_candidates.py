#!/usr/bin/env python3
"""
verify_subset_candidates.py — Download and verify top ATB candidate assemblies
for high-priority subset locus replacement (KL601, KL713, KL742).

For each candidate:
  1. Fetches assembly FASTA via ENA / NCBI API
  2. Extracts the locus region (LexicMap sstart..send) + a 10 kb downstream
     buffer to capture any genes that extend beyond the current truncated query
  3. Writes extracted FASTA to DB/subset_candidates/<KL>/<accession>.fasta
  4. Annotates with Prokka (if on PATH) and counts CDS
  5. Checks extracted region for 100-N spacers (KL713 only)
  6. Prints a verification summary table

Usage:
    python scripts/verify_subset_candidates.py

Requires: requests, Biopython
Optional: prokka on PATH (falls back to Biopython ORF count if absent)

Targets (top candidates from subset_candidates.tsv):
  KL601  SAMN12943639  contig00001  sstart=887872  send=900893
  KL713  SAMN07152258  contig00009  sstart=106678  send=107539
  KL742  SAMEA115390294 contig00001 sstart=2052256 send=2066952
"""

import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

REPO_ROOT   = Path(__file__).resolve().parents[1]
OUT_DIR     = REPO_ROOT / "DB" / "subset_candidates"
QUERY_DIR   = REPO_ROOT / "DB" / "subset_queries"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Buffer downstream of the match end (to capture extra genes beyond truncated query)
DOWNSTREAM_BUFFER = 10_000
# Buffer upstream (to be safe with annotation context)
UPSTREAM_BUFFER   = 2_000

# Top candidates identified from subset_candidates.tsv
# Coordinates (sstart/send) are from the ATB assembly — not used directly.
# The NCBI assembly for the same BioSample differs; locus is located via blastn.
TARGETS = [
    dict(
        locus   = "KL601",
        accession = "SAMN12943639",
        contig  = "SAMN12943639.contig00001",
        sstart  = 887872,
        send    = 900893,
        note    = "Single HSP covers full KL601 query from pos 1; verify 3 missing 3'-end genes",
        current_cds = 17,
    ),
    dict(
        locus   = "KL713",
        accession = "SAMN07152258",
        contig  = "SAMN07152258.contig00009",
        sstart  = 106678,
        send    = 107539,
        note    = "100% identity on single contig; verify no 100-N spacer in full locus",
        current_cds = 15,
    ),
    dict(
        locus   = "KL742",
        accession = "SAMN41228912",
        contig  = "SAMN41228912.contig00001",
        sstart  = 2048072,
        send    = 2062768,
        note    = "Full qcov from pos 1 (rank 2); verify 8 missing rhamnose-synthesis genes",
        current_cds = 15,
    ),
]


# ---------------------------------------------------------------------------
# Assembly download helpers
# ---------------------------------------------------------------------------

def fetch_assembly_url_ena(biosample: str) -> "str | None":
    """Query ENA portal for the submitted FASTA FTP URL for a BioSample."""
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={biosample}&result=assembly"
        "&fields=submitted_ftp,ftp&format=tsv"
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        lines = r.text.strip().splitlines()
        if len(lines) < 2:
            return None
        for line in lines[1:]:
            parts = line.split("\t")
            for part in parts:
                # Take first .fasta.gz or .fa.gz URL
                for ftp_url in part.split(";"):
                    ftp_url = ftp_url.strip()
                    if ftp_url.endswith((".fasta.gz", ".fa.gz", ".fna.gz")):
                        return ftp_url.replace("ftp://", "https://")
    except Exception as e:
        print(f"    ENA API error: {e}")
    return None


def fetch_assembly_url_ncbi(biosample: str) -> "str | None":
    """Query NCBI Entrez to get the FTP path for the genomic FASTA."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    try:
        # Step 1: accession string → numeric BioSample UID
        r0 = requests.get(
            f"{base}/esearch.fcgi",
            params={"db": "biosample", "term": biosample, "retmode": "json"},
            timeout=30,
        )
        r0.raise_for_status()
        uid_list = r0.json()["esearchresult"]["idlist"]
        if not uid_list:
            return None
        bs_uid = uid_list[0]
        time.sleep(0.4)

        # Step 2: BioSample UID → Assembly IDs
        r = requests.get(
            f"{base}/elink.fcgi",
            params={"dbfrom": "biosample", "db": "assembly",
                    "id": bs_uid, "retmode": "json"},
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
        ids = []
        for linkset in data.get("linksets", []):
            for ld in linkset.get("linksetdbs", []):
                ids.extend(ld.get("links", []))
        if not ids:
            return None
        # Assembly ID → FTP path
        time.sleep(0.4)  # NCBI rate limit
        r2 = requests.get(
            f"{base}/esummary.fcgi",
            params={"db": "assembly", "id": ids[0], "retmode": "json"},
            timeout=30,
        )
        r2.raise_for_status()
        doc = r2.json()["result"][str(ids[0])]
        ftp = doc.get("ftppath_genbank", "") or doc.get("ftppath_refseq", "")
        if not ftp:
            return None
        asm_name = ftp.rstrip("/").split("/")[-1]
        url = f"{ftp}/{asm_name}_genomic.fna.gz"
        return url.replace("ftp://", "https://")
    except Exception as e:
        print(f"    NCBI API error: {e}")
    return None


def download_assembly(accession: str, out_path: Path) -> bool:
    """Try ENA then NCBI to download the assembly FASTA."""
    if out_path.exists():
        print(f"    Assembly already downloaded: {out_path.name}")
        return True

    print(f"    Fetching URL via ENA...", end=" ")
    url = fetch_assembly_url_ena(accession)
    if not url:
        print("not found.")
        print(f"    Fetching URL via NCBI...", end=" ")
        url = fetch_assembly_url_ncbi(accession)
    if not url:
        print("not found.")
        return False
    print(f"OK\n    Downloading: {url}")
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    tmp = out_path.with_suffix(".tmp.gz")
    with open(tmp, "wb") as fh:
        for chunk in r.iter_content(chunk_size=1 << 20):
            fh.write(chunk)
    tmp.rename(out_path)
    print(f"    Saved: {out_path.name} ({out_path.stat().st_size // 1024} KB)")
    return True


# ---------------------------------------------------------------------------
# Locus extraction via BLAST
# ---------------------------------------------------------------------------

def extract_locus_blast(fasta_gz: Path, locus: str,
                        upstream_buf: int, downstream_buf: int,
                        out_fasta: Path, work_dir: Path) -> "tuple[int, bool]":
    """
    Locate the K-locus in the assembly using blastn against the query FASTA,
    then extract [best_hit_start - upstream_buf .. best_hit_end + downstream_buf].
    Returns (extracted_length, has_n_spacer).
    """
    import gzip
    import csv

    query_fasta = QUERY_DIR / f"{locus}_query.fasta"
    if not query_fasta.exists():
        print(f"    ERROR: query FASTA not found: {query_fasta}")
        return 0, False

    # Decompress assembly to tmp file for BLAST
    asm_fasta = work_dir / f"{fasta_gz.stem}.fna"
    if not asm_fasta.exists():
        print(f"    Decompressing assembly...")
        import gzip as gz
        with gz.open(fasta_gz, "rt") as fin, open(asm_fasta, "w") as fout:
            for line in fin:
                fout.write(line)

    # Make BLAST db
    db_path = work_dir / "blastdb"
    if not (work_dir / "blastdb.nsq").exists() and not (work_dir / "blastdb.nsi").exists():
        subprocess.run(
            ["makeblastdb", "-in", str(asm_fasta), "-dbtype", "nucl",
             "-out", str(db_path), "-parse_seqids"],
            capture_output=True, check=True,
        )

    # Run blastn
    blast_out = work_dir / "blast_result.tsv"
    subprocess.run(
        ["blastn", "-query", str(query_fasta), "-db", str(db_path),
         "-out", str(blast_out), "-outfmt", "6 sseqid sstart send sstrand pident qcovs",
         "-perc_identity", "90", "-max_target_seqs", "1", "-num_threads", "4"],
        capture_output=True, check=True,
    )

    # Parse best hit
    hits = []
    with open(blast_out) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row:
                continue
            sseqid, sstart, send, sstrand, pident, qcovs = row
            hits.append((sseqid, int(sstart), int(send), sstrand,
                         float(pident), float(qcovs)))
    if not hits:
        print(f"    WARNING: no BLAST hit ≥90% identity found for {locus}")
        return 0, False

    # Take hit with highest qcovs
    hits.sort(key=lambda x: x[5], reverse=True)
    sseqid, sstart, send, sstrand, pident, qcovs = hits[0]
    print(f"    BLAST hit: {sseqid} {sstart}-{send} ({sstrand}) "
          f"pident={pident:.1f}% qcov={qcovs:.1f}%")

    # Load assembly and extract
    # BLAST with -parse_seqids returns 'gb|ACC|' style IDs; normalise to bare accession
    bare_id = sseqid.split("|")[1] if "|" in sseqid else sseqid
    records = {r.id: r for r in SeqIO.parse(asm_fasta, "fasta")}
    rec = records.get(bare_id) or records.get(sseqid)
    if rec is None:
        print(f"    WARNING: hit contig {sseqid} (bare: {bare_id}) not found")
        return 0, False

    # Normalise orientation (blastn sstart/send: sstart may > send on minus strand)
    lo = min(sstart, send) - 1      # 0-based
    hi = max(sstart, send)
    start = max(0, lo - upstream_buf)
    end   = min(len(rec.seq), hi + downstream_buf)
    region_seq = rec.seq[start:end]

    has_n_spacer = bool(re.search(r"N{20,}", str(region_seq)))

    extracted = SeqRecord(
        region_seq,
        id=f"{sseqid}_{start+1}_{end}",
        description=f"BLAST-extracted {locus} region from {fasta_gz.name}",
    )
    with open(out_fasta, "w") as fh:
        SeqIO.write(extracted, fh, "fasta")

    return len(region_seq), has_n_spacer


# ---------------------------------------------------------------------------
# CDS counting
# ---------------------------------------------------------------------------

def count_cds_pyrodigal(fasta: Path) -> "int | None":
    """Run pyrodigal (Prodigal Python port) and return CDS count."""
    try:
        import pyrodigal
    except ImportError:
        return None
    rec = next(SeqIO.parse(fasta, "fasta"))
    orf_finder = pyrodigal.GeneFinder(meta=True)   # metagenome mode: no training needed
    genes = orf_finder.find_genes(bytes(rec.seq))
    return len(genes)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    try:
        import pyrodigal
        have_pyrodigal = True
    except ImportError:
        have_pyrodigal = False
    print(f"pyrodigal available: {have_pyrodigal}")
    print()

    results = []

    for t in TARGETS:
        locus     = t["locus"]
        accession = t["accession"]
        print(f"{'='*60}")
        print(f"{locus}  |  {accession}  |  {t['note']}")
        print(f"{'='*60}")

        locus_dir = OUT_DIR / locus
        locus_dir.mkdir(exist_ok=True)

        # 1. Download
        fasta_gz = locus_dir / f"{accession}_genomic.fna.gz"
        ok = download_assembly(accession, fasta_gz)
        if not ok:
            print(f"  SKIPPED — could not download {accession}\n")
            results.append(dict(locus=locus, accession=accession,
                                status="download_failed"))
            continue

        # 2. Extract locus region via BLAST
        region_fasta = locus_dir / f"{accession}_{locus}_region.fasta"
        print(f"  Locating {locus} in assembly via BLAST "
              f"(+{UPSTREAM_BUFFER}/{DOWNSTREAM_BUFFER} bp buffer)...")
        region_len, has_spacer = extract_locus_blast(
            fasta_gz, locus,
            UPSTREAM_BUFFER, DOWNSTREAM_BUFFER,
            region_fasta, locus_dir,
        )
        if region_len == 0:
            results.append(dict(locus=locus, accession=accession,
                                status="extraction_failed"))
            continue
        print(f"  Extracted region: {region_len:,} bp | N-spacer: {has_spacer}")

        # 3. Count CDS
        if have_pyrodigal:
            cds = count_cds_pyrodigal(region_fasta)
            method = "pyrodigal (Prodigal metagenome mode)"
        else:
            cds = None
            method = "unavailable (install pyrodigal)"

        print(f"  CDS count ({method}): {cds}  |  current rep: {t['current_cds']}")
        if cds is not None and cds > t["current_cds"]:
            verdict = "MORE GENES — GOOD CANDIDATE"
        elif locus == "KL713" and not has_spacer:
            verdict = "NO N-SPACER — GOOD CANDIDATE"
        else:
            verdict = "no improvement detected"
        print(f"  Verdict: {verdict}")
        print()

        results.append(dict(
            locus=locus, accession=accession,
            region_len=region_len, has_n_spacer=has_spacer,
            cds=cds, cds_method=method,
            current_cds=t["current_cds"], verdict=verdict,
        ))

    # Summary table
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"{'Locus':<8} {'Accession':<18} {'Region bp':>10} {'CDS':>5} "
          f"{'Curr':>5} {'N-spacer':<10} Verdict")
    print("-"*70)
    for r in results:
        if "cds" not in r:
            print(f"{r['locus']:<8} {r['accession']:<18}  {r.get('status','?')}")
            continue
        print(f"{r['locus']:<8} {r['accession']:<18} "
              f"{r.get('region_len',0):>10,} {str(r.get('cds','?')):>5} "
              f"{r.get('current_cds','?'):>5} "
              f"{'YES' if r.get('has_n_spacer') else 'no':<10} "
              f"{r.get('verdict','')}")

    print(f"\nExtracted FASTAs in: {OUT_DIR}/")
    print("Next: inspect FASTAs, confirm gene content, then replace GenBank records.")


if __name__ == "__main__":
    main()
