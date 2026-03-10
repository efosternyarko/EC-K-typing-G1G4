#!/usr/bin/env python3
"""
download_novel_rep_fastas.py

Download representative FASTA assemblies for all novel ATB loci (KL392–KL968)
listed in DB/novel_kl_summary.tsv.

For each accession, tries in order:
  1. ATB S3 public bucket (fastest for ATB assemblies)
  2. ENA portal API
  3. NCBI datasets CLI

Downloads to DB/novel_rep_fastas/{accession}.fa
Skips accessions already present.

Usage:
    python3 scripts/download_novel_rep_fastas.py [--threads N]
"""

import argparse
import gzip
import json
import shutil
import subprocess
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO_DIR   = Path(__file__).resolve().parent.parent
SUMMARY    = REPO_DIR / "DB" / "novel_kl_summary.tsv"
FASTA_DIR  = REPO_DIR / "DB" / "novel_rep_fastas"


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

def _download_url(url: str, dest: Path) -> bool:
    """Download url (plain or .gz) to dest as plain FASTA. Returns True on success."""
    tmp = dest.with_suffix(".tmp.gz" if url.endswith(".gz") else ".tmp.fa")
    http_url = url.replace("ftp://", "http://")
    try:
        urllib.request.urlretrieve(http_url, str(tmp))
    except Exception:
        tmp.unlink(missing_ok=True)
        return False
    if url.endswith(".gz"):
        try:
            with gzip.open(str(tmp), "rb") as fin, open(dest, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            tmp.unlink()
        except Exception:
            tmp.unlink(missing_ok=True)
            dest.unlink(missing_ok=True)
            return False
    else:
        tmp.rename(dest)
    return dest.exists() and dest.stat().st_size > 0


def try_s3(acc: str, dest: Path) -> bool:
    url = f"https://allthebacteria-assemblies.s3.amazonaws.com/{acc}.fa.gz"
    return _download_url(url, dest)


def try_ena(acc: str, dest: Path) -> bool:
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={acc}&result=assembly&fields=submitted_ftp,ftp&format=json"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read())
            if data:
                ftp = data[0].get("submitted_ftp") or data[0].get("ftp", "")
                for f in ftp.split(";"):
                    f = f.strip()
                    if f.endswith((".fasta.gz", ".fa.gz", ".fasta", ".fa")):
                        return _download_url(f, dest)
    except Exception:
        pass
    return False


def try_ncbi(acc: str, dest: Path) -> bool:
    try:
        result = subprocess.run(
            ["datasets", "summary", "genome", "accession", acc, "--as-json-lines"],
            capture_output=True, text=True, timeout=30
        )
        for line in result.stdout.splitlines():
            d = json.loads(line)
            ftp = d.get("assembly_info", {}).get("assembly_stats", {}).get("ftp_path", "")
            if ftp:
                name = ftp.rstrip("/").split("/")[-1]
                return _download_url(f"{ftp}/{name}_genomic.fna.gz", dest)
    except Exception:
        pass
    return False


def download_one(acc: str, kl: str, dest: Path) -> tuple:
    """Try all sources; return (acc, kl, success, source)."""
    if dest.exists() and dest.stat().st_size > 0:
        return acc, kl, True, "cached"
    for fn, label in [(try_s3, "S3"), (try_ena, "ENA"), (try_ncbi, "NCBI")]:
        if fn(acc, dest):
            return acc, kl, True, label
    return acc, kl, False, "all_failed"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=8,
                        help="Parallel download threads (default: 8)")
    args = parser.parse_args()

    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    # Load targets
    targets = []
    with open(SUMMARY) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                kl, acc = parts[0], parts[1]
                targets.append((acc, kl))

    print(f"Targets:  {len(targets)} novel loci")
    cached = sum(1 for acc, _ in targets
                 if (FASTA_DIR / f"{acc}.fa").exists()
                 and (FASTA_DIR / f"{acc}.fa").stat().st_size > 0)
    print(f"Cached:   {cached}")
    print(f"To fetch: {len(targets) - cached}")
    print(f"Threads:  {args.threads}")
    print()

    results = {"S3": 0, "ENA": 0, "NCBI": 0, "cached": 0, "all_failed": 0}
    failed = []

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(download_one, acc, kl, FASTA_DIR / f"{acc}.fa"): (acc, kl)
            for acc, kl in targets
        }
        done = 0
        for fut in as_completed(futures):
            acc, kl, ok, source = fut.result()
            done += 1
            results[source] += 1
            if not ok:
                failed.append((acc, kl))
            if done % 50 == 0 or done == len(targets):
                pct = 100 * done / len(targets)
                print(f"  [{done}/{len(targets)}  {pct:.0f}%]  "
                      f"S3:{results['S3']}  ENA:{results['ENA']}  "
                      f"NCBI:{results['NCBI']}  cached:{results['cached']}  "
                      f"failed:{results['all_failed']}")

    n_ok = len(targets) - results["all_failed"]
    print(f"\nDownloaded: {n_ok}/{len(targets)}")
    if failed:
        print(f"Failed ({len(failed)}):")
        for acc, kl in failed:
            print(f"  {kl}  {acc}")

    # Write a mapping file for use with type_normalized.py
    mapping_out = REPO_DIR / "DB" / "novel_rep_kl_map.tsv"
    with open(mapping_out, "w") as fh:
        fh.write("KL\tsource_assembly\tlength_bp\tgroup\n")
        with open(SUMMARY) as sf:
            next(sf)
            for line in sf:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    kl, acc, seq_len = parts[0], parts[1], parts[2]
                    fa = FASTA_DIR / f"{acc}.fa"
                    if fa.exists() and fa.stat().st_size > 0:
                        fh.write(f"{kl}\t{acc}\t{seq_len}\tG1/G4\n")
    n_mapped = sum(1 for _ in open(mapping_out)) - 1
    print(f"\nMapping written: {n_mapped} entries → {mapping_out}")


if __name__ == "__main__":
    main()
