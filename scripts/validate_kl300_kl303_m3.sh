#!/bin/bash
# =============================================================================
# validate_kl300_kl303_m3.sh
#
# Extract top candidate assemblies for KL300 and KL303 from ATB archives
# on M3 and run Kaptive v3 normalised scoring to confirm self-typing.
#
# Prerequisites:
#   scp scripts/validate_kl300_kl303_m3.sh ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/scripts/
#   scp scripts/normalise_novel_scores.py   ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/scripts/
#
# Submit:
#   sbatch scripts/validate_kl300_kl303_m3.sh
# =============================================================================

#SBATCH --job-name=val_kl300_303
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/logs/val_kl300_303_%j.log
#SBATCH --error=/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/logs/val_kl300_303_%j.log

set -euo pipefail

module load minimap2/2.28

REPO_DIR="/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4"
DB="${REPO_DIR}/DB/EC-K-typing_group1and4_v0.8.gbk"
ARCHIVES_DIR="/scratch2/js66/atb_archives"
FASTA_DIR="${REPO_DIR}/DB/kl300_303_fastas"
SCORES_TSV="${REPO_DIR}/DB/kaptive_scores_kl300_303.tsv"
MAPPING_FILE="${REPO_DIR}/DB/kl300_303_kl_map.tsv"
THREADS=8

mkdir -p "$FASTA_DIR" "${REPO_DIR}/logs"

export PATH="$HOME/.local/bin:$PATH"

echo "=================================================="
echo " KL300/KL303 candidate validation"
echo " $(date)"
echo "=================================================="

# ── Step 1: Ensure Kaptive v3 and dependencies ───────────────────────────────
echo "[1] Checking Kaptive and Python dependencies..."
if ! python3 -c "import kaptive; from kaptive._version import __version__; v=__version__; assert int(v.split('.')[0])>=3" 2>/dev/null; then
    pip install --user --quiet "kaptive>=3.0"
fi
pip install --user --quiet pandas
kaptive --version 2>&1 | head -1

# ── Step 2: Download target FASTAs from ENA/NCBI ─────────────────────────────
echo "[2] Downloading candidate FASTAs from ENA/NCBI..."

python3 << 'PYEOF'
import subprocess
import urllib.request
import urllib.error
import json
from pathlib import Path

FASTA_DIR    = Path("/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/DB/kl300_303_fastas")
MAPPING_FILE = Path("/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/DB/kl300_303_kl_map.tsv")
FASTA_DIR.mkdir(parents=True, exist_ok=True)

# Top 3 candidates per locus (99.8% qcov, 100% pident from LexicMap)
targets = {
    "SAMD00053151":   "KL300",
    "SAMEA103908486": "KL300",
    "SAMN12532904":   "KL300",
    "SAMEA6656333":   "KL303",
    "SAMEA6656328":   "KL303",
    "SAMN20959615":   "KL303",
}
print(f"    Target accessions: {len(targets)}")


def ena_get_assembly_ftp(biosample_acc):
    """Use ENA portal API to get FTP path for assembly FASTA from a BioSample."""
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={biosample_acc}&result=assembly"
           f"&fields=submitted_ftp,ftp&format=json")
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read())
            if data:
                ftp = data[0].get("submitted_ftp") or data[0].get("ftp", "")
                # May be semicolon-separated; pick the .fasta.gz or .fa.gz file
                for f in ftp.split(";"):
                    if f.endswith((".fasta.gz", ".fa.gz", ".fasta", ".fa")):
                        return f.strip()
    except Exception:
        pass
    return None


def ncbi_get_assembly_url(biosample_acc):
    """Use NCBI datasets to get FTP path for a BioSample."""
    try:
        result = subprocess.run(
            ["datasets", "summary", "genome", "accession", biosample_acc,
             "--as-json-lines"],
            capture_output=True, text=True, timeout=30
        )
        for line in result.stdout.splitlines():
            d = json.loads(line)
            ftp = (d.get("assembly_info", {})
                    .get("assembly_stats", {})
                    .get("ftp_path", ""))
            if ftp:
                name = ftp.rstrip("/").split("/")[-1]
                return f"{ftp}/{name}_genomic.fna.gz"
    except Exception:
        pass
    return None


def download_fasta(url, dest_fa):
    """Download (possibly gzipped) FASTA from HTTP/FTP URL, write as plain FA."""
    import gzip, shutil
    http_url = url.replace("ftp://", "http://")
    tmp = dest_fa.with_suffix(".tmp.gz" if url.endswith(".gz") else ".tmp.fa")
    print(f"      Downloading {http_url} ...")
    try:
        urllib.request.urlretrieve(http_url, str(tmp))
    except Exception as e:
        print(f"      ERROR: {e}")
        return False
    if url.endswith(".gz"):
        with gzip.open(str(tmp), 'rb') as fin, open(dest_fa, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        tmp.unlink()
    else:
        tmp.rename(dest_fa)
    return dest_fa.exists() and dest_fa.stat().st_size > 0


already = {p.stem: p for p in FASTA_DIR.glob("*.fa")}
downloaded = 0
failed = []

for acc, kl in targets.items():
    dest = FASTA_DIR / f"{acc}.fa"
    if acc in already:
        print(f"    {acc}: already exists, skipping")
        downloaded += 1
        continue

    print(f"    {acc} ({kl}):")

    # Try ENA first (covers SAMEA* and SAMD*)
    url = ena_get_assembly_ftp(acc)
    if url:
        print(f"      ENA URL found: {url}")
        if download_fasta(url, dest):
            print(f"      Downloaded: {dest.stat().st_size // 1024} KB")
            downloaded += 1
            continue
        else:
            print(f"      ENA download failed, trying NCBI...")

    # Try NCBI datasets (covers SAMN*)
    url = ncbi_get_assembly_url(acc)
    if url:
        print(f"      NCBI URL found: {url}")
        if download_fasta(url, dest):
            print(f"      Downloaded: {dest.stat().st_size // 1024} KB")
            downloaded += 1
            continue

    # Try ATB S3 bucket as last resort
    s3_url = f"https://allthebacteria-assemblies.s3.amazonaws.com/{acc}.fa.gz"
    print(f"      Trying ATB S3: {s3_url}")
    if download_fasta(s3_url, dest):
        print(f"      Downloaded from S3: {dest.stat().st_size // 1024} KB")
        downloaded += 1
        continue

    print(f"      FAILED: could not download {acc}")
    failed.append(acc)

print(f"\n    Downloaded: {downloaded}/{len(targets)}")
if failed:
    print(f"    Failed: {failed}")

# Write KL mapping for successfully downloaded assemblies
with open(MAPPING_FILE, 'w') as f:
    f.write("accession\tKL_type\tfasta_path\n")
    for acc, kl in targets.items():
        fa = FASTA_DIR / f"{acc}.fa"
        if fa.exists() and fa.stat().st_size > 0:
            f.write(f"{acc}\t{kl}\t{fa}\n")

n = sum(1 for line in open(MAPPING_FILE)) - 1
print(f"    Mapping written: {n} entries → {MAPPING_FILE}")
PYEOF

echo "    Extraction complete."
ls -lh "$FASTA_DIR"

# ── Step 3: Run Kaptive --scores ──────────────────────────────────────────────
echo "[3] Running Kaptive --scores on $(ls $FASTA_DIR/*.fa | wc -l) assemblies..."

# Remove any stale scores file to prevent appending
rm -f "$SCORES_TSV"

kaptive assembly "$DB" ${FASTA_DIR}/*.fa \
    --scores "$SCORES_TSV" \
    -t "$THREADS" \
    2>&1 | tail -5

echo "    Scores written: $SCORES_TSV"

# ── Step 4: Normalised scoring ────────────────────────────────────────────────
echo "[4] Running normalised scoring..."

python3 "${REPO_DIR}/scripts/normalise_novel_scores.py" \
    --scores  "$SCORES_TSV" \
    --db      "$DB" \
    --mapping "$MAPPING_FILE" \
    --out     "${REPO_DIR}/DB/kl300_303_validation_summary.tsv"

echo ""
echo "=================================================="
echo " Validation complete — $(date)"
echo " Results: ${REPO_DIR}/DB/kl300_303_validation_summary.tsv"
echo "=================================================="
