#!/bin/bash
# =============================================================================
# validate_novel_loci_m3.sh
#
# Self-typing validation for the 574 novel ATB loci (KL392-KL968) in v0.8.
#
# Steps:
#   1. Build accession→batch map from ATB screen results
#   2. Extract representative FASTAs from batch archives
#   3. Run Kaptive v3 --scores against v0.8 GenBank
#   4. Normalised scoring: check each assembly types to its correct KL
#
# Prerequisites:
#   scp DB/EC-K-typing_group1and4_v0.8.gbk ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/DB/
#   scp atb_clustering/novel_kl_summary.tsv  ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/DB/
#
# Submit:
#   sbatch scripts/validate_novel_loci_m3.sh
# =============================================================================

#SBATCH --job-name=val_novel
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/logs/val_novel_%j.log
#SBATCH --error=/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/logs/val_novel_%j.log

set -euo pipefail

module load minimap2/2.28

REPO_DIR="/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4"
DB="${REPO_DIR}/DB/EC-K-typing_group1and4_v0.8.gbk"
NOVEL_SUMMARY="${REPO_DIR}/DB/novel_kl_summary.tsv"
RESULTS_DIR="/home/ebenezef/js66_scratch/ebenn/atb_screen/results"
ARCHIVES_DIR="/scratch2/js66/atb_archives"
FASTA_DIR="${REPO_DIR}/DB/novel_rep_fastas"
SCORES_TSV="${REPO_DIR}/DB/kaptive_scores_novel_v0.8norm.tsv"
THREADS=16

mkdir -p "$FASTA_DIR"
mkdir -p "${REPO_DIR}/logs"

echo "=================================================="
echo " Novel loci self-typing validation — v0.8"
echo " $(date)"
echo "=================================================="

# ── Step 1: Install Kaptive v3 ───────────────────────────────────────────────
echo "[1] Checking Kaptive version..."
if python3 -c "import kaptive; from kaptive._version import __version__; v=__version__; assert int(v.split('.')[0])>=3" 2>/dev/null; then
    echo "    Kaptive v3+ already installed."
else
    echo "    Installing Kaptive v3..."
    pip install --user --quiet "kaptive>=3.0"
    echo "    Done."
fi
KAPTIVE=$(python3 -c "import os,sys; print(os.path.join(os.path.dirname(sys.executable),'..','bin','kaptive'))" 2>/dev/null || echo "kaptive")
# Try user bin if module kaptive is old
export PATH="$HOME/.local/bin:$PATH"
kaptive --version 2>&1 | head -1

# ── Step 2: Build accession→batch map and extract FASTAs ─────────────────────
echo "[2] Building accession→batch map and extracting FASTAs..."

python3 << 'PYEOF'
import subprocess
from pathlib import Path
from collections import defaultdict

CANDIDATES_FILE = Path("/home/ebenezef/js66_scratch/ebenn/atb_screen/batches_with_candidates.txt")
ARCHIVES_DIR    = Path("/scratch2/js66/atb_archives")
NOVEL_SUMMARY   = Path("/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/DB/novel_kl_summary.tsv")
FASTA_DIR       = Path("/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/DB/novel_rep_fastas")

# Load target accessions
targets = {}
with open(NOVEL_SUMMARY) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        kl, acc = parts[0], parts[1]
        targets[acc] = kl
print(f"    Target accessions: {len(targets)}")

# Read candidate batch numbers
with open(CANDIDATES_FILE) as f:
    candidate_batches = [line.strip() for line in f if line.strip()]

# Build accession→archive index by scanning archive contents
# (faster than hits files since representatives may not appear in hits)
print(f"    Scanning {len(candidate_batches)} candidate archives for target accessions...")
acc_to_archive = {}
remaining = set(targets.keys())

for batch_num in candidate_batches:
    if not remaining:
        break
    archive = ARCHIVES_DIR / f"atb.assembly.incr_release.202408.batch.{batch_num}.tar.xz"
    if not archive.exists():
        continue
    result = subprocess.run(
        ["tar", "tf", str(archive)],
        capture_output=True, text=True
    )
    for member in result.stdout.splitlines():
        acc = Path(member).stem  # e.g. SAMN30284532
        if acc in remaining:
            acc_to_archive[acc] = (archive, member)
            remaining.discard(acc)

print(f"    Found: {len(acc_to_archive)}/{len(targets)}")
if remaining:
    print(f"    Not found in any candidate batch: {len(remaining)} — {list(remaining)[:5]}")

# Extract FASTAs
already = {p.stem for p in FASTA_DIR.glob("*.fa")}
n_extracted = 0

# Group by archive for efficiency
archive_to_accs = defaultdict(list)
for acc, (archive, member) in acc_to_archive.items():
    if acc not in already:
        archive_to_accs[archive].append(member)

for archive, members in sorted(archive_to_accs.items()):
    cmd = ["tar", "xf", str(archive), "-C", str(FASTA_DIR),
           "--strip-components=1"] + members
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        print(f"    WARNING: tar error for {archive.name}: {result.stderr.decode()[:200]}")
    extracted = sum(1 for m in members if (FASTA_DIR / Path(m).name).exists())
    n_extracted += extracted
    print(f"    {archive.name}: {extracted}/{len(members)} extracted")

print(f"    Total FASTAs available: {n_extracted + len(already)} "
      f"({len(already)} pre-existing, {n_extracted} newly extracted)")

# Write accession→KL mapping
mapping_file = Path("/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/DB/novel_rep_kl_map.tsv")
with open(mapping_file, 'w') as f:
    f.write("accession\tKL_type\tfasta_path\n")
    for acc, kl in targets.items():
        fa = FASTA_DIR / f"{acc}.fa"
        if fa.exists():
            f.write(f"{acc}\t{kl}\t{fa}\n")
n_mapped = sum(1 for line in open(mapping_file)) - 1
print(f"    Mapping written: {n_mapped} entries → {mapping_file}")
PYEOF

echo "    Extraction complete."

# ── Step 3: Run Kaptive --scores ──────────────────────────────────────────────
echo "[3] Running Kaptive --scores on $(ls $FASTA_DIR/*.fa | wc -l) assemblies..."

kaptive assembly "$DB" ${FASTA_DIR}/*.fa \
    --scores "$SCORES_TSV" \
    -t "$THREADS" \
    2>&1 | tail -5

echo "    Scores written: $SCORES_TSV"

# ── Step 4: Normalised scoring ────────────────────────────────────────────────
echo "[4] Running normalised scoring..."

python3 "${REPO_DIR}/scripts/normalise_novel_scores.py" \
    --scores   "$SCORES_TSV" \
    --db       "$DB" \
    --mapping  "${REPO_DIR}/DB/novel_rep_kl_map.tsv" \
    --out      "${REPO_DIR}/DB/novel_validation_summary_v0.8norm.tsv"

echo ""
echo "=================================================="
echo " Validation complete — $(date)"
echo "=================================================="
