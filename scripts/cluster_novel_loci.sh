#!/bin/bash
# =============================================================================
# cluster_novel_loci.sh
#
# Cluster novel G1/G4 loci extracted from AllTheBacteria to identify new KL types.
#
# Steps:
#   1. Combine all novel loci FASTAs into one file
#   2. Cluster with MMseqs2 at 95% nucleotide identity, 80% coverage
#   3. Extract representative sequence per cluster
#   4. Assign new KL type numbers (KL344+)
#   5. Write summary TSV
#
# Run interactively or as a single SLURM job (not an array).
# Submit as:
#   sbatch scripts/cluster_novel_loci.sh
# =============================================================================

#SBATCH --job-name=cluster_novel
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=/home/ebenezef/js66_scratch/ebenn/atb_screen/logs/cluster_%j.log
#SBATCH --error=/home/ebenezef/js66_scratch/ebenn/atb_screen/logs/cluster_%j.log

set -eo pipefail

module load mmseqs2/17-b804f-cuda12
module load blast/2.15.0

NOVEL_DIR="/home/ebenezef/js66_scratch/ebenn/atb_screen/novel_loci"
OUT_DIR="/home/ebenezef/js66_scratch/ebenn/atb_screen/clustering"
REF_FASTA="/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/DB/reference_loci_v0.6.fasta"
EXISTING_KL_COUNT=93   # KL300-KL343 + K24 + K96; new types start at KL344

PIDENT=0.95   # 95% nucleotide identity
QCOV=0.80     # 80% query coverage
THREADS=16

mkdir -p "$OUT_DIR"

echo "=========================================="
echo " Novel G1/G4 loci clustering"
echo " $(date)"
echo "=========================================="

# ── Step 1: Combine all novel loci ──────────────────────────────────────────
echo "[1] Combining novel loci FASTAs..."
COMBINED="${OUT_DIR}/all_novel_loci.fasta"
cat "${NOVEL_DIR}"/*.fasta > "$COMBINED"
N_TOTAL=$(grep -c "^>" "$COMBINED" || echo 0)
echo "    Total sequences: ${N_TOTAL}"

if [ "$N_TOTAL" -eq 0 ]; then
    echo "No novel loci found. Exiting."
    exit 0
fi

# ── Step 2: Cluster with MMseqs2 ─────────────────────────────────────────────
echo "[2] Clustering with MMseqs2 (${PIDENT} id, ${QCOV} cov, ${THREADS} threads)..."
MMSEQS_DB="${OUT_DIR}/mmseqs_db"
MMSEQS_CLUDB="${OUT_DIR}/mmseqs_clu"
MMSEQS_TMP="${OUT_DIR}/mmseqs_tmp"

mkdir -p "$MMSEQS_TMP"

# Create MMseqs2 database
mmseqs createdb "$COMBINED" "$MMSEQS_DB"

# Cluster: --min-seq-id = identity threshold, -c = coverage, --cov-mode 0 = bidirectional
mmseqs cluster \
    "$MMSEQS_DB" \
    "$MMSEQS_CLUDB" \
    "$MMSEQS_TMP" \
    --min-seq-id "$PIDENT" \
    -c "$QCOV" \
    --cov-mode 0 \
    --cluster-mode 0 \
    --threads "$THREADS" \
    -v 2

# Extract cluster representatives
REPS_FA="${OUT_DIR}/cluster_representatives.fasta"
mmseqs result2repseq "$MMSEQS_DB" "$MMSEQS_CLUDB" "${OUT_DIR}/mmseqs_reps"
mmseqs convert2fasta "${OUT_DIR}/mmseqs_reps" "$REPS_FA"

N_CLUSTERS=$(grep -c "^>" "$REPS_FA" || echo 0)
echo "    Clusters found: ${N_CLUSTERS}"

# Extract cluster membership TSV
CLUSTER_TSV="${OUT_DIR}/cluster_membership.tsv"
mmseqs createtsv "$MMSEQS_DB" "$MMSEQS_DB" "$MMSEQS_CLUDB" "$CLUSTER_TSV"
echo "    Cluster membership written to: ${CLUSTER_TSV}"

# ── Step 3: Assign KL type numbers ───────────────────────────────────────────
echo "[3] Assigning KL type numbers starting at KL$((300 + EXISTING_KL_COUNT + 1))..."
python3 << 'PYEOF'
import os, re
from Bio import SeqIO

OUT_DIR      = "/home/ebenezef/js66_scratch/ebenn/atb_screen/clustering"
REPS_FA      = f"{OUT_DIR}/cluster_representatives.fasta"
CLU_TSV      = f"{OUT_DIR}/cluster_membership.tsv"
EXISTING     = 93    # existing KL types in DB
KL_START     = 344   # KL344 is first new type (KL300-KL343 = 44 types, + K24 + K96 = 46...
                     # but new numbering continues from last assigned)

# Read cluster sizes
from collections import defaultdict
clu_sizes = defaultdict(int)
with open(CLU_TSV) as f:
    for line in f:
        rep, member = line.strip().split('\t')
        clu_sizes[rep] += 1

# Read representatives and assign KL numbers
records = list(SeqIO.parse(REPS_FA, 'fasta'))

# Sort by cluster size (largest first) — biggest clusters get lowest KL numbers
rep_ids = [r.id for r in records]
rep_ids_sorted = sorted(rep_ids, key=lambda x: -clu_sizes.get(x, 0))

kl_map = {}
for i, rep_id in enumerate(rep_ids_sorted):
    kl_num = KL_START + i
    kl_map[rep_id] = f"KL{kl_num}"

# Write renamed representatives FASTA
out_fa = f"{OUT_DIR}/novel_kl_types.fasta"
out_tsv = f"{OUT_DIR}/novel_kl_summary.tsv"

with open(out_fa, 'w') as fa, open(out_tsv, 'w') as tsv:
    tsv.write("KL_type\trep_genome\tseq_len\tcluster_size\n")
    for rec in records:
        kl = kl_map.get(rec.id, "unknown")
        size = clu_sizes.get(rec.id, 1)
        fa.write(f">{kl} {rec.id} {len(rec.seq)}bp cluster_size={size}\n{str(rec.seq)}\n")
        tsv.write(f"{kl}\t{rec.id}\t{len(rec.seq)}\t{size}\n")

print(f"Novel KL types assigned: {len(records)}")
print(f"KL range: KL{KL_START} – KL{KL_START + len(records) - 1}")
print(f"Written: {out_fa}")
print(f"Written: {out_tsv}")
PYEOF

# ── Step 4: Summary ───────────────────────────────────────────────────────────
echo ""
echo "=========================================="
echo " Clustering complete"
echo " $(date)"
echo "=========================================="
echo "  Input sequences : ${N_TOTAL}"
echo "  Novel KL clusters: ${N_CLUSTERS}"
echo "  Representatives : ${OUT_DIR}/novel_kl_types.fasta"
echo "  Summary TSV     : ${OUT_DIR}/novel_kl_summary.tsv"
echo "  Cluster membership: ${OUT_DIR}/cluster_membership.tsv"
echo ""
echo "Next step: inspect novel_kl_summary.tsv, then"
echo "  incorporate into G1/G4 database rebuild."
