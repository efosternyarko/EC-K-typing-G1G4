#!/bin/bash
# 03_blastp_allvsall.sh
#
# BLASTp all-vs-all on all 12,735 CDS proteins from the v1.2 database.
# Used to build protein family presence/absence matrix.
#
# Thresholds (Kat Holt / Gladstone convention):
#   Identity  >= 80%
#   Coverage  >= 80% (query coverage used here; reciprocal enforced in step 04)
#
# Requires: blastp + makeblastdb from my_python.env
#
# Usage: bash scripts/03_blastp_allvsall.sh

set -eo pipefail

BLASTP="/Users/lshef4/lshef4_to_copy/anaconda3/envs/my_python.env/bin/blastp"
MAKEBLASTDB="/Users/lshef4/lshef4_to_copy/anaconda3/envs/my_python.env/bin/makeblastdb"

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
IS_DIR="$REPO_DIR/DB/is_analysis"
FAA="$IS_DIR/all_proteins.faa"
DB="$IS_DIR/all_proteins_db"
OUT="$IS_DIR/blastp_allvsall.tsv"

if [[ ! -f "$FAA" ]]; then
    echo "ERROR: $FAA not found. Run 01_extract_proteins.py first."
    exit 1
fi

echo "Building BLAST protein database..."
"$MAKEBLASTDB" -in "$FAA" -dbtype prot -out "$DB" -logfile /dev/null

echo "Running BLASTp all-vs-all (this may take 10-30 min for 12k proteins)..."
"$BLASTP" \
    -query "$FAA" \
    -db "$DB" \
    -outfmt "6 qseqid sseqid pident qcovs scovs length evalue bitscore" \
    -evalue 1e-5 \
    -num_threads 8 \
    -out "$OUT"

echo "Done. Results: $OUT"
wc -l "$OUT" | xargs echo "Total hits:"

# Clean up BLAST db files
for ext in .phr .pin .psq .pdb .pot .ptf .pto .pjs; do
    rm -f "${DB}${ext}"
done
