#!/bin/bash
# 02c_run_panaroo.sh
#
# Run Panaroo on all 623 GFF3 locus files to build the gene presence/absence matrix.
#
# Panaroo settings (following Gladstone et al. 2025 for E. coli K-loci):
#   --threshold 0.70        70% AA identity for gene family assignment
#   --family_threshold 0.70 70% family identity
#   --clean-mode strict     Remove genes not supported by graph topology
#   --remove-invalid-genes  Strip short/partial CDSs before graph construction
#   --aligner mafft         Use MAFFT for core genome alignment (needed for graph)
#
# Runtime: ~5-20 min for 623 loci on a laptop
#
# Usage: bash scripts/02c_run_panaroo.sh

set -eo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
GFF_DIR="$REPO_DIR/DB/panaroo_inputs"
OUT_DIR="$REPO_DIR/DB/panaroo_out"

if [[ ! -d "$GFF_DIR" ]] || [[ $(ls "$GFF_DIR"/*.gff3 2>/dev/null | wc -l) -eq 0 ]]; then
    echo "ERROR: No GFF3 files found in $GFF_DIR"
    echo "Run 02b_gbk_to_gff3.py first."
    exit 1
fi

N_GFF=$(ls "$GFF_DIR"/*.gff3 | wc -l)
echo "Running Panaroo on $N_GFF GFF3 files..."
echo "Output: $OUT_DIR"

mkdir -p "$OUT_DIR"

panaroo \
    -i "$GFF_DIR"/*.gff3 \
    -o "$OUT_DIR" \
    --clean-mode strict \
    --remove-invalid-genes \
    -c 0.70 \
    -f 0.70 \
    --len_dif_percent 0.95 \
    --aligner mafft \
    -t 8 \
    2>&1 | tee "$OUT_DIR/panaroo.log"

echo ""
echo "Panaroo complete."
echo "Key outputs:"
echo "  $OUT_DIR/gene_presence_absence.csv"
echo "  $OUT_DIR/pan_genome_reference.fa"
echo "  $OUT_DIR/summary_statistics.txt"

if [[ -f "$OUT_DIR/summary_statistics.txt" ]]; then
    cat "$OUT_DIR/summary_statistics.txt"
fi
