#!/bin/bash
# 02_detect_is_proteins.sh
#
# Scan all predicted proteins with the ISEScan pHMM profiles using hmmsearch.
# This identifies transposase-like proteins without needing ISEScan itself.
#
# Requires: hmmbuild / hmmsearch (HMMER 3.x) — install with:
#   conda install -c bioconda hmmer
#
# ISEScan pHMMs: /tmp/ISEScan/pHMMs/clusters.faa.hmm (cloned from GitHub)
#
# Usage: bash scripts/02_detect_is_proteins.sh

set -eo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
IS_DIR="$REPO_DIR/DB/is_analysis"
PHMMS="/tmp/ISEScan/pHMMs/clusters.faa.hmm"
FAA="$IS_DIR/all_proteins.faa"
OUT="$IS_DIR/hmmsearch_is_hits.tbl"

if [[ ! -f "$FAA" ]]; then
    echo "ERROR: $FAA not found. Run 01_extract_proteins.py first."
    exit 1
fi

if [[ ! -f "$PHMMS" ]]; then
    echo "ERROR: ISEScan pHMMs not found at $PHMMS"
    echo "Clone ISEScan: git clone --depth=1 https://github.com/xiezhq/ISEScan.git /tmp/ISEScan"
    exit 1
fi

echo "Running hmmsearch: $(wc -l < "$FAA") lines in FAA..."
hmmsearch \
    --tblout "$OUT" \
    --cpu 4 \
    -E 1e-5 \
    "$PHMMS" \
    "$FAA" \
    > "$IS_DIR/hmmsearch_is_hits.log"

echo "Done. Hits written to $OUT"
grep -v "^#" "$OUT" | wc -l | xargs echo "IS protein hits:"
