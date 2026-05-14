#!/bin/bash
# run_is_analysis.sh
#
# Master script: IS detection, variant finding, and v1.3 database build.
#
# Steps:
#   01  Extract all CDS proteins from v1.2 GenBank
#   02  Scan proteins with ISEScan pHMM profiles (hmmsearch)
#   03  BLASTp all-vs-all protein comparison
#   04  Build gene presence/absence matrix
#   05  Find IS-variant and deletion-variant groups
#   06  Build v1.3 database (remove variants, strip IS annotations)
#
# Prerequisites:
#   - HMMER (hmmbuild, hmmsearch): conda install -c bioconda hmmer
#   - ISEScan pHMMs: git clone --depth=1 https://github.com/xiezhq/ISEScan.git /tmp/ISEScan
#   - blastp/makeblastdb: available at the paths below
#   - Python 3 + Biopython
#
# Usage:
#   bash scripts/run_is_analysis.sh [--skip-blast] [--dry-run]
#
# --skip-blast  Skip step 03 if blastp_allvsall.tsv already exists
# --dry-run     Run step 06 in dry-run mode (preview without writing)

set -eo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPTS="$REPO_DIR/scripts"
IS_DIR="$REPO_DIR/DB/is_analysis"
SKIP_BLAST=false
DRY_RUN=""

for arg in "$@"; do
    case $arg in
        --skip-blast) SKIP_BLAST=true ;;
        --dry-run) DRY_RUN="--dry-run" ;;
    esac
done

mkdir -p "$IS_DIR"

echo "============================================================"
echo "  E. coli G1/G4 IS Analysis Pipeline"
echo "  Repo:    $REPO_DIR"
echo "  Output:  $IS_DIR"
echo "============================================================"

# ── Parse flags ──────────────────────────────────────────────────────────────
SKIP_PANAROO=false
USE_BLAST_ONLY=false
for arg in "$@"; do
    case $arg in
        --skip-panaroo) SKIP_PANAROO=true ;;
        --blast-only)   USE_BLAST_ONLY=true ;;
    esac
done

# ── Check prerequisites ───────────────────────────────────────────────────────
if ! command -v hmmsearch &>/dev/null; then
    echo "ERROR: hmmsearch not found. Install with: conda install -c bioconda hmmer"
    exit 1
fi

if [[ ! -f /tmp/ISEScan/pHMMs/clusters.faa.hmm ]]; then
    echo "ERROR: ISEScan pHMMs not found."
    echo "Run: git clone --depth=1 https://github.com/xiezhq/ISEScan.git /tmp/ISEScan"
    exit 1
fi

PANAROO_AVAILABLE=false
if command -v panaroo &>/dev/null; then
    PANAROO_AVAILABLE=true
fi

# ── Step 01: Extract proteins ─────────────────────────────────────────────────
echo ""
echo "[01] Extracting CDS proteins from v1.2 GenBank..."
python3 "$SCRIPTS/01_extract_proteins.py"

# ── Step 02a: IS detection via hmmsearch ─────────────────────────────────────
echo ""
echo "[02a] Running hmmsearch IS detection..."
bash "$SCRIPTS/02_detect_is_proteins.sh"

# ── Step 02b+c: Panaroo (preferred) or BLASTp (fallback) ─────────────────────
if [[ "$PANAROO_AVAILABLE" == true && "$USE_BLAST_ONLY" == false ]]; then
    echo ""
    if [[ "$SKIP_PANAROO" == true && -f "$REPO_DIR/DB/panaroo_out/gene_presence_absence.csv" ]]; then
        echo "[02b] Skipping GFF3 conversion (--skip-panaroo, output exists)"
        echo "[02c] Skipping Panaroo (--skip-panaroo, output exists)"
    else
        echo "[02b] Converting GenBank to GFF3 for Panaroo..."
        python3 "$SCRIPTS/02b_gbk_to_gff3.py"

        echo ""
        echo "[02c] Running Panaroo (gene presence/absence clustering)..."
        bash "$SCRIPTS/02c_run_panaroo.sh"
    fi

    # ── Step 03: Skip BLASTp (Panaroo is sufficient) ──────────────────────────
    echo ""
    echo "[03] Panaroo used for gene families — BLASTp all-vs-all skipped"

    # ── Step 04: Find variants from Panaroo output ────────────────────────────
    echo ""
    echo "[04] Finding IS-variant and deletion-variant groups (Panaroo)..."
    python3 "$SCRIPTS/04b_panaroo_find_variants.py"

    MERGE_TSV="$IS_DIR/panaroo_loci_to_merge.tsv"

else
    # Fallback: BLASTp all-vs-all path
    echo ""
    if [[ "$SKIP_BLAST" == true && -f "$IS_DIR/blastp_allvsall.tsv" ]]; then
        echo "[03] Skipping BLASTp (--skip-blast, file exists)"
    else
        echo "[03] Running BLASTp all-vs-all (~20-60 min)..."
        bash "$SCRIPTS/03_blastp_allvsall.sh"
    fi

    echo ""
    echo "[04] Building gene presence/absence matrix (BLASTp)..."
    python3 "$SCRIPTS/04_build_presence_absence.py"

    echo ""
    echo "[05] Finding IS-variant and deletion-variant groups (BLASTp)..."
    python3 "$SCRIPTS/05_find_variants.py"

    MERGE_TSV="$IS_DIR/loci_to_merge.tsv"
fi

# ── Step 06: Build v1.3 database ─────────────────────────────────────────────
echo ""
echo "[06] Building v1.3 database..."
# Pass the correct merge TSV to the build script via env var
PANAROO_MERGE="$MERGE_TSV" python3 "$SCRIPTS/06_build_v1.3.py" $DRY_RUN

echo ""
echo "============================================================"
echo "  Pipeline complete."
echo "  Review IS_DIR for full results:"
echo "    $IS_DIR/"
echo ""
echo "  Key output files:"
if [[ "$PANAROO_AVAILABLE" == true && "$USE_BLAST_ONLY" == false ]]; then
    echo "    panaroo_is_variant_groups.tsv       IS-variant groups"
    echo "    panaroo_deletion_variant_groups.tsv Deletion-variant groups"
    echo "    panaroo_loci_to_merge.tsv           All loci flagged for removal"
    echo "    panaroo_is_families.txt             IS gene family names"
else
    echo "    is_variant_groups.tsv       IS-variant groups found"
    echo "    deletion_variant_groups.tsv Deletion-variant groups found"
    echo "    loci_to_merge.tsv           All loci flagged for removal"
fi
echo ""
if [[ -z "$DRY_RUN" ]]; then
    echo "  Database:"
    echo "    DB/EC-K-typing_group1and4_v1.3.gbk"
else
    echo "  (Dry-run: v1.3 GenBank not written)"
fi
echo "============================================================"
