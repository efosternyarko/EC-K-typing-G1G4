#!/bin/bash
# =============================================================================
# setup_novelty_db.sh
#
# Run ONCE before submitting extract_novel_loci job.
# Builds a BLAST nucleotide database from the 93 G1/G4 reference loci so
# that novelty screening can run without re-building it per task.
#
# Usage:
#   bash scripts/setup_novelty_db.sh
# =============================================================================

set -eo pipefail

module load blast/2.15.0

REPO="/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4"
REF_FASTA="${REPO}/DB/reference_loci_v0.6.fasta"
DB_DIR="/home/ebenezef/js66_scratch/ebenn/atb_screen/novelty_db"
DB_PATH="${DB_DIR}/g1g4_refs"

mkdir -p "$DB_DIR"

echo "Building BLAST DB from 93 G1/G4 reference loci..."
makeblastdb \
    -in "$REF_FASTA" \
    -dbtype nucl \
    -out "$DB_PATH" \
    -title "G1G4_reference_loci_v0.6"

echo "Done. DB at: ${DB_PATH}"
echo "Verify with: blastdbcmd -db ${DB_PATH} -info"
