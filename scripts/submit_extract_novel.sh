#!/bin/bash
# =============================================================================
# submit_extract_novel.sh
#
# SLURM array job: extract novel G1/G4 loci from AllTheBacteria candidates.
#
# Prerequisites:
#   1. Screening complete (665 batch_*_hits.tsv files exist)
#   2. combine_atb_hits.py has been run (atb_g1g4_candidates.tsv exists)
#   3. setup_novelty_db.sh has been run (BLAST DB built from reference loci)
#
# Submit:
#   bash scripts/setup_novelty_db.sh        # run once first
#   sbatch scripts/submit_extract_novel.sh
#
# After completion:
#   python scripts/cluster_novel_loci.py
# =============================================================================

#SBATCH --job-name=atb_extract
#SBATCH --array=1-665%50
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/ebenezef/js66_scratch/ebenn/atb_screen/logs/extract_%A_%a.log
#SBATCH --error=/home/ebenezef/js66_scratch/ebenn/atb_screen/logs/extract_%A_%a.log

set -eo pipefail

module load blast/2.15.0

SCRIPT="/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/scripts/extract_novel_loci.py"

echo "Task ${SLURM_ARRAY_TASK_ID} starting on $(hostname) at $(date)"
python3 "$SCRIPT" "${SLURM_ARRAY_TASK_ID}"
echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
