#!/bin/bash
# =============================================================================
# submit_atb_screen.sh
#
# SLURM array job: screen all 888 AllTheBacteria batches for G1/G4 flanking
# genes (galF + gnd). Runs up to 50 tasks concurrently to avoid overloading
# the shared filesystem.
#
# Prerequisites on M3:
#   module load blast/2.15.0
#   git clone https://github.com/efosternyarko/EC-K-typing-G1G4.git \
#       /home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4
#
# Submit:
#   cd /home/ebenezef/js66_scratch/ebenn/atb_screen
#   sbatch submit_atb_screen.sh
#
# Monitor:
#   squeue -u ebenezef
#   tail -f logs/screen_<jobid>_<taskid>.log
#
# After completion:
#   python combine_atb_hits.py
# =============================================================================

#SBATCH --job-name=atb_g1g4
#SBATCH --array=1-888%50
#SBATCH --partition=comp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/ebenezef/js66_scratch/ebenn/atb_screen/logs/screen_%A_%a.log
#SBATCH --error=/home/ebenezef/js66_scratch/ebenn/atb_screen/logs/screen_%A_%a.log

set -euo pipefail

module load blast/2.15.0

SCRIPT="/home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/scripts/atb_screen_batch.py"

echo "Task ${SLURM_ARRAY_TASK_ID} starting on $(hostname) at $(date)"
python3 "$SCRIPT" "${SLURM_ARRAY_TASK_ID}"
echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
