#!/bin/bash
#SBATCH -J avb3_domain_steer
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:V100:1
#SBATCH -C V100-16GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/domain_steer_%j.log
#SBATCH --error=logs/avb3-conformers/domain_steer_%j.err
set -euo pipefail

# Domain-preserving steering for AVB3 integrin.
# Uses centroid angle torques, domain restraints, or CV biasing instead of
# brute-force head/tail pulling that destroys domain folding.
#
# Usage:
#   sbatch --export=ALL,STEERING_PRESET=gentle_open submit_domain_steering_slurm.sh
#   sbatch --export=ALL,STEERING_PRESET=restrained_pull submit_domain_steering_slurm.sh
#   sbatch --export=ALL,STEERING_PRESET=cv_distance_extend submit_domain_steering_slurm.sh

module load cuda anaconda3

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
ROYALMD_ROOT="${ROYALMD_ROOT:-$HOME/scratch/RoyalMD}"
STEERING_PRESET="${STEERING_PRESET:-gentle_open}"
PRODUCTION_TIME="${PRODUCTION_TIME:-1000}"  # ps
WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/domain_steering/${STEERING_PRESET}_${SLURM_JOB_ID}}"
INPUT_PDB="${INPUT_PDB:-$ROYALMD_ROOT/test_systems/AVB3_clean.pdb}"

mkdir -p "$WORK_DIR" "$CONFORMERS_ROOT/logs/avb3-conformers"

# Activate OpenMM conda environment
conda activate "$CONFORMERS_ROOT/venv_openmm"

echo "=== Domain-Preserving Steering ==="
echo "Preset: $STEERING_PRESET"
echo "Input: $INPUT_PDB"
echo "Output: $WORK_DIR"
echo "Production: ${PRODUCTION_TIME} ps"

# Run the domain steering wrapper
python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/run_domain_steering.py" \
    --input-pdb "$INPUT_PDB" \
    --output-dir "$WORK_DIR" \
    --royalmd-root "$ROYALMD_ROOT" \
    --steering-preset "$STEERING_PRESET" \
    --production-time "$PRODUCTION_TIME"

echo "Domain steering complete: $WORK_DIR"
