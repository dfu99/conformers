#!/bin/bash
#SBATCH -J a5b1_cofold_avispy
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/protenix-a5b1/cofold_avispy_%j.log
#SBATCH --error=logs/protenix-a5b1/cofold_avispy_%j.err
set -euo pipefail

# Co-fold α5-Avi + β1-SpyCatcher heterodimer via Protenix
# 1789 total residues — needs A100 (trying without 80GB constraint first)

module load cuda

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
PROTENIX_ROOT="${PROTENIX_ROOT:-$HOME/scratch/protenix}"
INPUT_JSON="${INPUT_JSON:-$CONFORMERS_ROOT/data/a5b1/inputs/a5b1_avi_spycatcher_cofold.json}"
OUTPUT_DIR="${OUTPUT_DIR:-$CONFORMERS_ROOT/data/runs/a5b1/cofold_avi_spy/protenix_${SLURM_JOB_ID}}"
SEED="${SEED:-101}"
NUM_SAMPLES="${NUM_SAMPLES:-5}"

mkdir -p "$OUTPUT_DIR" "$CONFORMERS_ROOT/logs/protenix-a5b1"

# Activate Protenix environment
source "$PROTENIX_ROOT/venv_protenix/bin/activate" 2>/dev/null || \
    source "$CONFORMERS_ROOT/venv_protenix/bin/activate" 2>/dev/null || true

echo "=== A5B1 Avi-SpyCatcher Co-fold ==="
echo "Input: $INPUT_JSON"
echo "Output: $OUTPUT_DIR"
echo "Seed: $SEED"
echo "Samples: $NUM_SAMPLES"
echo "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo 'unknown')"

protenix pred \
    --input "$INPUT_JSON" \
    --out_dir "$OUTPUT_DIR" \
    --seeds "$SEED" \
    --sample_n "$NUM_SAMPLES" \
    --use_msa true \
    --use_template false

echo "Co-fold complete: $OUTPUT_DIR"
