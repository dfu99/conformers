#!/bin/bash
#SBATCH --job-name=a5b1_conjugates_first
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/protenix-a5b1/a5b1_conjugates_first_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/protenix-a5b1/a5b1_conjugates_first_%j.err
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH --constraint=A100-80GB
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

set -euo pipefail

module load cuda
module load python/3.12

PROTENIX_VENV="${PROTENIX_VENV:-$HOME/scratch/venv_protenix}"
if [[ ! -f "$PROTENIX_VENV/bin/activate" ]]; then
  echo "ERROR: Protenix venv not found: $PROTENIX_VENV" >&2
  exit 1
fi
source "$PROTENIX_VENV/bin/activate"

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
RUN_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/run_conjugates_first_pipeline.sh"
if [[ ! -f "$RUN_SCRIPT" ]]; then
  echo "ERROR: run script not found: $RUN_SCRIPT" >&2
  exit 1
fi

echo "========================================="
echo "Job ID:           $SLURM_JOB_ID"
echo "Node:             $(hostname)"
echo "GPU:              $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo N/A)"
echo "CONFORMERS_ROOT:  $CONFORMERS_ROOT"
echo "PROTENIX_VENV:    $PROTENIX_VENV"
echo "RUN_SCRIPT:       $RUN_SCRIPT"
echo "========================================="

srun bash "$RUN_SCRIPT"
