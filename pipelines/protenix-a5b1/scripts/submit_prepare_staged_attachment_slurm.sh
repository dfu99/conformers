#!/bin/bash
#SBATCH --job-name=a5b1_stage_setup
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/protenix-a5b1/a5b1_stage_setup_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/protenix-a5b1/a5b1_stage_setup_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1
#SBATCH -n4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

set -euo pipefail

module load python/3.12

PROTENIX_VENV="${PROTENIX_VENV:-$HOME/scratch/venv_protenix}"
if [[ ! -f "$PROTENIX_VENV/bin/activate" ]]; then
  echo "ERROR: Protenix venv not found: $PROTENIX_VENV" >&2
  echo "Set PROTENIX_VENV to your cluster protenix environment." >&2
  exit 1
fi
source "$PROTENIX_VENV/bin/activate"

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
RUN_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/run_prepare_staged_attachment.sh"
if [[ ! -f "$RUN_SCRIPT" ]]; then
  echo "ERROR: run script not found: $RUN_SCRIPT" >&2
  exit 1
fi

mkdir -p "$CONFORMERS_ROOT/logs/protenix-a5b1"

echo "========================================="
echo "Job ID:           $SLURM_JOB_ID"
echo "Node:             $(hostname)"
echo "CONFORMERS_ROOT:  $CONFORMERS_ROOT"
echo "PROTENIX_VENV:    $PROTENIX_VENV"
echo "RUN_SCRIPT:       $RUN_SCRIPT"
echo "========================================="

srun bash "$RUN_SCRIPT"
