#!/bin/bash
#SBATCH -J avb3_relax
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-rtx6000
#SBATCH --gres=gpu:RTX_6000:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/relax_%j.log
#SBATCH --error=logs/avb3-conformers/relax_%j.err
set -euo pipefail

# AVB3 relaxation MD — run equilibrium simulation from a pulled frame.
# Submits one job per frame. Use submit_relax_batch.sh to launch all.

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
ROYALMD_ROOT="${ROYALMD_ROOT:-$HOME/scratch/RoyalMD}"

WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/conformers}"
FRAME_PDB="${FRAME_PDB:?Set FRAME_PDB to the path of the pulled frame PDB}"
SEED_NAME="${SEED_NAME:?Set SEED_NAME, e.g., seed_030}"
RELAX_PRODUCTION_TIME="${RELAX_PRODUCTION_TIME:-1000}"

RELAX_OUT="$WORK_DIR/relax/$SEED_NAME"

module load cuda

mkdir -p "$RELAX_OUT"
mkdir -p "$CONFORMERS_ROOT/logs/avb3-conformers"

source "$ROYALMD_ROOT/venv/bin/activate" 2>/dev/null || true

echo "=== AVB3 Relaxation: $SEED_NAME ==="
echo "Input: $FRAME_PDB"
echo "Output: $RELAX_OUT"
echo "Production: ${RELAX_PRODUCTION_TIME} ps"

python "$ROYALMD_ROOT/RoyalMD.py" "$FRAME_PDB" \
    --out-dir "$RELAX_OUT" \
    --no-pulling \
    --production-time "$RELAX_PRODUCTION_TIME"

echo "Relaxation complete: $SEED_NAME"
