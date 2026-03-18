#!/bin/bash
#SBATCH -J avb3_pull
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-rtx6000
#SBATCH --gres=gpu:RTX_6000:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/pull_%j.log
#SBATCH --error=logs/avb3-conformers/pull_%j.err
set -euo pipefail

# AVB3 steered MD — pulling force on head/tail to generate extended conformers.
# RTX 6000 is sufficient for OpenMM MD (not deep learning inference).

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
ROYALMD_ROOT="${ROYALMD_ROOT:-$HOME/scratch/RoyalMD}"

INPUT_PDB="${INPUT_PDB:-$ROYALMD_ROOT/test_systems/AVB3_clean.pdb}"
WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/conformers}"
PULL_DIR="$WORK_DIR/pull"
PULL_FORCE_PN="${PULL_FORCE_PN:-2.0}"
PULL_PRODUCTION_TIME="${PULL_PRODUCTION_TIME:-1000}"

module load cuda

mkdir -p "$PULL_DIR"
mkdir -p "$CONFORMERS_ROOT/logs/avb3-conformers"

# Activate RoyalMD environment
source "$ROYALMD_ROOT/venv/bin/activate" 2>/dev/null || true

echo "=== AVB3 Pull Simulation ==="
echo "Input: $INPUT_PDB"
echo "Output: $PULL_DIR"
echo "Force: ${PULL_FORCE_PN} pN"
echo "Production: ${PULL_PRODUCTION_TIME} ps"
echo ""

python "$ROYALMD_ROOT/RoyalMD.py" "$INPUT_PDB" \
    --out-dir "$PULL_DIR" \
    --pulling-enabled \
    --pulling-force-pn "$PULL_FORCE_PN" \
    --production-time "$PULL_PRODUCTION_TIME"

echo "Pull simulation complete."
echo "Trajectory: $PULL_DIR/production.nc"
