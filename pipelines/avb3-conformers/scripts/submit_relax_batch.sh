#!/bin/bash
set -euo pipefail
#
# Submit one SLURM relax job per selected frame.
#
# Usage:
#   RELAX_FRAMES=10,20,30,40,50 bash submit_relax_batch.sh
#

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/conformers}"
PULL_FRAMES_DIR="$WORK_DIR/pull_frames"
RELAX_FRAMES="${RELAX_FRAMES:-10,20,30,40,50}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

IFS=',' read -ra FRAMES <<< "$RELAX_FRAMES"

for frame_idx in "${FRAMES[@]}"; do
    # Find the frame PDB
    frame_pdb=$(ls "$PULL_FRAMES_DIR"/*"$(printf '%03d' "$frame_idx")"*.pdb 2>/dev/null | head -1 || true)
    if [[ -z "$frame_pdb" || ! -f "$frame_pdb" ]]; then
        echo "WARNING: Frame $frame_idx not found in $PULL_FRAMES_DIR, skipping."
        continue
    fi

    seed_name="seed_$(printf '%03d' "$frame_idx")"
    echo "Submitting relax job: $seed_name ($frame_pdb)"

    sbatch \
        --export="ALL,FRAME_PDB=$frame_pdb,SEED_NAME=$seed_name,WORK_DIR=$WORK_DIR" \
        --job-name="avb3_relax_${seed_name}" \
        "$SCRIPT_DIR/submit_relax_slurm.sh"
done
