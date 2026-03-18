#!/bin/bash
set -euo pipefail
#
# MSA-Subsampled Conformer Validation Pipeline
#
# Validates whether pulled RoyalMD conformers are evolutionarily realistic
# by checking if any MSA subsample produces a Protenix prediction that
# resembles the pulled structure.
#
# Stages:
#   1. PREPARE  — Extract extra frames, subsample MSAs, build Protenix inputs
#   2. SUBMIT   — Run Protenix at each MSA depth on PACE
#   3. SCORE    — Compare pulled frames vs predictions, output heatmap
#
# Usage:
#   bash run_msa_validation.sh prepare
#   bash run_msa_validation.sh submit
#   bash run_msa_validation.sh score
#

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}"
ROYALMD_ROOT="${ROYALMD_ROOT:-$(cd "$CONFORMERS_ROOT/../RoyalMD" && pwd)}"
PIPELINE_DIR="$CONFORMERS_ROOT/pipelines/avb3-conformers"

WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/msa_validation}"
PULLS_DIR="${PULLS_DIR:-$ROYALMD_ROOT/results-avb3_04}"

# Frame indices to validate (0-120 already split, need 150-300)
FRAME_INDICES="${FRAME_INDICES:-0,10,20,30,40,50,60,70,80,90,100,110,120,150,180,210,240,270,300}"

# MSA depth fractions
MSA_FRACTIONS="${MSA_FRACTIONS:-1.0,0.5,0.25,0.1,0.05,0.01}"

# ═══════════════════════════════════════════════════════════════════════════════
stage_prepare() {
    echo "=== Stage 1: PREPARE ==="
    mkdir -p "$WORK_DIR/frames"

    # Ensure we have all needed frames split from the trajectory
    local traj_pdb="$PULLS_DIR/production_trajectory.pdb"
    if [[ ! -f "$traj_pdb" ]]; then
        echo "Trajectory PDB not found. Converting from .nc first..."
        echo "Run in $PULLS_DIR:"
        echo "  cpptraj nc_to_pdb.in"
        echo "Then re-run this stage."
        return 1
    fi

    # Check which frames we need to extract
    local splits_dir="$PULLS_DIR/splits"
    local need_split=0

    IFS=',' read -ra FRAMES <<< "$FRAME_INDICES"
    for idx in "${FRAMES[@]}"; do
        local padded
        padded=$(printf "%03d" "$idx")
        if ! ls "$splits_dir"/*"${padded}"*.pdb >/dev/null 2>&1; then
            need_split=1
            echo "  Frame $idx not found in $splits_dir"
        fi
    done

    if [[ "$need_split" == "1" ]]; then
        echo "Splitting additional frames from trajectory..."
        # Find max frame needed
        local max_frame=0
        for idx in "${FRAMES[@]}"; do
            [[ "$idx" -gt "$max_frame" ]] && max_frame="$idx"
        done

        python3 "$ROYALMD_ROOT/scripts/split_pdb_trajectory.py" \
            "$traj_pdb" \
            --start 0 \
            --end "$max_frame" \
            --prefix "production_trajectory"

        echo "Split complete."
    fi

    # Copy selected frames to validation work directory
    echo "Copying selected frames to $WORK_DIR/frames/..."
    for idx in "${FRAMES[@]}"; do
        local padded
        padded=$(printf "%03d" "$idx")
        local src
        src=$(ls "$splits_dir"/*"${padded}"*.pdb 2>/dev/null | head -1 || true)
        if [[ -n "$src" && -f "$src" ]]; then
            cp "$src" "$WORK_DIR/frames/frame_${padded}.pdb"
        else
            echo "  WARNING: Frame $idx not found, skipping."
        fi
    done

    local n_frames
    n_frames=$(ls "$WORK_DIR/frames"/*.pdb 2>/dev/null | wc -l)
    echo "  $n_frames frames staged"

    # Build Protenix inputs (sequence-only for now, MSA paths added when available)
    local ref_pdb
    ref_pdb=$(ls "$WORK_DIR/frames"/frame_000.pdb 2>/dev/null || ls "$WORK_DIR/frames"/*.pdb 2>/dev/null | head -1)

    echo ""
    echo "Building Protenix sweep inputs..."
    python3 "$PIPELINE_DIR/scripts/build_msa_sweep_inputs.py" \
        --reference-pdb "$ref_pdb" \
        --output-dir "$WORK_DIR/inputs" \
        --depths "$MSA_FRACTIONS" \
        --no-msa-baseline

    echo ""
    echo "Prepare complete. Next: submit MSA sweep to PACE."
    echo "  sbatch $PIPELINE_DIR/scripts/submit_msa_sweep_slurm.sh"
}

# ═══════════════════════════════════════════════════════════════════════════════
stage_submit() {
    echo "=== Stage 2: SUBMIT ==="
    echo "Submit MSA sweep job to PACE:"
    echo ""
    echo "  sbatch --export=ALL,WORK_DIR=$WORK_DIR \\"
    echo "    $PIPELINE_DIR/scripts/submit_msa_sweep_slurm.sh"
    echo ""
    echo "Or if MSAs are not yet generated, run Protenix with use_msa=true"
    echo "and let it generate MSAs during the first run (depth_1.00), then"
    echo "subsample for subsequent runs."
}

# ═══════════════════════════════════════════════════════════════════════════════
stage_score() {
    echo "=== Stage 3: SCORE ==="

    if [[ ! -d "$WORK_DIR/predictions" ]]; then
        echo "ERROR: No predictions found at $WORK_DIR/predictions/" >&2
        echo "Run the PACE sweep first." >&2
        return 1
    fi

    python3 "$PIPELINE_DIR/scripts/score_conformers.py" \
        --frames-dir "$WORK_DIR/frames" \
        --predictions-dir "$WORK_DIR/predictions" \
        --output-dir "$WORK_DIR/scores" \
        --frame-indices "$FRAME_INDICES"

    echo ""
    echo "Scores saved to $WORK_DIR/scores/"
    echo "Key outputs:"
    echo "  conformer_scores.csv              — full score matrix"
    echo "  conformer_summary.json            — best score per frame"
    echo "  conformer_validation_heatmap.png  — TM-score & RMSD heatmaps"
    echo "  conformer_validity_curve.png      — validity vs frame index"
}

# ═══════════════════════════════════════════════════════════════════════════════
usage() {
    echo "Usage: $0 <stage>"
    echo ""
    echo "Stages:"
    echo "  prepare   Split extra frames, build Protenix inputs"
    echo "  submit    Print PACE submission command"
    echo "  score     Compare pulled frames vs predictions"
    echo "  status    Show pipeline state"
}

stage_status() {
    echo "=== MSA Validation Status ==="
    echo "Work dir: $WORK_DIR"
    echo ""
    echo "Frames:"
    echo "  $(ls "$WORK_DIR/frames"/*.pdb 2>/dev/null | wc -l) staged"
    echo ""
    echo "Inputs:"
    if [[ -f "$WORK_DIR/inputs/sweep_config.json" ]]; then
        echo "  Sweep config exists"
        python3 -c "import json; c=json.loads(open('$WORK_DIR/inputs/sweep_config.json').read()); print(f'  {len(c)} MSA depths configured')"
    else
        echo "  Not prepared"
    fi
    echo ""
    echo "Predictions:"
    if [[ -d "$WORK_DIR/predictions" ]]; then
        for d in "$WORK_DIR/predictions"/*/; do
            echo "  $(basename "$d"): $(find "$d" -name '*.cif' -o -name '*.pdb' 2>/dev/null | wc -l) files"
        done
    else
        echo "  Not run"
    fi
    echo ""
    echo "Scores:"
    if [[ -f "$WORK_DIR/scores/conformer_scores.csv" ]]; then
        echo "  Scores CSV exists"
        [[ -f "$WORK_DIR/scores/conformer_validation_heatmap.png" ]] && echo "  Heatmap generated"
    else
        echo "  Not scored"
    fi
}

STAGE="${1:-}"
case "$STAGE" in
    prepare)  stage_prepare ;;
    submit)   stage_submit ;;
    score)    stage_score ;;
    status)   stage_status ;;
    *)        usage; exit 1 ;;
esac
