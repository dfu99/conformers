#!/usr/bin/env bash
# afm-overlay pipeline — fit a real HS-AFM video against a conformer
# library and render the fitted PDB overlay on top of the real video.
#
# Stages:
#   1. ingest       — load conformer library (PDBs + library.json)
#   2. pseudo-afm   — render one pseudo-AFM image per conformer for matching
#   3. fit          — per-frame correlation + SO(3) refinement + head tracking
#   4. tail-flips   — head-anchored 180° flip resolution
#   5. temporal     — rolling-median coord smoothing + per-frame head re-anchor
#   6. kinetics     — state classification (BC/EC/EO/Intermediate) + plots
#   7. overlay      — render fitted PDB projected over real AFM frames
#
# Usage:
#   bash run_pipeline.sh <afm_gif> <library_dir> <output_dir>

set -euo pipefail
AFM_GIF=${1:?"afm_gif (.gif) required"}
LIBRARY_DIR=${2:?"library_dir (folder + library.json) required"}
OUTPUT_DIR=${3:?"output_dir required"}

# Parameters (defaults match validated v7 αVβ3 run).
TIP_RADIUS_NM=${TIP_RADIUS_NM:-2.0}
TIP_ANGLE_DEG=${TIP_ANGLE_DEG:-20}
SKIP_FRAMES=${SKIP_FRAMES:-30}
HEAD_RES_A=${HEAD_RES_A:-1-440}
HEAD_RES_B=${HEAD_RES_B:-1-350}
TEMPORAL_WINDOW=${TEMPORAL_WINDOW:-7}
TAIL_FLIP_WINDOW=${TAIL_FLIP_WINDOW:-21}
TAIL_FLIP_HYSTERESIS=${TAIL_FLIP_HYSTERESIS:-0.05}
TAIL_FLIP_ITERS=${TAIL_FLIP_ITERS:-20}
PROTEIN=${PROTEIN:-avb3}
AFMFOLD_ROOT=${AFMFOLD_ROOT:-/home2/Documents/code/afmfold}
DEVICE=${DEVICE:-cpu}

mkdir -p "$OUTPUT_DIR"
SCRIPTS=$(dirname "$0")/scripts
SIM_PIPE_SCRIPTS=$(dirname "$0")/../sim-afm-video/scripts

echo "=== Stage 1: ingest conformer library ==="
python "$SIM_PIPE_SCRIPTS/load_conformer_library.py" \
    --library "$LIBRARY_DIR" \
    --output "$OUTPUT_DIR/ingest"

# Re-export per-frame PDBs from ingest (for downstream tools that
# expect a frames-dir of PDBs).
echo "=== Stage 2: build pseudo-AFM library for correlation matching ==="
python "$SCRIPTS/process_frames_to_afm.py" \
    --frames-dir "$LIBRARY_DIR" \
    --output-dir "$OUTPUT_DIR/pseudo_afm" \
    --afmfold-root "$AFMFOLD_ROOT" \
    --protein-name "$PROTEIN" \
    --min-tip-radius "$TIP_RADIUS_NM" --max-tip-radius "$TIP_RADIUS_NM" \
    --min-tip-angle "$TIP_ANGLE_DEG" --max-tip-angle "$TIP_ANGLE_DEG" \
    --batch-size 1 \
    --device "$DEVICE"

echo "=== Stage 3: fit (correlation matching + SO(3) + head tracking) ==="
python "$SCRIPTS/fit_with_head_tracking.py" \
    --gif "$AFM_GIF" \
    --output-dir "$OUTPUT_DIR/fit" \
    --frames-dir "$LIBRARY_DIR" \
    --matched-indices "$OUTPUT_DIR/pseudo_afm/matched_indices.npy" \
    --afmfold-root "$AFMFOLD_ROOT" \
    --skip-frames "$SKIP_FRAMES" \
    --tip-radius "$TIP_RADIUS_NM"

echo "=== Stage 4: tail-flip resolution ==="
python "$SCRIPTS/resolve_tail_flips.py" \
    --fitted-dir "$OUTPUT_DIR/fit" \
    --head-residues-a "$HEAD_RES_A" \
    --head-residues-b "$HEAD_RES_B" \
    --window "$TAIL_FLIP_WINDOW" \
    --hysteresis "$TAIL_FLIP_HYSTERESIS" \
    --iterations "$TAIL_FLIP_ITERS"

echo "=== Stage 5: temporal smoothing (rolling median + head re-anchor) ==="
python "$SCRIPTS/smooth_coords_temporal.py" \
    --fitted-dir "$OUTPUT_DIR/fit" \
    --window "$TEMPORAL_WINDOW"
python "$SCRIPTS/reprocess_smooth_coords.py" \
    --fitted-dir "$OUTPUT_DIR/fit" \
    --head-residues-a "$HEAD_RES_A" \
    --head-residues-b "$HEAD_RES_B"

echo "=== Stage 6: state kinetics ==="
python "$SCRIPTS/analyze_state_kinetics.py" \
    --fitted-dir "$OUTPUT_DIR/fit" \
    --output-dir "$OUTPUT_DIR/kinetics"

echo "=== Stage 7: overlay rendering ==="
python "$SCRIPTS/render_projection_overlay.py" \
    --fitted-dir "$OUTPUT_DIR/fit" \
    --gif "$AFM_GIF" \
    --output-dir "$OUTPUT_DIR/overlay" \
    --skip-frames "$SKIP_FRAMES"

echo "=== Done ==="
echo "Overlay: $OUTPUT_DIR/overlay/pdb_projection_video.gif"
echo "State kinetics: $OUTPUT_DIR/kinetics/state_kinetics.png"
echo "Fit metadata: $OUTPUT_DIR/fit/fitting_metadata.json"
