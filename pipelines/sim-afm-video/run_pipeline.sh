#!/usr/bin/env bash
# Sim HS-AFM video pipeline — end-to-end from a conformer library to a
# realistic stylized GIF.
#
# Stages:
#   1. ingest      — load PDB library + metadata into trajectory array
#   2. stabilize   — PCA-flatten + side-lock + stepwise yaw + head anchor
#   3. render      — forward through afmfold imaging model with z-blur
#   4. stylize     — copper colormap + substrate noise + slant + jitter
#                    + flash streaks (single-canvas, no inset border)
#
# Usage:
#   bash run_pipeline.sh <library_dir> <output_dir>
#
# Required: library_dir contains library.json (see load_conformer_library.py)
# and at least one PDB. afmfold must be available at AFMFOLD_ROOT.

set -euo pipefail
LIBRARY_DIR=${1:?"library_dir required"}
OUTPUT_DIR=${2:?"output_dir required"}

# Pipeline parameters (defaults match the validated v15 αVβ3 render).
WIDTH=${WIDTH:-60}
HEIGHT=${HEIGHT:-60}
RES_NM=${RES_NM:-0.98}
TIP_RADIUS_NM=${TIP_RADIUS_NM:-3.0}
TIP_ANGLE_DEG=${TIP_ANGLE_DEG:-20}
Z_BLUR_NM=${Z_BLUR_NM:-0.8}
Z_NOISE_NM=${Z_NOISE_NM:-0.35}
HEAD_RES_A=${HEAD_RES_A:-1-440}
HEAD_RES_B=${HEAD_RES_B:-1-350}
SIDE_CHAIN=${SIDE_CHAIN:-A}
YAW_THRESHOLD_DEG=${YAW_THRESHOLD_DEG:-50}
YAW_DWELL_FRAMES=${YAW_DWELL_FRAMES:-20}
YAW_STEP_CAP_DEG=${YAW_STEP_CAP_DEG:-30}
HEAD_ANCHOR_SIGMA=${HEAD_ANCHOR_SIGMA:-8}
UPSCALE=${UPSCALE:-5}
AFMFOLD_ROOT=${AFMFOLD_ROOT:-/home2/Documents/code/afmfold}
DEVICE=${DEVICE:-cpu}

mkdir -p "$OUTPUT_DIR"
SCRIPTS=$(dirname "$0")/scripts

echo "=== Stage 1: ingest ==="
python "$SCRIPTS/load_conformer_library.py" \
    --library "$LIBRARY_DIR" \
    --output "$OUTPUT_DIR/ingest"

echo "=== Stage 2: stabilize orientation ==="
python "$SCRIPTS/stabilize_orientation.py" \
    --fitted-coords "$OUTPUT_DIR/ingest/coords.npy" \
    --topology "$OUTPUT_DIR/ingest/topology.pdb" \
    --output "$OUTPUT_DIR/stable/coords_stable.npy" \
    --head-residues-a "$HEAD_RES_A" \
    --head-residues-b "$HEAD_RES_B" \
    --side-lock-chain "$SIDE_CHAIN" \
    --yaw-step-threshold "$YAW_THRESHOLD_DEG" \
    --yaw-step-cap "$YAW_STEP_CAP_DEG" \
    --yaw-dwell-frames "$YAW_DWELL_FRAMES" \
    --head-anchor-sigma "$HEAD_ANCHOR_SIGMA"

echo "=== Stage 3: forward-render through afmfold imaging model ==="
mkdir -p "$OUTPUT_DIR/render"
cp "$OUTPUT_DIR/ingest/topology.pdb" "$OUTPUT_DIR/render/"
cp "$OUTPUT_DIR/stable/coords_stable.npy" "$OUTPUT_DIR/render/"
# simulate_afm_video.py expects fitted_coords + head_positions_px +
# fitting_metadata.json. For a library-driven render (no real AFM
# comparison), we need a minimal head_positions file and metadata.
python -c "
import numpy as np
import json
coords = np.load('$OUTPUT_DIR/render/coords_stable.npy')
n = coords.shape[0]
np.save('$OUTPUT_DIR/render/head_positions_px.npy',
        np.full((n, 2), $WIDTH // 2, dtype=np.float32))
with open('$OUTPUT_DIR/render/fitting_metadata.json', 'w') as f:
    json.dump({'skip_frames': 0}, f)
"
# Touch a placeholder gif so the script doesn't fail on the comparison
# step. The output sim_images.npy is what we need downstream.
python -c "
from PIL import Image
import numpy as np
arr = np.zeros(($HEIGHT, $WIDTH), dtype=np.uint8)
Image.fromarray(arr, mode='L').save('$OUTPUT_DIR/render/placeholder.gif')
"
python "$SCRIPTS/simulate_afm_video.py" \
    --fitted-dir "$OUTPUT_DIR/render" \
    --coord-file coords_stable.npy \
    --gif "$OUTPUT_DIR/render/placeholder.gif" \
    --output-dir "$OUTPUT_DIR/render" \
    --width "$WIDTH" --height "$HEIGHT" --resolution-nm "$RES_NM" \
    --tip-radius "$TIP_RADIUS_NM" --tip-angle "$TIP_ANGLE_DEG" \
    --z-blur-nm "$Z_BLUR_NM" --noise-nm "$Z_NOISE_NM" \
    --afmfold-root "$AFMFOLD_ROOT" --device "$DEVICE"

echo "=== Stage 4: stylize (single-canvas, no inset border) ==="
python "$SCRIPTS/stylize_sim_afm.py" \
    --sim-images "$OUTPUT_DIR/render/sim_images.npy" \
    --output "$OUTPUT_DIR/sim_afm.gif" \
    --upscale "$UPSCALE"

echo "=== Done ==="
echo "Final GIF: $OUTPUT_DIR/sim_afm.gif"
