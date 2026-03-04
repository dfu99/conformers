#!/bin/bash
set -euo pipefail

# Minimal Protenix CLI inference wrapper (template-enabled).
# Mirrors Protenix repository example style:
#   protenix pred -i <input.json> -o <out_dir> -s 101 -n protenix_base_default_v1.0.0 \
#     --use_template true --use_default_params true

INPUT_JSON="${1:-data/runs/avb3/protenix_template/simple_pipeline/inputs/seed_090_frame_000_template.json}"
OUT_DIR="${2:-data/runs/avb3/protenix_template/simple_pipeline/outputs/output_base_v1}"
SEED="${3:-101}"
MODEL_NAME="${4:-protenix_base_default_v1.0.0}"
USE_TEMPLATE="${5:-true}"

if ! command -v protenix >/dev/null 2>&1; then
  echo "ERROR: protenix CLI not found in PATH. Activate your Protenix environment first."
  exit 1
fi

if [ "$USE_TEMPLATE" = "true" ]; then
  if ! command -v kalign >/dev/null 2>&1 && ! command -v kalign3 >/dev/null 2>&1; then
    echo "ERROR: template-enabled inference requires kalign (or kalign3) in PATH."
    exit 1
  fi
fi

if [ -z "${PROTENIX_ROOT_DIR:-}" ]; then
  echo "WARNING: PROTENIX_ROOT_DIR is not set."
  echo "         Set it if your Protenix runtime needs an explicit data root."
fi
if [ -z "${CUTLASS_PATH:-}" ]; then
  echo "WARNING: CUTLASS_PATH is not set."
  echo "         Set it if your environment requires it for deepspeed kernels."
fi

echo "Running protenix pred with:"
echo "  input_json:    $INPUT_JSON"
echo "  output_dir:    $OUT_DIR"
echo "  seed:          $SEED"
echo "  model_name:    $MODEL_NAME"
echo "  use_template:  $USE_TEMPLATE"

protenix pred \
  -i "$INPUT_JSON" \
  -o "$OUT_DIR" \
  -s "$SEED" \
  -n "$MODEL_NAME" \
  --use_template "$USE_TEMPLATE" \
  --use_default_params true
