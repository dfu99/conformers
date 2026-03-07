#!/bin/bash
set -euo pipefail

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
WORKFLOW_DIR="${WORKFLOW_DIR:-$CONFORMERS_ROOT/data/runs/a5b1/staged_attachment}"

SEQUENCE_FILE="${SEQUENCE_FILE:-$CONFORMERS_ROOT/data/a5b1/sequences/sequences_updated}"
HETERODIMER_PREDICTIONS_DIR="${HETERODIMER_PREDICTIONS_DIR:-$CONFORMERS_ROOT/data/a5b1/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/seed_101/predictions}"
HETERODIMER_CIF="${HETERODIMER_CIF:-}"

MODEL_NAME="${MODEL_NAME:-protenix_base_default_v1.0.0}"
SEED="${SEED:-101}"
DTYPE="${DTYPE:-bf16}"
USE_MSA="${USE_MSA:-true}"
USE_TEMPLATE="${USE_TEMPLATE:-false}"
USE_DEFAULT_PARAMS="${USE_DEFAULT_PARAMS:-true}"
TRIATT_KERNEL="${TRIATT_KERNEL:-torch}"
TRIMUL_KERNEL="${TRIMUL_KERNEL:-torch}"
FORCE_RERUN="${FORCE_RERUN:-0}"

A5B1_SETUP_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/setup_staged_attachment_workflow.py"
A5B1_INPUT_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/build_staged_protenix_inputs.py"
A5B1_MERGE_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/merge_staged_tagged_complex.py"

for f in "$A5B1_SETUP_SCRIPT" "$A5B1_INPUT_SCRIPT" "$A5B1_MERGE_SCRIPT"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: required script not found: $f" >&2
    exit 1
  fi
done
if [[ ! -f "$SEQUENCE_FILE" ]]; then
  echo "ERROR: sequence file not found: $SEQUENCE_FILE" >&2
  exit 1
fi
if [[ -z "$HETERODIMER_CIF" && ! -d "$HETERODIMER_PREDICTIONS_DIR" ]]; then
  echo "ERROR: heterodimer predictions dir not found: $HETERODIMER_PREDICTIONS_DIR" >&2
  echo "Set HETERODIMER_CIF to bypass auto-selection from predictions." >&2
  exit 1
fi

STAGE1_OUT="$WORKFLOW_DIR/outputs/protenix_stage1_spytag"
STAGE2_OUT="$WORKFLOW_DIR/outputs/protenix_stage2_streptavidin"
FINAL_OUT_DIR="$WORKFLOW_DIR/outputs/final"
mkdir -p "$STAGE1_OUT" "$STAGE2_OUT" "$FINAL_OUT_DIR"

echo "========================================="
echo "A5B1 complete tagged pipeline"
echo "CONFORMERS_ROOT=$CONFORMERS_ROOT"
echo "WORKFLOW_DIR=$WORKFLOW_DIR"
echo "SEQUENCE_FILE=$SEQUENCE_FILE"
echo "HETERODIMER_PREDICTIONS_DIR=$HETERODIMER_PREDICTIONS_DIR"
echo "HETERODIMER_CIF=${HETERODIMER_CIF:-<auto>}"
echo "MODEL_NAME=$MODEL_NAME"
echo "SEED=$SEED"
echo "DTYPE=$DTYPE"
echo "USE_MSA=$USE_MSA"
echo "USE_TEMPLATE=$USE_TEMPLATE"
echo "USE_DEFAULT_PARAMS=$USE_DEFAULT_PARAMS"
echo "FORCE_RERUN=$FORCE_RERUN"
echo "========================================="

SETUP_ARGS=(
  --sequence-file "$SEQUENCE_FILE"
  --predictions-dir "$HETERODIMER_PREDICTIONS_DIR"
  --outdir "$WORKFLOW_DIR"
)
if [[ -n "$HETERODIMER_CIF" ]]; then
  SETUP_ARGS+=(--heterodimer-cif "$HETERODIMER_CIF")
fi
python3 "$A5B1_SETUP_SCRIPT" "${SETUP_ARGS[@]}"

python3 "$A5B1_INPUT_SCRIPT" \
  --sequence-file "$SEQUENCE_FILE" \
  --workflow-dir "$WORKFLOW_DIR"

STAGE1_JSON="$WORKFLOW_DIR/inputs/protenix/stage1_spytag_input.json"
STAGE2_JSON="$WORKFLOW_DIR/inputs/protenix/stage2_streptavidin_input.json"

if [[ ! -f "$STAGE1_JSON" || ! -f "$STAGE2_JSON" ]]; then
  echo "ERROR: staged Protenix input JSONs not generated under $WORKFLOW_DIR/inputs/protenix" >&2
  exit 1
fi

get_job_name() {
  local path="$1"
  python3 - "$path" <<'PY'
import json
import sys
from pathlib import Path
p = Path(sys.argv[1])
obj = json.loads(p.read_text(encoding='utf-8'))
print(obj[0]['name'])
PY
}

STAGE1_JOB_NAME="$(get_job_name "$STAGE1_JSON")"
STAGE2_JOB_NAME="$(get_job_name "$STAGE2_JSON")"

STAGE1_PRED_DIR="$STAGE1_OUT/$STAGE1_JOB_NAME/seed_${SEED}/predictions"
STAGE2_PRED_DIR="$STAGE2_OUT/$STAGE2_JOB_NAME/seed_${SEED}/predictions"

run_stage_if_needed() {
  local stage_name="$1"
  local input_json="$2"
  local output_dir="$3"
  local pred_dir="$4"

  if [[ "$FORCE_RERUN" != "1" && -d "$pred_dir" ]] && ls "$pred_dir"/*_summary_confidence_sample_*.json >/dev/null 2>&1; then
    echo "[$stage_name] Found existing predictions at $pred_dir; skipping (FORCE_RERUN=$FORCE_RERUN)."
    return
  fi

  echo "[$stage_name] Running protenix pred"
  protenix pred \
    -i "$input_json" \
    -o "$output_dir" \
    -s "$SEED" \
    -n "$MODEL_NAME" \
    --dtype "$DTYPE" \
    --use_msa "$USE_MSA" \
    --use_template "$USE_TEMPLATE" \
    --use_default_params "$USE_DEFAULT_PARAMS" \
    --triatt_kernel "$TRIATT_KERNEL" \
    --trimul_kernel "$TRIMUL_KERNEL"

  if [[ ! -d "$pred_dir" ]] || ! ls "$pred_dir"/*_summary_confidence_sample_*.json >/dev/null 2>&1; then
    echo "ERROR: $stage_name predictions not found at expected path: $pred_dir" >&2
    exit 1
  fi
}

run_stage_if_needed "stage1_spytag" "$STAGE1_JSON" "$STAGE1_OUT" "$STAGE1_PRED_DIR"
run_stage_if_needed "stage2_streptavidin" "$STAGE2_JSON" "$STAGE2_OUT" "$STAGE2_PRED_DIR"

BASE_CIF="${BASE_CIF:-$(python3 - "$WORKFLOW_DIR/inputs/components.json" <<'PY'
import json
import sys
from pathlib import Path
p = Path(sys.argv[1])
obj = json.loads(p.read_text(encoding='utf-8'))
print(obj['heterodimer_cif'])
PY
)}"

FINAL_CIF="$FINAL_OUT_DIR/a5b1_tagged_complete.cif"
FINAL_PDB="$FINAL_OUT_DIR/a5b1_tagged_complete.pdb"

python3 "$A5B1_MERGE_SCRIPT" \
  --base-cif "$BASE_CIF" \
  --stage1-predictions-dir "$STAGE1_PRED_DIR" \
  --stage2-predictions-dir "$STAGE2_PRED_DIR" \
  --receptor-chains A,B \
  --stage1-ligand-chain C \
  --stage2-ligand-chain C \
  --out-stage1-chain C \
  --out-stage2-chain D \
  --out-cif "$FINAL_CIF" \
  --out-pdb "$FINAL_PDB"

echo "Pipeline complete."
echo "  Final CIF: $FINAL_CIF"
echo "  Final PDB: $FINAL_PDB"
