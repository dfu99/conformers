#!/bin/bash
set -euo pipefail

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
WORKFLOW_DIR="${WORKFLOW_DIR:-$CONFORMERS_ROOT/data/runs/a5b1/conjugates_first}"
SEQUENCE_FILE="${SEQUENCE_FILE:-$CONFORMERS_ROOT/data/a5b1/sequences/sequences_updated}"

MODEL_NAME="${MODEL_NAME:-protenix_base_default_v1.0.0}"
SEED="${SEED:-101}"
DTYPE="${DTYPE:-bf16}"
USE_MSA="${USE_MSA:-true}"
USE_TEMPLATE="${USE_TEMPLATE:-false}"
USE_DEFAULT_PARAMS="${USE_DEFAULT_PARAMS:-true}"
TRIATT_KERNEL="${TRIATT_KERNEL:-torch}"
TRIMUL_KERNEL="${TRIMUL_KERNEL:-torch}"
FORCE_RERUN="${FORCE_RERUN:-0}"
MERGE_SELECTION_MODE="${MERGE_SELECTION_MODE:-ranking}"

INPUT_BUILDER="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/build_conjugates_first_protenix_inputs.py"
MERGE_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/merge_staged_tagged_complex.py"
for f in "$INPUT_BUILDER" "$MERGE_SCRIPT"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: required script not found: $f" >&2
    exit 1
  fi
done
if [[ ! -f "$SEQUENCE_FILE" ]]; then
  echo "ERROR: sequence file not found: $SEQUENCE_FILE" >&2
  exit 1
fi

python3 "$INPUT_BUILDER" \
  --sequence-file "$SEQUENCE_FILE" \
  --workflow-dir "$WORKFLOW_DIR"

STAGE1_JSON="$WORKFLOW_DIR/inputs/protenix/stage1_alpha_streptavidin_input.json"
STAGE2_JSON="$WORKFLOW_DIR/inputs/protenix/stage2_beta_spytag_input.json"
STAGE3_JSON="$WORKFLOW_DIR/inputs/protenix/stage3_a5b1_docking_input.json"

for f in "$STAGE1_JSON" "$STAGE2_JSON" "$STAGE3_JSON"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: expected stage input missing: $f" >&2
    exit 1
  fi
done

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

stage_output_dir() {
  local stage_slug="$1"
  echo "$WORKFLOW_DIR/outputs/protenix_${stage_slug}"
}

run_stage_if_needed() {
  local stage_name="$1"
  local stage_slug="$2"
  local input_json="$3"

  local out_dir
  out_dir="$(stage_output_dir "$stage_slug")"
  local job_name
  job_name="$(get_job_name "$input_json")"
  local pred_dir="$out_dir/$job_name/seed_${SEED}/predictions"

  mkdir -p "$out_dir"

  if [[ "$FORCE_RERUN" != "1" && -d "$pred_dir" ]] && ls "$pred_dir"/*_summary_confidence_sample_*.json >/dev/null 2>&1; then
    echo "[$stage_name] Found existing predictions at $pred_dir; skipping (FORCE_RERUN=$FORCE_RERUN)."
    return
  fi

  echo "[$stage_name] Running protenix pred"
  protenix pred \
    -i "$input_json" \
    -o "$out_dir" \
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

best_summary_json() {
  local pred_dir="$1"
  python3 - "$pred_dir" <<'PY'
import json
import sys
from pathlib import Path
pred_dir = Path(sys.argv[1])
best = None
for s in sorted(pred_dir.glob("*_summary_confidence_sample_*.json")):
    score = float(json.loads(s.read_text(encoding='utf-8')).get("ranking_score", float("-inf")))
    cif = s.with_name(s.name.replace("_summary_confidence", "")).with_suffix(".cif")
    if not cif.exists():
        continue
    if best is None or score > best[0]:
        best = (score, s, cif)
if best is None:
    raise SystemExit(f"No valid summary/cif pairs in {pred_dir}")
score, summary, cif = best
print(json.dumps({"ranking_score": score, "summary_json": str(summary), "best_cif": str(cif)}))
PY
}

run_stage_if_needed "stage1_alpha_streptavidin" "stage1_alpha_streptavidin" "$STAGE1_JSON"
run_stage_if_needed "stage2_beta_spytag" "stage2_beta_spytag" "$STAGE2_JSON"
run_stage_if_needed "stage3_a5b1_docking" "stage3_a5b1_docking" "$STAGE3_JSON"

STAGE1_JOB="$(get_job_name "$STAGE1_JSON")"
STAGE2_JOB="$(get_job_name "$STAGE2_JSON")"
STAGE3_JOB="$(get_job_name "$STAGE3_JSON")"

STAGE1_PRED_DIR="$(stage_output_dir "stage1_alpha_streptavidin")/$STAGE1_JOB/seed_${SEED}/predictions"
STAGE2_PRED_DIR="$(stage_output_dir "stage2_beta_spytag")/$STAGE2_JOB/seed_${SEED}/predictions"
STAGE3_PRED_DIR="$(stage_output_dir "stage3_a5b1_docking")/$STAGE3_JOB/seed_${SEED}/predictions"

STAGE1_BEST_JSON="$(best_summary_json "$STAGE1_PRED_DIR")"
STAGE2_BEST_JSON="$(best_summary_json "$STAGE2_PRED_DIR")"
STAGE3_BEST_JSON="$(best_summary_json "$STAGE3_PRED_DIR")"

json_field() {
  local payload="$1"
  local key="$2"
  python3 - "$payload" "$key" <<'PY'
import json
import sys
obj = json.loads(sys.argv[1])
print(obj[sys.argv[2]])
PY
}

STAGE1_BEST_CIF="$(json_field "$STAGE1_BEST_JSON" "best_cif")"
STAGE2_BEST_CIF="$(json_field "$STAGE2_BEST_JSON" "best_cif")"
STAGE3_BEST_CIF="$(json_field "$STAGE3_BEST_JSON" "best_cif")"

FINAL_DIR="$WORKFLOW_DIR/outputs/final"
mkdir -p "$FINAL_DIR"

MERGED_CIF="$FINAL_DIR/a5b1_conjugates_first_combined.cif"
MERGED_PDB="$FINAL_DIR/a5b1_conjugates_first_combined.pdb"

python3 "$MERGE_SCRIPT" \
  --base-cif "$STAGE3_BEST_CIF" \
  --stage1-predictions-dir "$STAGE1_PRED_DIR" \
  --stage2-predictions-dir "$STAGE2_PRED_DIR" \
  --stage1-receptor-chain-map A:A \
  --stage2-receptor-chain-map B:A \
  --stage1-ligand-chain B \
  --stage2-ligand-chain B \
  --out-stage1-chain C \
  --out-stage2-chain D \
  --selection-mode "$MERGE_SELECTION_MODE" \
  --out-cif "$MERGED_CIF" \
  --out-pdb "$MERGED_PDB"

MERGE_SUMMARY_JSON="${MERGED_CIF%.cif}.merge_summary.json"
SUMMARY_JSON="$FINAL_DIR/conjugates_first_summary.json"

python3 - "$SUMMARY_JSON" \
  "$STAGE1_BEST_JSON" \
  "$STAGE2_BEST_JSON" \
  "$STAGE3_BEST_JSON" \
  "$MERGED_CIF" \
  "$MERGED_PDB" \
  "$MERGE_SUMMARY_JSON" <<'PY'
import json
import sys
from pathlib import Path
out = Path(sys.argv[1])
stage1 = json.loads(sys.argv[2])
stage2 = json.loads(sys.argv[3])
stage3 = json.loads(sys.argv[4])
merged_cif = sys.argv[5]
merged_pdb = sys.argv[6]
merge_summary = sys.argv[7]
payload = {
    "note": "Experimental conjugates-first branch. Stage3 provides heterodimer base frame and stage1/2 ligands are merged by receptor alignment.",
    "stage1_alpha_streptavidin": stage1,
    "stage2_beta_spytag": stage2,
    "stage3_a5b1_docking": stage3,
    "combined_output": {
        "merged_cif": merged_cif,
        "merged_pdb": merged_pdb,
        "merge_summary_json": merge_summary,
    },
}
out.write_text(json.dumps(payload, indent=2), encoding='utf-8')
print(f"Wrote {out}")
PY

echo "Conjugates-first pipeline complete."
echo "  Summary: $SUMMARY_JSON"
echo "  Combined CIF: $MERGED_CIF"
echo "  Combined PDB: $MERGED_PDB"
