#!/bin/bash
set -euo pipefail

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
WORKFLOW_DIR="${WORKFLOW_DIR:-$CONFORMERS_ROOT/data/runs/a5b1/staged_attachment}"
PRIMARY_HETERODIMER_PREDICTIONS_DIR="$CONFORMERS_ROOT/data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/seed_101/predictions"
LEGACY_HETERODIMER_PREDICTIONS_DIR="$CONFORMERS_ROOT/data/a5b1/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/seed_101/predictions"
HETERODIMER_PREDICTIONS_DIR="${HETERODIMER_PREDICTIONS_DIR:-$PRIMARY_HETERODIMER_PREDICTIONS_DIR}"

PIPELINE_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/run_complete_tagged_pipeline.sh"
CHECK_SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/check_tagged_structure_quality.py"

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
  echo "ERROR: pipeline script not found: $PIPELINE_SCRIPT" >&2
  exit 1
fi
if [[ ! -f "$CHECK_SCRIPT" ]]; then
  echo "ERROR: quality check script not found: $CHECK_SCRIPT" >&2
  exit 1
fi

SEED_LIST_CSV="${SEED_LIST_CSV:-101,202,303,404,505,606,707,808,909}"
BASE_SAMPLE_LIST_CSV="${BASE_SAMPLE_LIST_CSV:-0,1,2,3,4}"
BASE_CIF_GLOB="${BASE_CIF_GLOB:-}"
STAGE1_ANCHOR_LIST_CSV="${STAGE1_ANCHOR_LIST_CSV:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}"
STAGE2_ANCHOR_LIST_CSV="${STAGE2_ANCHOR_LIST_CSV:-1}"
SELECTION_MODE="${SELECTION_MODE:-tail_distance}"
FORCE_RERUN="${FORCE_RERUN:-0}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-}"

# Stricter geometry defaults (can be overridden via env).
MAX_ALPHA_TAIL_DISTANCE="${MAX_ALPHA_TAIL_DISTANCE:-20}"
MAX_BETA_TAIL_DISTANCE="${MAX_BETA_TAIL_DISTANCE:-20}"
ALPHA_ATTACHMENT_RESIDUE="${ALPHA_ATTACHMENT_RESIDUE:-966}"
ALPHA_ATTACHMENT_ATOM="${ALPHA_ATTACHMENT_ATOM:-NZ}"
MAX_ALPHA_ATTACHMENT_DISTANCE="${MAX_ALPHA_ATTACHMENT_DISTANCE:-6}"
BETA_ATTACHMENT_RESIDUE="${BETA_ATTACHMENT_RESIDUE:-735}"
BETA_ATTACHMENT_ATOM="${BETA_ATTACHMENT_ATOM:-NZ}"
BETA_PARTNER_RESIDUE="${BETA_PARTNER_RESIDUE:-10}"
BETA_PARTNER_ATOMS="${BETA_PARTNER_ATOMS:-OD1,OD2,CG}"
MAX_BETA_ATTACHMENT_DISTANCE="${MAX_BETA_ATTACHMENT_DISTANCE:-6}"
ALPHA_DISALLOWED_RESIDUES="${ALPHA_DISALLOWED_RESIDUES:-780-820}"
MIN_ALPHA_PARTNER_DISTANCE_TO_ALPHA_DISALLOWED="${MIN_ALPHA_PARTNER_DISTANCE_TO_ALPHA_DISALLOWED:-10}"
BETA_DISALLOWED_RESIDUES="${BETA_DISALLOWED_RESIDUES:-}"
MIN_BETA_PARTNER_DISTANCE_TO_BETA_DISALLOWED="${MIN_BETA_PARTNER_DISTANCE_TO_BETA_DISALLOWED:-}"

if [[ ! -d "$HETERODIMER_PREDICTIONS_DIR" && -d "$LEGACY_HETERODIMER_PREDICTIONS_DIR" ]]; then
  HETERODIMER_PREDICTIONS_DIR="$LEGACY_HETERODIMER_PREDICTIONS_DIR"
fi
if [[ ! -d "$HETERODIMER_PREDICTIONS_DIR" ]]; then
  echo "ERROR: heterodimer predictions dir not found: $HETERODIMER_PREDICTIONS_DIR" >&2
  exit 1
fi

FINAL_DIR="$WORKFLOW_DIR/outputs/final"
FINAL_CIF="$FINAL_DIR/a5b1_tagged_complete.cif"
FINAL_PDB="$FINAL_DIR/a5b1_tagged_complete.pdb"
FINAL_MERGE_SUMMARY="$FINAL_DIR/a5b1_tagged_complete.merge_summary.json"

TRIAL_DIR="$FINAL_DIR/trials"
TRIAL_REPORT_DIR="$TRIAL_DIR/reports"
mkdir -p "$TRIAL_DIR" "$TRIAL_REPORT_DIR"

IFS=',' read -r -a SEEDS <<< "$SEED_LIST_CSV"
IFS=',' read -r -a BASE_SAMPLES <<< "$BASE_SAMPLE_LIST_CSV"
IFS=',' read -r -a STAGE1_ANCHORS <<< "$STAGE1_ANCHOR_LIST_CSV"
IFS=',' read -r -a ANCHORS <<< "$STAGE2_ANCHOR_LIST_CSV"

echo "========================================="
echo "A5B1 auto-search until pass"
echo "SEED_LIST_CSV=$SEED_LIST_CSV"
echo "BASE_SAMPLE_LIST_CSV=$BASE_SAMPLE_LIST_CSV"
echo "BASE_CIF_GLOB=${BASE_CIF_GLOB:-<none>}"
echo "STAGE1_ANCHOR_LIST_CSV=$STAGE1_ANCHOR_LIST_CSV"
echo "STAGE2_ANCHOR_LIST_CSV=$STAGE2_ANCHOR_LIST_CSV"
echo "SELECTION_MODE=$SELECTION_MODE"
echo "FORCE_RERUN=$FORCE_RERUN"
echo "MAX_ATTEMPTS=${MAX_ATTEMPTS:-<none>}"
echo "MAX_ALPHA_TAIL_DISTANCE=$MAX_ALPHA_TAIL_DISTANCE"
echo "MAX_BETA_TAIL_DISTANCE=$MAX_BETA_TAIL_DISTANCE"
echo "ALPHA_ATTACHMENT_RESIDUE=$ALPHA_ATTACHMENT_RESIDUE"
echo "ALPHA_ATTACHMENT_ATOM=$ALPHA_ATTACHMENT_ATOM"
echo "MAX_ALPHA_ATTACHMENT_DISTANCE=$MAX_ALPHA_ATTACHMENT_DISTANCE"
echo "BETA_ATTACHMENT_RESIDUE=$BETA_ATTACHMENT_RESIDUE"
echo "BETA_ATTACHMENT_ATOM=$BETA_ATTACHMENT_ATOM"
echo "BETA_PARTNER_RESIDUE=$BETA_PARTNER_RESIDUE"
echo "BETA_PARTNER_ATOMS=$BETA_PARTNER_ATOMS"
echo "MAX_BETA_ATTACHMENT_DISTANCE=$MAX_BETA_ATTACHMENT_DISTANCE"
echo "ALPHA_DISALLOWED_RESIDUES=$ALPHA_DISALLOWED_RESIDUES"
echo "MIN_ALPHA_PARTNER_DISTANCE_TO_ALPHA_DISALLOWED=$MIN_ALPHA_PARTNER_DISTANCE_TO_ALPHA_DISALLOWED"
echo "BETA_DISALLOWED_RESIDUES=$BETA_DISALLOWED_RESIDUES"
echo "MIN_BETA_PARTNER_DISTANCE_TO_BETA_DISALLOWED=${MIN_BETA_PARTNER_DISTANCE_TO_BETA_DISALLOWED:-<none>}"
echo "HETERODIMER_PREDICTIONS_DIR=$HETERODIMER_PREDICTIONS_DIR"
echo "WORKFLOW_DIR=$WORKFLOW_DIR"
echo "========================================="

attempt=0
stop_due_to_max=0

build_base_label() {
  local path="$1"
  local name
  name="$(basename "$path")"
  name="${name%.cif}"
  echo "$name" | tr -cs '[:alnum:]' '_'
}

declare -a BASE_CIFS=()
if [[ -n "$BASE_CIF_GLOB" ]]; then
  # Intentional word splitting for shell glob expansion supplied by caller.
  # shellcheck disable=SC2086
  for base_cif in $BASE_CIF_GLOB; do
    [[ -f "$base_cif" ]] || continue
    BASE_CIFS+=("$base_cif")
  done
  if [[ ${#BASE_CIFS[@]} -eq 0 ]]; then
    echo "ERROR: no base CIFs matched BASE_CIF_GLOB=$BASE_CIF_GLOB" >&2
    exit 1
  fi
else
  for base_sample in "${BASE_SAMPLES[@]}"; do
    base_sample="$(echo "$base_sample" | xargs)"
    [[ -z "$base_sample" ]] && continue
    base_cif="$HETERODIMER_PREDICTIONS_DIR/integrin_alpha5_beta1_sample_${base_sample}.cif"
    if [[ ! -f "$base_cif" ]]; then
      echo "[base_sample=$base_sample] Missing base CIF: $base_cif; skipping."
      continue
    fi
    BASE_CIFS+=("$base_cif")
  done
fi

for base_cif in "${BASE_CIFS[@]}"; do
  base_label="$(build_base_label "$base_cif")"
  for seed in "${SEEDS[@]}"; do
    seed="$(echo "$seed" | xargs)"
    [[ -z "$seed" ]] && continue
    for stage1_anchor in "${STAGE1_ANCHORS[@]}"; do
      stage1_anchor="$(echo "$stage1_anchor" | xargs)"
      [[ -z "$stage1_anchor" ]] && continue
      for stage2_anchor in "${ANCHORS[@]}"; do
        stage2_anchor="$(echo "$stage2_anchor" | xargs)"
        [[ -z "$stage2_anchor" ]] && continue

        if [[ -n "$MAX_ATTEMPTS" && "$attempt" -ge "$MAX_ATTEMPTS" ]]; then
          stop_due_to_max=1
          break
        fi

        attempt=$((attempt + 1))
        tag="attempt_${attempt}_base_${base_label}_seed_${seed}_stage1_anchor_${stage1_anchor}_stage2_anchor_${stage2_anchor}"

        echo "-----------------------------------------"
        echo "[$tag] Running pipeline"

        if ! CONFORMERS_ROOT="$CONFORMERS_ROOT" \
          WORKFLOW_DIR="$WORKFLOW_DIR" \
          SEED="$seed" \
          HETERODIMER_CIF="$base_cif" \
          SELECTION_MODE="$SELECTION_MODE" \
          STAGE1_LIGAND_ANCHOR_RESIDUE="$stage1_anchor" \
          STAGE2_LIGAND_ANCHOR_RESIDUE="$stage2_anchor" \
          FORCE_RERUN="$FORCE_RERUN" \
          bash "$PIPELINE_SCRIPT"; then
          echo "[$tag] Pipeline run failed; continuing."
          continue
        fi

        trial_cif="$TRIAL_DIR/a5b1_tagged_complete_${tag}.cif"
        trial_pdb="$TRIAL_DIR/a5b1_tagged_complete_${tag}.pdb"
        trial_merge="$TRIAL_DIR/a5b1_tagged_complete_${tag}.merge_summary.json"
        trial_report="$TRIAL_REPORT_DIR/${tag}.quality.json"

        if [[ ! -f "$FINAL_CIF" || ! -f "$FINAL_PDB" || ! -f "$FINAL_MERGE_SUMMARY" ]]; then
          echo "[$tag] Missing final outputs after pipeline run; continuing."
          continue
        fi

        cp -f "$FINAL_CIF" "$trial_cif"
        cp -f "$FINAL_PDB" "$trial_pdb"
        cp -f "$FINAL_MERGE_SUMMARY" "$trial_merge"

        echo "[$tag] Checking quality"
        CHECK_ARGS=(
          --pdb "$trial_pdb"
          --alpha-chain A
          --beta-chain B
          --alpha-partner-chain D
          --beta-partner-chain C
          --expected-chain-count 4
          --max-alpha-tail-distance "$MAX_ALPHA_TAIL_DISTANCE"
          --max-beta-tail-distance "$MAX_BETA_TAIL_DISTANCE"
          --alpha-attachment-residue "$ALPHA_ATTACHMENT_RESIDUE"
          --alpha-attachment-atom "$ALPHA_ATTACHMENT_ATOM"
          --max-alpha-attachment-distance "$MAX_ALPHA_ATTACHMENT_DISTANCE"
          --beta-attachment-residue "$BETA_ATTACHMENT_RESIDUE"
          --beta-attachment-atom "$BETA_ATTACHMENT_ATOM"
          --beta-partner-residue "$BETA_PARTNER_RESIDUE"
          --beta-partner-atoms "$BETA_PARTNER_ATOMS"
          --max-beta-attachment-distance "$MAX_BETA_ATTACHMENT_DISTANCE"
          --report-json "$trial_report"
        )
        if [[ -n "$ALPHA_DISALLOWED_RESIDUES" ]]; then
          CHECK_ARGS+=(--alpha-disallowed-residues "$ALPHA_DISALLOWED_RESIDUES")
        fi
        if [[ -n "$MIN_ALPHA_PARTNER_DISTANCE_TO_ALPHA_DISALLOWED" ]]; then
          CHECK_ARGS+=(
            --min-alpha-partner-distance-to-alpha-disallowed
            "$MIN_ALPHA_PARTNER_DISTANCE_TO_ALPHA_DISALLOWED"
          )
        fi
        if [[ -n "$BETA_DISALLOWED_RESIDUES" ]]; then
          CHECK_ARGS+=(--beta-disallowed-residues "$BETA_DISALLOWED_RESIDUES")
        fi
        if [[ -n "$MIN_BETA_PARTNER_DISTANCE_TO_BETA_DISALLOWED" ]]; then
          CHECK_ARGS+=(
            --min-beta-partner-distance-to-beta-disallowed
            "$MIN_BETA_PARTNER_DISTANCE_TO_BETA_DISALLOWED"
          )
        fi

        if python3 "$CHECK_SCRIPT" "${CHECK_ARGS[@]}"; then
          echo "[$tag] PASS"
          cp -f "$trial_cif" "$FINAL_CIF"
          cp -f "$trial_pdb" "$FINAL_PDB"
          cp -f "$trial_merge" "$FINAL_MERGE_SUMMARY"
          cp -f "$trial_report" "$FINAL_DIR/a5b1_tagged_complete.quality.json"
          exit 0
        else
          echo "[$tag] FAIL"
        fi
      done
      if [[ "$stop_due_to_max" == "1" ]]; then
        break
      fi
    done
    if [[ "$stop_due_to_max" == "1" ]]; then
      break
    fi
  done
  if [[ "$stop_due_to_max" == "1" ]]; then
    break
  fi
done

if [[ "$stop_due_to_max" == "1" ]]; then
  echo "Stopped after MAX_ATTEMPTS=$MAX_ATTEMPTS without a passing structure." >&2
  exit 2
fi

echo "No passing structure found across all attempts." >&2
exit 2
