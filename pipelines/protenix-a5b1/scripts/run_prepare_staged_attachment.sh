#!/bin/bash
set -euo pipefail

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
SCRIPT="$CONFORMERS_ROOT/pipelines/protenix-a5b1/scripts/setup_staged_attachment_workflow.py"

SEQUENCE_FILE="${SEQUENCE_FILE:-$CONFORMERS_ROOT/data/a5b1/sequences/sequences_updated}"
PREDICTIONS_DIR="${PREDICTIONS_DIR:-$CONFORMERS_ROOT/data/a5b1/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/seed_101/predictions}"
HETERODIMER_CIF="${HETERODIMER_CIF:-}"
OUTDIR="${OUTDIR:-$CONFORMERS_ROOT/data/runs/a5b1/staged_attachment}"

ALPHA_CHAIN="${ALPHA_CHAIN:-A}"
BETA_CHAIN="${BETA_CHAIN:-B}"
ALPHA_NAME="${ALPHA_NAME:-Integrin alpha5-Avi}"
BETA_NAME="${BETA_NAME:-Integrin beta1-spycatcher}"
SPYTAG_NAME="${SPYTAG_NAME:-Spytag}"
STREPTAVIDIN_NAME="${STREPTAVIDIN_NAME:-Streptavidin}"

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: script not found: $SCRIPT" >&2
  exit 1
fi
if [[ ! -f "$SEQUENCE_FILE" ]]; then
  echo "ERROR: sequence file not found: $SEQUENCE_FILE" >&2
  exit 1
fi
if [[ -z "$HETERODIMER_CIF" && ! -d "$PREDICTIONS_DIR" ]]; then
  echo "ERROR: predictions dir not found: $PREDICTIONS_DIR" >&2
  echo "Set HETERODIMER_CIF to bypass auto-selection by ranking_score." >&2
  exit 1
fi
if [[ -n "$HETERODIMER_CIF" && ! -f "$HETERODIMER_CIF" ]]; then
  echo "ERROR: heterodimer cif not found: $HETERODIMER_CIF" >&2
  exit 1
fi

ARGS=(
  --sequence-file "$SEQUENCE_FILE"
  --predictions-dir "$PREDICTIONS_DIR"
  --outdir "$OUTDIR"
  --alpha-name "$ALPHA_NAME"
  --beta-name "$BETA_NAME"
  --spytag-name "$SPYTAG_NAME"
  --streptavidin-name "$STREPTAVIDIN_NAME"
  --alpha-chain "$ALPHA_CHAIN"
  --beta-chain "$BETA_CHAIN"
)
if [[ -n "$HETERODIMER_CIF" ]]; then
  ARGS+=(--heterodimer-cif "$HETERODIMER_CIF")
fi

echo "Preparing staged A5B1 attachment workflow"
echo "  CONFORMERS_ROOT=$CONFORMERS_ROOT"
echo "  SEQUENCE_FILE=$SEQUENCE_FILE"
echo "  PREDICTIONS_DIR=$PREDICTIONS_DIR"
echo "  HETERODIMER_CIF=${HETERODIMER_CIF:-<auto>}"
echo "  OUTDIR=$OUTDIR"

python3 "$SCRIPT" "${ARGS[@]}"

echo "Done. Stage manifests written under: $OUTDIR/inputs/manifests"
