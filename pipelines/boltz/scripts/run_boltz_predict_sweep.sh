#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_boltz_predict_sweep.sh \
    --sequence-file <fasta_like.txt> \
    --outdir <dir> \
    [--template-cif <path.cif>] \
    [--chain-a-msa <path.a3m>] [--chain-b-msa <path.a3m>] \
    [--force-thresholds <csv>] [--extended-only] \
    [--diffusion-samples 12] [--recycles 6]
EOF
}

SEQUENCE_FILE=""
OUTDIR=""
TEMPLATE_CIF=""
CHAIN_A_MSA=""
CHAIN_B_MSA=""
FORCE_THRESHOLDS="0.60,0.75,0.90"
EXTENDED_ONLY=0
DIFFUSION_SAMPLES=12
RECYCLES=6
NAME_A="Integrin alpha5-Avi"
NAME_B="Integrin beta1-spycatcher"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sequence-file) SEQUENCE_FILE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --template-cif) TEMPLATE_CIF="$2"; shift 2 ;;
    --chain-a-msa) CHAIN_A_MSA="$2"; shift 2 ;;
    --chain-b-msa) CHAIN_B_MSA="$2"; shift 2 ;;
    --force-thresholds) FORCE_THRESHOLDS="$2"; shift 2 ;;
    --extended-only) EXTENDED_ONLY=1; shift ;;
    --diffusion-samples) DIFFUSION_SAMPLES="$2"; shift 2 ;;
    --recycles) RECYCLES="$2"; shift 2 ;;
    --name-a) NAME_A="$2"; shift 2 ;;
    --name-b) NAME_B="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$SEQUENCE_FILE" || -z "$OUTDIR" ]]; then
  echo "Missing required args." >&2
  usage
  exit 1
fi
if [[ ! -f "$SEQUENCE_FILE" ]]; then
  echo "Sequence file not found: $SEQUENCE_FILE" >&2
  exit 1
fi
if [[ -n "$TEMPLATE_CIF" && ! -f "$TEMPLATE_CIF" ]]; then
  echo "Template cif not found: $TEMPLATE_CIF" >&2
  exit 1
fi
if [[ "$EXTENDED_ONLY" -eq 1 && -z "$TEMPLATE_CIF" ]]; then
  echo "--extended-only requires --template-cif" >&2
  exit 1
fi
if ! command -v boltz >/dev/null 2>&1; then
  echo "ERROR: boltz command not found in PATH." >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"
SEQ_DIR="$OUTDIR/sequences"
JOB_DIR="$OUTDIR/jobs"
BOLTZ_OUT="$OUTDIR/boltz_outputs"

mkdir -p "$SEQ_DIR" "$JOB_DIR" "$BOLTZ_OUT"

echo "[1/4] Extracting target integrin sequences"
python3 "$SCRIPT_DIR/extract_integrin_sequences.py" \
  --sequence-file "$SEQUENCE_FILE" \
  --outdir "$SEQ_DIR" \
  --name-a "$NAME_A" \
  --name-b "$NAME_B"

echo "[2/4] Building Boltz sweep jobs"
SWEEP_ARGS=(
  --chain-a-seq-file "$SEQ_DIR/chain_A.seq"
  --chain-b-seq-file "$SEQ_DIR/chain_B.seq"
  --outdir "$JOB_DIR"
  --force-thresholds "$FORCE_THRESHOLDS"
)
if [[ -n "$CHAIN_A_MSA" ]]; then
  SWEEP_ARGS+=(--chain-a-msa "$CHAIN_A_MSA")
fi
if [[ -n "$CHAIN_B_MSA" ]]; then
  SWEEP_ARGS+=(--chain-b-msa "$CHAIN_B_MSA")
fi
if [[ -n "$TEMPLATE_CIF" ]]; then
  SWEEP_ARGS+=(--template-cif "$TEMPLATE_CIF")
fi
if [[ "$EXTENDED_ONLY" -eq 1 ]]; then
  SWEEP_ARGS+=(--extended-only)
fi

python3 "$SCRIPT_DIR/build_boltz_predict_sweep.py" "${SWEEP_ARGS[@]}"

echo "[3/4] Running Boltz predict over sweep"
boltz predict "$JOB_DIR" \
  --out_dir "$BOLTZ_OUT" \
  --diffusion_samples "$DIFFUSION_SAMPLES" \
  --recycling_steps "$RECYCLES" \
  --override

echo "[4/4] Ranking for extension proxy"
python3 "$SCRIPT_DIR/rank_extended.py" \
  --pred-root "$BOLTZ_OUT" \
  --out "$OUTDIR/extended_rank.tsv"

echo "Done: $OUTDIR"
