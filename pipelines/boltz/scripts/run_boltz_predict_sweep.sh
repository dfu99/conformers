#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF2'
Usage:
  run_boltz_predict_sweep.sh \
    --outdir <dir> \
    [--chain-a-seq-file <path>] [--chain-b-seq-file <path>] \
    [--sequence-file <fasta_like.txt> --name-a <header> --name-b <header>] \
    [--template-cif <path.cif>] \
    [--chain-a-msa <path.a3m>] [--chain-b-msa <path.a3m>] \
    [--force-thresholds <csv>] [--extended-only] \
    [--diffusion-samples 12] [--recycles 6]

Notes:
- Provide either:
  1) --chain-a-seq-file + --chain-b-seq-file
  2) --sequence-file (then sequences are extracted by header name)
EOF2
}

SEQUENCE_FILE=""
CHAIN_A_SEQ_FILE=""
CHAIN_B_SEQ_FILE=""
OUTDIR=""
TEMPLATE_CIF=""
CHAIN_A_MSA=""
CHAIN_B_MSA=""
FORCE_THRESHOLDS="0.60,0.75,0.90"
EXTENDED_ONLY=0
DIFFUSION_SAMPLES=12
RECYCLES=6
NAME_A="Integrin alphaV"
NAME_B="Integrin beta3"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sequence-file) SEQUENCE_FILE="$2"; shift 2 ;;
    --chain-a-seq-file) CHAIN_A_SEQ_FILE="$2"; shift 2 ;;
    --chain-b-seq-file) CHAIN_B_SEQ_FILE="$2"; shift 2 ;;
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

if [[ -z "$OUTDIR" ]]; then
  echo "Missing required --outdir." >&2
  usage
  exit 1
fi

if [[ -n "$CHAIN_A_SEQ_FILE" || -n "$CHAIN_B_SEQ_FILE" ]]; then
  if [[ -z "$CHAIN_A_SEQ_FILE" || -z "$CHAIN_B_SEQ_FILE" ]]; then
    echo "Provide both --chain-a-seq-file and --chain-b-seq-file." >&2
    exit 1
  fi
  if [[ ! -f "$CHAIN_A_SEQ_FILE" ]]; then
    echo "chain A seq file not found: $CHAIN_A_SEQ_FILE" >&2
    exit 1
  fi
  if [[ ! -f "$CHAIN_B_SEQ_FILE" ]]; then
    echo "chain B seq file not found: $CHAIN_B_SEQ_FILE" >&2
    exit 1
  fi
else
  if [[ -z "$SEQUENCE_FILE" ]]; then
    echo "Provide --sequence-file or both chain seq files." >&2
    usage
    exit 1
  fi
  if [[ ! -f "$SEQUENCE_FILE" ]]; then
    echo "Sequence file not found: $SEQUENCE_FILE" >&2
    exit 1
  fi
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

if [[ -n "$CHAIN_A_SEQ_FILE" && -n "$CHAIN_B_SEQ_FILE" ]]; then
  echo "[1/4] Using provided chain sequence files"
  cp "$CHAIN_A_SEQ_FILE" "$SEQ_DIR/chain_A.seq"
  cp "$CHAIN_B_SEQ_FILE" "$SEQ_DIR/chain_B.seq"
else
  echo "[1/4] Extracting target chain sequences from FASTA-like input"
  python3 "$SCRIPT_DIR/extract_integrin_sequences.py" \
    --sequence-file "$SEQUENCE_FILE" \
    --outdir "$SEQ_DIR" \
    --name-a "$NAME_A" \
    --name-b "$NAME_B"
fi

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
