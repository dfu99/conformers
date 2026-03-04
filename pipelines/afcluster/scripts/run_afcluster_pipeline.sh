#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_afcluster_pipeline.sh \
    --chain-a-seq-file <path> \
    --chain-b-seq-file <path> \
    --chain-a-msa <path.a3m> \
    --chain-b-msa <path.a3m> \
    --outdir <dir> \
    [--template-cif <path.cif>] \
    [--top-a 8] [--top-b 8] \
    [--diffusion-samples 8] [--recycles 3]
EOF
}

CHAIN_A_SEQ_FILE=""
CHAIN_B_SEQ_FILE=""
CHAIN_A_MSA=""
CHAIN_B_MSA=""
OUTDIR=""
TEMPLATE_CIF=""
TOP_A=8
TOP_B=8
DIFFUSION_SAMPLES=8
RECYCLES=3

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chain-a-seq-file) CHAIN_A_SEQ_FILE="$2"; shift 2 ;;
    --chain-b-seq-file) CHAIN_B_SEQ_FILE="$2"; shift 2 ;;
    --chain-a-msa) CHAIN_A_MSA="$2"; shift 2 ;;
    --chain-b-msa) CHAIN_B_MSA="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --template-cif) TEMPLATE_CIF="$2"; shift 2 ;;
    --top-a) TOP_A="$2"; shift 2 ;;
    --top-b) TOP_B="$2"; shift 2 ;;
    --diffusion-samples) DIFFUSION_SAMPLES="$2"; shift 2 ;;
    --recycles) RECYCLES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

for req in "$CHAIN_A_SEQ_FILE" "$CHAIN_B_SEQ_FILE" "$CHAIN_A_MSA" "$CHAIN_B_MSA"; do
  if [[ -z "$req" ]]; then
    echo "Missing required args." >&2
    usage
    exit 1
  fi
done
if [[ -z "$OUTDIR" ]]; then
  echo "Missing --outdir" >&2
  usage
  exit 1
fi

for f in "$CHAIN_A_SEQ_FILE" "$CHAIN_B_SEQ_FILE" "$CHAIN_A_MSA" "$CHAIN_B_MSA"; do
  if [[ ! -f "$f" ]]; then
    echo "File not found: $f" >&2
    exit 1
  fi
done
if [[ -n "$TEMPLATE_CIF" && ! -f "$TEMPLATE_CIF" ]]; then
  echo "Template cif not found: $TEMPLATE_CIF" >&2
  exit 1
fi

if ! command -v boltz >/dev/null 2>&1; then
  echo "ERROR: boltz command not found in PATH." >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"
CLUSTER_A_DIR="$OUTDIR/clusters/A"
CLUSTER_B_DIR="$OUTDIR/clusters/B"
JOB_DIR="$OUTDIR/boltz_jobs"
BOLTZ_OUT="$OUTDIR/boltz_outputs"

mkdir -p "$CLUSTER_A_DIR" "$CLUSTER_B_DIR" "$JOB_DIR" "$BOLTZ_OUT"

echo "[1/4] Clustering chain A MSA"
python3 "$SCRIPT_DIR/cluster_chain_msa.py" \
  --input-a3m "$CHAIN_A_MSA" \
  --outdir "$CLUSTER_A_DIR" \
  --keep-top "$TOP_A"

echo "[2/4] Clustering chain B MSA"
python3 "$SCRIPT_DIR/cluster_chain_msa.py" \
  --input-a3m "$CHAIN_B_MSA" \
  --outdir "$CLUSTER_B_DIR" \
  --keep-top "$TOP_B"

echo "[3/4] Building Boltz YAML job set"
JOB_ARGS=(
  --chain-a-seq-file "$CHAIN_A_SEQ_FILE"
  --chain-b-seq-file "$CHAIN_B_SEQ_FILE"
  --chain-a-cluster-dir "$CLUSTER_A_DIR"
  --chain-b-cluster-dir "$CLUSTER_B_DIR"
  --top-a "$TOP_A"
  --top-b "$TOP_B"
  --outdir "$JOB_DIR"
  --job-prefix "integrin_ab_afcluster"
  --include-empty-msa-control
)
if [[ -n "$TEMPLATE_CIF" ]]; then
  JOB_ARGS+=(--template-cif "$TEMPLATE_CIF")
fi

python3 "$SCRIPT_DIR/make_boltz_jobs_from_clusters.py" "${JOB_ARGS[@]}"

echo "[4/4] Running Boltz on generated YAML jobs"
boltz predict "$JOB_DIR" \
  --out_dir "$BOLTZ_OUT" \
  --diffusion_samples "$DIFFUSION_SAMPLES" \
  --recycling_steps "$RECYCLES" \
  --override

echo "[post] Ranking structures for extension proxy"
python3 "$SCRIPT_DIR/rank_extended.py" \
  --pred-root "$BOLTZ_OUT" \
  --out "$OUTDIR/extended_rank.tsv"

echo "Done: $OUTDIR"
