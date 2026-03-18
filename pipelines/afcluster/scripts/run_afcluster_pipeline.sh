#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF2'
Usage:
  run_afcluster_pipeline.sh \
    --chain-a-seq-file <path> \
    --chain-b-seq-file <path> \
    --chain-a-msa <path.a3m> \
    --chain-b-msa <path.a3m> \
    --outdir <dir> \
    [--template-cif <path.cif>] \
    [--top-a 8] [--top-b 8] \
    [--backend boltzgen|boltz] \
    [--protocol protein-anything] [--num-designs 10000] [--budget 200] \
    [--diffusion-samples 8] [--recycles 3]

Notes:
- AFCluster always clusters the per-chain MSA first.
- Backend `boltzgen` uses: boltzgen run <yaml> --output ... --protocol ... --num_designs ... --budget ...
- Backend `boltz` keeps legacy boltz predict behavior.
EOF2
}

CHAIN_A_SEQ_FILE=""
CHAIN_B_SEQ_FILE=""
CHAIN_A_MSA=""
CHAIN_B_MSA=""
OUTDIR=""
TEMPLATE_CIF=""
TOP_A=8
TOP_B=8

BACKEND="boltzgen"
PROTOCOL="protein-anything"
NUM_DESIGNS=10000
BUDGET=200
BOLTZGEN_BIN="boltzgen"
BOLTZGEN_EXTRA=""

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
    --backend) BACKEND="$2"; shift 2 ;;
    --protocol) PROTOCOL="$2"; shift 2 ;;
    --num-designs) NUM_DESIGNS="$2"; shift 2 ;;
    --budget) BUDGET="$2"; shift 2 ;;
    --boltzgen-bin) BOLTZGEN_BIN="$2"; shift 2 ;;
    --boltzgen-extra) BOLTZGEN_EXTRA="$2"; shift 2 ;;
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

if [[ "$BACKEND" != "boltzgen" && "$BACKEND" != "boltz" ]]; then
  echo "Unsupported --backend: $BACKEND (expected boltzgen or boltz)" >&2
  exit 1
fi

# BoltzGen uses huggingface_hub for checkpoints and inference data. Keep those
# downloads on scratch by default to avoid exhausting the home-directory quota.
HF_CACHE_ROOT="${HF_HOME:-$HOME/scratch/.cache/huggingface}"
export HF_HOME="$HF_CACHE_ROOT"
export HF_HUB_CACHE="${HF_HUB_CACHE:-$HF_HOME/hub}"
export HF_XET_CACHE="${HF_XET_CACHE:-$HF_HOME/xet}"
mkdir -p "$HF_HUB_CACHE" "$HF_XET_CACHE"

if [[ "$BACKEND" == "boltzgen" ]]; then
  if ! command -v "$BOLTZGEN_BIN" >/dev/null 2>&1; then
    echo "ERROR: $BOLTZGEN_BIN not found in PATH." >&2
    exit 1
  fi
else
  if ! command -v boltz >/dev/null 2>&1; then
    echo "ERROR: boltz command not found in PATH." >&2
    exit 1
  fi
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"
CLUSTER_A_DIR="$OUTDIR/clusters/A"
CLUSTER_B_DIR="$OUTDIR/clusters/B"
JOB_DIR="$OUTDIR/boltz_jobs"
BOLTZ_OUT="$OUTDIR/boltz_outputs"

mkdir -p "$CLUSTER_A_DIR" "$CLUSTER_B_DIR" "$JOB_DIR" "$BOLTZ_OUT"
echo "HF_HOME=$HF_HOME"

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

echo "[3/4] Building YAML job set"
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

if [[ "$BACKEND" == "boltzgen" ]]; then
  echo "[4/4] Running BoltzGen over generated YAML jobs"
  mapfile -t SPECS < <(find "$JOB_DIR" -maxdepth 1 -type f -name '*.yaml' | sort)
  if [[ ${#SPECS[@]} -eq 0 ]]; then
    echo "ERROR: no YAML specs found in $JOB_DIR" >&2
    exit 1
  fi

  for spec in "${SPECS[@]}"; do
    job_name="$(basename "${spec%.yaml}")"
    job_out="$BOLTZ_OUT/$job_name"
    mkdir -p "$job_out"

    cmd=(
      "$BOLTZGEN_BIN" run "$spec"
      --output "$job_out"
      --protocol "$PROTOCOL"
      --num_designs "$NUM_DESIGNS"
      --budget "$BUDGET"
    )
    if [[ -n "$BOLTZGEN_EXTRA" ]]; then
      # shellcheck disable=SC2206
      extra_parts=($BOLTZGEN_EXTRA)
      cmd+=("${extra_parts[@]}")
    fi

    echo "Running BoltzGen job: $job_name"
    "${cmd[@]}"
  done
else
  echo "[4/4] Running Boltz predict over generated YAML jobs"
  boltz predict "$JOB_DIR" \
    --out_dir "$BOLTZ_OUT" \
    --diffusion_samples "$DIFFUSION_SAMPLES" \
    --recycling_steps "$RECYCLES" \
    --override
fi

echo "[post] Ranking structures for extension proxy"
python3 "$SCRIPT_DIR/rank_extended.py" \
  --pred-root "$BOLTZ_OUT" \
  --out "$OUTDIR/extended_rank.tsv"

echo "Done: $OUTDIR"
