#!/bin/bash
#SBATCH --job-name=afcluster_boltz
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/afcluster/afcluster_boltz_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/afcluster/afcluster_boltz_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

set -euo pipefail

# AVB3 stream defaults for extended-conformation search.
CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
ROOT="${ROOT:-$CONFORMERS_ROOT/pipelines/afcluster}"

SEED_PDB="${SEED_PDB:-$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb}"
CHAIN_A_ID="${CHAIN_A_ID:-A}"
CHAIN_B_ID="${CHAIN_B_ID:-B}"

SEQ_DIR_DEFAULT="$CONFORMERS_ROOT/data/runs/afcluster/avb3/seq"
SEQ_A="${SEQ_A:-$SEQ_DIR_DEFAULT/chain_A.seq}"
SEQ_B="${SEQ_B:-$SEQ_DIR_DEFAULT/chain_B.seq}"

MSA_A="${MSA_A:-$CONFORMERS_ROOT/data/avb3/template_example/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$CONFORMERS_ROOT/data/avb3/template_example/msa/1/non_pairing.a3m}"
TEMPLATE_CIF="${TEMPLATE_CIF:-$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.cif}"

OUTDIR="${OUTDIR:-$CONFORMERS_ROOT/data/runs/afcluster/avb3/slurm_afcluster_boltz_${SLURM_JOB_ID}}"
TOP_A="${TOP_A:-8}"
TOP_B="${TOP_B:-8}"

# Backend can be boltzgen (default) or boltz.
BACKEND="${BACKEND:-boltzgen}"
PROTOCOL="${PROTOCOL:-protein-anything}"
NUM_DESIGNS="${NUM_DESIGNS:-10000}"
BUDGET="${BUDGET:-200}"
BOLTZGEN_BIN="${BOLTZGEN_BIN:-boltzgen}"
BOLTZGEN_EXTRA="${BOLTZGEN_EXTRA:-}"

# Legacy boltz-only params.
DIFFUSION_SAMPLES="${DIFFUSION_SAMPLES:-8}"
RECYCLES="${RECYCLES:-3}"

# Environment split: AFCluster/BoltzGen pipeline uses boltz venv.
BOLTZ_VENV="${BOLTZ_VENV:-$HOME/scratch/venv_boltz}"
if [[ ! -f "$BOLTZ_VENV/bin/activate" ]]; then
  echo "ERROR: Boltz venv not found at $BOLTZ_VENV" >&2
  echo "Set BOLTZ_VENV to your boltz environment (e.g., ~/scratch/venv_boltz)." >&2
  exit 1
fi
source "$BOLTZ_VENV/bin/activate"

echo "CONFORMERS_ROOT=$CONFORMERS_ROOT"
echo "ROOT=$ROOT"
echo "SEED_PDB=$SEED_PDB"
echo "CHAIN_A_ID=$CHAIN_A_ID"
echo "CHAIN_B_ID=$CHAIN_B_ID"
echo "OUTDIR=$OUTDIR"
echo "TOP_A=$TOP_A"
echo "TOP_B=$TOP_B"
echo "BACKEND=$BACKEND"
echo "PROTOCOL=$PROTOCOL"
echo "NUM_DESIGNS=$NUM_DESIGNS"
echo "BUDGET=$BUDGET"
echo "BOLTZGEN_BIN=$BOLTZGEN_BIN"
echo "BOLTZ_VENV=$BOLTZ_VENV"
python - <<'PY'
import gemmi
print(f"gemmi={gemmi.__version__}")
PY

if [[ ! -f "$SEQ_A" || ! -f "$SEQ_B" ]]; then
  EXTRACT_SCRIPT="$CONFORMERS_ROOT/pipelines/afcluster/scripts/extract_chain_sequences_from_pdb.py"
  if [[ ! -f "$EXTRACT_SCRIPT" ]]; then
    echo "ERROR: sequence extractor not found: $EXTRACT_SCRIPT" >&2
    exit 1
  fi
  if [[ ! -f "$SEED_PDB" ]]; then
    echo "ERROR: seed pdb not found for seq generation: $SEED_PDB" >&2
    exit 1
  fi
  mkdir -p "$(dirname "$SEQ_A")"
  python3 "$EXTRACT_SCRIPT" \
    --pdb "$SEED_PDB" \
    --outdir "$(dirname "$SEQ_A")" \
    --chain-a "$CHAIN_A_ID" \
    --chain-b "$CHAIN_B_ID"
fi

cd "$ROOT"
RUN_ARGS=(
  --chain-a-seq-file "$SEQ_A"
  --chain-b-seq-file "$SEQ_B"
  --chain-a-msa "$MSA_A"
  --chain-b-msa "$MSA_B"
  --template-cif "$TEMPLATE_CIF"
  --outdir "$OUTDIR"
  --top-a "$TOP_A"
  --top-b "$TOP_B"
  --backend "$BACKEND"
)

if [[ "$BACKEND" == "boltzgen" ]]; then
  RUN_ARGS+=(
    --protocol "$PROTOCOL"
    --num-designs "$NUM_DESIGNS"
    --budget "$BUDGET"
    --boltzgen-bin "$BOLTZGEN_BIN"
  )
  if [[ -n "$BOLTZGEN_EXTRA" ]]; then
    RUN_ARGS+=(--boltzgen-extra "$BOLTZGEN_EXTRA")
  fi
else
  RUN_ARGS+=(
    --diffusion-samples "$DIFFUSION_SAMPLES"
    --recycles "$RECYCLES"
  )
fi

bash scripts/run_afcluster_pipeline.sh "${RUN_ARGS[@]}"
