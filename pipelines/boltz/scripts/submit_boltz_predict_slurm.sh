#!/bin/bash
#SBATCH --job-name=boltz_sweep
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/boltz/boltz_sweep_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/boltz/boltz_sweep_%j.err
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
ROOT="${ROOT:-$CONFORMERS_ROOT/pipelines/boltz}"
SEED_PDB="${SEED_PDB:-$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb}"
CHAIN_A_ID="${CHAIN_A_ID:-A}"
CHAIN_B_ID="${CHAIN_B_ID:-B}"

SEQ_DIR_DEFAULT="$CONFORMERS_ROOT/data/runs/boltz/avb3/seq"
SEQ_A="${SEQ_A:-$SEQ_DIR_DEFAULT/chain_A.seq}"
SEQ_B="${SEQ_B:-$SEQ_DIR_DEFAULT/chain_B.seq}"

MSA_A="${MSA_A:-$CONFORMERS_ROOT/data/avb3/template_example/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$CONFORMERS_ROOT/data/avb3/template_example/msa/1/non_pairing.a3m}"
TEMPLATE_CIF="${TEMPLATE_CIF:-$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.cif}"

OUTDIR="${OUTDIR:-$CONFORMERS_ROOT/data/runs/boltz/avb3/slurm_boltz_${SLURM_JOB_ID}}"

# Environment split: use Boltz venv, not Protenix venv.
BOLTZ_VENV="${BOLTZ_VENV:-$HOME/scratch/venv_boltz}"
if [[ ! -f "$BOLTZ_VENV/bin/activate" ]]; then
  echo "ERROR: Boltz venv not found at $BOLTZ_VENV" >&2
  echo "Set BOLTZ_VENV to your boltz environment (e.g., ~/scratch/venv_boltz)." >&2
  exit 1
fi
source "$BOLTZ_VENV/bin/activate"

# Extended-state focused defaults.
EXTENDED_ONLY="${EXTENDED_ONLY:-1}"
FORCE_THRESHOLDS="${FORCE_THRESHOLDS:-0.85,0.90,0.95}"
DIFFUSION_SAMPLES="${DIFFUSION_SAMPLES:-12}"
RECYCLES="${RECYCLES:-6}"

echo "CONFORMERS_ROOT=$CONFORMERS_ROOT"
echo "ROOT=$ROOT"
echo "SEED_PDB=$SEED_PDB"
echo "CHAIN_A_ID=$CHAIN_A_ID"
echo "CHAIN_B_ID=$CHAIN_B_ID"
echo "OUTDIR=$OUTDIR"
echo "BOLTZ_VENV=$BOLTZ_VENV"
echo "EXTENDED_ONLY=$EXTENDED_ONLY"
echo "FORCE_THRESHOLDS=$FORCE_THRESHOLDS"
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
    echo "ERROR: seed pdb not found: $SEED_PDB" >&2
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
  --force-thresholds "$FORCE_THRESHOLDS"
  --diffusion-samples "$DIFFUSION_SAMPLES"
  --recycles "$RECYCLES"
  --outdir "$OUTDIR"
)

if [[ -n "$TEMPLATE_CIF" ]]; then
  RUN_ARGS+=(--template-cif "$TEMPLATE_CIF")
fi

if [[ "$EXTENDED_ONLY" = "1" || "$EXTENDED_ONLY" = "true" ]]; then
  RUN_ARGS+=(--extended-only)
fi

bash scripts/run_boltz_predict_sweep.sh "${RUN_ARGS[@]}"
