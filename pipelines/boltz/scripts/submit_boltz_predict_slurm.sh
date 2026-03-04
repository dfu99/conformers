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

# A5B1 stream only: data/a5b1 + data/runs/a5b1
CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
ROOT="${ROOT:-$CONFORMERS_ROOT/pipelines/boltz}"
SEQUENCE_FILE="${SEQUENCE_FILE:-$CONFORMERS_ROOT/data/a5b1/sequences/sequences_updated}"

# Use heterodimer-only MSA stream from prior Protenix run.
MSA_A="${MSA_A:-$CONFORMERS_ROOT/data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$CONFORMERS_ROOT/data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/1/non_pairing.a3m}"

# Default template from A5B1 stream.
TEMPLATE_CIF="${TEMPLATE_CIF:-$CONFORMERS_ROOT/data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/seed_101/predictions/integrin_alpha5_beta1_sample_0.cif}"

# Boltz stream output root.
OUTDIR="${OUTDIR:-$CONFORMERS_ROOT/data/runs/boltz/slurm_boltz_${SLURM_JOB_ID}}"

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
echo "SEQUENCE_FILE=$SEQUENCE_FILE"
echo "OUTDIR=$OUTDIR"
echo "BOLTZ_VENV=$BOLTZ_VENV"
echo "EXTENDED_ONLY=$EXTENDED_ONLY"
echo "FORCE_THRESHOLDS=$FORCE_THRESHOLDS"
python - <<'PY'
import gemmi
print(f"gemmi={gemmi.__version__}")
PY

cd "$ROOT"
RUN_ARGS=(
  --sequence-file "$SEQUENCE_FILE"
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
