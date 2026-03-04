#!/bin/bash
#SBATCH --job-name=boltz_sweep
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/boltz_sweep_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/boltz_sweep_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

set -euo pipefail

# Defaults aligned with scratch/conformers layout used by Protenix jobs.
ROOT="${ROOT:-$HOME/scratch/conformers/codex-Boltz}"
SEQUENCE_FILE="${SEQUENCE_FILE:-$HOME/scratch/conformers/data/sequences_updated}"
TEMPLATE_CIF="${TEMPLATE_CIF:-$HOME/scratch/conformers/data/template_example/seed_090_frame_000.cif}"
MSA_A="${MSA_A:-$HOME/scratch/conformers/data/seed_090_frame_000/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$HOME/scratch/conformers/data/seed_090_frame_000/msa/1/non_pairing.a3m}"
OUTDIR="${OUTDIR:-$HOME/scratch/conformers/data/boltz_outputs/slurm_boltz_${SLURM_JOB_ID}}"

# Extended-state focused defaults.
EXTENDED_ONLY="${EXTENDED_ONLY:-1}"
FORCE_THRESHOLDS="${FORCE_THRESHOLDS:-0.85,0.90,0.95}"
DIFFUSION_SAMPLES="${DIFFUSION_SAMPLES:-12}"
RECYCLES="${RECYCLES:-6}"

echo "ROOT=$ROOT"
echo "OUTDIR=$OUTDIR"
echo "EXTENDED_ONLY=$EXTENDED_ONLY"
echo "FORCE_THRESHOLDS=$FORCE_THRESHOLDS"

cd "$ROOT"
RUN_ARGS=(
  --sequence-file "$SEQUENCE_FILE"
  --template-cif "$TEMPLATE_CIF"
  --chain-a-msa "$MSA_A"
  --chain-b-msa "$MSA_B"
  --force-thresholds "$FORCE_THRESHOLDS"
  --diffusion-samples "$DIFFUSION_SAMPLES"
  --recycles "$RECYCLES"
  --outdir "$OUTDIR"
)

if [[ "$EXTENDED_ONLY" = "1" || "$EXTENDED_ONLY" = "true" ]]; then
  RUN_ARGS+=(--extended-only)
fi

bash scripts/run_boltz_predict_sweep.sh "${RUN_ARGS[@]}"
