#!/bin/bash
#SBATCH --job-name=boltz_sweep
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

# Edit these paths for your cluster.
ROOT="${ROOT:-$HOME/scratch/Protenix/codex-Boltz}"
SEQUENCE_FILE="${SEQUENCE_FILE:-$ROOT/../data/sequences_updated}"
TEMPLATE_CIF="${TEMPLATE_CIF:-$ROOT/../data/template_example/seed_090_frame_000.cif}"
MSA_A="${MSA_A:-$ROOT/../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$ROOT/../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/1/non_pairing.a3m}"
OUTDIR="${OUTDIR:-$ROOT/runs/slurm_boltz_${SLURM_JOB_ID}}"

echo "ROOT=$ROOT"
echo "OUTDIR=$OUTDIR"

cd "$ROOT"
bash scripts/run_boltz_predict_sweep.sh \
  --sequence-file "$SEQUENCE_FILE" \
  --template-cif "$TEMPLATE_CIF" \
  --chain-a-msa "$MSA_A" \
  --chain-b-msa "$MSA_B" \
  --outdir "$OUTDIR"
