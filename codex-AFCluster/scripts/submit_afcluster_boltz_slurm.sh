#!/bin/bash
#SBATCH --job-name=afcluster_boltz
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

# Edit these paths for your cluster.
ROOT="${ROOT:-$HOME/scratch/Protenix/codex-AFCluster}"
SEQ_A="${SEQ_A:-$ROOT/../data/template_example/workflow_outputs/inputs/chain_A.seq}"
SEQ_B="${SEQ_B:-$ROOT/../data/template_example/workflow_outputs/inputs/chain_B.seq}"
MSA_A="${MSA_A:-$ROOT/../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$ROOT/../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/1/non_pairing.a3m}"
TEMPLATE_CIF="${TEMPLATE_CIF:-$ROOT/../data/template_example/seed_090_frame_000.cif}"
OUTDIR="${OUTDIR:-$ROOT/runs/slurm_afcluster_boltz_${SLURM_JOB_ID}}"

echo "ROOT=$ROOT"
echo "OUTDIR=$OUTDIR"

cd "$ROOT"
bash scripts/run_afcluster_pipeline.sh \
  --chain-a-seq-file "$SEQ_A" \
  --chain-b-seq-file "$SEQ_B" \
  --chain-a-msa "$MSA_A" \
  --chain-b-msa "$MSA_B" \
  --template-cif "$TEMPLATE_CIF" \
  --outdir "$OUTDIR"
