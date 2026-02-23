#!/bin/bash
#SBATCH --job-name=afcluster_boltz
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/afcluster_boltz_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/afcluster_boltz_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

set -euo pipefail

# Defaults aligned with scratch/conformers layout used by Protenix jobs.
ROOT="${ROOT:-$HOME/scratch/conformers/codex-AFCluster}"
SEQ_A="${SEQ_A:-$HOME/scratch/conformers/data/template_example/workflow_outputs/inputs/chain_A.seq}"
SEQ_B="${SEQ_B:-$HOME/scratch/conformers/data/template_example/workflow_outputs/inputs/chain_B.seq}"
MSA_A="${MSA_A:-$HOME/scratch/conformers/data/seed_090_frame_000/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$HOME/scratch/conformers/data/seed_090_frame_000/msa/1/non_pairing.a3m}"
TEMPLATE_CIF="${TEMPLATE_CIF:-$HOME/scratch/conformers/data/template_example/seed_090_frame_000.cif}"
OUTDIR="${OUTDIR:-$HOME/scratch/conformers/data/afcluster_boltz_outputs/slurm_afcluster_boltz_${SLURM_JOB_ID}}"

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
