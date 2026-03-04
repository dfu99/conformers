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

# A5B1 stream only: data/a5b1 + data/runs/a5b1
CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
ROOT="${ROOT:-$CONFORMERS_ROOT/pipelines/afcluster}"

SEQUENCE_FILE="${SEQUENCE_FILE:-$CONFORMERS_ROOT/data/a5b1/sequences/sequences_updated}"
SEQ_DIR_DEFAULT="$CONFORMERS_ROOT/data/runs/afcluster/seq"
SEQ_A="${SEQ_A:-$SEQ_DIR_DEFAULT/chain_A.seq}"
SEQ_B="${SEQ_B:-$SEQ_DIR_DEFAULT/chain_B.seq}"

# Use heterodimer-only MSA stream from prior Protenix run.
MSA_A="${MSA_A:-$CONFORMERS_ROOT/data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/0/non_pairing.a3m}"
MSA_B="${MSA_B:-$CONFORMERS_ROOT/data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/1/non_pairing.a3m}"

# Default template from A5B1 stream.
TEMPLATE_CIF="${TEMPLATE_CIF:-$CONFORMERS_ROOT/data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/seed_101/predictions/integrin_alpha5_beta1_sample_0.cif}"

# AFCluster stream output root.
OUTDIR="${OUTDIR:-$CONFORMERS_ROOT/data/runs/afcluster/slurm_afcluster_boltz_${SLURM_JOB_ID}}"

# Environment split: AFCluster+Boltz pipeline uses boltz venv.
BOLTZ_VENV="${BOLTZ_VENV:-$HOME/scratch/venv_boltz}"
if [[ ! -f "$BOLTZ_VENV/bin/activate" ]]; then
  echo "ERROR: Boltz venv not found at $BOLTZ_VENV" >&2
  echo "Set BOLTZ_VENV to your boltz environment (e.g., ~/scratch/venv_boltz)." >&2
  exit 1
fi
source "$BOLTZ_VENV/bin/activate"

echo "CONFORMERS_ROOT=$CONFORMERS_ROOT"
echo "ROOT=$ROOT"
echo "SEQUENCE_FILE=$SEQUENCE_FILE"
echo "OUTDIR=$OUTDIR"
echo "BOLTZ_VENV=$BOLTZ_VENV"
python - <<'PY'
import gemmi
print(f"gemmi={gemmi.__version__}")
PY

# Auto-generate chain seq files from sequences_updated if missing.
if [[ ! -f "$SEQ_A" || ! -f "$SEQ_B" ]]; then
  EXTRACT_SCRIPT="$CONFORMERS_ROOT/pipelines/boltz/scripts/extract_integrin_sequences.py"
  if [[ ! -f "$EXTRACT_SCRIPT" ]]; then
    echo "ERROR: sequence extractor not found: $EXTRACT_SCRIPT" >&2
    exit 1
  fi
  if [[ ! -f "$SEQUENCE_FILE" ]]; then
    echo "ERROR: sequence file not found for seq generation: $SEQUENCE_FILE" >&2
    exit 1
  fi
  mkdir -p "$(dirname "$SEQ_A")"
  python3 "$EXTRACT_SCRIPT" \
    --sequence-file "$SEQUENCE_FILE" \
    --outdir "$(dirname "$SEQ_A")"
fi

cd "$ROOT"
bash scripts/run_afcluster_pipeline.sh \
  --chain-a-seq-file "$SEQ_A" \
  --chain-b-seq-file "$SEQ_B" \
  --chain-a-msa "$MSA_A" \
  --chain-b-msa "$MSA_B" \
  --template-cif "$TEMPLATE_CIF" \
  --outdir "$OUTDIR"
