#!/bin/bash
#SBATCH -J avb3_proteinttt
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/proteinttt_%j.log
#SBATCH --error=logs/avb3-conformers/proteinttt_%j.err
set -euo pipefail

# ProteinTTT: test-time training of ESMFold on AVB3 integrin sequence.
# Fine-tunes ESMFold's language model on the target sequence via masked
# language modeling, then predicts structure with improved pLDDT.
# Can also score pulled conformers by comparing TTT predictions to them.

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
PROTEINTTT_ROOT="${PROTEINTTT_ROOT:-$HOME/scratch/ProteinTTT}"
WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/proteinttt}"

# MSA path for AVB3 (if available, improves TTT quality)
MSA_PATH="${MSA_PATH:-}"

# Reference PDB to evaluate against (optional)
REFERENCE_PDB="${REFERENCE_PDB:-}"

# TTT hyperparameters
TTT_STEPS="${TTT_STEPS:-30}"
TTT_LR="${TTT_LR:-4e-4}"
TTT_LORA_RANK="${TTT_LORA_RANK:-8}"
TTT_BATCH_SIZE="${TTT_BATCH_SIZE:-4}"
TTT_CROP_SIZE="${TTT_CROP_SIZE:-1024}"

module load cuda

mkdir -p "$WORK_DIR" "$CONFORMERS_ROOT/logs/avb3-conformers"

# Set up ProteinTTT environment
# Use Python 3.12 from the spack module (3.9 default is too old for fair-esm)
module load python/3.12 2>/dev/null || true

VENV_DIR="$HOME/scratch/venv_proteinttt"
SETUP_MARKER="$VENV_DIR/.setup_complete"
if [[ ! -f "$SETUP_MARKER" ]]; then
    echo "Creating/repairing ProteinTTT venv..."
    rm -rf "$VENV_DIR"
    python3 -m venv "$VENV_DIR"
    source "$VENV_DIR/bin/activate"
    pip install --upgrade pip
    pip install torch --index-url https://download.pytorch.org/whl/cu121
    pip install fair-esm omegaconf pandas biopython requests
    pip install -e "$PROTEINTTT_ROOT"
    touch "$SETUP_MARKER"
else
    source "$VENV_DIR/bin/activate"
fi

echo "=== ProteinTTT on AVB3 ==="
echo "Work dir: $WORK_DIR"
echo "Steps: $TTT_STEPS, LR: $TTT_LR, LoRA rank: $TTT_LORA_RANK"
echo ""

python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/run_proteinttt.py" \
    --conformers-root "$CONFORMERS_ROOT" \
    --proteinttt-root "$PROTEINTTT_ROOT" \
    --work-dir "$WORK_DIR" \
    --ttt-steps "$TTT_STEPS" \
    --ttt-lr "$TTT_LR" \
    --ttt-lora-rank "$TTT_LORA_RANK" \
    --ttt-batch-size "$TTT_BATCH_SIZE" \
    --ttt-crop-size "$TTT_CROP_SIZE" \
    ${MSA_PATH:+--msa-path "$MSA_PATH"} \
    ${REFERENCE_PDB:+--reference-pdb "$REFERENCE_PDB"}

echo "ProteinTTT complete. Results in $WORK_DIR"
