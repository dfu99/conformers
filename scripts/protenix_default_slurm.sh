#!/bin/bash

#SBATCH --job-name=protenix_template
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/protenix_template_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/protenix_template_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

# Load necessary modules
module load cuda
module load python/3.12

# Activate virtual environment (adjust if you have a protenix-specific venv)
source ~/scratch/venv_protenix/bin/activate

# --- Paths (adjust as needed) ---
PROTENIX_ROOT="$HOME/scratch/Protenix"
MMCIF_SOURCE_DIR="${PROTENIX_MMCIF_DIR:-$HOME/scratch/conformers/mmcif}"
MMCIF_TARGET_DIR="$PROTENIX_ROOT/mmcif"

# Protenix requires template mmcif DB under: $PROTENIX_ROOT_DIR/mmcif
export PROTENIX_ROOT_DIR="$PROTENIX_ROOT"

if [ ! -d "$MMCIF_SOURCE_DIR" ]; then
    echo "ERROR: mmcif directory not found: $MMCIF_SOURCE_DIR"
    echo "Set PROTENIX_MMCIF_DIR to an existing mmcif path before running sbatch."
    exit 1
fi

if [ "$MMCIF_SOURCE_DIR" != "$MMCIF_TARGET_DIR" ]; then
    if [ -e "$MMCIF_TARGET_DIR" ] && [ ! -L "$MMCIF_TARGET_DIR" ]; then
        echo "ERROR: $MMCIF_TARGET_DIR exists and is not a symlink."
        echo "Rename/remove it, or point PROTENIX_MMCIF_DIR to $MMCIF_TARGET_DIR."
        exit 1
    fi
    ln -sfn "$MMCIF_SOURCE_DIR" "$MMCIF_TARGET_DIR"
fi

if [ ! -d "$MMCIF_TARGET_DIR" ]; then
    echo "ERROR: mmcif target directory unavailable: $MMCIF_TARGET_DIR"
    exit 1
fi

cd "$PROTENIX_ROOT"
echo "PROTENIX_ROOT_DIR: $PROTENIX_ROOT_DIR"
echo "mmcif source:      $MMCIF_SOURCE_DIR"
echo "mmcif target:      $MMCIF_TARGET_DIR"

N_sample=5
N_step=200
N_cycle=10
seed=103
input_json_path="./examples/examples_with_template/example_9fm7.json"
dump_dir="./test_outputs/sh/output_m_9fm7"
model_name="protenix_base_default_v1.0.0"

srun python3 runner/inference.py \
    --model_name ${model_name} \
    --seeds ${seed} \
    --dump_dir ${dump_dir} \
    --input_json_path ${input_json_path} \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step} \
    --triangle_attention "cuequivariance" \
    --use_seeds_in_json true \
    --triangle_multiplicative "cuequivariance" \
    --use_template true
