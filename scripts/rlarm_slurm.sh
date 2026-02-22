#!/bin/bash
# This slurm script is for running RL-Arm training on a SLURM cluster.

#SBATCH --job-name=rlarm_train
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/logs/rlarm_train_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/logs/rlarm_train_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:RTX_6000:1
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

# Load necessary modules
module load cuda
module load python/3.12

# Activate the Python virtual environment
source ~/scratch/venv_planar/bin/activate

# Assign arguments to variables
NUM_LINKS=$1
STEPS=$2
JERK_PENALTY=$3
JERK_THRESHOLD=$4
MODEL_DIM=$5
OUTPUT_DIR=$6

# Run the Python script with the passed arguments
srun python demo_v7.1.py \
    --train \
    --wandb \
    --steps "$STEPS" \
    --num-links "$NUM_LINKS" \
    --output-dir "$OUTPUT_DIR" \
    --jerk-penalty "$JERK_PENALTY" \
    --jerk-threshold "$JERK_THRESHOLD" \
    --model-d-model "$MODEL_DIM"
