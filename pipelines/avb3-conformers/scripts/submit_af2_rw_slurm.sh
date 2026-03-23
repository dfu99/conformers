#!/bin/bash
#SBATCH -J avb3_af2_rw
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/af2_rw_%j.log
#SBATCH --error=logs/avb3-conformers/af2_rw_%j.err
set -euo pipefail

# AF2 RandomWalk: weight perturbation for conformational diversity.
# Runs AF2 multiple times with Gaussian noise added to model weights.

module purge
module load alphafold/2.3.2
export PYTHONNOUSERSITE=1

CONFORMERS_ROOT="$HOME/scratch/conformers"
FASTA="$CONFORMERS_ROOT/data/runs/avb3/af2_msa_validation/fasta/avb3_multimer.fasta"
OUTPUT_DIR="$CONFORMERS_ROOT/data/runs/avb3/af2_random_walk_${SLURM_JOB_ID}"
DATA_DIR="/storage/coda1/d-pace_community/0/alphafold/alphafold_2.3.2_data"

N_PERTURBATIONS="${N_PERTURBATIONS:-5}"
NOISE_SCALE="${NOISE_SCALE:-0.01}"

mkdir -p "$OUTPUT_DIR" "$CONFORMERS_ROOT/logs/avb3-conformers"

echo "=== AF2 RandomWalk ==="
echo "Perturbations: $N_PERTURBATIONS, Noise scale: $NOISE_SCALE"
echo "FASTA: $FASTA"
echo "Output: $OUTPUT_DIR"

# Generate FASTA if it doesn't exist
if [[ ! -f "$FASTA" ]]; then
    mkdir -p "$(dirname "$FASTA")"
    python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/extract_avb3_fasta.py" \
        "$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb" \
        "$(dirname "$FASTA")"
fi

python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/af2_random_walk.py" \
    --fasta "$FASTA" \
    --output-dir "$OUTPUT_DIR" \
    --data-dir "$DATA_DIR" \
    --n-perturbations "$N_PERTURBATIONS" \
    --noise-scale "$NOISE_SCALE" \
    --skip-unperturbed \
    --db-preset reduced_dbs

echo ""
echo "=== AF2 RandomWalk Complete ==="
find "$OUTPUT_DIR" -name "ranked_0.pdb" | wc -l
echo "successful runs"
