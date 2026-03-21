#!/bin/bash
#SBATCH -J avb3_af2_msa_test
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/af2_msa_test_%j.log
#SBATCH --error=logs/avb3-conformers/af2_msa_test_%j.err
set -euo pipefail

# AlphaFold2 MSA-subsampled conformer validation for AVB3.
# Runs AF2 at full MSA and shallow MSA (5%) to compare predicted
# conformations, testing whether MSA depth diversifies AF2 output
# (unlike Protenix which was MSA-depth-invariant).

module purge
module load alphafold/2.3.2

export DOWNLOAD_DIR=/storage/coda1/d-pace_community/0/alphafold/alphafold_2.3.2_data
export NSLOTS=${SLURM_CPUS_PER_TASK:-8}

CONFORMERS_ROOT="$HOME/scratch/conformers"
WORK_DIR="$CONFORMERS_ROOT/data/runs/avb3/af2_msa_validation"
REF_PDB="$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb"

mkdir -p "$WORK_DIR/fasta" "$WORK_DIR/predictions" "$CONFORMERS_ROOT/logs/avb3-conformers"

echo "=== AF2 MSA-Subsampled Validation Test ==="

# Step 1: Extract AVB3 sequence and write FASTA
python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/extract_avb3_fasta.py" \
    "$REF_PDB" "$WORK_DIR/fasta"

# Step 2: Run AF2 with full database (default)
echo ""
echo "=== Running AF2 with full MSA ==="
OUTPUT_FULL="$WORK_DIR/predictions/full_msa"
mkdir -p "$OUTPUT_FULL"

if ls "$OUTPUT_FULL"/ranked_*.pdb >/dev/null 2>&1; then
    echo "  Full MSA predictions already exist, skipping."
else
    export OUTPUT_DIR="$OUTPUT_FULL"
    alphafold \
        --fasta_paths="$WORK_DIR/fasta/avb3_multimer.fasta" \
        --max_template_date=2020-05-14 \
        --model_preset=multimer \
        --db_preset=full \
        --output_dir="$OUTPUT_FULL" \
        --data_dir="$DOWNLOAD_DIR" \
        --use_gpu_relax=true \
        || echo "WARNING: AF2 full MSA run failed"
fi

# Step 3: Run AF2 with reduced database (proxy for shallow MSA)
echo ""
echo "=== Running AF2 with reduced MSA ==="
OUTPUT_REDUCED="$WORK_DIR/predictions/reduced_msa"
mkdir -p "$OUTPUT_REDUCED"

if ls "$OUTPUT_REDUCED"/ranked_*.pdb >/dev/null 2>&1; then
    echo "  Reduced MSA predictions already exist, skipping."
else
    export OUTPUT_DIR="$OUTPUT_REDUCED"
    alphafold \
        --fasta_paths="$WORK_DIR/fasta/avb3_multimer.fasta" \
        --max_template_date=2020-05-14 \
        --model_preset=multimer \
        --db_preset=reduced_dbs \
        --output_dir="$OUTPUT_REDUCED" \
        --data_dir="$DOWNLOAD_DIR" \
        --use_gpu_relax=true \
        || echo "WARNING: AF2 reduced MSA run failed"
fi

echo ""
echo "=== AF2 MSA Test Complete ==="
echo "Full MSA predictions:"
ls "$OUTPUT_FULL"/ranked_*.pdb 2>/dev/null || echo "  (none)"
echo "Reduced MSA predictions:"
ls "$OUTPUT_REDUCED"/ranked_*.pdb 2>/dev/null || echo "  (none)"
