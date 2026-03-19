#!/bin/bash
#SBATCH -J avb3_msa_test
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
#SBATCH --output=logs/avb3-conformers/msa_test_%j.log
#SBATCH --error=logs/avb3-conformers/msa_test_%j.err
set -euo pipefail

# Minimal 2-depth MSA validation test: full MSA vs 5% subsample.
# Tests whether Protenix predictions differ by MSA depth for AVB3.

module load cuda
source ~/scratch/venv_protenix/bin/activate

CONFORMERS_ROOT="$HOME/scratch/conformers"
WORK_DIR="$CONFORMERS_ROOT/data/runs/avb3/msa_validation"
MSA_A="$CONFORMERS_ROOT/data/avb3/template_example/msa/0/non_pairing.a3m"
MSA_B="$CONFORMERS_ROOT/data/avb3/template_example/msa/1/non_pairing.a3m"
REF_PDB="$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb"

mkdir -p "$WORK_DIR/msa_subsamples/chain_A" \
         "$WORK_DIR/msa_subsamples/chain_B" \
         "$WORK_DIR/inputs_test" \
         "$WORK_DIR/predictions"

echo "=== MSA Validation Test (2 depths: full + 5%) ==="

# Step 1: Subsample MSAs
echo "--- Subsampling MSAs ---"
python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/subsample_msa.py" \
    --input "$MSA_A" \
    --output-dir "$WORK_DIR/msa_subsamples/chain_A" \
    --fractions "1.0,0.05"

python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/subsample_msa.py" \
    --input "$MSA_B" \
    --output-dir "$WORK_DIR/msa_subsamples/chain_B" \
    --fractions "1.0,0.05"

# Step 2: Build input JSONs and run Protenix for each depth
for DEPTH in "1.00" "0.05"; do
    LABEL="depth_${DEPTH}"
    PRED_DIR="$WORK_DIR/predictions/${LABEL}"
    INPUT_JSON="$WORK_DIR/inputs_test/${LABEL}_input.json"
    MSA_A_SUB="$WORK_DIR/msa_subsamples/chain_A/${LABEL}.a3m"
    MSA_B_SUB="$WORK_DIR/msa_subsamples/chain_B/${LABEL}.a3m"

    mkdir -p "$PRED_DIR"

    echo ""
    echo "=== ${LABEL} ==="
    echo "  MSA A: $(grep -c '^>' "$MSA_A_SUB") sequences"
    echo "  MSA B: $(grep -c '^>' "$MSA_B_SUB") sequences"

    # Build Protenix input JSON
    python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/build_msa_test_input.py" \
        --ref-pdb "$REF_PDB" \
        --msa-a "$MSA_A_SUB" \
        --msa-b "$MSA_B_SUB" \
        --label "$LABEL" \
        --output "$INPUT_JSON"

    # Run Protenix if predictions don't already exist
    if ls "$PRED_DIR"/*/seed_*/predictions/*sample_*.cif >/dev/null 2>&1; then
        echo "  Predictions already exist, skipping."
    else
        echo "  Running protenix pred..."
        protenix pred \
            -i "$INPUT_JSON" \
            -o "$PRED_DIR" \
            -s 101 \
            -n protenix_base_default_v1.0.0 \
            --dtype bf16 \
            --use_msa true \
            --use_template false
    fi

    echo "  Done: ${LABEL}"
done

echo ""
echo "=== MSA Test Complete ==="
echo "Predictions:"
find "$WORK_DIR/predictions/" -name "*sample_*.cif" 2>/dev/null
