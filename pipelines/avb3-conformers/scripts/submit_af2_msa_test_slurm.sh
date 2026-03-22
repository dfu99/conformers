#!/bin/bash
#SBATCH -J avb3_af2_msa_test
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
#SBATCH --output=logs/avb3-conformers/af2_msa_test_%j.log
#SBATCH --error=logs/avb3-conformers/af2_msa_test_%j.err
set -euo pipefail

# AlphaFold2 MSA-subsampled conformer validation for AVB3.

module purge
module load alphafold/2.3.2

export DOWNLOAD_DIR=/storage/coda1/d-pace_community/0/alphafold/alphafold_2.3.2_data
export NSLOTS=${SLURM_CPUS_PER_TASK:-8}
export PYTHONNOUSERSITE=1

CONFORMERS_ROOT="$HOME/scratch/conformers"
WORK_DIR="$CONFORMERS_ROOT/data/runs/avb3/af2_msa_validation"
REF_PDB="$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb"

mkdir -p "$WORK_DIR/fasta" "$WORK_DIR/predictions" "$CONFORMERS_ROOT/logs/avb3-conformers"

echo "=== AF2 MSA-Subsampled Validation Test ==="

# Step 1: Extract AVB3 sequence and write FASTA
python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/extract_avb3_fasta.py" \
    "$REF_PDB" "$WORK_DIR/fasta"

# Common AF2 database paths
DB_COMMON=(
    --uniref90_database_path="$DOWNLOAD_DIR/uniref90/uniref90.fasta"
    --mgnify_database_path="$DOWNLOAD_DIR/mgnify/mgy_clusters_2022_05.fa"
    --template_mmcif_dir="$DOWNLOAD_DIR/pdb_mmcif/mmcif_files"
    --obsolete_pdbs_path="$DOWNLOAD_DIR/pdb_mmcif/obsolete.dat"
    --uniprot_database_path="$DOWNLOAD_DIR/uniprot/uniprot.fasta"
    --pdb_seqres_database_path="$DOWNLOAD_DIR/pdb_seqres/pdb_seqres.txt"
)

# Step 2: Run AF2 multimer with reduced_dbs (faster, smaller MSA)
echo ""
echo "=== Running AF2 multimer with reduced_dbs ==="
OUTPUT_REDUCED="$WORK_DIR/predictions/reduced_dbs"
mkdir -p "$OUTPUT_REDUCED"

if ls "$OUTPUT_REDUCED"/ranked_*.pdb >/dev/null 2>&1; then
    echo "  Reduced MSA predictions already exist, skipping."
else
    alphafold \
        --fasta_paths="$WORK_DIR/fasta/avb3_multimer.fasta" \
        --max_template_date=2020-05-14 \
        --model_preset=multimer \
        --db_preset=reduced_dbs \
        --output_dir="$OUTPUT_REDUCED" \
        --data_dir="$DOWNLOAD_DIR" \
        --small_bfd_database_path="$DOWNLOAD_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" \
        "${DB_COMMON[@]}" \
        --use_gpu_relax=true \
        2>&1 | tee "$WORK_DIR/af2_reduced_dbs.log"
fi

# Step 3: Run AF2 multimer with full_dbs (full MSA, much deeper)
echo ""
echo "=== Running AF2 multimer with full_dbs ==="
OUTPUT_FULL="$WORK_DIR/predictions/full_dbs"
mkdir -p "$OUTPUT_FULL"

if ls "$OUTPUT_FULL"/ranked_*.pdb >/dev/null 2>&1; then
    echo "  Full MSA predictions already exist, skipping."
else
    alphafold \
        --fasta_paths="$WORK_DIR/fasta/avb3_multimer.fasta" \
        --max_template_date=2020-05-14 \
        --model_preset=multimer \
        --db_preset=full_dbs \
        --output_dir="$OUTPUT_FULL" \
        --data_dir="$DOWNLOAD_DIR" \
        --bfd_database_path="$DOWNLOAD_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt" \
        --uniref30_database_path="$DOWNLOAD_DIR/uniref30/UniRef30_2021_03" \
        "${DB_COMMON[@]}" \
        --use_gpu_relax=true \
        2>&1 | tee "$WORK_DIR/af2_full_dbs.log"
fi

echo ""
echo "=== AF2 MSA Test Complete ==="
echo "Reduced MSA predictions:"
ls "$OUTPUT_REDUCED"/ranked_*.pdb 2>/dev/null || echo "  (none)"
echo "Full MSA predictions:"
ls "$OUTPUT_FULL"/ranked_*.pdb 2>/dev/null || echo "  (none)"
