#!/bin/bash
# SLURM script for AF-Cluster integrin sampling via LocalColabFold (AF2).
# Clusters the MSA by gap pattern and runs ColabFold on each subset.
#
# Prerequisite: install LocalColabFold
#   https://github.com/YoshitakaMo/localcolabfold
#   bash install_colabfold.sh ~/scratch/localcolabfold

#SBATCH --job-name=afcluster_integrin
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/afcluster_integrin_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/afcluster_integrin_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

module load cuda
module load python/3.12

# Activate venv (needs scikit-learn: pip install scikit-learn)
source ~/scratch/venv_protenix/bin/activate

# LocalColabFold bin directory on PATH
export PATH="$HOME/scratch/localcolabfold/colabfold-conda/bin:$PATH"

# --- Paths ---
PROTENIX_ROOT="$HOME/scratch/conformers"
SCRIPT="${PROTENIX_ROOT}/claude-AFCluster/afcluster_workflow.py"
INPUT_PDB="${PROTENIX_ROOT}/data/template_example/seed_090_frame_000.pdb"
# MSA root with per-chain subdirs (0/non_pairing.a3m, 1/non_pairing.a3m, ...)
# Compatible with Protenix-generated MSA or ColabFold --save-all output.
MSA_ROOT="${PROTENIX_ROOT}/data/seed_090_frame_000/msa"
WORKFLOW_DIR="${PROTENIX_ROOT}/data/afcluster_outputs"

# --- Sanity checks ---
for f in "$SCRIPT" "$INPUT_PDB"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"; exit 1
    fi
done
if [ ! -d "$MSA_ROOT" ]; then
    echo "ERROR: MSA root not found: $MSA_ROOT"; exit 1
fi
if ! command -v colabfold_batch &>/dev/null; then
    echo "ERROR: colabfold_batch not in PATH."
    echo "Install LocalColabFold: https://github.com/YoshitakaMo/localcolabfold"
    exit 1
fi

echo "========================================="
echo "Job ID:          $SLURM_JOB_ID"
echo "Node:            $(hostname)"
echo "GPU:             $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo N/A)"
echo "colabfold_batch: $(which colabfold_batch)"
echo "Script:          $SCRIPT"
echo "MSA root:        $MSA_ROOT"
echo "Workflow dir:    $WORKFLOW_DIR"
echo "========================================="

srun python "$SCRIPT" \
    --input_pdb "$INPUT_PDB" \
    --msa_root "$MSA_ROOT" \
    --workflow_dir "$WORKFLOW_DIR" \
    --n_clusters 8 \
    --min_cluster_size 5 \
    --cluster_on_chain 0 \
    --num_seeds 3 \
    --num_recycles 3 \
    --num_models 5 \
    --num_relax 0 \
    --run

RC=$?
[ $RC -ne 0 ] && { echo "ERROR: workflow exited $RC"; exit $RC; }

COUNT=$(find "$WORKFLOW_DIR/outputs" -name '*.pdb' -o -name '*.cif' 2>/dev/null | wc -l)
echo "Done. $COUNT structure files under $WORKFLOW_DIR/outputs"
