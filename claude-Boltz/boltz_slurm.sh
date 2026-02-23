#!/bin/bash
# SLURM script for Boltz-1 conformation sampling of integrin.
# Fully independent of Protenix — install boltz, point at PDB, run.
#
# Install: pip install boltz

#SBATCH --job-name=boltz_integrin
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/boltz_integrin_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/conformers/logs/boltz_integrin_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu

module load cuda
module load python/3.12

# Activate venv with boltz installed (pip install boltz)
source ~/scratch/venv_boltz/bin/activate

# --- Paths ---
PROTENIX_ROOT="$HOME/scratch/conformers"
SCRIPT="${PROTENIX_ROOT}/claude-Boltz/boltz_workflow.py"
INPUT_PDB="${PROTENIX_ROOT}/data/template_example/seed_090_frame_000.pdb"
WORKFLOW_DIR="${PROTENIX_ROOT}/data/boltz_outputs"

# --- Sanity checks ---
if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: $SCRIPT not found"; exit 1
fi
if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: $INPUT_PDB not found"; exit 1
fi
if ! command -v boltz &>/dev/null; then
    echo "ERROR: boltz not in PATH. Run: pip install boltz"
    exit 1
fi

echo "========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $(hostname)"
echo "GPU:          $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo N/A)"
echo "boltz:        $(which boltz)"
echo "Input PDB:    $INPUT_PDB"
echo "Workflow dir: $WORKFLOW_DIR"
echo "========================================="

srun python "$SCRIPT" \
    --input_pdb "$INPUT_PDB" \
    --workflow_dir "$WORKFLOW_DIR" \
    --seeds 101,202,303,404,505,606,707,808,909 \
    --diffusion_samples 5 \
    --recycling_steps 3 \
    --sampling_steps 200 \
    --accelerator gpu \
    --run

RC=$?
[ $RC -ne 0 ] && { echo "ERROR: boltz_workflow.py exited $RC"; exit $RC; }

CIF_COUNT=$(find "$WORKFLOW_DIR/outputs" -name '*.cif' 2>/dev/null | wc -l)
echo "Done. $CIF_COUNT .cif files under $WORKFLOW_DIR/outputs"
