#!/bin/bash
# SLURM script for running template_self_sampling_workflow.py on PACE cluster.
# Runs Protenix inference with a self-template to bias toward extended integrin conformation.

#SBATCH --job-name=protenix_template
#SBATCH --output=/storage/home/hcoda1/6/dfu71/scratch/logs/protenix_template_%j.out
#SBATCH --error=/storage/home/hcoda1/6/dfu71/scratch/logs/protenix_template_%j.err
#SBATCH -A gts-yke8
#SBATCH -N1 --gres=gpu:RTX_6000:1
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
SCRIPT="${PROTENIX_ROOT}/scripts/template_self_sampling_workflow.py"
INPUT_PDB="${PROTENIX_ROOT}/data/template_example/seed_090_frame_000.pdb"
MSA_ROOT="$HOME/scratch/Protenix/data/seed_090_frame_000/msa"
WORKFLOW_DIR="${PROTENIX_ROOT}/data/template_example/workflow_outputs"
TEMPLATE_CIF="${PROTENIX_ROOT}/data/template_example/seed_090_frame_000.cif"
# Expected chain order for integrin alpha/beta input.
# Override at submit time if needed, e.g.:
#   CHAIN_ORDER=A sbatch scripts/protenix_template_slurm.sh
CHAIN_ORDER="${CHAIN_ORDER:-A,B}"

# Protenix requires template mmcif DB under: $PROTENIX_ROOT_DIR/mmcif
export PROTENIX_ROOT_DIR="$PROTENIX_ROOT"
# If your mmcif is stored elsewhere on cluster, set this env var before sbatch:
#   export PROTENIX_MMCIF_DIR=/path/to/existing/mmcif
MMCIF_SOURCE_DIR="${PROTENIX_MMCIF_DIR:-$PROTENIX_ROOT/mmcif}"
MMCIF_TARGET_DIR="$PROTENIX_ROOT/mmcif"

# Add kalign and hmmer to PATH (built from source)
export PATH="$HOME/scratch/kalign/build/src:$HOME/scratch/hmmer/bin:$PATH"

# kalign must be in PATH or specify explicitly
KALIGN_BIN=$(which kalign 2>/dev/null || which kalign3 2>/dev/null)
if [ -z "$KALIGN_BIN" ]; then
    echo "ERROR: kalign not found in PATH. Install with: conda install -c bioconda kalign3"
    exit 1
fi

# Verify critical files exist
for f in "$SCRIPT" "$INPUT_PDB"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done
if [ ! -d "$MSA_ROOT" ]; then
    echo "ERROR: MSA directory not found: $MSA_ROOT"
    exit 1
fi
if [ ! -d "$MMCIF_SOURCE_DIR" ]; then
    echo "ERROR: mmcif directory not found: $MMCIF_SOURCE_DIR"
    echo "Protenix template mode requires mmcif at \$PROTENIX_ROOT_DIR/mmcif."
    echo "Set PROTENIX_MMCIF_DIR to your existing cluster mmcif path, or create $MMCIF_TARGET_DIR."
    exit 1
fi

# Preflight: verify expected chains exist in INPUT_PDB and print CA counts.
python - "$INPUT_PDB" "$CHAIN_ORDER" <<'PY'
import sys
from collections import defaultdict
from pathlib import Path

pdb_path = Path(sys.argv[1])
expected = [x.strip() for x in sys.argv[2].split(",") if x.strip()]

counts = defaultdict(int)
seen = defaultdict(set)
with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
    for line in handle:
        if not line.startswith("ATOM  "):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain_id = line[21].strip() or "_"
        resseq = line[22:26].strip()
        icode = line[26].strip()
        key = (resseq, icode)
        if key in seen[chain_id]:
            continue
        seen[chain_id].add(key)
        counts[chain_id] += 1

if not counts:
    print(f"ERROR: No protein CA atoms found in {pdb_path}")
    sys.exit(1)

print("Preflight chain CA counts:")
for chain_id in sorted(counts):
    print(f"  chain {chain_id}: length={counts[chain_id]}")

missing = [c for c in expected if c not in counts]
if missing:
    available = ",".join(sorted(counts))
    print(
        "ERROR: Missing expected chain(s): "
        + ",".join(missing)
        + f" in {pdb_path}. Available chains: {available}"
    )
    sys.exit(1)
PY
if [ $? -ne 0 ]; then
    exit 1
fi

# Ensure $PROTENIX_ROOT/mmcif points to the configured source location.
if [ "$MMCIF_SOURCE_DIR" != "$MMCIF_TARGET_DIR" ]; then
    mkdir -p "$PROTENIX_ROOT"
    ln -sfn "$MMCIF_SOURCE_DIR" "$MMCIF_TARGET_DIR"
fi

echo "========================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $(hostname)"
echo "GPU:           $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Script:        $SCRIPT"
echo "Input PDB:     $INPUT_PDB"
echo "Chain order:   $CHAIN_ORDER"
echo "MSA root:      $MSA_ROOT"
echo "Workflow dir:  $WORKFLOW_DIR"
echo "PROTENIX_ROOT_DIR: $PROTENIX_ROOT_DIR"
echo "mmcif source:  $MMCIF_SOURCE_DIR"
echo "mmcif target:  $MMCIF_TARGET_DIR"
echo "kalign:        $KALIGN_BIN"
echo "========================================="

srun python "$SCRIPT" \
    --input_pdb "$INPUT_PDB" \
    --chain_order "$CHAIN_ORDER" \
    --msa_root "$MSA_ROOT" \
    --workflow_dir "$WORKFLOW_DIR" \
    --template_cif "$TEMPLATE_CIF" \
    --template_json_mode templatesPath \
    --template_entry_id s090 \
    --register_template_mmcif \
    --template_mmcif_dir "$MMCIF_TARGET_DIR" \
    --seeds 101,202,303,404,505,606,707,808,909 \
    --samples_per_seed 5 \
    --dtype bf16 \
    --kalign_binary_path "$KALIGN_BIN" \
    --run --run_template_only

RC=$?
if [ $RC -ne 0 ]; then
    echo "ERROR: Workflow exited with non-zero status: $RC"
    exit $RC
fi

# Protenix may log dataloader/template failures without propagating a non-zero code.
# Guard on expected output artifact presence.
PRED_ROOT="${WORKFLOW_DIR}/outputs/msa_self_template"
if ! find "$PRED_ROOT" -type f -name '*.cif' -print -quit | grep -q .; then
    echo "ERROR: No prediction .cif files found under: $PRED_ROOT"
    echo "Check stderr for template/mmcif dependency failures."
    exit 2
fi

echo "Done. Exit code: $RC"
