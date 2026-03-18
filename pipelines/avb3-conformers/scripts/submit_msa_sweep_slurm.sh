#!/bin/bash
#SBATCH -J avb3_msa_sweep
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/msa_sweep_%j.log
#SBATCH --error=logs/avb3-conformers/msa_sweep_%j.err
set -euo pipefail

# MSA-subsampled Protenix sweep for AVB3 conformer validation.
# Runs Protenix at each MSA depth to produce reference predictions.
# After this job, run score_conformers.py locally to compare against pulled frames.

CONFORMERS_ROOT="${CONFORMERS_ROOT:-$HOME/scratch/conformers}"
WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/msa_validation}"
SWEEP_CONFIG="${SWEEP_CONFIG:-$WORK_DIR/inputs/sweep_config.json}"

SEED="${SEED:-101}"
MODEL_NAME="${MODEL_NAME:-protenix_base_default_v1.0.0}"
DTYPE="${DTYPE:-bf16}"

module load cuda

mkdir -p "$CONFORMERS_ROOT/logs/avb3-conformers"

# Activate Protenix environment
source "$CONFORMERS_ROOT/venv_protenix/bin/activate" 2>/dev/null || true

echo "=== AVB3 MSA-Subsampled Validation Sweep ==="
echo "Sweep config: $SWEEP_CONFIG"
echo "Seed: $SEED"
echo ""

if [[ ! -f "$SWEEP_CONFIG" ]]; then
    echo "ERROR: sweep config not found at $SWEEP_CONFIG" >&2
    exit 1
fi

# Read sweep config and run Protenix for each MSA depth
python3 - "$SWEEP_CONFIG" "$WORK_DIR" "$SEED" "$MODEL_NAME" "$DTYPE" << 'PYEOF'
import json
import subprocess
import sys
from pathlib import Path

sweep_config = json.loads(Path(sys.argv[1]).read_text())
work_dir = Path(sys.argv[2])
seed = sys.argv[3]
model_name = sys.argv[4]
dtype = sys.argv[5]

for entry in sweep_config:
    label = entry["label"]
    input_json = entry["input_json"]
    use_msa = entry.get("use_msa", True)
    out_dir = work_dir / "predictions" / label

    print(f"\n{'='*60}")
    print(f"Running {label} (MSA: {use_msa})")
    print(f"  Input: {input_json}")
    print(f"  Output: {out_dir}")
    print(f"{'='*60}\n")

    out_dir.mkdir(parents=True, exist_ok=True)

    # Check if predictions already exist
    pred_dir = list(out_dir.glob("*/seed_*/predictions"))
    if pred_dir:
        print(f"  Predictions already exist at {pred_dir[0]}, skipping.")
        continue

    cmd = [
        "protenix", "pred",
        "-i", str(input_json),
        "-o", str(out_dir),
        "-s", seed,
        "-n", model_name,
        "--dtype", dtype,
        "--use_msa", str(use_msa).lower(),
        "--use_template", "false",
    ]

    print(f"  CMD: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"  WARNING: {label} failed with exit code {result.returncode}")
    else:
        print(f"  {label} completed successfully.")

print("\n=== Sweep complete ===")
PYEOF

echo "MSA sweep complete. Run score_conformers.py to compare against pulled frames."
