#!/bin/bash
#SBATCH -J proteinttt_install
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH -C A100-80GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/proteinttt_install_%j.log
#SBATCH --error=logs/avb3-conformers/proteinttt_install_%j.err
set -euo pipefail

# Install OpenFold + ProteinTTT dependencies on a GPU node (needs nvcc).
# After this completes, submit submit_proteinttt_slurm.sh for actual runs.

module load cuda
module load python/3.12.5

PROTEINTTT_ROOT="$HOME/scratch/ProteinTTT"
VENV_DIR="$HOME/scratch/venv_proteinttt"

echo "=== Installing OpenFold + ProteinTTT ==="
echo "Python: $(python3 --version)"
echo "nvcc: $(nvcc --version | tail -1)"

# Rebuild venv from scratch with Python 3.12
rm -rf "$VENV_DIR"
python3 -m venv --upgrade-deps "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# Verify pip works
python3 -m pip --version
python3 -m pip install setuptools wheel numpy

echo "--- Installing PyTorch ---"
python3 -m pip install torch --index-url https://download.pytorch.org/whl/cu121

echo "--- Installing ninja (for faster OpenFold build) ---"
python3 -m pip install ninja

echo "--- Installing OpenFold ---"
python3 -m pip install --no-build-isolation git+https://github.com/aqlaboratory/openfold.git

echo "--- Installing ESM + ProteinTTT deps ---"
python3 -m pip install fair-esm omegaconf pandas biopython requests tqdm ml_collections dm-tree

echo "--- Installing ProteinTTT ---"
python3 -m pip install -e "$PROTEINTTT_ROOT"

# Mark setup complete
touch "$VENV_DIR/.setup_complete"

echo ""
echo "=== Verifying imports ==="
python3 -c "
import torch; print(f'torch {torch.__version__}, CUDA: {torch.cuda.is_available()}')
import openfold; print('openfold OK')
import esm; print('esm OK')
from esm.esmfold.v1.esmfold import ESMFold; print('ESMFold import OK')
from proteinttt.base import TTTConfig; print('ProteinTTT import OK')
from proteinttt.models.esmfold import ESMFoldTTT; print('ESMFoldTTT import OK')
print('All imports successful!')
"

echo ""
echo "=== Install Complete ==="
echo "Now submit the actual ProteinTTT run:"
echo "  sbatch pipelines/avb3-conformers/scripts/submit_proteinttt_slurm.sh"
