#!/bin/bash
#SBATCH -J avb3_domain_map
#SBATCH -A gts-yke8
#SBATCH --partition=gpu-rtx6000
#SBATCH --gres=gpu:RTX_6000:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH --output=logs/avb3-conformers/domain_map_%j.log
#SBATCH --error=logs/avb3-conformers/domain_map_%j.err
set -euo pipefail

# Domain-by-domain mapping of αIIbβ3 pathway structures onto αVβ3.
# Aligns head, calf, and leg/tail domains independently.

CONFORMERS_ROOT="$HOME/scratch/conformers"
AIIB3_DIR="$HOME/scratch/principalcurve_integrin_structures/data/image_structures_relaxed"
AVB3_PDB="$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb"
OUTPUT_DIR="$CONFORMERS_ROOT/data/runs/avb3/aiib3_domain_mapped"

mkdir -p "$OUTPUT_DIR" "$CONFORMERS_ROOT/logs/avb3-conformers"

echo "=== Domain-by-domain mapping: αIIbβ3 → αVβ3 ==="
python3 "$CONFORMERS_ROOT/pipelines/avb3-conformers/scripts/map_aiib3_to_avb3_domains.py" \
    --aiib3-dir "$AIIB3_DIR" \
    --avb3-pdb "$AVB3_PDB" \
    --output-dir "$OUTPUT_DIR"

echo ""
echo "=== Done ==="
ls "$OUTPUT_DIR"/*.pdb 2>/dev/null | wc -l
echo "domain-mapped PDB files produced"
