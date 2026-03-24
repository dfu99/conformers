#!/bin/bash
#SBATCH -J avb3_afcluster
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
#SBATCH --output=logs/avb3-conformers/afcluster_avb3_%j.log
#SBATCH --error=logs/avb3-conformers/afcluster_avb3_%j.err
set -euo pipefail

# AVB3 AFCluster sweep: cluster MSAs per chain, then run BoltzGen on
# individual chains (full heterodimer crashes BoltzGen with cuBLAS error).
# This produces diverse single-chain conformers as comparison baselines.

module load cuda

CONFORMERS_ROOT="$HOME/scratch/conformers"
BOLTZ_VENV="$HOME/scratch/venv_boltz"
source "$BOLTZ_VENV/bin/activate"

export HF_HOME="$HOME/scratch/.cache/huggingface"
export HF_HUB_CACHE="$HOME/scratch/.cache/huggingface/hub"
export HF_XET_CACHE="$HOME/scratch/.cache/huggingface/xet"

WORK_DIR="$CONFORMERS_ROOT/data/runs/avb3/afcluster_sweep_${SLURM_JOB_ID}"
MSA_A="$CONFORMERS_ROOT/data/avb3/template_example/msa/0/non_pairing.a3m"
MSA_B="$CONFORMERS_ROOT/data/avb3/template_example/msa/1/non_pairing.a3m"
SEED_PDB="$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.pdb"
TEMPLATE_CIF="$CONFORMERS_ROOT/data/avb3/template_example/seed_090_frame_000.cif"

TOP_A="${TOP_A:-3}"
TOP_B="${TOP_B:-2}"
BUDGET="${BUDGET:-100}"
NUM_DESIGNS="${NUM_DESIGNS:-5000}"

mkdir -p "$WORK_DIR/clusters/A" "$WORK_DIR/clusters/B" "$WORK_DIR/boltz_outputs" \
         "$CONFORMERS_ROOT/logs/avb3-conformers"

echo "=== AVB3 AFCluster Sweep ==="
echo "Work dir: $WORK_DIR"
echo "Top A clusters: $TOP_A, Top B clusters: $TOP_B"
echo "BoltzGen budget: $BUDGET, designs: $NUM_DESIGNS"

# Step 1: Cluster MSAs
echo ""
echo "--- Clustering chain A MSA ---"
python3 "$CONFORMERS_ROOT/pipelines/afcluster/scripts/cluster_chain_msa.py" \
    --input-a3m "$MSA_A" \
    --outdir "$WORK_DIR/clusters/A" \
    2>&1 || echo "WARNING: Chain A clustering failed"

echo ""
echo "--- Clustering chain B MSA ---"
python3 "$CONFORMERS_ROOT/pipelines/afcluster/scripts/cluster_chain_msa.py" \
    --input-a3m "$MSA_B" \
    --outdir "$WORK_DIR/clusters/B" \
    2>&1 || echo "WARNING: Chain B clustering failed"

echo ""
echo "--- Cluster results ---"
echo "Chain A clusters:"
cat "$WORK_DIR/clusters/A/clusters.tsv" 2>/dev/null || echo "  (none)"
echo "Chain B clusters:"
cat "$WORK_DIR/clusters/B/clusters.tsv" 2>/dev/null || echo "  (none)"

# Step 2: Generate BoltzGen specs for individual chains
# (Full heterodimer crashes BoltzGen — run per-chain instead)
echo ""
echo "--- Generating per-chain BoltzGen specs ---"
python3 << 'PYEOF'
import json, sys
from pathlib import Path

work_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
seed_pdb = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("seed.pdb")
template_cif = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("template.cif")
top_a = int(sys.argv[4]) if len(sys.argv) > 4 else 3
top_b = int(sys.argv[5]) if len(sys.argv) > 5 else 2

specs_dir = work_dir / "specs"
specs_dir.mkdir(parents=True, exist_ok=True)

# Read cluster manifests
specs = []
for chain, chain_id, top_n in [("A", "A", top_a), ("B", "B", top_b)]:
    clusters_tsv = work_dir / "clusters" / chain / "clusters.tsv"
    if not clusters_tsv.exists():
        print(f"  No clusters for chain {chain}, skipping")
        continue

    cluster_files = []
    for line in clusters_tsv.read_text().splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            cluster_files.append(parts[1])

    for i, cluster_a3m in enumerate(cluster_files[:top_n]):
        spec_name = f"avb3_chain{chain}_cluster{i}"
        spec = {
            "entities": [{
                "file": {
                    "path": str(template_cif),
                    "include": [{
                        "chain": {
                            "id": chain_id,
                            "msa": cluster_a3m,
                        }
                    }],
                    "structure_groups": "all",
                }
            }]
        }
        spec_path = specs_dir / f"{spec_name}.yaml"

        import yaml
        spec_path.write_text(yaml.dump(spec, default_flow_style=False))
        specs.append({"name": spec_name, "path": str(spec_path), "chain": chain, "cluster": i})
        print(f"  {spec_name} -> {spec_path}")

manifest = specs_dir / "manifest.json"
manifest.write_text(json.dumps(specs, indent=2))
print(f"\nManifest: {manifest} ({len(specs)} specs)")
PYEOF
python3 /dev/stdin "$WORK_DIR" "$SEED_PDB" "$TEMPLATE_CIF" "$TOP_A" "$TOP_B" 2>&1 \
    || echo "WARNING: Spec generation failed, trying without yaml"

# Fallback: generate specs without pyyaml (write YAML manually)
if [[ ! -f "$WORK_DIR/specs/manifest.json" ]]; then
    echo "Generating specs with manual YAML..."
    python3 - "$WORK_DIR" "$TEMPLATE_CIF" "$TOP_A" "$TOP_B" << 'PYEOF2'
import json, sys
from pathlib import Path

work_dir = Path(sys.argv[1])
template_cif = sys.argv[2]
top_a = int(sys.argv[3])
top_b = int(sys.argv[4])

specs_dir = work_dir / "specs"
specs_dir.mkdir(parents=True, exist_ok=True)
specs = []

for chain, chain_id, top_n in [("A", "A", top_a), ("B", "B", top_b)]:
    clusters_tsv = work_dir / "clusters" / chain / "clusters.tsv"
    if not clusters_tsv.exists():
        continue
    cluster_files = []
    for line in clusters_tsv.read_text().splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            cluster_files.append(parts[1])

    for i, cluster_a3m in enumerate(cluster_files[:top_n]):
        spec_name = f"avb3_chain{chain}_cluster{i}"
        yaml_content = f"""entities:
- file:
    path: {template_cif}
    include:
    - chain:
        id: {chain_id}
        msa: {cluster_a3m}
    structure_groups: all
"""
        spec_path = specs_dir / f"{spec_name}.yaml"
        spec_path.write_text(yaml_content)
        specs.append({"name": spec_name, "path": str(spec_path), "chain": chain, "cluster": i})
        print(f"  {spec_name} -> {spec_path}")

(specs_dir / "manifest.json").write_text(json.dumps(specs, indent=2))
print(f"\n{len(specs)} specs generated")
PYEOF2
fi

# Step 3: Run BoltzGen on each per-chain spec
echo ""
echo "--- Running BoltzGen per-chain predictions ---"
if [[ -f "$WORK_DIR/specs/manifest.json" ]]; then
    python3 - "$WORK_DIR" "$BUDGET" "$NUM_DESIGNS" << 'PYEOF3'
import json, subprocess, sys
from pathlib import Path

work_dir = Path(sys.argv[1])
budget = sys.argv[2]
num_designs = sys.argv[3]

manifest = json.loads((work_dir / "specs/manifest.json").read_text())
for spec in manifest:
    name = spec["name"]
    spec_path = spec["path"]
    out_dir = work_dir / "boltz_outputs" / name

    print(f"\n=== {name} ===")
    if list(out_dir.glob("final_ranked_designs/*.cif")):
        print(f"  Already has results, skipping")
        continue

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "boltzgen", "run", spec_path,
        "--output", str(out_dir),
        "--protocol", "protein-anything",
        "--num_designs", num_designs,
        "--budget", budget,
    ]
    print(f"  CMD: {' '.join(cmd)}")
    result = subprocess.run(cmd, timeout=3600)
    if result.returncode != 0:
        print(f"  FAILED (exit {result.returncode})")
    else:
        print(f"  SUCCESS")
        designs = list(out_dir.glob("final_ranked_designs/*.cif"))
        print(f"  Produced {len(designs)} designs")
PYEOF3
else
    echo "  No specs manifest found, skipping BoltzGen"
fi

echo ""
echo "=== AFCluster Sweep Complete ==="
echo "Results: $WORK_DIR"
find "$WORK_DIR/boltz_outputs" -name "*.cif" 2>/dev/null | wc -l
echo "CIF files produced"
