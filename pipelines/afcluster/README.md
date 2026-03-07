# afcluster

AFCluster-style conformational sampling for AVB3 with selectable backend execution.

## Data stream
- Seed structure: `data/avb3/template_example/seed_090_frame_000.pdb`
- MSA/template source: `data/avb3/template_example/...`
- Outputs: `data/runs/afcluster/avb3/...`

## Main entrypoints
- `scripts/run_afcluster_pipeline.sh`
- `scripts/submit_afcluster_boltz_slurm.sh`

## Backend modes
- `boltzgen` (default): runs each generated YAML spec via
  `boltzgen run <spec.yaml> --output <dir> --protocol <name> --num_designs <N> --budget <K>`
- `boltz` (legacy): runs `boltz predict` over generated YAML sweep.

## SLURM submit (cluster)
```bash
sbatch pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh
```

Useful overrides:
```bash
BACKEND=boltzgen PROTOCOL=protein-anything NUM_DESIGNS=10000 BUDGET=200 \
sbatch pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh
```

The submit script activates `venv_boltz` by default:
- `BOLTZ_VENV=${BOLTZ_VENV:-$HOME/scratch/venv_boltz}`

`chain_A.seq` and `chain_B.seq` are auto-generated from `seed_090_frame_000.pdb` if missing.
