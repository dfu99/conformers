# afcluster

AFCluster-style conformational sampling for AVB3 using Boltz as the structure backend.

## Data stream
- Seed structure: `data/avb3/template_example/seed_090_frame_000.pdb`
- MSA/template source: `data/avb3/template_example/...`
- Outputs: `data/runs/afcluster/avb3/...`

## Main entrypoints
- `scripts/run_afcluster_pipeline.sh`
- `scripts/submit_afcluster_boltz_slurm.sh`

## SLURM submit (cluster)
```bash
sbatch pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh
```

The submit script activates `venv_boltz` by default:
- `BOLTZ_VENV=${BOLTZ_VENV:-$HOME/scratch/venv_boltz}`

`chain_A.seq` and `chain_B.seq` are auto-generated from `seed_090_frame_000.pdb` if missing.
