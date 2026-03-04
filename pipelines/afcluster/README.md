# afcluster

AFCluster-style conformational sampling for A5B1 using Boltz as the structure backend.

## Data stream
- Input sequence: `data/a5b1/sequences/sequences_updated`
- Prior MSA/template source: `data/runs/a5b1/protenix/...`
- Outputs: `data/runs/afcluster/...`

## Main entrypoints
- `scripts/run_afcluster_pipeline.sh`
- `scripts/submit_afcluster_boltz_slurm.sh`

## SLURM submit (cluster)
```bash
sbatch pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh
```

The submit script activates `venv_boltz` by default:
- `BOLTZ_VENV=${BOLTZ_VENV:-$HOME/scratch/venv_boltz}`

`chain_A.seq` and `chain_B.seq` are auto-generated from `data/a5b1/sequences/sequences_updated` if missing.
