# boltz

Boltz-focused conformational sampling pipeline for finding rare extended AVB3 conformers.

## Data stream
- Seed structure: `data/avb3/template_example/seed_090_frame_000.pdb`
- MSA/template source: `data/avb3/template_example/...`
- Outputs: `data/runs/boltz/avb3/...`

## Main entrypoints
- `scripts/run_boltz_predict_sweep.sh`
- `scripts/submit_boltz_predict_slurm.sh`

## SLURM submit (cluster)
```bash
sbatch pipelines/boltz/scripts/submit_boltz_predict_slurm.sh
```

The submit script activates `venv_boltz` by default:
- `BOLTZ_VENV=${BOLTZ_VENV:-$HOME/scratch/venv_boltz}`

Override defaults at submit time as needed, for example:
```bash
SEED_PDB=$HOME/scratch/conformers/data/avb3/template_example/seed_090_frame_000.pdb \
OUTDIR=$HOME/scratch/conformers/data/runs/boltz/avb3/custom_run \
sbatch pipelines/boltz/scripts/submit_boltz_predict_slurm.sh
```
