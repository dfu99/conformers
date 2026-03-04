# boltz

Boltz-focused conformational sampling pipeline for finding extended A5B1 states.

## Data stream
- Input sequence: `data/a5b1/sequences/sequences_updated`
- Prior MSA/template source: `data/runs/a5b1/protenix/...`
- Outputs: `data/runs/boltz/...`

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
SEQUENCE_FILE=$HOME/scratch/conformers/data/a5b1/sequences/sequences_updated \
OUTDIR=$HOME/scratch/conformers/data/runs/boltz/custom_run \
sbatch pipelines/boltz/scripts/submit_boltz_predict_slurm.sh
```
