# protenix-avb3-template

AVB3 template registration and template-enabled Protenix inference pipeline.

## Data stream
- Seed/template inputs: `data/avb3/template_example/...`
- Run outputs: `data/runs/avb3/protenix_template/...`
- Logs: `logs/protenix-avb3/...`

## Main entrypoints
- `scripts/setup_protenix_template_inputs.py`
- `scripts/register_protenix_template_entry.py`
- `scripts/run_protenix_pred_simple.sh`
- `scripts/submit_slurm.sh`

## SLURM submit (cluster)
```bash
sbatch pipelines/protenix-avb3-template/scripts/submit_slurm.sh
```

The submit script activates `venv_protenix` by default:
- `PROTENIX_VENV=${PROTENIX_VENV:-$HOME/scratch/venv_protenix}`
