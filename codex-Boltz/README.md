# codex-Boltz

Boltz-focused conformational sampling branch for extended integrin states.

This branch has two tracks:

1. Boltz `predict` sweep (implemented now)
- Generates many YAML jobs varying template force and MSA usage.
- Runs Boltz inference across the generated job set.
- Ranks outputs by extension proxy.

2. BoltzGen track (starter wrapper)
- Provides a launch wrapper you can point to a BoltzGen design spec once finalized.

## Implemented files

- `scripts/extract_integrin_sequences.py`
- `scripts/build_boltz_predict_sweep.py`
- `scripts/run_boltz_predict_sweep.sh`
- `scripts/run_boltzgen_attempt.sh`
- `scripts/rank_extended.py`

## Requirements

- Python: `pyyaml`, `gemmi`
- CLI: `boltz` (and optionally `boltzgen`)

```bash
pip install pyyaml gemmi
```

## Quick start (predict sweep)

```bash
cd codex-Boltz

./scripts/run_boltz_predict_sweep.sh \
  --sequence-file ../data/sequences_updated \
  --template-cif ../data/template_example/seed_090_frame_000.cif \
  --outdir ./runs/integrin_ab_boltz \
  --chain-a-msa ../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/0/non_pairing.a3m \
  --chain-b-msa ../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/1/non_pairing.a3m
```

Then inspect:

- `runs/integrin_ab_boltz/jobs/`
- `runs/integrin_ab_boltz/boltz_outputs/`
- `runs/integrin_ab_boltz/extended_rank.tsv`

## BoltzGen wrapper

```bash
./scripts/run_boltzgen_attempt.sh \
  --design-spec /path/to/your_boltzgen_design_spec.yaml \
  --outdir ./runs/boltzgen_try_01
```

This wrapper is intentionally thin; BoltzGen spec details are project-specific.
