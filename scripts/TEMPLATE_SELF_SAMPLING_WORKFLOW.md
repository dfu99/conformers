# Template Self-Sampling Workflow (Protenix)

This workflow starts from:

- `data/template_example/seed_090_frame_000.pdb`

and generates two Protenix jobs:

1. `msa_only` (MSA on, template off)
2. `msa_self_template` (MSA on, template on, using the same structure as template)

It is designed to test whether Protenix can sample alternative conformations even when guided by an extended-state template.

## 1. Install runtime deps (Colab)

```bash
pip install -q -U protenix gemmi
apt-get update -y && apt-get install -y kalign dssp
```

## 2. Generate inputs only (no inference)

```bash
python scripts/template_self_sampling_workflow.py \
  --input_pdb data/template_example/seed_090_frame_000.pdb \
  --msa_root /content/drive/MyDrive/colab_cache/afmfold-data/AVB3/seed_090_frame_000/msa \
  --workflow_dir data/template_example/workflow_outputs
```

This writes:

- `data/template_example/workflow_outputs/inputs/seed_090_frame_000_msa_only.json`
- `data/template_example/workflow_outputs/inputs/seed_090_frame_000_msa_self_template.json`
- `data/template_example/seed_090_frame_000.cif` (template converted from PDB)
- `data/template_example/workflow_outputs/inputs/template_search_results/*.a3m`

## 3. Run both jobs (recommended)

```bash
python scripts/template_self_sampling_workflow.py \
  --input_pdb data/template_example/seed_090_frame_000.pdb \
  --msa_root /content/drive/MyDrive/colab_cache/afmfold-data/AVB3/seed_090_frame_000/msa \
  --workflow_dir data/template_example/workflow_outputs \
  --template_converter auto \
  --template_json_mode templatesPath \
  --template_entry_id s090 \
  --register_template_mmcif \
  --template_mmcif_dir "$PROTENIX_ROOT_DIR/mmcif" \
  --seeds 101,202,303,404,505,606,707,808,909 \
  --samples_per_seed 5 \
  --run
```

Outputs:

- `data/template_example/workflow_outputs/outputs/msa_only/...`
- `data/template_example/workflow_outputs/outputs/msa_self_template/...`

## 4. Run one mode only

MSA-only:

```bash
python scripts/template_self_sampling_workflow.py --run --run_msa_only
```

Self-template only:

```bash
python scripts/template_self_sampling_workflow.py --run --run_template_only
```

## Notes

- Chain mapping is explicit and correct: chain `A` sequence uses template chain `A`, chain `B` uses `B`.
- Default output uses Protenix v1 template schema: `templatesPath` pointing to per-chain `.a3m`.
- Template hit headers include both `ENTRY_CHAIN` and `ENTRY_CHAIN/start-end` variants for parser compatibility.
- When `--register_template_mmcif` is set, the script installs your self-template as
  `mmcif/<entry[1:3]>/<entry>.cif.gz` so Protenix can resolve it from `.a3m`.
- Use `--template_json_mode legacy_templates` only if you need the old inline
  `templates: [{mmcif_path, chain_id}]` format.
- Template guidance is still soft; sampling multiple seeds is required to observe alternative conformations.
- `--template_converter auto` prefers `mkdssp` for PDB->mmCIF conversion and falls back to gemmi.
- Template-enabled runs require `kalign` in `PATH`; optionally pass `--kalign_binary_path /path/to/kalign`.
- If present, `mkdssp` is used automatically; you can also pass `--mkdssp_binary_path /path/to/mkdssp`.
