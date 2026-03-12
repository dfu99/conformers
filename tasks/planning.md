# Planning

## Purpose
Keep A5B1 and AVB3 efforts decoupled at script, data, and output levels while iterating toward extended conformations.

## Implemented Refactor (2026-03-04)

### 1) Canonical pipeline roots under `pipelines/`
- `pipelines/protenix-a5b1`
- `pipelines/protenix-avb3-template`
- `pipelines/boltz`
- `pipelines/afcluster`

### 2) Canonical data streams (clean switch)
- A5B1 sequence source: `data/a5b1/sequences/sequences_updated`
- AVB3 template source: `data/avb3/template_example/...`
- A5B1 Protenix outputs: `data/runs/a5b1/protenix/...`
- AVB3 Protenix outputs: `data/runs/avb3/...`
- Boltz outputs: `data/runs/boltz/...`
- AFCluster outputs: `data/runs/afcluster/...`

### 3) Legacy aliases removed
- Removed `codex-AFCluster` and `codex-Boltz` aliases.
- Removed legacy `scripts/*` pipeline symlink entrypoints.
- Reduced `scripts/` to reference-only files.

### 4) Environment split enforced in submit scripts
- Boltz/AFCluster submit scripts activate `venv_boltz`.
- Protenix AVB3 submit script activates `venv_protenix`.

## Active Workstreams

### Workstream A: A5B1 heterodimer + staged partner assembly
1. Predict and rank heterodimer-only core (A5-Avi + B1-SpyCatcher).
2. Run staged docking/refinement for SpyTag and Streptavidin separately.
3. Merge via overlapping heterodimer residues and validate linker geometry.

### Workstream B: AVB3 extended-state discovery
1. Keep Protenix template registration/inference isolated in AVB3 pipeline.
2. Run Boltz and AFCluster sweeps as independent branches.
3. Rank by extension proxy + interface plausibility and compare across branches.

## Implemented: A5B1 Complete Tagged Pipeline (2026-03-06)

Added end-to-end staged pipeline under `pipelines/protenix-a5b1/scripts`:
- `run_complete_tagged_pipeline.sh`
- `submit_complete_tagged_pipeline_slurm.sh`
- `build_staged_protenix_inputs.py`
- `merge_staged_tagged_complex.py`
- existing `setup_staged_attachment_workflow.py`

Execution flow:
1. Select accepted heterodimer CIF from `data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/...`.
2. Build stage-1 (`A5B1 + SpyTag`) and stage-2 (`A5B1 + Streptavidin`) Protenix input JSONs.
3. Run `protenix pred` for both stages.
4. Select best sample per stage by `ranking_score`.
5. Receptor-align stage outputs to the accepted heterodimer and write one merged final tagged CIF/PDB.

Final outputs:
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.cif`
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.pdb`

## Next Implementation Tasks
- [x] Add dedicated A5B1 staged pipeline runner + sbatch entrypoint.
- [x] Add deterministic stage-merging utility for final tagged complex export.
- [ ] Add post-merge geometry checks (tail distance/interface sanity) before ranking final structure.
- [ ] Add one unified ranking table script across Protenix/Boltz/AFCluster outputs.

## Implemented: AFCluster BoltzGen CLI Alignment (2026-03-06)
- Updated `pipelines/afcluster/scripts/run_afcluster_pipeline.sh` to support backend selection and BoltzGen-native execution.
- Default backend is now `boltzgen` with parameters: `--protocol`, `--num-designs`, `--budget`.
- Updated `pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh` to pass BoltzGen parameters by default while retaining legacy `boltz` mode.
- Corrected helper launcher `pipelines/boltz/scripts/run_boltzgen_attempt.sh` to use `boltzgen run ... --output ... --protocol ... --num_designs ... --budget ...`.

## Implemented: PACE Minimal Remote Control + Queue-Test GPU Target (2026-03-09)
- Added `pipelines/protenix-a5b1/scripts/pace_minimal.sh` with:
  - `check` (connectivity + SLURM command availability check)
  - `submit` (remote `git pull` + `sbatch`)
  - `watch` (remote `squeue`/`sacct` poll loop)
  - `fetch` (job logs + result sync back to local workspace; supports `.log`, `.err`, and `.out`)
  - `smoke` (tiny submit/watch/fetch validation job with configurable poll/sleep; default sleep 60s)
- Temporarily targeted `RTX_6000` for queue-friendly smoke validation during remote automation bring-up.
- Documented minimal usage in `pipelines/protenix-a5b1/README.md`.
- Verified against `dfu71@login-phoenix.pace.gatech.edu` with completed smoke jobs (`4734209`, `4734230`), successful `squeue -u dfu71` status checks, and local fetch of `.log`, `.err`, and smoke result artifacts.

## Updated: A100 Runtime Requirement for Protenix A5B1 (2026-03-10)
- Reverted A5B1 Protenix submit scripts to A100 80GB requests:
  - `pipelines/protenix-a5b1/scripts/submit_complete_tagged_pipeline_slurm.sh`
  - `pipelines/protenix-a5b1/scripts/submit_conjugates_first_pipeline_slurm.sh`
- Added legacy heterodimer predictions fallback in `run_complete_tagged_pipeline.sh`:
  - primary: `data/runs/a5b1/protenix/.../predictions`
  - fallback: `data/a5b1/outputs/.../predictions`
- Goal: avoid immediate path-failure for tail-aware run and avoid RTX6000 OOM behavior by scheduling on A100 80GB.

## Implemented: Claude Skill for PACE Job Ops (2026-03-10)
- Added `.claude/skills/pace-slurm-ops/SKILL.md` with explicit trigger conditions and deterministic command flow for:
  - SSH preflight/auth checks
  - smoke submit/watch/fetch (`hello world`, 60s sleep)
  - real pipeline submit/watch/fetch
  - queue/status verification (`squeue -u dfu71`, `sacct`)
  - common failure handling (VPN, key auth, account requirement)
- Added `CLAUDE.md` pointer to this skill so Claude agents can discover and apply it directly instead of re-deriving workflow details.

## Implemented: Tail-Aware Selection + Conjugates-First Branch (2026-03-10)
- Updated `pipelines/protenix-a5b1/scripts/merge_staged_tagged_complex.py` to support per-stage prediction selection modes:
  - `ranking`
  - `tail_distance`
  - `hybrid` (ranking among near-tail candidates, fallback by distance)
- Added explicit tail/anchor arguments and recorded candidate tables + selected sample metadata in merge summary JSON.
- Updated `pipelines/protenix-a5b1/scripts/run_complete_tagged_pipeline.sh` to expose selection env controls:
  - `SELECTION_MODE`
  - `STAGE1_LIGAND_ANCHOR_RESIDUE`, `STAGE2_LIGAND_ANCHOR_RESIDUE`
  - `STAGE1_MAX_TAIL_DISTANCE`, `STAGE2_MAX_TAIL_DISTANCE`
- Updated `pipelines/protenix-a5b1/scripts/setup_staged_attachment_workflow.py` to default SpyTag reactive residue to D10 (from collaborator design reference).
- Added experimental conjugates-first branch:
  - `pipelines/protenix-a5b1/scripts/build_conjugates_first_protenix_inputs.py`
  - `pipelines/protenix-a5b1/scripts/run_conjugates_first_pipeline.sh`
- Updated `pipelines/protenix-a5b1/README.md` with tail-aware mode usage and conjugates-first workflow docs.

## Updated: Conjugates-First Combined Complex Output (2026-03-11)
- Extended `merge_staged_tagged_complex.py` with stage-specific receptor chain mappings:
  - `--stage1-receptor-chain-map`
  - `--stage2-receptor-chain-map`
- Updated `run_conjugates_first_pipeline.sh` to produce merged final outputs from best stage predictions:
  - `data/runs/a5b1/conjugates_first/outputs/final/a5b1_conjugates_first_combined.cif`
  - `data/runs/a5b1/conjugates_first/outputs/final/a5b1_conjugates_first_combined.pdb`
  - `data/runs/a5b1/conjugates_first/outputs/final/a5b1_conjugates_first_combined.merge_summary.json`
- `conjugates_first_summary.json` now includes `combined_output` paths.
