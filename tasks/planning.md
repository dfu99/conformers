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

## Next Priority
1. **Post-merge geometry checks** — Add tail distance and interface sanity validation before ranking the final merged structure.
2. **Unified ranking table** — Build one script that collects and ranks outputs across Protenix/Boltz/AFCluster pipelines.
3. **Run real A5B1 pipeline on PACE** — Use pace_minimal.sh to submit the full staged tagged pipeline and fetch results.

## Recently Completed
- [x] Add dedicated A5B1 staged pipeline runner + sbatch entrypoint. (obj-001, 2026-03-06)
- [x] Add deterministic stage-merging utility for final tagged complex export. (obj-002, 2026-03-06)
- [x] Build PACE minimal remote control script with smoke test validation. (obj-003, 2026-03-09)
- [x] Add Claude skill for deterministic PACE job operations. (obj-004, 2026-03-10)

## Implemented: AFCluster BoltzGen CLI Alignment (2026-03-06)
- Updated `pipelines/afcluster/scripts/run_afcluster_pipeline.sh` to support backend selection and BoltzGen-native execution.
- Default backend is now `boltzgen` with parameters: `--protocol`, `--num-designs`, `--budget`.
- Updated `pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh` to pass BoltzGen parameters by default while retaining legacy `boltz` mode.
- Corrected helper launcher `pipelines/boltz/scripts/run_boltzgen_attempt.sh` to use `boltzgen run ... --output ... --protocol ... --num_designs ... --budget ...`.

## Implemented: PACE Minimal Remote Control + RTX_6000 Test Target (2026-03-09)
- Added `pipelines/protenix-a5b1/scripts/pace_minimal.sh` with:
  - `check` (connectivity + SLURM command availability check)
  - `submit` (remote `git pull` + `sbatch`)
  - `watch` (remote `squeue`/`sacct` poll loop)
  - `fetch` (job logs + result sync back to local workspace; supports `.log`, `.err`, and `.out`)
  - `smoke` (tiny submit/watch/fetch validation job with configurable poll/sleep; default sleep 60s)
- Updated `pipelines/protenix-a5b1/scripts/submit_complete_tagged_pipeline_slurm.sh` GPU request from `A100` to `RTX_6000` and removed A100-only constraint for queue availability.
- Documented minimal usage in `pipelines/protenix-a5b1/README.md`.
- Verified against `dfu71@login-phoenix.pace.gatech.edu` with completed smoke jobs (`4734209`, `4734230`), successful `squeue -u dfu71` status checks, and local fetch of `.log`, `.err`, and smoke result artifacts.

## Implemented: Claude Skill for PACE Job Ops (2026-03-10)
- Added `.claude/skills/pace-slurm-ops/SKILL.md` with explicit trigger conditions and deterministic command flow for:
  - SSH preflight/auth checks
  - smoke submit/watch/fetch (`hello world`, 60s sleep)
  - real pipeline submit/watch/fetch
  - queue/status verification (`squeue -u dfu71`, `sacct`)
  - common failure handling (VPN, key auth, account requirement)
- Added `CLAUDE.md` pointer to this skill so Claude agents can discover and apply it directly instead of re-deriving workflow details.
