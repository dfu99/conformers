# CLAUDE.md

This file is an implementation contract for Claude-assisted refactors in this repository.

## 1) Canonical Structure (Do Not Revert)
As of 2026-03-04, runnable pipelines were consolidated into:

- `pipelines/protenix-a5b1`
- `pipelines/protenix-avb3-template`
- `pipelines/boltz`
- `pipelines/afcluster`

Do not reintroduce legacy runnable roots like:

- `codex-Boltz/`
- `codex-AFCluster/`
- legacy pipeline launchers under `scripts/` (except explicit reference-only files)

## 2) Canonical Data Streams
Use these paths only:

- A5B1 inputs: `data/a5b1/...`
- AVB3 template inputs: `data/avb3/...`
- A5B1 Protenix outputs: `data/runs/a5b1/protenix/...`
- AVB3 Protenix outputs: `data/runs/avb3/...`
- Boltz outputs: `data/runs/boltz/...`
- AFCluster outputs: `data/runs/afcluster/...`

Do not write new code against old paths like:

- `data/sequences_updated`
- `data/template_example`
- `data/outputs`

## 3) Runtime Entry Points
Use only these submit scripts for cluster runs:

- `pipelines/protenix-avb3-template/scripts/submit_slurm.sh`
- `pipelines/boltz/scripts/submit_boltz_predict_slurm.sh`
- `pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh`

## 4) Environment Separation (Hard Requirement)
Dependency conflicts are real (e.g., gemmi mismatch). Keep strict env split:

- Protenix pipelines: `venv_protenix`
- Boltz + AFCluster pipelines: `venv_boltz`

Never collapse these into one shared environment.

## 5) Logs Are Intentionally Tracked
Logs are part of debugging sync between local and cluster and must stay sync-friendly.

Current scoped log layout:

- `logs/protenix-avb3/`
- `logs/boltz/`
- `logs/afcluster/`

Do not disable log tracking as part of “cleanup”.

## 6) Refactor Guardrails
When refactoring Claude-owned directories (e.g., `claude-AFCluster/`, `claude-Boltz/`):

1. Align implementations to canonical pipeline roots under `pipelines/`.
2. Preserve path contracts in existing submit scripts unless intentionally versioned.
3. Keep pipeline boundaries isolated; do not couple A5B1 and AVB3 streams.
4. Avoid introducing compatibility shims unless explicitly requested.
5. Prefer small, stratified commits:
   - pipeline/code moves
   - data/log path handling
   - docs/planning updates

## 7) Memory + Planning Updates
Before major changes, read:

- `AGENTS.md`
- `tasks/planning.md`
- `tasks/lessons.md`

After major refactors or repeated failure fixes, update `tasks/planning.md` and/or `tasks/lessons.md`.

## 8) Quick Validation Checklist
After refactor work:

1. `bash -n` on touched shell submit scripts.
2. `python -m py_compile` on touched Python scripts.
3. Verify no new references to legacy paths (`codex-*`, `data/outputs`, `data/template_example`, `data/sequences_updated`).
4. Verify submit scripts still map to correct venv (`venv_protenix` vs `venv_boltz`).
