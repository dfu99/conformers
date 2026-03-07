# Lessons

## Purpose
Capture recurring operational/debugging lessons for existing pipelines.

## Current Lessons

### Protenix Template Detection
- Symptom: logs show `Found 0 templates for sequence` despite `use_template=true`.
- Common causes:
1. Template entry not registered in Protenix mmcif DB layout (`mmcif/<id[1:3]>/<id>.cif.gz`).
2. Template entry missing from `common/release_date_cache.json`.
3. a3m hit headers not matching expected format (`ENTRY_CHAIN`).
4. Chain IDs in template mmCIF not matching query chain IDs.

### Protenix OOM on Smaller GPUs
- Symptom: CUDA OOM on ~24 GB GPUs for full integrin complexes.
- Observed behavior: Protenix can auto-fallback to fp32/torch kernels on older compute capability GPUs, increasing memory demand.
- Practical implication: A100 80GB class resources are preferred for large template-enabled runs.

### mkdssp Availability
- mkdssp is not strictly mandatory in your current workflow if gemmi fallback is active.
- Missing mkdssp alone does not explain zero-template detection when registration/parsing is wrong.

### Environment Isolation
- Symptom: conflicting dependency requirements across Protenix and Boltz (e.g., gemmi version mismatch).
- Resolution: enforce separate virtual environments per pipeline (`venv_protenix` vs `venv_boltz`) in submit scripts.
- Practical implication: avoid sharing one venv between Protenix and Boltz/AFCluster workflows.

### Protenix-Dock Scope
- Symptom: trying to use Protenix-Dock for SpyTag/Streptavidin attachment to A5B1.
- Likely cause: tool mismatch. Protenix-Dock CLI expects protein + ligand SDF style inputs (protein-ligand docking), not full protein-protein docking between integrin and protein tags.
- Action: use staged multichain Protenix inference for each attachment branch, then merge branches by receptor alignment.

### A5B1 Tagged Assembly Strategy
- Symptom: one-shot 4-chain inference can place tags in unrealistic internal positions.
- Likely cause: unconstrained global co-folding of all chains at once.
- Action: split into two stage predictions (`A5B1+SpyTag`, `A5B1+Streptavidin`) and merge ligand placements onto an accepted heterodimer receptor frame.

## Runbook Reminder
Before diagnosing new failures, verify in order:
1. Input JSON points to expected MSA/template files.
2. Template mmcif entry exists at registered path.
3. Release date cache contains custom entry ID.
4. a3m headers match entry and chain IDs.
5. Runtime logs confirm template parsing and non-zero template hits.

### BoltzGen CLI Contract
- Symptom: scripts invoked `boltz predict` or `boltzgen run` with `--out_dir`, causing mismatch with BoltzGen docs.
- Likely cause: confusion between Boltz and BoltzGen CLIs.
- Action: for BoltzGen use `boltzgen run <spec.yaml> --output <dir> --protocol <name> --num_designs <N> --budget <K>`.
- Practical implication: keep backend-explicit scripts (`boltzgen` vs `boltz`) and avoid mixing command/flag conventions.

### AFCluster API Drift (`max_dist` Constructor Error)
- Symptom: `TypeError: AFCluster.__init__() got an unexpected keyword argument 'max_dist'`.
- Likely cause: installed `afcluster` version uses a different constructor signature than expected by script.
- Action: introspect AFCluster constructor/cluster signatures and pass only supported kwargs at runtime.
- Practical implication: keep `cluster_chain_msa.py` API-compatible across AFCluster releases.
