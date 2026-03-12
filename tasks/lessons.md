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

### PACE SSH Timeout During Non-Interactive Automation
- Command context: `pipelines/protenix-a5b1/scripts/pace_minimal.sh smoke 5`.
- Symptom: `ssh: connect to host pace port 22: Operation timed out`.
- Likely cause: VPN/session routing not active to PACE login host, or SSH alias/host not reachable from current network.
- Action: run `pace_minimal.sh check` first, verify VPN is connected, and confirm `PACE_HOST` resolves to the correct PACE login endpoint.

### PACE Non-Interactive Auth Requires Authorized SSH Key
- Command context: `PACE_HOST='dfu71@login-phoenix.pace.gatech.edu' pipelines/protenix-a5b1/scripts/pace_minimal.sh check`.
- Symptom: `Permission denied (publickey,gssapi-keyex,gssapi-with-mic,password,hostbased)` in BatchMode.
- Likely cause: local private key exists but matching public key is not installed in remote `~/.ssh/authorized_keys`.
- Action: perform one-time key enrollment (`ssh-copy-id` or manual append) and then re-run `pace_minimal.sh check`.

### PACE Smoke Submit Requires Explicit Account
- Command context: `pace_minimal.sh smoke 5 60`.
- Symptom: `sbatch: error: --account option required`.
- Likely cause: cluster policy requires account on ad hoc smoke jobs.
- Action: pass `-A` in smoke submit (`PACE_ACCOUNT`, default `gts-yke8`).

### SSH Quoting Bug in Remote Watch Loop
- Command context: `pace_minimal.sh smoke` right after successful submit.
- Symptom: local shell error `line 164: JOB_ID: unbound variable` before watch/fetch.
- Likely cause: remote `$JOB_ID`/`$POLL_SECONDS` expanded locally under `set -u`.
- Action: escape remote variable references (`\$JOB_ID`, `\$POLL_SECONDS`) in the watch command string.

### A5B1 Staged Constraints Were Prepared But Not Enforced
- Command context: `run_complete_tagged_pipeline.sh` with staged attachment inputs.
- Symptom: tags can land near integrin head despite tail-intent manifests.
- Likely cause: tail constraints are written to JSON/manifests but Protenix stage inputs keep `covalent_bonds: []`, and runner calls `protenix pred` without consuming constraint files.
- Action: enforce tail-aware sample selection during merge (`tail_distance`/`hybrid`) and expose anchor/tail controls in runner.

### Ranking-Only Selection Can Oppose Tail Objective
- Command context: tail-distance diagnostics over synced stage outputs (`seed_101`, samples 0-4).
- Symptom: highest ranking sample was `sample_0` for both stages, while closest-to-tail sample was `sample_4`.
- Likely cause: ranking objective does not encode tail proximity directly.
- Action: when tail placement matters, use `SELECTION_MODE=tail_distance` or `SELECTION_MODE=hybrid` with max tail-distance thresholds.

### RTX6000 Incompatibility for A5B1 Protenix Production Runs
- Command context: `a5b1_conjugates_first` on `gpu-rtx6000` (job `4738472`).
- Symptom: runtime forced `dtype: fp32` on Compute Capability 7.x and OOMed in pairformer/template attention (`torch.OutOfMemoryError`).
- Likely cause: Protenix compatibility path disables bf16/deepspeed kernels on RTX 6000, increasing peak memory beyond ~24 GB.
- Action: run A5B1 production jobs on A100 80GB (`--partition=gpu-a100 --gres=gpu:A100:1 --constraint=A100-80GB`).

### Heterodimer Predictions Path Drift Between Local/Remote Layouts
- Command context: `a5b1_tagged_full` tail-aware submit (job `4738744`).
- Symptom: runner exited immediately with missing predictions dir under `data/runs/a5b1/protenix/...`.
- Likely cause: remote accepted heterodimer outputs still lived under legacy `data/a5b1/outputs/...` path.
- Action: add runner fallback to legacy predictions path (or set `HETERODIMER_CIF` explicitly).

### Conjugates-First Initially Produced No Combined Complex Artifact
- Command context: `run_conjugates_first_pipeline.sh` output inspection after successful A100 job.
- Symptom: only `conjugates_first_summary.json` existed; no merged heterodimer+tags CIF/PDB.
- Likely cause: workflow stopped after selecting best stage1/2/3 predictions and never executed a merge.
- Action: add explicit merge step using `merge_staged_tagged_complex.py` with stage-specific receptor mapping (`A:A` for stage1, `B:A` for stage2) and emit combined CIF/PDB + merge summary.
