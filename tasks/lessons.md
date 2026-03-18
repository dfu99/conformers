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

### `sbatch --export` Comma Pitfall for CSV Environment Lists
- Command context: auto-search submit with `--export=ALL,SEED_LIST_CSV=101,202,...,STAGE2_ANCHOR_LIST_CSV=1,10,...`.
- Symptom: only first value of each list reached the job (`SEED_LIST_CSV=101`, `STAGE2_ANCHOR_LIST_CSV=1`), causing a single-attempt run.
- Likely cause: `--export` itself is comma-delimited; unescaped commas split variable assignments.
- Action: for comma-delimited values, set them in script defaults or pass via env file/alternate delimiter; avoid raw CSV commas directly inside `--export`.

### Current A5B1 Complete Merge Still Fails Tail-Proximity Checks (Seed 101)
- Command context: `check_tagged_structure_quality.py` on merged final output.
- Symptom: 4 chains are present, but distances remain large (`A tail->D ~92-101 A`, `B tail->C ~50 A`) across seed-101 anchor sweeps.
- Likely cause: stage-2 placements are still head-biased and stage-1 alignment to base frame is poor for this seed set.
- Action: run multi-seed auto-search (`run_complete_tagged_until_pass.sh`) and require pass/fail gating on tail distances.

### Stage-1 Anchor Range Must Match Ligand Chain Residues
- Command context: `run_complete_tagged_until_pass.sh` with `SELECTION_MODE=tail_distance`.
- Symptom: repeated merge failures like `ValueError: Residue 24 with CA not found in chain 'C'`.
- Likely cause: stage-1 ligand chain `C` has CA residues only for `1..16`; broader anchor lists include invalid residue IDs.
- Action: keep stage-1 anchor sweep constrained to valid residue IDs (`1..16`) and stage-2 sweep targeted unless evidence shows sensitivity.

### Large A100 Seed Sweeps Should Use Per-Seed `sbatch --export`
- Command context: launching many new seeds after cached-seed sweeps failed strict thresholds.
- Symptom: comma-delimited seed lists are awkward in `sbatch --export` and easy to misconfigure.
- Likely cause: `--export` comma parsing conflicts with CSV payloads.
- Action: submit one seed per job (`SEED_LIST_CSV=<single_seed>`) and let each job sweep stage anchors locally.

### Additional Seeds (1001/1111/1222/1333/1444/1555) Did Not Improve Tail Geometry
- Command context: six A100 jobs (`4901896`..`4901901`) with `AUTO_SEARCH_UNTIL_PASS=1`, stage1 anchors `1..16`, stage2 anchor `1`.
- Symptom: all jobs exited `FAILED` with script exit `2` (no strict pass found).
- Observed best outcomes: still from cached-seed runs near `A~26.44 A`, `B~26.12 A` (`seed 404`), while new seeds were substantially worse.
- Action: next search space must change model inputs/objective (not just seed) if strict `<=25 A` is required.

### Seed-404 Stage2 Anchor Sweep Appears Geometry-Invariant
- Command context: targeted cached-prediction sweep for `seed 404`, `stage1 anchor=1`, stage2 anchors `1..N`.
- Symptom: repeated reports with identical distances (`A=26.4368 A`, `B=26.1193 A`) across tested stage2 anchors.
- Likely cause: stage2 sample choice is unchanged under tail-distance scoring for this seed/configuration.
- Action: deprioritize wider stage2-anchor brute force for current input setup; focus on changing base heterodimer/template constraints or merge objective.

### Base Heterodimer Frame Is the Dominant Lever for Tail Proximity
- Command context: post-alignment rescoring over cached stage predictions against alternate heterodimer base samples (`integrin_alpha5_beta1_sample_0..4`).
- Symptom: with base `sample_0`, best achievable max tail distance stayed at `26.44 A`; with base `sample_2`, best achievable max dropped to `16.06 A`.
- Likely cause: receptor frame choice shifts aligned ligand placement enough to determine whether strict tail-distance gates are feasible.
- Action: include base heterodimer sample/frame sweep in the automation before launching expensive new-seed inference.

### Strict Passing 4-Chain Structure Achieved via Base-2 + Seed-404 Best Samples
- Command context: forced merge of stage1 `seed_404 sample_2` + stage2 `seed_404 sample_4` onto base `integrin_alpha5_beta1_sample_2`.
- Outcome: 4 chains present, strict check passed (`alpha=16.06 A`, `beta=6.73 A`, threshold `<=25 A`).
- Output: `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete_base2_seed404_forced.{cif,pdb}` and copied to canonical `a5b1_tagged_complete.{cif,pdb}` with matching quality JSON.

### CA-Only Tail Checks Can Accept Geometrically Wrong AviTag/Streptavidin Placement
- Command context: manual inspection of canonical final PDB after `<=25 A` tail gating.
- Symptom: chain D remained close to alpha-leg body residues (around `A786/A798`) despite tail-distance pass.
- Likely cause: CA tail proximity alone does not penalize partner proximity to disallowed receptor domains.
- Action: add stricter chemistry-aware and exclusion-aware checks:
  - `A966:NZ -> chain D` max distance
  - `B735:NZ -> C10:[OD1,OD2,CG]` max distance
  - `D` minimum distance from alpha disallowed range (`A780-820`) must exceed threshold.

### Seed-404 + Base-Sweep Still Fails Under Stricter Exclusion Criteria
- Command context: `run_complete_tagged_until_pass.sh` with `seed=404`, `base=0..4`, strict gates (`tail<=20 A`, attachment<=`6 A`, `A780-820` exclusion>=`10 A`), `FORCE_RERUN=0`.
- Symptom: all 80 attempts failed strict gating.
- Likely cause: cached stage predictions for that seed cannot satisfy both tail/chemistry and alpha-body exclusion simultaneously.
- Action: continue with broader multi-seed cached sweep and in parallel submit fresh A100 reruns with stricter gate as acceptance criterion.

### Use Separate Workflow Roots for Fresh Strict Reruns
- Command context: launching strict A100 reruns while a cached strict sweep is active.
- Symptom risk: concurrent runs can overwrite shared `outputs/final` and trial artifacts when using one `WORKFLOW_DIR`.
- Action: submit reruns with per-seed isolated workflow directories (e.g., `.../staged_attachment_strict_rerun_<seed>`) to avoid artifact collisions.

### Full-Day A100 Strict Sweeps Are a Poor Default for A5B1 Tagged Search
- Command context: strict rerun jobs `4925696`, `4925697`, `4925698` on A100 80GB, each with `--time=24:00:00`, `FORCE_RERUN=1`.
- Symptom: all three jobs hit walltime after only ~36 attempts and still produced zero strict passes.
- Likely cause: per-attempt stage reruns are too expensive for broad brute-force search, and the core issue is likely model/geometry mismatch rather than insufficient seed count.
- Action: cap future Protenix search jobs with `MAX_ATTEMPTS` and prefer cheaper heterodimer base-frame exploration before spending more A100 hours.

### A5B1 Heterodimer Base Search Should Reuse Cached Tag Branches
- Command context: after stricter gates rejected all Protenix cached/rerun attempts.
- Symptom: repeated failures suggest the dominant unknown is the heterodimer base frame, not the downstream merge mechanics alone.
- Likely cause: alternative A/B backbone conformations are needed, but rerunning full staged tag inference for every trial is too expensive.
- Action: sample base heterodimers with cheaper tools (`BoltzGen`, `AFCluster -> BoltzGen`) and then reuse cached stage predictions for strict merge rescoring.

### Sequential RTX6000 Jobs Are Preferable to Concurrent A100 Jobs for Exploratory BoltzGen Runs
- Command context: A5B1 budget-constrained heterodimer exploration on March 14, 2026.
- Symptom/risk: concurrent 24h A100 jobs consumed budget too quickly relative to information gained.
- Action: submit short `gpu-rtx6000` jobs sequentially with dependencies (one active GPU at a time), using very small search settings first (`top_a/top_b`, `num_designs`, `budget`).

### Current BoltzGen Helpers Emit Legacy Boltz YAML, Not Installed BoltzGen Schema
- Command context: direct BoltzGen job `4951450` and AFCluster -> BoltzGen job `4951452` on March 14, 2026.
- Symptom: both jobs failed quickly during `boltzgen run` config validation with `ValueError: Found invalid keys in yaml file: {'version', 'sequences', ...}`.
- Likely cause: `pipelines/boltz/scripts/build_boltz_predict_sweep.py` and `pipelines/afcluster/scripts/make_boltz_jobs_from_clusters.py` still write legacy Boltz-style YAML (`version`, `sequences`, nested `templates`) while the installed `boltzgen` expects the newer schema rooted in `entities`.
- Important nuance: AFCluster itself worked for A5B1 under the tiny test and produced clustered A3Ms plus a manifest; the failure happened only at the BoltzGen spec-validation boundary.
- Action: fix the YAML generators to emit the installed BoltzGen schema before submitting any more BoltzGen or AFCluster->BoltzGen jobs.

### `boltzgen check` Is the Right Cheap Preflight for YAML Spec Changes
- Command context: PACE validation on March 15, 2026 after patching the YAML emitters.
- Symptom/risk: schema changes can look plausible by inspection but still fail later when the installed parser resolves files and canonical residue data.
- Action: run `boltzgen check <spec.yaml>` on representative direct and AFCluster-generated specs before submitting any GPU jobs.
- Practical implication: this catches parser/schema regressions cheaply on the login node and avoids wasting GPU allocation on immediate config failures.

### AFCluster Size-Based Ranking Over-Prioritizes `cluster_outliers.a3m`
- Command context: inspecting `clusters.tsv` from the tiny A5B1 AFCluster run (`afcluster_boltzgen_24776`) during the March 15, 2026 BoltzGen fix.
- Symptom: chain A `clusters.tsv` ranked `cluster_outliers.a3m` ahead of `cluster_0.a3m`, so the nominal “top” AFCluster BoltzGen job could be driven by outlier sequences.
- Likely cause: `cluster_chain_msa.py` ranked cluster A3Ms purely by file size, and `make_boltz_jobs_from_clusters.py` consumed `clusters.tsv` in that order.
- Action: rank numbered clusters ahead of `cluster_outliers.a3m`, and treat the outlier bucket as fallback only when no numbered cluster exists.
- Practical implication: AFCluster budget runs now spend their first attempts on bona fide clusters instead of the outlier bucket.

### BoltzGen On PACE Must Not Use Home-Directory Hugging Face Cache
- Command context: patched BoltzGen reruns on March 15, 2026 (`4957942`, `4957943`) after schema validation was fixed.
- Symptom: both jobs passed `boltzgen check`, started model staging, and then failed with `Disk quota exceeded` while `huggingface_hub` tried to write into `~/.cache/huggingface`.
- Observed state: home quota was exactly full (`20480M / 20480M`); the cached BoltzGen model tree under `~/.cache/huggingface/hub/models--boltzgen--boltzgen-1` was already about `3.7G`.
- Action: default `HF_HOME`, `HF_HUB_CACHE`, and `HF_XET_CACHE` to scratch in the BoltzGen launchers.
- Practical implication: on PACE, any BoltzGen job should use scratch-backed HF cache unless there is confirmed free home quota.

### `boltzgen check` Writes Visualization CIFs To The Current Working Directory Unless `--output` Is Set
- Command context: manual schema validation on March 15, 2026 from the PACE home directory.
- Symptom: `a5b1_schema_check_*.cif` and `a5b1_afc_schema_check_*.cif` were written into `~` instead of under the scratch run folders.
- Likely cause: `boltzgen check` defaults to writing the visualization mmCIF into the current working directory when no `--output` directory is provided.
- Action: when validating on PACE, either `cd ~/scratch/...` first or pass `--output ~/scratch/...` explicitly.
- Practical implication: avoid cluttering or exhausting home quota with validation artifacts.

### MSA Subsampling Can Validate Conformer Realism (Literature)
- Source: Wayment-Steele et al., Nature Comms 2024; AFsample2
- Key insight: subsampling MSA depth causes AF2 to explore different conformational states; predicted frequency correlates with experimental populations (>80% accuracy vs NMR).
- Implication for pulled conformers: a valid extended state should be "reachable" by some MSA subsample with decent pLDDT. An over-stretched invalid state should have no MSA depth that agrees.
- Action: use progressive MSA subsampling (100% → 5%) as a cheap conformer validation filter before investing in full pipeline runs.

### ProteinTTT Can Bypass HuggingFace Residue Limits Locally
- Source: anton-bushuiev/ProteinTTT (ICLR 2026)
- Symptom: HuggingFace Space rejects AVB3 (>400 residues).
- Action: install locally and run on PACE A100 80GB — no hard residue limit in the code, only GPU memory.
- Note: ProteinTTT fine-tunes ESMFold per-protein via test-time training, reportedly 2x pLDDT improvement on hard targets.

### BoltzGen Design Step OOMs on RTX 6000 for Full A5B1 Integrin
- Command context: BoltzGen jobs 4963854 and 4963855 on RTX 6000.
- Symptom: all 64 diffusion batches hit "ran out of memory" warnings. Zero designs produced. Inverse folding then fails with "No designs found".
- Likely cause: A5B1 heterodimer (~1600+ residues) exceeds RTX 6000 VRAM (~24GB) for BoltzGen's diffusion step, and compute capability 7.5 disables optimized kernels.
- Action: BoltzGen for full-size integrins must use A100 80GB. Alternatively, run on individual chains/domains.
