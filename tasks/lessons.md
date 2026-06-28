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

### Temporal Flip Resolution Must Anchor to Global Targets, Not Previous Frames
- Command context: `resolve_flips()` in `fit_with_head_tracking.py` for HS-AFM overlay post-processing.
- Symptom: head positions in smoothed coordinates drifted up to 10.4nm from where SO(3) fitting placed them. 399/1266 frames drifted >1nm. Visible as overlay misalignment at frames 174, 606, 678, 750.
- Likely cause: comparing each frame to the previous frame propagates small errors, creating cumulative drift.
- Action: use `resolve_flips_head_anchored()` which compares to the tracked AFM head position (a global anchor), then re-centers. Result: zero drift.

### Pip-Installed OpenMM May Lack CUDA Platform
- Command context: running OpenMM steering MD on RunPod A4500.
- Symptom: `Platform.getNumPlatforms()` returned only Reference and CPU — no CUDA despite `nvidia-smi` showing a working GPU. Simulation silently fell back to CPU (would take months).
- Likely cause: `pip install openmm` gives the CPU-only build; CUDA kernels are a separate package.
- Action: always `pip install openmm-cuda` after installing OpenMM on GPU machines. Verify with `Platform.getPlatformByName("CUDA")` before launching long runs.

### Scipy NetCDF Writer Has 2GB File Size Limit
- Command context: mdtraj NetCDFReporter writing `production.nc` during OpenMM steering MD.
- Symptom: process died silently at exactly step 199,500 (same point on two independent runs). The production.nc file was 1.9 GB when it crashed.
- Likely cause: scipy's NetCDF implementation uses int32 file offsets, causing a crash near the 2 GB limit.
- Action: `pip install netCDF4` before running MD. The mdtraj warning about missing netCDF4 is not cosmetic — scipy's fallback breaks on large trajectories.

### Conformer Library Must Be Validated Against Target CV Range Before Fitting
- Command context: HS-AFM overlay showing bent conformers at frames where AFM shows extension.
- Symptom: at frame 1110, correlation matching selected a compact conformer (CV0=53A) with r=0.72. AFM clearly shows extended features.
- Likely cause: the pseudo-AFM library used only 66/264 conformers, all below CV0=79A. No extended states existed to match against.
- Action: before running the full pipeline on new AFM data, verify that the conformer library spans the expected CV range. Use the CV distribution histogram to check coverage.

### RunPod Workspace Disk Quotas Kill Long MD Runs Silently
- Command context: 3ns OpenMM steering MD on RunPod A4500 with `--save-frames`.
- Symptom: process died silently after ~40 PDB frames (5.2 GB), with "Disk quota exceeded" on subsequent writes. No error message in the log, no OOM signal.
- Likely cause: each solvated PDB frame is 131 MB (430K atoms). 40 frames = 5.2 GB. Full 3ns = 400 frames = 52 GB — far exceeds RunPod per-pod volume limits (~20 GB network storage).
- Action: never use `--save-frames` for solvated systems on RunPod. Rely on `production.nc` trajectory (compressed, ~5 GB for 3ns) and extract protein-only frames post-hoc with mdtraj. Add disk monitoring to cron watchdog.

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

### BoltzGen cuBLAS Crash on A100 for Full Integrin — Not Just an OOM Issue
- Command context: BoltzGen A100 jobs 5043221 and 5043229.
- Symptom: `RuntimeError: CUDA error: CUBLAS_STATUS_INVALID_VALUE` in `boltzgen/model/modules/masker.py:58` at `torch.bmm()`. Crashes on the very first batch, not OOM.
- Likely cause: full A5B1 integrin (~1600+ residues) produces token/pair tensors that exceed cuBLAS matrix dimension limits for bf16 batched matrix multiply. A100 uses optimized kernels (capability 8.0) and batch_size=10, but the tensor dimensions are too large.
- Key observation: RTX 6000 failed with OOM (didn't even attempt computation); A100 got further (kernels=True, batch_size=10) but hit a cuBLAS numerical limit.
- Action: BoltzGen cannot handle full-size integrins in its current form. Options: (1) run on individual chains/domains, (2) check if BoltzGen has a `--max_tokens` or chunking option, (3) file a bug report with BoltzGen team.

### Protenix Is MSA-Depth-Invariant for AVB3 Conformational Diversity
- Command context: MSA validation test job 5073313, full MSA (7425+2068 seqs) vs 5% subsample (371+103 seqs).
- Symptom: full and 5% MSA produce identical predictions — TM-score between any frame and either depth's prediction differs by <0.001.
- Likely cause: Protenix's architecture may collapse MSA information more aggressively than AF2's evoformer, or the AVB3 bent conformation is so dominant in training data that MSA depth doesn't matter.
- Action: use AF2 (which showed MSA-depth sensitivity in Wayment-Steele et al.) instead of Protenix for conformational diversity via MSA subsampling.
- Practical implication: Protenix TM-score still works as a conformer validity filter (monotonic 0.96→0.06), just not as a conformational diversity generator.

### OpenMM pip CUDA Requires Matched openmm + openmm-cuda + numpy<2
- Command context: RunPod A5000, pip-installed OpenMM 8.5.0 with numpy 2.4.4.
- Symptom: `openmm-cuda` package installs OpenMM 8.1.1 but system has numpy 2.x, causing `ImportError: numpy.core.multiarray failed to import`.
- Likely cause: `openmm-cuda` pins to older OpenMM which was compiled against numpy 1.x.
- Action: for pip-based OpenMM with CUDA, install `openmm-cuda` first (which pulls matching openmm), then `pip install "numpy<2"`. Verify with `Platform.getPlatformByName('CUDA')`.
- Conda avoids this entirely but isn't always available on RunPod.

### OpenMM ionicStrength Must Use Unit-Bearing Quantity
- Command context: domain steering solvation step on PACE (jobs 6195166-6195169).
- Symptom: `AttributeError: 'float' object has no attribute 'value_in_unit'` in `modeller.addSolvent()`.
- Likely cause: `ionicStrength=0.15 * kJ/mol / kJ/mol` evaluates to a bare float (0.15), not a Quantity with molar units.
- Action: use `ionicStrength=0.15 * molar` and import `molar` from `openmm.unit`.

### OpenMM Position Restraints + Centroid Pulling Is Numerically Unstable
- Command context: domain steering `restrained_pull` preset on RunPod A5000.
- Symptom: `OpenMMException: Particle coordinate is NaN` during equilibration, even after reducing k from 500→200 and pull force from 1.0→0.5 pN.
- Likely cause: harmonic position restraints on all CA atoms create competing forces with inter-domain pulling, causing the integrator to diverge.
- Action: this steering method needs gradual force ramping during equilibration, or should be replaced by the CV distance bias method which is stable.

### CV Distance Bias Is the Effective Steering Method for Integrin Opening
- Command context: comparison of 4 domain steering presets on AVB3 (RunPod A5000, 2026-04-07).
- Key result: gentle_open and moderate_open (centroid angle torques, k=50-200) produced negligible change (<5° angles, <4Å distances). cv_distance_extend (flat-bottom distance biases, k=200, targets 16-20nm) produced +88° leg opening and +113Å head-tail extension.
- Action: use `cv_distance_extend` preset for all future steering runs. Angle torques at these force constants are too weak to overcome integrin hinge stiffness in 1ns.

### AFMFold Pipeline: Steered MD Frames Can Replace Protenix Conditional Generation
- Command context: building AFMFold CNN pipeline for AVB3 (2026-04-08).
- Key insight: AFMFold's step 1 (candidate generation) uses Protenix guided inference to populate a CV grid. Our steered MD trajectory already covers the CV space continuously.
- Action: feed steering PDB frames directly into AFMFold's image generation (step 2) and CNN training (step 3), skipping the expensive Protenix generation step entirely.
- Bridge script: `process_frames_to_afm.py` already connects steering outputs to afmfold's `generate_images()`.

### OpenFold Requires GPU Node Install (nvcc)
- Command context: ProteinTTT venv setup on PACE login node.
- Symptom: `pip install openfold` fails with `FileNotFoundError: No such file or directory: '/usr/local/cuda/bin/nvcc'`.
- Action: install OpenFold from a SLURM GPU job with `module load cuda`, not from the login node.

### Protenix components.cif Download Can Corrupt on PACE
- Command context: MSA test job 5060745.
- Symptom: Protenix tried to download 490MB `components.cif` from Chinese CDN; got only 12MB due to network interruption. Corrupted file then blocked subsequent runs.
- Action: if Protenix fails with `ContentTooShortError`, delete `~/common/components.cif` and retry.

### AF2 Also Locked to Bent Conformation for AVB3 — No MSA-Depth Diversity
- Command context: AF2 2.3.2 multimer with reduced_dbs on AVB3 (job 5357284), 25 ranked predictions.
- Symptom: all 25 AF2 predictions have pairwise RMSD 0.1-2.9Å. TM-score vs pulled frames: 0.99 (frame 0, bent) → 0.06 (frame 300, impossible). Profile nearly identical to Protenix.
- Likely cause: AVB3 integrin bent conformation is overwhelmingly dominant in PDB training data. Unlike the proteins in Wayment-Steele et al. (which had comparable populations of open/closed states), AVB3's extended state is extremely rare in crystallography, so neither AF2 nor Protenix can sample it via MSA subsampling.
- Action: MSA subsampling for conformational diversity is unlikely to work for AVB3 specifically. Need physics-based approaches (steered MD) or test-time training (ProteinTTT) instead.
- Practical implication: the TM-score scoring methodology works perfectly as a conformer validity filter, but MSA-based conformer generation is not viable for mechanobiologically sensitive proteins with rare extended states.

### Standard AF2/Protenix Cannot Predict Unseen Conformations — Use Perturbation Methods
- Context: Literature review for alternative conformation prediction (March 2026).
- Key insight: AF2 reproduces PDB training data conformations. For rare states (like extended integrin), need specialized methods:
  1. AlphaFold-RandomWalk (AF-RW): noise injection into model weights, published 2026 Scripps Research.
  2. AFsample3: random MSA column masking for AF3, published Jan 2026.
  3. Energetic frustration + AF2: identifies where to "push" AF2 (PNAS 2024).
  4. Finite temperature string method: Dasetty et al. (bioRxiv 2025) solved bent→extended for αIIbβ3 integrin specifically.
- Action: prioritize AF-RW and the αIIbβ3 string method paper as most directly applicable to AVB3.

### Additional Methods for Integrin Extended Conformations (Literature Review March 2026)

#### Targeted MD / Metadynamics for Integrins
- **Funnel Metadynamics on αIIbβ3** (PMC 2024): FM applied to αIIbβ3 + RGD peptide. Does not require a priori information about interactions. Could be applied to αVβ3 for enhanced sampling of activation pathway.
- **Force-Regulated Conformational Changes of αVβ3** (ACS Nano 2024): SMD simulations on αVβ3 found hexa-stable intermediate states with ~7 hydrogen bonds disrupted consecutively during unbending. Energy landscape mapped with intermediate barriers. Directly applicable — uses the same αVβ3 we're studying.
- **Conformational response of αIIbβ3 and αVβ3 to force** (Biophysj 2023): Equilibrium and enhanced MD simulations comparing both integrins. Provides reference force profiles.
- Action: the ACS Nano 2024 paper's intermediate states could serve as targets for our steered MD frames — validating which pulled frames correspond to known intermediates.

#### Coarse-Grained Models
- **GōMartini 3** (Nature Comms 2025): Combines Go-model structure-based potentials with MARTINI 3 coarse-grained force field. Demonstrated for large conformational changes, protein-membrane binding, and AFM force profile calculations. Directly relevant for integrin bent→extended transitions.
- **Switching Gō-Martini** (JCTC 2024): Simulates large-scale protein conformational transitions between different states. Uses switching between Go-model contacts for two endpoint structures. Could model bent↔extended integrin transition.
- **SBMOpenMM** (JCIM 2021): Structure-based model builder for OpenMM. Could implement Go-like models for integrin with our existing OpenMM installation on PACE.
- Action: GōMartini 3 is the most promising CG approach. It would require GROMACS (not OpenMM) but could sample the bent→extended transition orders of magnitude faster than all-atom MD.

#### String Method Implementations
- **Ferg-Lab/principalcurve_integrin_structures** (already cloned on PACE): The αIIbβ3 string method uses OpenMM + SEEKR2. Their code is available but requires significant adaptation for αVβ3 (different membrane setup, force field parameterization, collective variables).
- **SBMOpenMM** could be used to build a simplified string method for integrin using structure-based models, avoiding the full membrane embedding.
- Action: adapting the Ferg-Lab string method to αVβ3 is a major project (weeks). For near-term, continue with homology mapping from their αIIbβ3 results.

#### Feasibility Summary
| Method | Feasibility | Timeline | Value |
|--------|------------|----------|-------|
| αIIbβ3 homology mapping | Done | Immediate | Medium — different α subunit |
| AF2-RW weight perturbation | Submitted | Days | Unknown — may not escape bent |
| GōMartini 3 CG MD | Needs GROMACS | 1-2 weeks | High — can sample transition |
| Metadynamics on αVβ3 | Needs setup | 1-2 weeks | High — enhanced sampling |
| Full string method for αVβ3 | Major project | Weeks | Very high — MFE pathway |
| ProteinTTT | Install pending | Days | Unknown — single-sequence |

### AFMFold Pseudo-AFM Tip Radius Must Match Experimental HS-AFM
- Command context: CNN training with `min_tip_radius=6.0, max_tip_radius=12.0` (pixels at 0.98 nm/px = 5.9–11.8 nm).
- Symptom: CNN predicted near-constant CVs (std=0.04 Å) extrapolating outside training range (82 Å vs max 79 Å).
- Root cause: experimental HS-AFM tip is 1–2 nm; our pseudo-AFM tip was 3–6× too large, producing featureless blobs.
- Fix: set `min_tip_radius=1.0, max_tip_radius=2.0` (≈1–2 nm). afmfold defaults (1–3 px) were actually correct.
- Result: corrected-tip CNN predicts within training range (mean 69 Å) but still near-constant (std=0.2 Å). Correlation matching with corrected-tip pseudo-AFM library is the better inference method (std=8.5 Å).
- Action: always match tip parameters to the experimental setup. For Linz HS-AFM: tip 1–2 nm.

### CNN Domain Gap Persists Even With Correct Tip — Use Correlation Matching
- Command context: corrected-tip CNN (val_loss=2.69 Å) inference on Linz HS-AFM GIFs.
- Symptom: predictions within training range but near-constant; no conformational discrimination.
- Likely cause: sim→real domain gap in noise, contrast, background. Pseudo-AFM images are clean height maps; real HS-AFM has substrate noise, scanning artifacts, variable contrast.
- Workaround: direct image correlation matching against the pseudo-AFM library captures real conformational dynamics (r up to 0.967).
- Next step to close gap: histogram matching during training, noise augmentation matching real HS-AFM characteristics, or fine-tuning with labeled real AFM data.

### 50 Pre-Saved SO(3) Rotations Are Too Sparse for Orientation Recovery
- Command context: rendering atomistic PDB overlays using 50 epoch rotations from pseudo-AFM library training.
- Symptom: persistent head-tail orientation flipping between adjacent frames, producing visually jarring GIFs.
- Root cause: 50 rotations cover ~2.4% of SO(3) uniformly — insufficient to find the orientation that best matches each frame. AFMFold's paper uses exhaustive sampling (2048+ per frame).
- Action: use `fit_rigid_body.py` with `--n-rotations 2048` for per-frame SO(3) search using AFMFold's correlation-based fitting.

### AFMFold `summarize_results` Bug With rot_batch > 1
- Command context: investigating `afmfold/src/afmfold/rigid_body_fitting.py` line 283.
- Symptom: `if i == best_rot_idx` compares step index `i` against flat rotation index `best_rot_idx` (which spans all steps × rot_batch).
- Likely cause: code assumes rot_batch=1 where step index equals flat rotation index.
- Action: wrote custom `fit_orientations()` in `fit_rigid_body.py` that tracks best rotations manually per frame instead of using the buggy `RigidBodyFitting.summarize_results()`.

### Piping Long-Running Python Through `head` Causes Silent SIGPIPE Kill
- Command context: launching `fit_rigid_body.py 2>&1 | head -30` as a background task.
- Symptom: output directory empty, process silently killed after producing >30 lines of output.
- Root cause: `head -30` closes the pipe after 30 lines, sending SIGPIPE to the Python process.
- Action: never pipe long-running GPU processes through `head`. Use `2>&1` only, or redirect to a log file.

### COM-Based Centering Offsets PDB From AFM Head
- Command context: `fit_rigid_body.py` v3 overlays.
- Symptom: PDB overlay visually offset from AFM features despite high per-frame correlations. Head domain displaced from brightest AFM feature.
- Likely cause: center-of-mass includes the entire molecule (head + legs). Extended legs shift the COM away from the globular head, causing 2+ nm offset.
- Action: track the head position in AFM frames (brightest topographic feature) and align PDB head domain to that position instead of full-structure COM. This improved correlations from 0.69→0.97 (video1) and 0.80→0.94 (video2).

### Optimal AFM Tip Radius is 2-3 nm
- Command context: `fit_with_head_tracking.py --compare-tips` on both Linz GIFs.
- Symptom: needed to determine which simulated tip radius best matches experimental data.
- Result: 2nm best for video2 (0.944), 3nm best for video1 (0.976). 1nm and 5nm both worse but close. All between 0.93-0.98 — tip radius is a secondary factor compared to centering.
- Action: default to 2-3nm tip for future rigid-body fitting against these Linz GIFs.

### First 30 HS-AFM Frames Are Scan Window Adjustment
- Command context: Linz AVB3 HS-AFM GIFs.
- Symptom: early frames show imaging artifacts as the scan window stabilizes on the specimen.
- Action: always skip the first 30 frames when fitting or analyzing these GIFs. Use `--skip-frames 30`.

### generate_images batch_size Truncation Bug
- Command context: `afmfold.images.generate_images()` via `process_frames_to_afm.py`.
- Symptom: pseudo-AFM library labels had narrow CV range [78.6, 80.8] Å despite conformer library spanning [52.9, 85.0] Å. Only first ~31 conformers were sampled.
- Root cause: with `batch_size=16` and `dataset_size=500`, each step generates `len(xyz) × batch_size = 309 × 16 = 4944` images, then truncates to 500. Since frame order is sequential, only conformers 0-30 survive truncation.
- Action: use `batch_size=1` so each step generates 309 images (one rotation per conformer). `ceil(500/309) = 2` steps per epoch ensures all conformers are covered. Also dramatically faster (1.1s vs 26s/epoch) and uses minimal memory.

### RunPod Disk Quota Silently Truncates np.save
- Command context: `fit_with_head_tracking.py` running on RunPod workspace.
- Symptom: script completed fitting and flip resolution but crashed without error during file save. `fitted_coords.npy` truncated to exactly 256 MiB (of expected 365 MiB).
- Root cause: RunPod workspace disk quota exceeded. Old steering backup files (12.5 GB) consumed the pod's quota. `np.save` gets an `OSError: Disk quota exceeded` during write but the error wasn't visible in the nohup log because Python crashed.
- Action: always clean up old data on RunPod before running jobs that produce large output. Check `du -sh /workspace/conformers/*` before starting.

### RunPod GPU vs CPU for SO(3) Fitting
- Command context: `fit_with_head_tracking.py` with 1266 frames × 2048 rotations.
- Result: A4500 GPU completed in 10 min vs ~4 hours on RunPod CPU (same machine). Always use `--device cuda` when GPU is available for fitting.

### Steering Forces Cause Minimization Hang
- Command context: `run_domain_steering.py` on RunPod A4500 / RTX 2000 Ada.
- Symptom: `simulation.minimizeEnergy()` hangs 30+ minutes at 100% CPU with no progress when steering forces are applied before minimization.
- Root cause: LBFGS minimization has to resolve the strong steering forces (k=200, targets 4-5nm from current state) simultaneously with the standard forcefield. The two force systems create overdetermined constraints.
- Action: always minimize the system BEFORE adding steering forces. Create initial system with just forcefield, run `minimizeEnergy()`, then add CV-bias forces and reinitialize the Simulation with the updated System. See `run_domain_steering.py` commit 266f3c6.

### pyparsing Breaks mdtraj.reporters Silently
- Command context: fresh RunPod pod, `pip install openmm-cuda mdtraj netCDF4`.
- Symptom: MD simulation runs for hours producing no `production.nc` trajectory file. Script doesn't error because it catches `ImportError` silently.
- Root cause: Ubuntu's system pyparsing (in /usr/lib/python3/dist-packages) is too old for mdtraj. mdtraj.reporters fails with `cannot import name 'OpAssoc' from 'pyparsing'`.
- Action: on fresh RunPod, always run `pip install --upgrade pyparsing` along with `pip install 'numpy<2' openmm-cuda mdtraj netCDF4`. Test with `python3 -c 'import mdtraj.reporters; print("OK")'` before launching MD.

### generate_images batch_size Truncates Conformer Sampling
- Already documented above, re-emphasizing: with `dataset_size=500, batch_size=16, N_frames=309`, `ceil(500/(309*16))=1` step produces 309*16=4944 images truncated to first 500 — only the first ~31 conformers are represented. Use `batch_size=1` and any `dataset_size >= N_frames` to get full coverage.

### RunPod Cleanup: Check Both Locations
- Command context: cleaning up after a project's work on a shared RunPod node.
- Symptom: PI audit caught 23 GB still on the pod after I reported it cleared.
- Root cause: I cleared `/workspace/<proj>` (where I had launched MD) but `mc runpod sync` had also pushed the entire repo to `/root/projects/<proj>` earlier. Two separate locations.
- Action: when cleaning up, check both `/workspace/<proj>` AND `/root/projects/<proj>`. `find / -maxdepth 6 -type d -path '*<proj>*'` catches anything else that snuck in.

### Sim AFM Realism Requires Multiple Low-Pass Stages
- Command context: forward-rendering fitted PDB trajectories to simulated HS-AFM.
- Symptom: sim AFM looked too sharp and too detailed compared to real Linz GIFs even with `idilation` (tip dilation) applied.
- Root cause: `idilation` is morphological-max — preserves edges sharply. Real HS-AFM has additional low-passes from cantilever bandwidth, feedback PID time constant, and tip-sample force response.
- Action: post-dilation anisotropic Gaussian blur (σ_y ≈ 0.7, σ_x ≈ 1.2 — wider along scan direction). Plus Gaussian sensor noise σ ~0.08. Plus soft-clip at ~0.92 (real AFM never saturates fully white). All three combine in the sim_afm_v11 renderer.

### Surface-Adsorbed Molecules Should Lie Flat (PCA-Constrained)
- Command context: HS-AFM specimens are physisorbed on mica/glass — they LIE on a side, head/knee/tail all in surface contact.
- Symptom: per-frame SO(3)-fitted PDBs sometimes had the molecule oriented with one domain up and another on the surface.
- Action: per-frame PCA of CA coords. Rotate so smallest principal axis (thinnest direction) is +Z. Sign-align PCA axes between consecutive frames (or you get 180° flips when SVD picks opposite signs). The remaining DoF is yaw around Z.

### Step-wise Yaw Beats Continuous for Surface-Bound Molecules
- Command context: real HS-AFM shows molecules holding an orientation for seconds, then snapping to a new one — never spinning smoothly.
- Action: hysteresis-based step detector on the raw yaw signal. Threshold ~50°, min dwell 20 frames, cap each step magnitude at ~30°. Gives "prefers an orientation, occasionally re-orients" behavior. Smooth (Gaussian) yaw produces visually-wrong continuous tumbling.

### Steering Presets Need Explicit Pair Lists
- Command context: `cv_distance_headopen` preset on `domain_steering.py`.
- Symptom: 790 ps run, target name says "head-head opening", but CV2 (α-head ↔ β-head separation) stayed at 34 ± 0.4 Å throughout.
- Root cause: the preset's `target_values` list was applied to AVB3_HINGE_DISTANCES (the DEFAULT pair list) which contains only head-tail pairs — no head-head pair. The preset's *name* claimed head-head opening but the underlying pairs didn't include head-head.
- Action: define a parallel pair list (`AVB3_HEADOPEN_DISTANCES`) with the (α-head, β-head) pair at index 0. Update `apply_steering_preset` to route preset name → correct pair list.

### Single-Canvas Stylization: No Inset Paste Between Mismatched Sizes
- Command context: forward-rendering pseudo-AFM with substrate noise + molecule render.
- Symptom: v11 produced visible rectangular border around the inner molecule whenever the molecule extended past its inset edge — clearly an internal cutoff artifact in supposedly natural-looking output.
- Root cause: noise was generated on a larger canvas and the molecule render was pasted into a smaller inner inset. Pasted boundaries always show.
- Action: do EVERY stylization step (zoom, substrate noise, anisotropic blur, slant, row jitter, flash streaks, soft-clip, colormap) on ONE canvas at ONE resolution. See `pipelines/sim-afm-video/scripts/stylize_sim_afm.py`. Never composite via paste-into-larger-canvas.
- Practical implication: any future image stylization that has both noise field + signal must keep them on the same array from the start.

### Matplotlib Image Subplots Need Explicit Aspect Lock
- Command context: 3-panel figure with figsize=(13, 4.2) showing AFM frame + fitted PDB overlay + CV trajectory at width_ratios=[1, 1, 1.5].
- Symptom: AFM and overlay panels rendered horizontally squished — visibly non-square despite displaying square arrays.
- Root cause: matplotlib's default sizing fills the available axes box per subplot, so the first two panels get stretched horizontally to match the figure's aspect.
- Action: call `ax.set_aspect('equal', adjustable='box')` on every panel that displays an image. The 'box' adjustable shrinks the axes box to match the data aspect; without it the data gets stretched to fill the box.
- Practical implication: any image-display panel inside a wide figure needs the explicit aspect lock — the default behavior is wrong for this case.

### Sim AFM Z-Axis Needs Explicit Blur + Noise to Match Real HS-AFM
- Command context: forward-render of fitted PDB through `simulate_afm_video.py`.
- Symptom: sim AFM had crisp, precise z-resolution that no real HS-AFM has — molecule features looked like a digital height map rather than a measurement.
- Root cause: pipeline was applying only xy blur (~0.1 nm) and minimal z-noise (0.05 nm). Real HS-AFM has cantilever-oscillation averaging plus 0.3-0.5 nm RMS z-noise from feedback loop dynamics.
- Action: apply Gaussian σ in nm to the dilated height map BEFORE per-frame normalize (new `--z-blur-nm` flag) AND increase z-noise floor to match the real instrument (`--noise-nm 0.35`). Validated config: `--z-blur-nm 0.8 --noise-nm 0.35`. Documented in `pipelines/sim-afm-video/README.md`.
- Practical implication: forward-render fidelity depends on modelling z-axis bandwidth limits, not just xy tip dilation.

### Anisotropic Blur Is a Workaround for Unstable Yaw — Not a Default
- Command context: sim AFM stylization v11 used σx=1.2 / σy=0.7 to compensate for free-tumbling renders.
- Symptom: after PCA-flatten + side-lock + stepwise yaw landed (v13+), the long axis is consistently along x; anisotropic blur visibly squashed the molecule horizontally.
- Root cause: the asymmetric blur was tuned for unstabilized frames where molecule could lie at any angle. Once orientation is locked, the asymmetry becomes a directional artifact.
- Action: default to isotropic σx=σy=1.0 once orientation is locked. The earlier defaults are only justified for non-locked rendering.
- Practical implication: parameters tuned for one stage of a pipeline can become bugs after a downstream stage changes assumptions; revisit defaults after each meaningful upstream change.

### RunPod A4500 Pod Quota Hits at ~1.6 GB Cumulative
- Command context: 3 ns OpenMM MD on solvated αVβ3 with `--report-interval 1500`.
- Symptom: process dies silently with no error after ~1.5 hours of production. production.nc plateaus at ~1.4 GB. Total project dir ~1.6 GB.
- Action: (a) use `--report-interval 5000` or higher (3× smaller .nc), (b) deploy a watchdog that deletes equilibration.nc and minimized.pdb once production starts (saves ~800 MB), (c) run a local backup loop that pulls production.nc to WD_BLACK every 30 min — as insurance, not as the primary fix.


### Published αVβ3 Ectodomain PDBs Are All Bent (Route D Eliminated)
- Command context: obj-041, scoring all 5 published full-ectodomain αVβ3 crystal structures (1JV2, 1L5G, 4G1E, 4G1M, 4MMX) on the v7 domain CVs.
- Symptom: All 5 cluster tightly in CV0 ≈ 51-52 Å (BC band), CV1 ≈ 26-28 Å, CV2 ≈ 36 Å. Even cilengitide-bound 1L5G (open headpiece internally) crystallizes bent overall.
- Root cause: full-ectodomain αVβ3 only crystallizes bent because the legs collapse onto the head in solution; cryo-EM of EO state exists but at lower resolution and not all deposited as PDBs.
- Action: do NOT propose "download an EO PDB" as a route to EO endpoints for αVβ3. Enhanced-sampling MD (string method, metadynamics, REMD, Switching Gō-Martini) is the only path. Documented as `docs/eo_coverage_strategy.md` with the recommendation revised from "D + A" to "A primary + E parallel backup" after this finding.
- Practical implication: empirically test "free shortcut" routes before betting on them in strategy docs. obj-041 took 1 day of $0 compute and saved a multi-week dead-end.

### V1↔V2 Per-Residue RMSF Pearson r=0.998 Validates Pipeline Reproducibility
- Command context: obj-042, head-aligned per-residue RMSF over V1 (379 frames) and V2 (1266 frames) fitted-trajectory CAs.
- Symptom: Per-video means agree to 4 sig figs (V1=19.29, V2=19.30 Å), per-residue Pearson r = 0.998. 500-bootstrap CI width median 2.41 Å.
- Action: this cross-video reproducibility level is much higher than typical pipeline noise — it cannot be a coding artefact since V1 and V2 are processed independently from raw HS-AFM streams. Use this as the headline reproducibility number in the paper. The mechanical-sensitivity ranking (obj-028 composite v2) is now quantitatively defended against finite-sample noise.
- Practical implication: when two independent runs of the same pipeline give nearly-identical per-residue scalars, that's strong evidence for *both* the pipeline being deterministic and the underlying biology being reproducible across HS-AFM datasets.

### BLOSUM62 αV-αIIb Pairwise Alignment: Per-Position Lookup Mandatory
- Command context: route-A kickoff doc, αV(1JV2 chain A) ↔ αIIb(3FCS chain A) global alignment.
- Symptom: 38.6 % identity, 949-column alignment. Naive expectation: a single global Δ residue offset would let us port αIIbβ3 CV definitions to αVβ3.
- Reality: the offset varies through the structure: −2 to +0 at the N-terminus (W1-W2 propeller blades), +12 to +16 in the mid-propeller (W4-W7), +13 across genu + calf-1, +6 in calf-2, +4 in membrane-proximal. αV D218 (RGD-Arg pocket) maps to αIIb F231 — Δ +13 *and* non-conservative substitution.
- Action: never use a single global Δ for cross-integrin homology porting. Always run a real pairwise alignment, save the per-position lookup table (results/route_a/av_aiib_alignment.json), and feed it through a remap script (pipelines/route_a/scripts/remap_cvs.py).
- Practical implication: this is generally true for cross-paralog porting in protein families. Apparently-similar sequence lengths can hide structurally-meaningful insertion/deletion patterns that break naive offset assumptions.

### Bayesian Rate Bound Salvages "We Don't Have Data There" Reviewer Concerns
- Command context: obj-043 audit-2026-05-05 §13. Reviewer A asks "what happens at CV0 ≥ 90 Å?" but the v7 fitted set has 42/1645 frames at CV0 ≥ 85 Å and zero at CV2 ≥ 50 Å (the EO definition).
- Naive response: "we don't sample EO; the answer requires enhanced-sampling MD" — concedes the reviewer's point fully.
- Better response: report a **Jeffreys (Beta(0.5,0.5)) two-sided 95 % upper bound** on the unobserved-rate parameter. With n=1645 and k=42, P(CV0 ≥ 85 Å) ≤ 3.40 % at 95 % confidence → ΔG_EO ≥ −kT·log(0.034) = 2.02 kcal/mol relative to the Intermediate minimum. The reviewer's geometric question is still unanswered, but the *thermodynamic* claim is bounded.
- Action: when a reviewer flags missing data in a region, look for a Bayesian rate bound or a one-sided concentration bound before retreating to "we need more compute". Often the existing data already constrains the unobserved region enough to make a publishable claim.
- Practical implication: Jeffreys priors give well-calibrated bounds for Bernoulli-type tail-occupancy questions even when the count is zero or near-zero. Use scipy.stats.beta.ppf(1-α/2, 0.5+k, 0.5+n-k) — no asymptotic-normal approximation, valid for any (k, n).

### Audit Document Integrity: Verify Numerical Claims Against On-Disk Data
- Command context: audit-2026-05-05 §13.1 deepening pass v5. Cross-checked every figure, JSON, and metric cited in §1-§12 against on-disk reality.
- Symptom: 23/23 cited artifacts existed but **1 numerical drift**: docs/route_a_kickoff.md and audit §12.2 claimed "Alignment score: 167" while the actual av_aiib_alignment.json has alignment_score=1507. The 167 figure was implausible (typical BLOSUM62 score for 949-column 38.6 %-identity alignment is 1000-1500), but had been silently propagated to two documents.
- Root cause: a metric was likely derived incorrectly during initial doc drafting (wrong sum or mis-scaled value) and never re-verified against the saved JSON.
- Action: any audit document that makes specific numerical claims should cycle through the source artifacts at least once and re-verify each metric. Add a per-pass integrity check as a routine step of audit deepening; it took ~10 minutes here and caught one real drift.
- Practical implication: audit/decision documents are most useful when they're trusted; one wrong number in a high-stakes doc poisons the rest. Build a metric-extraction lineage so each claim can point to a specific JSON key and value.

### HS-AFM State Populations Are Non-Stationary Across Acquisition Time
- Command context: obj-044 audit §14. Windowed (50-frame block) re-analysis of obj-043's pooled state populations in V1 (379 frames) and V2 (1266 frames) of fitted CV0.
- Symptom: under the stationary-Bernoulli null with Bonferroni-correction across V1+V2 (128 block-state pairs), V1 has 7/28 and V2 has 27/100 pairs significantly outside the pooled mean. Maximum |z| = 6.70 (V1) and 7.50 (V2).
- Root cause: the surface-bound HS-AFM ensemble evolves on the experiment timescale — at 1 fps acquisition this is seconds-to-minutes. The pooled population fractions (BC 25 % / Inter 46 % / EC 24 % / EO* 5 %) are an ensemble *average* across this drift, not a steady-state distribution.
- Action: when reporting fractional state populations from HS-AFM (or any time-series molecule sampling), always (a) check stationarity by windowed analysis, (b) include both the pooled average and a per-window breakdown in the supplementary, (c) note explicitly in the discussion that the ensemble is non-stationary if so. Pooling alone hides drift.
- Practical implication: this also explains the V1↔V2 KS p = 5.5e-4 from obj-043 — V1 and V2 sample different parts of an evolving ensemble, not different ensembles. The cross-video reproducibility at the *ranking* level (RMSF Pearson r = 0.998) is preserved; only the marginal *distribution* differs in the tails.

### Sim + Overlay Must Project the Same Coord Variant That Fed the Sim Render
- Command context: obj-071. The v7 pipeline at `results/afm_pipeline/v7_smoothed_final/video<N>/` stores at least 8 named coord variants per video (smooth, stable, smooth_pretemporal, smooth_tail, median_smooth, median_reanchored, temporal, raw). They are all 1266×N_atoms×3 same-shape arrays and differ only in which post-processing has been applied.
- Symptom: when I rendered a v16 sim AFM overlay using `render_projection_overlay.py` (defaults to `fitted_coords_smooth.npy`) and compared against the sim AFM image (which was rendered by `simulate_afm_video.py` reading `fitted_coords_stable.npy`), the projected atoms sat at a visible 30° tilt off the AFM blob. The bug was *invisible to the script*: r=0.91 in the overlay caption because the script computes r vs. the AFM image's intensity correlation, which is rotation-tolerant enough at 35-px scale.
- Root cause: `simulate_afm_video.py` consumes `--coord-file coords_stable.npy` (default in `run_pipeline.sh`), but `render_projection_overlay.py` defaults to `fitted_coords_smooth.npy`. The two differ by per-frame cosmetic rotations applied by `stabilize_orientation.py` (PCA-flatten + side-lock + sign-align + stepwise yaw + head xy anchor smoothing). At frame 100 the head→tail axis differs by 66° between smooth and stable.
- Action: when rendering an overlay onto a sim AFM, always trace which coord variant the sim was built from and project the **same** variant. If the overlay script doesn't expose a `--coord-file` flag, build a shadow fitted-dir where `fitted_coords_smooth.npy` is a symlink to the variant you actually want; the script will read through transparently.
- Practical implication: any future "overlay on sim" rendering must include a one-line provenance print at top showing the absolute path of the coord file actually loaded. The render script's silence about this is the trap — add a `print(f"Loaded coords from: {smooth_path.resolve()}")` after the load to make the choice loud. Future agents grep that line first when alignment looks off.

### Always Verify a PDB File Actually Matches Its Claimed Identity Before Use
- Command context: obj-070. Planning priority #6 directed me to extend `multi_integrin_first_principles.py` to "6WOV α5β1 bent cryo-EM"; `data/multi_integrin/raw_pdbs/6WOV.cif` was already in the repo.
- Symptom: `grep _struct.title` on 6WOV.cif returned "Cryo-EM structure of recombinant mouse Ryanodine Receptor type 2 wild type in complex with FKBP12.6" — a 565 kDa Ca²⁺ channel tetramer, not an integrin. The PDB ID had been mis-cited in `docs/integrin_heterodimer_plan.md`.
- Root cause: someone (possibly during the initial integrin-heterodimer survey) recorded a wrong PDB ID for α5β1 full ectodomain cryo-EM. Reasonable guess: confusion between PDB IDs starting with "6W" in 2020 deposits. The correct ID is **7NXD** (Schumacher 2021, "human integrin alpha5beta1 in the half-bent conformation").
- Action: before running any per-PDB pipeline on a downloaded structure, do a one-line sanity check on `_struct.title` (CIF) or `HEADER`/`TITLE` (PDB) and confirm the PDB ID actually matches the protein you think it is. Especially for: (a) cryo-EM structures with C2/C4 symmetry (an "asymmetric unit shows 4 copies" complaint is often a symptom of a wrong PDB), (b) historical references in planning docs that haven't been re-verified, (c) bulk-downloaded sets pulled by someone else.
- Practical implication: a single 30-second grep saves hours of wrong-protein analysis. Add `_struct.title` to the validation checklist for any new PDB added to a registry.

### When Re-rendering Pipeline Output, Re-render ALL Cases
- Command context: obj-069. The May-3 overlay redesign (commit c1e7b0f, locked-square AFM + stacked CV-trajectory line plot with state bands + playhead) was applied to the script but only re-rendered for v2. v1 was left at the old style for 10 days until the PI asked for the pair over Slack.
- Symptom: cousin outputs drift out of sync. "Final" outputs in `results/afm_pipeline/v7_smoothed_final/video{1,2}/` were stylistically inconsistent — paper figure can't pair them.
- Root cause: when iterating on a multi-input pipeline, the engineer often re-renders against the case they're staring at (here, v2 because it's the harder, longer movie) and forgets the rest.
- Action: when a rendering / plotting / analysis script is materially modified, immediately regenerate output for ALL canonical input cases (the v7_smoothed_final/video{1,2}/ pair, both BC and EC reference frames, both V1 and V2 trajectories, etc.). Treat "the easy / quick case" as part of the same commit. If the cluster cost makes that prohibitive, gate it with a comment in the script's docstring listing which cases still need re-rendering.
- Practical implication: a "final" subfolder is only as final as the most-recent file in it. Add modification-time consistency to the validation checklist for canonical output directories.

### Frame Compute Asks to the PI as Staged Early-Stop Gates, Not a Flat Total
- Command context: 2026-06-27 Slack exchange. The route-A docs requested a flat "4-week PACE A100-80GB block (~$800)" for the αIIbβ3 string-method port. The PI pushed back: *"Why do we need $800 of A100 wall time? What are you running? I thought we established that no model is actually giving us the conformational diversity we want. So is it something else?"*
- Symptom: a single-number ask ($800) with no decision-gate structure reads as an all-or-nothing gamble, and the PI (managing ~20 projects) had lost the thread that this is *enhanced-sampling MD (the string method)*, not another structure predictor. The flat ask invited a flat "no / why?".
- Root cause: the planning docs (a) buried the "this is MD, not a model" distinction, and (b) presented the budget as one lump instead of the gated structure that already existed in `docs/route_a_risk_register.md` (4 gates, each a cheap early-stop). The compute *plan* was staged; the *ask* was not.
- Action: when requesting non-trivial compute/$$ from the PI, always (a) lead with what the run *is and is not* in plain language (here: enhanced-sampling MD that yields a ΔG curve — explicitly NOT the predictors we already showed collapse to bent), and (b) present a staged budget with a cheap first decision-gate, recommending approval of *only* the first stage. Convert an $800 gamble into a $200 test. State the honest joint success odds (~40 % here) so the gates' value is obvious.
- Practical implication: this is the codified form of the PI's pushback (autonomy rule: corrections belong in docs/lessons, not just memory). `docs/route_a_kickoff.md` §4–§5 rewritten 2026-06-27 into a 3-stage gate table (Stage 1 ~$200 / ~1 wk topology+CV-check, Stages 2–3 ~$600 gated on Stage 1). Reuse this template for any future RunPod/PACE allocation request.

### Do Not Inherit GPU VRAM Specs Across Job Classes — MD ≠ ML Inference
- Command context: 2026-06-27 Slack. PI asked "do you really need an A100-80GB?" for the route-A string-method run. Answer: no. The 80 GB was copy-inherited from this repo's structure-prediction jobs (Protenix/Boltz/AF2), where it is genuinely needed for attention over the full ~1654-residue sequence.
- Symptom: an over-specced (and over-priced) GPU request that drew justified PI pushback. A100-80GB is the 6th-cheapest of PACE's 9 card tiers; the job actually fits on the cheapest (V100-16GB).
- Root cause: molecular dynamics and deep-learning inference are different memory regimes. MD (OpenMM/GROMACS/AMBER) holds positions + forces + neighbor lists — roughly a few GB even for ~1 M-atom membrane-embedded explicit-solvent systems (well under 16 GB). ML inference holds large weight/activation tensors that scale with sequence length and need 40–80 GB. Defaulting the MD job to the prediction job's card conflated the two.
- Action: pick the GPU from the *job class*, not from whatever the last job used. For MD: VRAM need ≈ a few GB; size to the system, and remember MD cost is driven by *throughput* (ns/day per $/hr), NOT by VRAM. The cheapest card that fits is usually right; benchmark ns/day on a cheap card first, then lock the production tier from real numbers. The string method parallelizes 1 image/GPU, so prefer "several modest cards" over "one large card."
- Practical implication: per CLAUDE.md, PACE cards are ordered cheapest-first (V100-16GB … A100-80GB … H200). For any MD allocation, start at the cheap end and only move up if a ns/day benchmark justifies it. `docs/route_a_kickoff.md` §4 now specifies V100-16GB (Stage 1) and A100-40GB/L40S (production), not A100-80GB.

### `df` on a Shared Network Filesystem Reports the Whole Cluster, Not Your Allocation
- Command context: 2026-06-28. Inspecting the RunPod pod, `df -h /workspace` returned `mfs#euro.runpod.net:9421  2.3P  1.7P  591T  75% /workspace`. I reported this to the PI as "we have a 2.3 PB persistent volume" and wrote it into `docs/route_a_kickoff.md` + `tasks/planning.md`. PI: "there is no such thing." Correct — it was a hallucinated capacity.
- Symptom: a wildly wrong storage figure (2.3 *petabytes*) presented as our allocation, propagated into committed planning docs and three Slack messages before being caught.
- Root cause: `/workspace` is a *shared* MooseFS network mount (`mfs#...:9421`, fuse). `df` on a distributed/network filesystem reports the **entire backing cluster's** total/used/free, which is shared across all tenants — it says nothing about the calling user's quota or allocation. The real per-pod truth: container/root overlay was 40 GB; actual data under /workspace was only 6.4 GB (`du -sh`).
- Action: never infer "our available space" from `df` on a network/shared mount (MooseFS/NFS/Lustre/GlusterFS/fuse). To size storage: (a) `du -sh` the actual directory for current usage, (b) check the provider's volume config (RunPod volume size set at deploy, or `RUNPOD_*` env / dashboard) for the real allocation, (c) treat any `df` total in the PB/hundreds-of-TB range on a `mfs#`/`nfs`/`:`-style mount as a *cluster* number, not yours. Sanity-check absurd magnitudes against the plausible (a hobby/research pod does not have a 2.3 PB private volume).
- Practical implication: faithful reporting > confident reporting. A storage/region/quota figure that feeds a PI decision must come from the provider's allocation config, not a filesystem df. When a number looks too big to be true, it is — verify before stating it as fact.
