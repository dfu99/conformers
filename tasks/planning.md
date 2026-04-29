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
1. **Matsumoto 2008 switch-residue overlay** — quick ~1h CPU job, mechanical-sensitivity validation against the literature. Direct overlay of our top-10 hinges against their ENM switch-residue list.
2. **Rotation-corrected RMSF** — align each frame internally to the head before RMSF. Small code change in `residue_rmsf.py`. ~1h.
3. **Multi-integrin pipeline** — αIIbβ3 (3FCS) + α5β1 (3VI4) + αVβ6 (4UM9) + αVβ8 (6OM2) etc. First-principles bent/extended distribution prediction. Plan in `docs/integrin_heterodimer_plan.md`. αIIbβ3 carries dual value: also yields EO templates (string-method structures from Ferg-Lab).
4. **AF2-Multimer ablation** — reviewer B push. ~50 GPU-hours; can run on A40 once free.
5. **Third independent HS-AFM dataset** — transferability test; falsifying if corr < 0.7.
6. **CNN retraining with corrected tip size (1-2 nm)** — Correlation matching (std=8.5 Å) remains the better inference method. CNN needs real-AFM fine-tuning.
7. **A5B1 Protenix co-fold on RunPod A100** — Blocked by PACE billing limits.

## Confirmed Negative Results
- **EC vs EO SMD steering (obj-025, 2026-04-29)** — Two iterations (k=250, k=1000) on RunPod A40 with corrected (α-head, β-head) pair list. Even k=1000 only opens CV2 by 0.9 Å in 620 ps (~0.07 Å/ps). Reaching the 60 Å EO target would need ~320 ns wall-clock, an order of magnitude beyond the 3-ns budget. Both runs disk-quota-bounded at ~620 ps regardless of force constant. **EO coverage now flagged as a hard pipeline limitation in `intuition.md` falsifiability section.** To recover EO state space, alternative routes are (a) αIIbβ3 string-method structures, (b) metadynamics with CV2 as collective variable, (c) replica exchange. Any of those is at least 2× the current MD budget.

## Recently Completed
- [x] Fix `cv_distance_headopen` preset (2026-04-29, obj-024) — first run (v3) showed CV2 flat at 34 Å despite the preset name. Cause: preset used default AVB3_HINGE_DISTANCES (head-tail pairs only). Added AVB3_HEADOPEN_DISTANCES with explicit (α-head, β-head) pair, target 6 nm. Routed via preset name in apply_steering_preset. v4 launched with disk-pressure mitigations: --report-interval 5000, watchdog cleanup of equilibration files, 30-min WD_BLACK backup loop.
- [x] Sim HS-AFM realism iteration v3-v11 (2026-04-23 to 2026-04-25, obj-022) — 11-version refinement loop on the forward-renderer with PI feedback at each step. Final settings: copper colormap, zoom-out 35→49 px, substrate noise std 0.022 baseline 0.15, tail-direction surface slant ±4.5%, row-jitter σ=0.008, partial-width flash streaks (2% freq), localized horizontal distortions on molecule rows, PCA-aligned flat-on-surface orientation, step-wise yaw via hysteresis (50° threshold, 20-frame dwell, 30° cap), Gaussian σ=16 head smoothing, post-dilation anisotropic blur σx=1.2/σy=0.7, soft-clip 0.92. Tip R=3 nm.
- [x] RunPod yield + cleanup (2026-04-23 to 2026-04-25, obj-023) — switched from busy-poll to event-driven `mc runpod subscribe` after new infra rolled out. Yielded to halulujah. Cleaned up both `/workspace/conformers/` (12 MB) and `/root/projects/conformers/` (23 GB) on the A40 after PI audit caught the second one. WD_BLACK mirror discipline maintained (15 GB results + 561 MB figures + 5.8 GB MD trajectories).
- [x] AFK overnight: sim-HS-AFM + mechanical sensitivity analysis (2026-04-23) — forward-rendering the v7 fitted trajectory gives sim-vs-real corr 0.82 (V1), 0.72 (V2), validating pipeline self-consistency. Three flexibility metrics (RMSF, cross-conformer CA std, CA-CA-CA angle σ) combined into composite figure. Top hinge match: β-knee B:353 (σ=25.4°) pre-registered at ~352. Top triple-agreement hotspots: C-terminal coils B:689, A:864, A:958, A:861, B:652. Scripts: simulate_afm_video.py, residue_rmsf.py, cross_conformer_rmsd.py, hinge_angles.py. Report: results/mechanical_sensitivity_report_v1.md.
- [x] Rolling-median coordinate smoothing (2026-04-22/23) — 77% jitter reduction (V1 45→10Å, V2 38→9Å) vs per-frame fit. Plus per-frame head re-anchor to fix V2 drift. Saved to results/afm_pipeline/v7_smoothed_final/{v1,v2}/.
- [x] v7 pipeline with bent steering library (2026-04-20) — ran `cv_distance_bent` steering on RunPod A4500 (3ns, 5.47 ns/day), extracted 306 bent protein frames, merged with existing 309 extend frames → 615-frame library spanning CV0 [47.3, 85.0]Å. v7 fits: video1 BC 12.4%→43.5%, video2 BC 12.5%→18.6%. Library is discriminatingly sensitive to the data. Also: AFK tasks completed — updated intuition.md, generated library coverage, 1JV2 comparison, tip calibration, steering manifold figures; wrote docs/pipeline_rationale.md and results/paper_draft_v1.md.
- [x] Rebuild pseudo-AFM library + rerun overlay pipeline (v6, 2026-04-16) — Fixed critical batch_size=16 sampling bug (only first 31/309 conformers sampled). Rebuilt with batch_size=1 covering full CV range [52.9, 85.0] Å. SO(3) fitting on RunPod GPU (10 min vs 4h CPU). Video1: 379 frames, corr=0.966, 67 extended. Video2: 1266 frames, corr=0.939, 179 extended.
- [x] Fix head-alignment drift + diagnose conformer coverage gap (obj-014/015, 2026-04-14) — `resolve_flips_head_anchored()` eliminates all position drift (was up to 10.4nm). Diagnosed conformer library gap: 66 conformers span CV0 53-79A, no extended states. Launched 3ns steering extension on RunPod A4500.
- [x] Head-tracked fitting + tip comparison (2026-04-12) — `fit_with_head_tracking.py`: head-based centering (not COM), skip first 30 frames, tip comparison (2-3nm optimal). Correlation: video1 0.69→0.97, video2 0.80→0.94.
- [x] Exhaustive SO(3) rigid-body fitting + rendering (2026-04-12) — `fit_rigid_body.py` (2048 rotations/frame), temporal flip resolution (38%→8% flips), updated `render_atomistic_overlay.py`. Both videos fitted and overlays rendered. Video2 corr=0.80, Video1 corr=0.69.
- [x] Corrected-tip pseudo-AFM correlation matching (2026-04-11) — Regenerated 13,200 pseudo-AFM images with tip 1-2 nm (was 6-12 nm). Correlation matching on both GIFs shows improved state detection: video1 extended 11%→47%, video2 extended 29%→60%.
- [x] Tip size investigation (2026-04-11) — Identified that pseudo-AFM tip radius was 6-12 nm (3-6x too large vs experimental 1-2 nm). Generated comparison grids showing structural detail at each tip size.
- [x] Initial HS-AFM pipeline application (2026-04-10) — Applied AFMFold pipeline to two Linz HS-AFM GIFs (409 + 1296 frames). CNN predictions near-constant due to domain gap, but correlation matching found real conformational dynamics.
- [x] AFMFold CNN pipeline built (obj-010, 2026-04-08) — End-to-end pipeline: domain steering → simulated AFM images → CNN training → GIF inference. CV distance extend is the effective steering method (+88° leg opening, +113Å extension). Pipeline scripts: train_afmfold_cnn.py, predict_from_afm_gif.py.
- [x] Domain steering experiments on RunPod A5000 (obj-010, 2026-04-08) — 3/4 presets completed. CV distance extend dramatically outperformed angle torques. Gentle/moderate had negligible effect in 1ns.
- [x] AF2 pLDDT vs displacement analysis (obj-009, 2026-03-23) — AF2 can't score arbitrary PDBs but pLDDT identifies flexible leg/tail domains (82-84). TM-score is the practical conformer filter.
- [x] AF2 reduced_dbs conformer validation (obj-007, 2026-03-23) — AF2 also locked to bent conformation. 25 predictions have pairwise RMSD 0.1-2.9Å.
- [x] MSA-subsampled Protenix conformer validation (obj-006, 2026-03-20) — TM-score validates frame realism but Protenix is MSA-depth-invariant for AVB3.
- [x] Build AVB3 conformer + pseudo-AFM image pipeline (obj-005, 2026-03-18)
- [x] Add dedicated A5B1 staged pipeline runner + sbatch entrypoint. (obj-001, 2026-03-06)
- [x] Add deterministic stage-merging utility for final tagged complex export. (obj-002, 2026-03-06)
- [x] Build PACE minimal remote control script with smoke test validation. (obj-003, 2026-03-09)
- [x] Add Claude skill for deterministic PACE job operations. (obj-004, 2026-03-10)

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

## Updated: A5B1 Auto-Search Until Pass (2026-03-12)
- Added merge-quality validator:
  - `pipelines/protenix-a5b1/scripts/check_tagged_structure_quality.py`
  - Checks: chain count (expect 4), alpha-tail to AviTag-branch distance, beta-tail to SpyCatcher-branch distance.
- Added automatic seed/anchor sweep runner:
  - `pipelines/protenix-a5b1/scripts/run_complete_tagged_until_pass.sh`
  - Iterates over `SEED_LIST_CSV` and `STAGE2_ANCHOR_LIST_CSV`, runs complete pipeline, and hard-fails if no attempt passes thresholds.
- Updated submit entrypoint:
  - `pipelines/protenix-a5b1/scripts/submit_complete_tagged_pipeline_slurm.sh`
  - Supports `AUTO_SEARCH_UNTIL_PASS=1` to run the new sweep workflow.

## Updated: A5B1 Anchor Sweep + Seed Expansion (2026-03-13)
- Refined anchor sweep defaults in `run_complete_tagged_until_pass.sh`:
  - `STAGE1_ANCHOR_LIST_CSV=1..16` (valid stage-1 chain-C residue range)
  - `STAGE2_ANCHOR_LIST_CSV=1` (targeted default while stage-2 anchor sensitivity remains low)
- Confirmed cached-seed (`101..909`) strict sweep still has no pass at `<=25 A` with best near-pass around `26.44/26.12 A`.
- Began fresh-seed A100 expansion via per-seed submits (`SEED_LIST_CSV=<single_seed>`) to avoid `sbatch --export` CSV pitfalls and continue autonomous pass search.
- Completed first fresh-seed batch (`1001,1111,1222,1333,1444,1555`; jobs `4901896`..`4901901`): all failed strict pass gating; best structure remained the same near-pass from seed `404`.
- Tested targeted stage2-anchor expansion for `seed 404` (`stage1=1`, stage2 incremental sweep): observed invariant geometry across tested anchors, so broad stage2 brute force is currently low-yield.
- Ran post-merge rescoring against alternate heterodimer base frames (`sample_0..4`) and identified that base `sample_2` unlocks a strict pass with existing stage predictions.
- Produced strict-passing final artifact using base `sample_2` + seed `404` best stage samples and promoted it to canonical final outputs:
  - `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.cif`
  - `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.pdb`
  - `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.quality.json`

## Updated: Stricter Chemistry + Exclusion Gating (2026-03-13)
- Extended `pipelines/protenix-a5b1/scripts/check_tagged_structure_quality.py` beyond CA-tail checks:
  - explicit attachment distances (`A966:NZ -> chain D`, `B735:NZ -> C10:[OD1,OD2,CG]`)
  - optional receptor exclusion zones (currently `A780-820` for D-branch rejection)
  - richer JSON diagnostics with nearest atoms and failing constraints.
- Extended `pipelines/protenix-a5b1/scripts/run_complete_tagged_until_pass.sh`:
  - base heterodimer frame sweep via `BASE_SAMPLE_LIST_CSV` (defaults `0..4`)
  - strict defaults (`20 A` tail gates, `6 A` attachment gates, alpha-body exclusion min-distance)
  - direct `HETERODIMER_CIF` override per attempt so base-frame search is in-loop.
- Synced stricter scripts to PACE and started:
  1. cached strict scan for `seed=404`, `base=0..4`, `stage1=1..16` (completed, no pass).
  2. broader cached strict background scan across seeds `101..1555` (in progress, log:
     `logs/protenix-a5b1/a5b1_strict_cached_sweep_20260313_160120.log`).
  3. strict A100 rerun submits for new seeds (queued):
     - `4925696` (`seed=1666`)
     - `4925697` (`seed=1777`)
     - `4925698` (`seed=1888`)

## Updated: Budget-Constrained Search Strategy (2026-03-14)
- Added bounded search controls to `pipelines/protenix-a5b1/scripts/run_complete_tagged_until_pass.sh`:
  - `MAX_ATTEMPTS` to hard-cap one sweep chunk
  - `BASE_CIF_GLOB` to test arbitrary heterodimer base CIF sets (not just Protenix sample `0..4`)
- Exposed `TOP_A` / `TOP_B` in `pipelines/afcluster/scripts/submit_afcluster_boltz_slurm.sh` so AFCluster search breadth can be explicitly kept small.
- Strategy change:
  1. Stop long unconstrained A100 Protenix sweeps as the default.
  2. Use cheap heterodimer-frame exploration first (`BoltzGen`, `AFCluster -> BoltzGen`) on RTX 6000 with short walltimes.
  3. Reuse cached tag-attachment stage predictions and rescore merged outputs under the stricter chemistry/exclusion gate.
- Launched budget-conscious sequential RTX 6000 jobs for A5B1 heterodimer exploration:
  - direct BoltzGen single-spec run: job `4951450`
  - AFCluster -> BoltzGen tiny clustered run (`top_a=2`, `top_b=1`) dependent on `4951450`: job `4951452`
- Result: both jobs failed quickly at BoltzGen config validation because the current YAML generators still emit legacy Boltz schema rather than the installed BoltzGen `entities` schema.
- Immediate next implementation task: update the Boltz/AFCluster YAML emitters before spending more cluster time on this branch.

## Updated: BoltzGen Schema Fix + AFCluster Outlier Ordering (2026-03-15)
- Patched BoltzGen YAML emitters to match the installed `entities` schema:
  - `pipelines/boltz/scripts/build_boltz_predict_sweep.py`
  - `pipelines/afcluster/scripts/make_boltz_jobs_from_clusters.py`
- Template-guided jobs now emit one `file`-context spec instead of legacy per-chain `templates` blocks, because the installed parser does not accept the old nested keys.
- Validated the emitted direct and AFCluster specs on PACE with `boltzgen check`, using real A5B1 sequence/MSA/template inputs. The schema-validation failure is resolved.
- Patched AFCluster cluster ordering to avoid spending budget on `cluster_outliers.a3m` before numbered clusters:
  - `pipelines/afcluster/scripts/cluster_chain_msa.py`
  - `pipelines/afcluster/scripts/make_boltz_jobs_from_clusters.py`
- Relaunched budget sequential RTX 6000 jobs after the fix:
  - direct BoltzGen template-context run: job `4957942`
  - AFCluster -> BoltzGen tiny run (`top_a=2`, `top_b=1`) dependent on `4957942`: job `4957943`

## Updated: BoltzGen Cache Redirect To Scratch (2026-03-16)
- March 15 relaunches (`4957942`, `4957943`) proved the schema fix worked, but both jobs then failed during Hugging Face checkpoint staging with `Disk quota exceeded`.
- Root cause: PACE home quota was fully exhausted (`20480M / 20480M`), and BoltzGen was defaulting to `~/.cache/huggingface`.
- Patched launchers to default Hugging Face caches to scratch:
  - `pipelines/boltz/scripts/run_boltzgen_attempt.sh`
  - `pipelines/afcluster/scripts/run_afcluster_pipeline.sh`
- Relaunched the same tiny sequential RTX 6000 jobs with scratch-backed HF cache:
  - direct BoltzGen cache-fixed run: job `4963854`
  - AFCluster -> BoltzGen cache-fixed run (`top_a=2`, `top_b=1`) dependent on `4963854`: job `4963855`
