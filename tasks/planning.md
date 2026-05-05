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
1. **αIIbβ3 steering MD on RunPod** — First-principles prediction (obj-029) confirms αIIbβ3 (3FCS) is the right second integrin to validate the AFMFold pipeline. 14% more head-leg buried SASA and 13% smaller head-tail than αVβ3, which agrees with Li 2024 / Chen 2023 published activation hierarchy. Run `cv_distance_bent` + `cv_distance_extend` steering on 3FCS once a RunPod window opens. **Pipeline path now locked**: feed 3FCS conformers + library.json directly into `pipelines/sim-afm-video/run_pipeline.sh` (override head residues + side chain). For overlay against any future αIIbβ3 HS-AFM data, use `pipelines/afm-overlay/run_pipeline.sh`. Domain definitions need re-mapping for αIIbβ3 (use `pipelines/avb3-conformers/scripts/map_aiib3_to_avb3_domains.py` as starting point).
2. **AF2-Multimer ablation** — reviewer B push. ~50 GPU-hours; can run on A40 once free.
3. **Third independent HS-AFM dataset** — transferability test; falsifying if corr < 0.7.
4. **Multi-integrin first-principles extension** — extend to 6WOV (full α5β1 cryo-EM bent — CIF-only, needs PDB→CIF parser switch). Skipped this round because mdtraj's CIF support flagged the asymmetric unit aggregation as 4 copies. Could also retry once a clean ectodomain extraction script is written.
5. **CNN retraining with corrected tip size (1-2 nm)** — Correlation matching (std=8.5 Å) remains the better inference method. CNN needs real-AFM fine-tuning.
6. **A5B1 Protenix co-fold on RunPod A100** — Blocked by PACE billing limits.

## Confirmed Negative Results
- **EC vs EO SMD steering (obj-025, 2026-04-29)** — Two iterations (k=250, k=1000) on RunPod A40 with corrected (α-head, β-head) pair list. Even k=1000 only opens CV2 by 0.9 Å in 620 ps (~0.07 Å/ps). Reaching the 60 Å EO target would need ~320 ns wall-clock, an order of magnitude beyond the 3-ns budget. Both runs disk-quota-bounded at ~620 ps regardless of force constant. **EO coverage now flagged as a hard pipeline limitation in `intuition.md` falsifiability section.** To recover EO state space, alternative routes are (a) αIIbβ3 string-method structures, (b) metadynamics with CV2 as collective variable, (c) replica exchange. Any of those is at least 2× the current MD budget.

## Recently Completed
- [x] Switching Gō-Martini parallel track kickoff (2026-05-05) — `docs/go_martini_kickoff.md` + concept figure `figures/go_martini_eo_trajectory.png`. Cross-validation track for route A: independent line of evidence using CG MD with switching Gō-contacts. Gō-A = 1JV2 bent αVβ3, Gō-B = literature EO cryo-EM CV target (CV0 ≈ 92, CV2 ≈ 58 — since obj-041 confirmed no published EO PDB). Software stack: Martinize3 + GROMACS 2024 in isolated `venv_gomartini`. 14-day plan, ~$200 GPU on A40. Acceptance: monotonic trajectory through (CV0, CV2) space, endpoint CV0 ≥ 80 + CV2 ≥ 50, MIDAS SASA *increases* on extension (reverses obj-039 negative result). Joint pass ~50 %.
- [x] Route-A kickoff one-pager (2026-05-05) — `docs/route_a_kickoff.md`. Real BLOSUM62 alignment of αV (1JV2 chain A, 927 CAs) vs αIIb (3FCS chain A, 913 CAs): 38.6 % identity, 949-column alignment. Per-position residue mapping saved as `results/route_a/av_aiib_alignment.json`. Key landmark: αV D218 (RGD-Arg pocket) maps to αIIb F231 (Δ +13) — non-conservative substitution requiring custom CV remap. Offset structure: −2 to +0 at N-term, +12 to +16 in propeller blades 3-7, +6 in calf-2, +4 in membrane-proximal. **Single global Δ does NOT work** — per-position lookup table is mandatory. 7 force-field gotchas documented (force-field family, lipid composition, glycosylation, Ca²⁺ binding loops, MIDAS shared, αIIb heavy/light chain split, β3 PSI disulfide shared). Day-1 deliverables: remap_cvs.py + build_av_topology.py (queued). PI sign-off requested for a 4-week PACE A100-80GB block from 2026-05-12. Companion docs: `docs/route_a_risk_register.md` (5 risks, ~40 % joint pass probability after gating).
- [x] Route-A risk register (2026-05-05) — `docs/route_a_risk_register.md`. 5 highest-probability failure modes for αIIbβ3→αVβ3 port: residue mismatch (P_fail 30 %), membrane lipid composition (25 %), string optimization convergence (30 %), αV glycosylation (15 %), force-field choice (10 %). Each with mitigation + gate to next stage. Cumulative ~28 % independent / ~40 % correlated joint pass.
- [x] Library coverage v3 — route D null result (2026-05-05, obj-041) — Audit §10.9 P=3 follow-up. Empirically tested EO-strategy route D ("download published EO PDBs as endpoints") by scoring all 5 published full-ectodomain αVβ3 crystal structures (1JV2, 1L5G, 4G1E, 4G1M, 4MMX). All 5 sit in the BC band: CV0 ≈ 51-52 Å, CV2 ≈ 36 Å (vs EO threshold ≥ 50 Å). Even cilengitide-bound 1L5G crystallizes bent. Route D eliminated; EO sampling now requires enhanced-sampling MD only. Strategy doc revised: "A primary + E parallel backup" replaces "D + A". Script `pipelines/conformer-library/scripts/library_coverage_v3.py`. Figure: `figures/library_coverage_v3.png`.
- [x] F4 Hertzian contact-mechanics control (2026-05-05, obj-040) — Reviewer C #4 blocker closed. R_eff = 0.857 nm, E* = 1.09 GPa (E_protein = 1 GPa, E_tip = 130 GPa, ν = 0.3). At HS-AFM imaging forces 50-200 pN, indentation δ = 0.11-0.28 nm (5-14 % of a 2 nm probe sphere) — all below the 1 nm vertical noise floor (Ando 2013). Hard-sphere pseudo-AFM is a defensible first-order model; Hertzian correction is a quantifiable known systematic. Optional follow-up: `--hertz-force` flag on `simulate_afm_video.py`. Script `pipelines/sim-afm-video/scripts/contact_mechanics_control.py`. Figure: `figures/contact_mechanics_control.png`.
- [x] EO-coverage strategy one-pager (2026-05-05) — `docs/eo_coverage_strategy.md`. Audit-2026-05-05 §10.9 P=1 follow-up. Compares 5 routes: (A) αIIbβ3 string method via Ferg-Lab `principalcurve_integrin_structures` (already cloned on PACE); (B) Funnel metadynamics with CV2 + Plumed; (C) replica exchange; (D) published EO PDBs (1L5G, 2VDR) as endpoints; (E) Switching Gō-Martini CG. **Recommendation: D + A in sequence** — D ships in 1 day for $0 (closes Reviewer A geometric concern); A produces the EC→EO trajectory needed for ΔG (closes Reviewer D + E). PI sign-off requested by 2026-05-08; approve D for tomorrow + a 4-week PACE A100-80GB block (~$800) for A.
- [x] Audit deepening pass evening (2026-05-05) — F1 library coverage v2 (figures/library_coverage_v2.png), F2 calibration controls (calibration_controls_v1.png), F3 ΔG bootstrap (free_energy_profile_v2.png), composite refresh (audit-2026-05-05.png). Reviewer tally moved 5/6/12 → 8/4/10 (3 concerns closed: Reviewer C tip+random-baseline, Reviewer D ΔG bootstrap, Reviewer E RGD PoC negative). EO blocker sharpened by obj-039.
- [x] RGD-binding-pocket accessibility — bent vs extended (2026-05-05, obj-039) — Reviewer-E follow-up. Negative result with diagnostic value: legs retreat 52 Å (60→112 Å midas-leg), but MIDAS-pocket SASA *drops* 35 % (459→298 Å²) and steric occlusion rises 11 %. Interpretation: v7 "extended" frames are EC (extended-closed); EO headpiece-opening never occurred. Real Vina docking would just rediscover the EO-coverage gap from obj-025. RGD-affinity comparison gated on generating EO frames (αIIbβ3 string method / metadynamics / REMD). Script `pipelines/conformer-library/scripts/rgd_pocket_accessibility.py`. Figures: `figures/rgd_docking_v1.png`, `results/afm_pipeline/rgd_docking/`.
- [x] Drop-list cleanup: shadow scripts deprecated (2026-05-05, audit follow-up) — `pipelines/avb3-conformers/scripts/render_projection_overlay.py` (byte-identical) and `stylize_sim_afm.py` (older σx=1.2/σy=0.7 variant) removed. Verified zero live importers (no Python `import`, no shell launcher reference; only README/changelog mentions). Canonical paths: `pipelines/afm-overlay/scripts/render_projection_overlay.py` and `pipelines/sim-afm-video/scripts/stylize_sim_afm.py`. Updated `results/afm_pipeline/v7_smoothed_final/README.md` regeneration block.
- [x] Experimental ΔG(CV0) from V1+V2 unbiased fitted trajectories (2026-05-05, obj-038) — Computed −kT log P(CV0) at 300 K from 1645 unbiased fitted-trajectory frames (V1=379 + V2=1266). Free-energy minimum sits in Intermediate range (CV0 ≈ 70-75 Å); BC and EC barriers ≈ 4-5 kcal/mol. EO region (>85 Å) hatched to make missing library coverage visually unambiguous. Addresses Reviewer D's #1 concern ("Where is the free-energy profile along your CV distance?") without requiring additional MD. New script `pipelines/afm-overlay/scripts/free_energy_profile.py`. Figure: `figures/free_energy_profile_v1.png`.
- [x] Audit 2026-05-05 (2026-05-05) — `tasks/audit-2026-05-05.md` + `figures/audit-2026-05-05.png`. Findings: 35 objectives shipped, but reviewer panel is 5 addressed / 6 partial / 12 open. Single biggest blocker is EO state coverage (SMD cannot open headpiece in 3 ns at any k). Drop list: three-predictor framing, 6WOV α5β1 cryo-EM analysis, A5B1 tagged-assembly auto-search. Queue refilled with 6 follow-ups; top of stack was the free-energy profile (now done).
- [x] Redesign overlay panel: locked-square AFM + CV-trajectory line plot (2026-05-03, obj-037) — `pipelines/afm-overlay/scripts/render_projection_overlay.py` now uses `set_aspect('equal', adjustable='box')` on AFM/overlay panels (fixes horizontal squish) and a stacked 3-row CV trajectory line plot with `axhspan` state bands + red `axvline` playhead + scatter dot in place of the moving bar chart. Output: `results/afm_pipeline/v7_smoothed_final/video2/overlay_v2/pdb_projection_video.gif` (1233×508, 150 frames, 14 MB).
- [x] Make sim AFM stylize blur isotropic σx=σy=1.0 (2026-05-03, obj-036) — earlier σx=1.2/σy=0.7 was tuned for non-orientation-locked frames; after PCA-flatten + side-lock the long axis is consistently along x and anisotropic blur visibly squashed the molecule horizontally. New defaults match the orientation-locked render. `figures/sim_afm_v16_isotropic_blur.png`.
- [x] Promote afm-overlay to PRIMARY pipeline (2026-04-30, obj-035) — built `pipelines/afm-overlay/` with README + 7-stage `run_pipeline.sh` (ingest → pseudo-AFM library → SO(3) fit with head tracking → tail-flip resolve → temporal median → kinetics → overlay render). Reference fits at `results/afm_pipeline/v7_smoothed_final/{video1,video2}/`. Three PRIMARY pipelines now form complete forward + inverse AFMFold workflow.
- [x] Promote sim-afm-video + conformer-library to PRIMARY pipelines (2026-04-30, obj-034) — created `pipelines/sim-afm-video/` and `pipelines/conformer-library/` with READMEs documenting input contract (folder of PDBs + `library.json` schema with conformers/cv_definitions/states/trajectory_order/adjacency), validated defaults, and provenance. Wrote top-level `pipelines.md` indexing all three. Updated memory entry `project_primary_pipelines.md`.
- [x] Add z-axis blur + realistic z-noise to sim AFM (2026-04-30, obj-033) — added `--z-blur-nm` flag to `simulate_afm_video.py` applying Gaussian σ in nm to dilated height map BEFORE per-frame normalize. Validated `--z-blur-nm 0.8 --noise-nm 0.35`. Sim z-axis now reads as blob-like rather than precision height map. `figures/sim_afm_v15_zblur.png`.
- [x] Enlarge sim AFM canvas to 60 px ~59 nm (2026-04-30, obj-032) — V2 stable coords have per-frame xy span up to 26 nm; 35-px canvas (34 nm) was clipping extended frames. Re-rendered video2 stable_big at 60 px. `figures/sim_afm_v14_big_canvas.png`.
- [x] Restore head + side anchoring + stepwise yaw via stabilize_orientation.py (2026-04-30, obj-031) — wrote `pipelines/sim-afm-video/scripts/stabilize_orientation.py` with 5-stage orientation lock: PCA-flatten → side-lock (chain A side stays down) → long-axis sign-align → stepwise yaw with hysteresis (50° threshold, 20-frame dwell, 30° cap) → head xy Gaussian smoothing (σ=8). V2 yaw on video2 committed 3 discrete 30° steps over 1266 frames. `figures/sim_afm_v13_stable.png`.
- [x] Fix sim AFM v11 inset-border crop (2026-04-30, obj-030) — diagnosed v11 generated noise on a larger canvas and pasted the molecule render into a smaller inset, producing a visible rectangular border whenever molecule extended past inset edge. Built `pipelines/sim-afm-video/scripts/stylize_sim_afm.py` doing all 8 stages on a single matched canvas. Border eliminated. `figures/sim_afm_v12_border_fix.png`.
- [x] Multi-integrin first-principles bent prediction (2026-04-29, obj-029) — Downloaded 7 PDBs (1JV2, 3FCS, 3VI4, 4UM9, 5FFG, 6OM2, 6UJC). αIIbβ3 (3FCS) predicted *more bent* than αVβ3 (1JV2) by both head-leg buried SASA (15,298 vs 13,442 Å²) and head-tail distance (37.2 vs 42.6 Å). Pre-registered first-principles prediction agrees with Li 2024 / Chen 2023 activation hierarchy. Headpiece-only structures excluded from bent-fraction analysis but preserved for future CV2 cross-integrin reference. Report: `results/multi_integrin_first_principles_v1.md`.
- [x] Mechanical-sensitivity composite v2 (2026-04-29, obj-028) — `figures/mechanical_sensitivity_composite_v2.png` rebuilt with rotation-corrected RMSF + cross-conformer std + angular σ. v1 top-5 hotspots all in v2 top-10. β-coil cluster (B:649, B:652) and α-coil-1 doublet (A:842, A:843) sharpen to top. Matsumoto residues overlaid by classification. Report: `results/mechanical_sensitivity_report_v2.md`.
- [x] Rotation-corrected RMSF (2026-04-29, obj-027) — added `--head-align` flag to `residue_rmsf.py` doing per-frame Kabsch alignment on 790 headpiece CAs (αV:1-440 + β3:1-350). Reduced V1 mean RMSF 53.5→19.3 Å (64%) and V2 71.2→19.3 Å (73%). Headpiece now near-rigid (7-10 Å) as expected; legs and coils retain real 20-47 Å motion. C-terminal coil ranking preserved → robust pipeline finding. Two-video consistency improves dramatically (was 53/71 Å, now 19.29/19.30 Å).
- [x] Matsumoto 2008 switch-residue overlay (2026-04-29, obj-026) — 5/5 backbone-hinge-relevant Matsumoto residues at ≥80th percentile in our v7 angular-σ. Cys374 (rank 11) + Leu375 (rank 47) direct hits; primary snap β3:Arg633 at 82.7th percentile. The 6 misses are Interaction-B partners or α-constraint Ser305 that stabilize via non-bonded networks rather than backbone hinging — a meaningful methodological distinction. Cross-method concordance with foundational ENM-NMA paper at residue level. Report: `results/matsumoto_overlay_v1.md`.
- [x] EC vs EO SMD steering — confirmed negative (2026-04-29, obj-025) — Two iterations (k=250 v4, k=1000 v5) on RunPod A40. Even at k=1000 only 0.9 Å of CV2 opening in 620 ps (~0.07 Å/ps). Reaching the 60 Å EO target would need ~320 ns. Both runs disk-quota-bounded at ~620 ps regardless of force constant. **EO coverage is now a confirmed pipeline limitation**; alternative routes are αIIbβ3 string-method structures, metadynamics, REMD.
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
