# Paper Draft v1 — Claim-to-Evidence Mapping

## Title (working)
Conformer-library template matching recovers single-molecule αVβ3 integrin
bent↔extended dynamics from high-speed AFM without CNN training

## Primary claim

A 615-frame MD-derived αVβ3 conformer library, projected through a
corrected 1–2 nm pseudo-AFM imaging model, fits two independent HS-AFM
recordings with per-frame correlation > 0.93 and recovers state
occupancies that are sensitive to library completeness — demonstrating a
structure-to-image template-matching framework for single-molecule
conformational dynamics.

## Evidence map

### Claim 1: Head-tracked SO(3) fitting achieves corr > 0.93 on both videos
- *Evidence*: `results/afm_pipeline/video1_v7/fitting_metadata.json`
  (corr=0.965), `results/afm_pipeline/video2_v7/fitting_metadata.json`
  (corr=0.939).
- *Figures*: `results/afm_pipeline/video1_v7/pdb_projection_video_opt.gif`,
  `results/afm_pipeline/video2_v7/pdb_projection_video_opt.gif`.
- *Status*: supported.

### Claim 2: Corrected tip size (1–2 nm) drives the correlation jump
- *Evidence*: pre-correction baseline corr=0.69 (video1) and 0.80 (video2)
  with 6–12 nm tip; after correction jumps to 0.97 and 0.94.
- *Figures*: earlier `tip_size_comparison.png` (to be re-rendered with v7 data);
  `figures/tip_calibration_sphere.png` for tip mechanics sanity check.
- *Status*: supported — but reviewer C will want Hertzian vs hard-sphere
  contact-mechanics validation. Pending.

### Claim 3: Library completeness drives state recovery
- *Evidence*: v6 library (309 frames, zero bent) → video1 BC=12.4%,
  video2 BC=12.5%. v7 library (+306 bent frames) → video1 BC=*43.5%*,
  video2 BC=18.6%. The divergence (video1 shifted dramatically, video2
  shifted modestly) shows the pipeline reads the data, not the library.
- *Figures*: `figures/library_coverage_v6_v7.png` (the discriminating
  figure), `results/afm_pipeline/video1_v7/state_kinetics.png`,
  `results/afm_pipeline/video2_v7/state_kinetics.png`.
- *Status*: supported — this is the key discriminating experiment.

### Claim 4: Steered MD produces conformers covering CV0 47.3 – 85.0 Å
- *Evidence*: extend steering CV0 range [52.9, 85.0] Å; bent steering
  CV0 range [47.3, 63.7] Å. Combined library 615 frames. Bent library
  reaches *below* PDB 1JV2 (CV0=51.4Å) — reviewer A's concern that we
  don't reach the bent crystal is not quite right: we actually over-shoot
  the crystal.
- *Figures*: `figures/steering_cv_trajectory.png`,
  `figures/pdb_1jv2_vs_library.png`.
- *Status*: supported.

### Claim 5: Head-anchored flip resolution eliminates position drift
- *Evidence*: previous-frame anchoring produced up to 10.4 nm cumulative
  drift; head-anchored zero-flip resolution is stable (0/1266 flips
  required in video2 v7). Tail-flip smoothing additionally corrects
  18–28% of frames.
- *Figures*: earlier `head_drift_comparison.png` (obj-014).
- *Status*: supported.

### Claim 6: CNN inference fails due to sim-real domain gap
- *Evidence*: CNN trained on simulated pseudo-AFM gives near-constant CV
  predictions on real Linz data. Documented in obj-010.
- *Status*: supported — but reviewer B will push for deeper analysis of
  *why* the gap exists. Pending.

## Claims that are currently NOT supported

### Pending 1: EC vs EO differentiation
Our library has CV2 (head-head) = 34.3–35.5 Å across *all* 615 conformers.
The headpiece never opens. This means we cannot call any state "EO" —
only BC, Intermediate, and EC (closed). Needs CV2-targeted steering or
αIIbβ3 string-method structures.

### Pending 2: Ablation against AlphaFold-Multimer baseline
No direct comparison yet. Estimated ~50 GPU-hours of AF2 runs needed.
See `docs/pipeline_rationale.md`.

### Pending 3: Random-conformer baseline
Pre-registered falsification: corr < 0.4 for uniform random conformer
assignment. Not yet run — would take a few hours locally.

### Pending 4: Third-dataset generalization
Entire claim is based on two GIFs from one source (Linz). A second
independent HS-AFM dataset would falsify or strengthen transferability.

### Pending 5: Contact-mechanics validation
Current pseudo-AFM uses hard-sphere dilation. Reviewer C will ask for
Hertzian or JKR mechanics comparison. Not yet tested — see
`figures/tip_calibration_sphere.png` for initial geometric sanity only.

## Key figures for the paper

1. *Fig 1*: Overlay video grids showing PDB fits on AFM frames (v7).
2. *Fig 2*: Library coverage v6 vs v7 — the discriminating experiment.
3. *Fig 3*: State kinetics (strip + time-course + pie) for both videos.
4. *Fig 4*: Steering CV manifolds + 1JV2 reference.
5. *Fig 5 (pending)*: AF2 ablation, random baseline, EC/EO sweep.

## Word count
~500 words. Full paper will expand methods + discussion to ~6000.
