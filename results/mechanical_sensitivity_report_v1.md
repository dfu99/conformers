# Mechanical Sensitivity Report v1 — 2026-04-23

## Summary

The 615-frame αVβ3 conformer library (extended + bent steering) plus
v7 fitted HS-AFM trajectories provide three independent, per-residue
flexibility metrics. The triple-metric composite
(`figures/mechanical_sensitivity_composite.png`) identifies C-terminal
coil residues as the most flexible, while bond-angle variance alone
correctly recovers the β-knee (B:353) as the #1 classical hinge.

## Methods

1. **Per-residue RMSF** (`figures/rmsf_v7/{v1,v2}/`) — RMS fluctuation
   of each CA across the v7 fitted trajectory after head-anchored
   alignment. Rotation-dominated because the head centroid is locked
   but orientation varies across frames.
2. **Cross-conformer CA std** (`figures/cross_conformer_v7/`) —
   Kabsch-align all 615 library conformers to conformer 0, compute the
   std of each CA's position across the library. Captures intrinsic
   conformational diversity of the library.
3. **CA-CA-CA angle σ** (`figures/hinges_v7/`) — for each residue i,
   the bond angle ∠(CA_{i-1}, CA_i, CA_{i+1}) across the 615 library
   conformers; its std is the hinge-flexibility metric.
4. **Triple-metric agreement**
   (`figures/mechanical_sensitivity_composite.png`) — product of the
   three metrics' z-scores (rectified), highlighting residues robust
   across all three.
5. **Forward-rendering validation**
   (`results/afm_pipeline/sim_afm/{v1,v2}/`) — v7 fitted trajectories
   run back through the afmfold imaging pipeline. Sim-vs-real cosine
   correlation: V1 0.82, V2 0.72.

## Headline findings

### Finding 1 — β-knee validated by hinge-angle metric
The #1 angular-variance hinge is *B:353* (σ=25.4°), which is literally
the β-head/β-tail boundary — i.e., the β-knee. Pre-registered prediction
matched at the exact residue. α-genu at A:420 ranks 13 (σ=15.6°), within
the expected α-head-thigh/α-calf boundary region.

### Finding 2 — C-terminal coils dominate intrinsic flexibility
Top-20 cross-conformer variable residues are all in:
- Chain A: 956-962 (α-coil, membrane-proximal)
- Chain A: 757-767 (α-calf-2 / α-coil junction)
- Chain B: 643-692 (β-coil, C-terminus)
All with per-residue std 36-40 Å.

These dominate the triple-agreement score. Top-5 most-robust
mechanical hotspots: B:689, A:864, A:958, A:861, B:652 — every one
is a C-terminal coil residue.

### Finding 3 — α-calf flexibility cluster (potentially new)
The α-calf domain shows unexpected flexibility clusters at A:592-595
(σ 19-22°) and A:732-741 (σ 18-21°). These are not classical hinge
residues in the Xiong 2001 / Luo 2007 review. Possible interpretations:
1. Our CV-distance steering drove the α-head-thigh ↔ α-calf distance
   and may have recruited motion at these intra-calf positions as a
   mechanical consequence of the forced reaction coordinate.
2. They represent genuine allosteric coupling within the calf domain
   not captured by classical NMA.
3. Artifact of the fitting — bears cross-checking against Matsumoto
   2008's ENM switch residue list.

### Finding 4 — pipeline self-consistency at 0.72-0.82 corr
Forward-rendering the fitted v7 trajectories reproduces the real Linz
AFM at per-frame cosine correlation 0.82 (V1) and 0.72 (V2). Both
well above random baseline (0.65, 0.43 resp). Gap to ideal (1.0) is
consistent with unmodeled contact mechanics and scan artifacts.

## Relation to prior work

Matsumoto et al. 2008 (ENM NMA on αVβ3) identified key switch residues
at the interdomain interfaces. Their work treated the bent ectodomain
as an elastic network — low-frequency modes reproduce the bent→extended
transition. Their hotspot list was dominated by electrostatic switch
residues at the genu. Direct comparison (not yet done): overlay our
top-10 hinges against their switch list. This should be AFK-NEXT-1.

Chen et al. 2023 (ACS Nano) showed αVβ3 is bistable even without
force — consistent with our v7 finding that video1 shows ~50/50
bent/extended occupancy.

The C-terminal coil dominance in our analysis is *not* a feature
emphasized in classical ENM work, which focuses on interdomain
rigid-body motions. HS-AFM specifically senses vertical height, so
C-terminal coils that flop around in Z would show up prominently. This
may be a unique strength of the imaging-informed flexibility approach.

## Open questions / next steps

1. **Overlay with Matsumoto 2008 switch residues** — explicit agreement
   table between our top hinges and their list. Validates classical
   literature overlap.
2. **Separate α-calf cluster** — is A:592-595 / A:732-741 flexibility
   real or a library-construction artifact? Need a targeted MD run
   that holds CV0/CV1 constant and varies only the α-calf internal
   angles to test.
3. **Head rotation confound** — the RMSF metric is rotation-dominated.
   To get rotation-corrected RMSF we should align each frame internally
   (e.g., to the head domain only) before computing RMSF. Small code
   change.
4. **EC vs EO gap persists** — CV2 constant at 34.3-35.5 Å across library.
   Headpiece opening not captured. Same caveat as v7.
5. **Publish readiness** — composite figure is figure 4-candidate for
   the paper. Claim-to-evidence map updated.

## Evidence artifacts

| Claim | Artifact |
|-------|----------|
| β-knee validated | `figures/hinges_v7/hinge_candidates.md` rank 1 |
| C-terminal coil dominance | `figures/cross_conformer_v7/top20_flexible_residues.md` |
| α-calf cluster | `figures/hinges_v7/hinge_candidates.md` rank 2-12 |
| Triple-agreement hotspots | `results/mech_sensitivity_top5_agreement.md` |
| Self-consistency | `results/afm_pipeline/sim_afm/README.md` |
| Lit context | `docs/literature_mechanical_sensitivity.md` |
| Composite figure | `figures/mechanical_sensitivity_composite.png` |
