# Mechanical Sensitivity Report v2 — 2026-04-29

## What changed from v1

Two pipeline upgrades feed into this revision:

1. **Rotation-corrected RMSF** (`obj-027`) — replaces the
   rotation-dominated RMSF panel. V1 mean RMSF dropped 53.5→19.3 Å
   (64%), V2 71.2→19.3 Å (73%). Headpiece is now near-rigid (7-10 Å)
   as expected; absolute RMSF values are now interpretable instead of
   dominated by global head orientation.
2. **Matsumoto 2008 cross-validation** (`obj-026`) — overlaid 12
   ENM-NMA switch residues on the per-residue σ map. 5/5 backbone-
   hinge-relevant Matsumoto residues sit at ≥80th percentile in our
   angular-σ. The 6 misses are non-bonded interaction partners that
   stabilize the bent state without backbone hinge activity.

The composite figure has been rebuilt as
`figures/mechanical_sensitivity_composite_v2.png` with rotation-
corrected RMSF + cross-conformer std + bond-angle σ, plus Matsumoto
residues overlaid by classification color and top-10 hotspots circled.

## Updated top-10 triple-agreement hotspots

Composite score = product of rectified z-scores across the three
metrics. Higher = robust across all three.

| rank | chain | resSeq | composite | z(RMSF) | z(cross) | z(angle) | location |
|---:|:---:|---:|---:|---:|---:|---:|---|
| 1 | B | 689 | 10.47 | 2.26 | 1.84 | 2.52 | β-coil C-terminus |
| 2 | B | 652 | 10.14 | 1.75 | 1.77 | 3.28 | β-coil C-terminus |
| 3 | A | 842 | 10.03 | 1.69 | 1.62 | 3.65 | α-coil 1 |
| 4 | B | 649 |  9.92 | 1.85 | 1.74 | 3.08 | β-coil C-terminus |
| 5 | A | 843 |  9.85 | 1.62 | 1.54 | 3.95 | α-coil 1 |
| 6 | A | 864 |  9.46 | 2.01 | 1.65 | 2.85 | α-coil 2 |
| 7 | A | 958 |  9.25 | 2.36 | 2.06 | 1.90 | α-coil deep C-term |
| 8 | A | 777 |  9.24 | 1.83 | 1.56 | 3.24 | α-calf-2/α-coil junction |
| 9 | A | 861 |  8.56 | 1.80 | 1.59 | 3.00 | α-coil 2 |
| 10 | A | 760 |  7.89 | 2.35 | 2.03 | 1.65 | α-calf-2/α-coil junction |

The v1 top-5 was {B:689, A:864, A:958, A:861, B:652} — all five remain
in the v2 top-10. Rotation correction promotes the β-coil cluster
(B:649, B:652) and the α-coil-1 doublet (A:842, A:843) into the very
top because their backbone bond-angle σ now stands out without RMSF
saturation hiding the signal.

## Key findings (revised)

### Finding 1 — β-knee + Matsumoto cross-validation
B:353 still ranks #1 by angular σ (25.4°). Adding the Matsumoto
overlay: Cys374 (rank 11, σ=18.0°) and Leu375 (rank 47, σ=12.2°)
sit immediately adjacent — the snap-sandwich neighborhood that
Matsumoto identified as locking Arg633. Cross-method concordance
strong.

### Finding 2 — C-terminal coil dominance survives rotation correction
Rotation-corrected RMSF preserves the C-terminal coil ranking; the
top-10 hotspots are dominated by α-coil (A:842-843, A:864, A:861,
A:958, A:777, A:760) and β-coil (B:649, B:652, B:689) residues. With
the RMSF artifact removed, this is now a clean, interpretable
finding: the C-terminal coils that flop in vertical (Z) are the
same ones that show high cross-conformer std and high backbone
bond-angle σ across the 615-frame library.

### Finding 3 — α-calf cluster persists, narrower
The α-calf clusters at A:592-595 (σ 19-22°) and A:732-741 (σ 18-21°)
are still in the top-20 angular-σ list. Neither overlaps with
Matsumoto's switch residues. The composite score's α-calf
contribution is at A:760-777 (the calf-2/coil junction), not at
the deeper calf cluster — suggesting calf-internal flexibility may
be a steering-induced motion that decouples from the bend/extend
hinge proper.

### Finding 4 — Interaction-network residues invisible to σ
Matsumoto's Asp550, Cys549, Leu548 (Interaction B) are at 38-57th
percentile in our σ — clearly *not* hinge residues by backbone
flexibility. They are still mechanistically important per Matsumoto
because they stabilize via salt-bridge / hydrophobic networks. The
v2 composite figure flags these as green (miss) overlays — making
clear that backbone-σ and interaction-energy methods are
complementary, not duplicative.

### Finding 5 — pipeline self-consistency at 0.72-0.82 corr (unchanged)
Forward-rendering the v7 fitted trajectories reproduces the real
Linz AFM at per-frame cosine correlation 0.82 (V1) and 0.72 (V2),
above random baseline (0.65, 0.43). This is a function of the
imaging model, not the flexibility analysis, so unchanged from v1.

## Open questions / next steps (revised)

1. **Multi-integrin pipeline** — αIIbβ3 + α5β1 + αVβ6 + αVβ8 first-
   principles bent/extended distribution. Plan in
   `docs/integrin_heterodimer_plan.md`.
2. **AF2-Multimer ablation** — reviewer push, 50 GPU-hours, can run
   on next A40 window.
3. **EO state coverage** — confirmed unreachable by SMD on 3-ns
   timescale (`obj-025`). Needs αIIbβ3 string-method, metadynamics,
   or REMD as alternative.
4. ~~Rotation-corrected RMSF~~ — done (`obj-027`).
5. ~~Matsumoto switch-residue overlay~~ — done (`obj-026`), 5/5
   backbone-hinge residues at ≥80th percentile.

## Provenance

| Claim | Artifact |
|-------|----------|
| Composite v2 figure | `figures/mechanical_sensitivity_composite_v2.png` |
| Composite v2 numerics | `figures/mechanical_sensitivity_composite_v2.json` |
| Composite v2 script | `pipelines/avb3-conformers/scripts/mechanical_sensitivity_composite.py` |
| Rotation-corrected RMSF | `figures/rmsf_v7_corrected/{video1,video2}/` |
| Matsumoto overlay | `figures/matsumoto_overlay/matsumoto_overlay.png` |
| Earlier report | `results/mechanical_sensitivity_report_v1.md` |
| Self-consistency | `results/afm_pipeline/sim_afm/README.md` |
