# Rotation-Corrected RMSF — v1 (2026-04-29)

## Why

The original `figures/rmsf_v7/{v1,v2}/` per-residue RMSF was dominated
by global rotation: the head centroid was anchored, but orientation
varied freely across frames, so even rigid headpiece residues showed
RMSF of 40-55 Å. This blocked any clean attribution of flexibility to
real backbone motion vs. coordinate-frame rotation.

## What changed

Added `--head-align` flag to `pipelines/avb3-conformers/scripts/residue_rmsf.py`.
Before computing RMSF, every frame is Kabsch-aligned to frame 0 using
only the rigid headpiece CAs (αV residues 1–440 + β3 residues 1–350,
790 CAs total). The full atomic frame is then rotated and translated
to bring the head into reference coincidence, leaving genuine
non-rotational motion in the residual.

## Result

| Domain | V1 orig (Å) | V1 corrected (Å) | V1 ratio | V2 orig (Å) | V2 corrected (Å) | V2 ratio |
|---|---:|---:|---:|---:|---:|---:|
| α-head-thigh | 43.7 | **8.3** | 0.19 | 53.5 | **7.4** | 0.14 |
| β-head | 42.7 | **9.7** | 0.23 | 56.2 | **8.5** | 0.15 |
| α-calf | 45.4 | 20.5 | 0.45 | 66.3 | 20.6 | 0.31 |
| β-tail | 62.2 | 25.7 | 0.41 | 82.7 | 26.4 | 0.32 |
| α-coil | 87.1 | 44.5 | 0.51 | 118.9 | 46.9 | 0.39 |

Mean RMSF: V1 53.5 → 19.3 Å (64% reduction); V2 71.2 → 19.3 Å (73%).

Rotation-corrected RMSF passes the consistency test:

1. **Headpiece is near-rigid** — RMSF 7-10 Å, consistent with internal
   side-chain plus loop motion only. The 0.14-0.23 reduction ratio
   confirms the residual after rotation removal is small at the head.
2. **Legs and coils retain real motion** — α-calf (20.5/20.6 Å) and
   β-tail (25.7/26.4 Å) still show meaningful flexibility because
   the headpiece-anchored leg sweep is *real* internal motion, not
   global rotation.
3. **C-terminal coils dominate (V1 top-20)** — A:758-766, A:909-912,
   A:956-962, B:691-692. Same ranking as the original cross-conformer
   analysis (`figures/cross_conformer_v7/top20_flexible_residues.md`).
   Confirms the C-terminal coil dominance finding survives rotation
   correction.

## Two videos consistent

V1 and V2 give nearly identical rotation-corrected mean RMSF
(19.29 vs 19.30 Å). Originally they differed (53 vs 71 Å) because
V2's longer trajectory (1266 vs 379 frames) accumulated more
rotational drift in the head-anchored frame. Rotation correction
brings them to the same flexibility baseline, as expected for the
same molecular system.

## Provenance

| Item | Path |
|------|------|
| Updated script | `pipelines/avb3-conformers/scripts/residue_rmsf.py` |
| Comparison util | `pipelines/avb3-conformers/scripts/rmsf_compare.py` |
| V1 corrected RMSF | `figures/rmsf_v7_corrected/video1/` |
| V2 corrected RMSF | `figures/rmsf_v7_corrected/video2/` |
| V1 comparison plot | `figures/rmsf_v7_corrected/comparison_video1.png` |
| V2 comparison plot | `figures/rmsf_v7_corrected/comparison_video2.png` |

## Implications for the paper

The mechanical-sensitivity composite figure
(`figures/mechanical_sensitivity_composite.png`) currently uses the
rotation-dominated RMSF as one of three metrics. The rotation-
corrected version should replace it before submission so the
composite represents true backbone flexibility rather than mostly
global head rotation. This is a one-line script change; the
top-residue ranking is preserved, but absolute values become
interpretable.
