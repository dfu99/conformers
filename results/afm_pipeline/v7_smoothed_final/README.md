# v7 Smoothed Final — 2026-04-23

Canonical final overlays with:
- 615-frame bent+extended conformer library (CV0 range 47.3-85.0 Å)
- Head-anchored SO(3) fit (A4500 GPU, ~10 min per video)
- Tail-flip smoothing (V1 28.5%, V2 18.1% corrected)
- Rolling-median temporal smoothing (window=7, jitter V1 45→10.5Å, V2 38→9.9Å)
- Per-frame head re-anchor (overlays stay locked to tracked AFM head)

## Results
- video1/  — 379 frames, corr 0.965, BC 43.5% / Intermediate 6.3% / EC 50.1%
- video2/  — 1266 frames, corr 0.939, BC 18.6% / Intermediate 2.8% / EC 78.6%

## Regeneration
```
python3 pipelines/avb3-conformers/scripts/fit_with_head_tracking.py ... --device cuda
python3 pipelines/avb3-conformers/scripts/resolve_tail_flips.py --fitted-dir ... --window 21
python3 pipelines/avb3-conformers/scripts/smooth_coords_temporal.py --fitted-dir ... --overwrite
python3 pipelines/avb3-conformers/scripts/render_projection_overlay.py ...
```
