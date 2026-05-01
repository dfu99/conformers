# afm-overlay — PRIMARY PIPELINE

**Status: validated 2026-04-22 (V1+V2 αVβ3 with v7 conformer library, 615 frames)**

End-to-end inference pipeline that fits a real HS-AFM video against a
conformer library and renders the fitted PDB overlay on top of the
real video. This is the AFMFold *template-matching* arm of the
framework — the structural-state inference output that pairs with
`pipelines/sim-afm-video` (the forward-modelling arm).

## Why this is a primary pipeline

The validated v7 fits at `results/afm_pipeline/v7_smoothed_final/{video1,video2}/`
gave the discriminating result that established the project: video1
shifted from 12.4% bent-conformer (BC) → 43.5% BC after the bent-state
library was added, while video2 stayed at ~18% BC. The library is
sensitive to the data, not the other way around. That conclusion
depends on this pipeline's correctness end-to-end. Future heterodimers
(αIIbβ3, α5β1, etc.) will reuse the same pipeline with adapted domain
definitions.

## Scope

In: a real HS-AFM video (`.gif`) + a conformer library (folder of
PDBs + `library.json` — the same input contract as `sim-afm-video`).
Out: a PDB-overlay video showing the fitted conformer projected on
each AFM frame, plus state-kinetics charts, CV trajectories, and
all intermediate npy artifacts for audit.

Not in scope: producing the conformer library or generating sim-AFM
forward-renders. Use `pipelines/conformer-library` and
`pipelines/sim-afm-video` for those.

## Input Contract

Identical to `pipelines/sim-afm-video`. See its README § Input
Contract for the `library.json` schema. The afm-overlay pipeline
additionally requires:

- `--afm-gif`: a single HS-AFM video (35×35 px or larger, will be
  cropped square + normalized).

## Pipeline Stages

1. **Ingest** — load PDBs in `library.json` order
2. **Pseudo-AFM library** (`process_frames_to_afm.py`) — render one
   pseudo-AFM image per conformer using the fixed tip parameters
   used for matching (default 2 nm tip, 20° angle, batch_size=1 to
   avoid the truncation bug, see `tasks/lessons.md`)
3. **Per-frame fit** (`fit_with_head_tracking.py`):
   - Skip first 30 frames (scan window stabilization)
   - Track AFM head position via thresholded centroid
   - Best-conformer correlation matching (cosine similarity)
   - SO(3) refinement (2048 random rotations per frame)
   - Anchor fitted PDB head to tracked AFM head
4. **Tail-flip resolution** (`resolve_tail_flips.py`):
   - Detect frames where head is stable but tail rotates 180°
   - Hysteresis-based fix with 21-frame rolling window
   - Up to 20 iterations of flip-relabeling for convergence
5. **Temporal smoothing** (`smooth_coords_temporal.py` +
   `reprocess_smooth_coords.py`):
   - Rolling-median (window=7) on `fitted_coords` directly
   - Per-frame head re-anchor to pre-temporal position to fix
     V2-length drift
6. **State kinetics** (`analyze_state_kinetics.py`):
   - Classify each frame as BC / EC / EO / Intermediate by CV
     thresholds
   - Plot occupancy time series + state-residence histograms
7. **Overlay render** (`render_projection_overlay.py`):
   - Project fitted PDB CAs onto real-AFM coordinate frame
   - Per-frame composite: AFM frame + PDB projection (color by
     chain or domain)
   - Save as GIF

## Validated Defaults

| Parameter | Value | Reason |
|---|---|---|
| `TIP_RADIUS_NM` | 2.0 | Calibrated to Linz HS-AFM (was 6-12 nm too large; see obj-014) |
| `TIP_ANGLE_DEG` | 20 | Standard AFM tip angle |
| `SKIP_FRAMES` | 30 | Scan window stabilization period in Linz data |
| `HEAD_RES_A` / `HEAD_RES_B` | 1-440 / 1-350 | αVβ3 headpiece definition |
| `TEMPORAL_WINDOW` | 7 | Rolling median window in frames; balances jitter vs CV reactivity |
| `TAIL_FLIP_WINDOW` | 21 | Wider window stabilizes tail-flip detection |
| `TAIL_FLIP_HYSTERESIS` | 0.05 | Lower hysteresis catches more flips without thrashing |
| `TAIL_FLIP_ITERS` | 20 | Most fits converge in <10 iterations |
| pseudo-AFM `batch_size` | 1 | Critical: batch>1 truncates sampling (see lessons.md) |

All overridable via env var when invoking `run_pipeline.sh`.

## Run

```bash
bash run_pipeline.sh \\
    inbox/linz_avb3_video2.gif \\
    data/runs/avb3/conformers/all_frames_bent_extended \\
    results/afm_pipeline/my_overlay
```

Outputs:
- `results/.../overlay/pdb_projection_video.gif` — overlay video
- `results/.../kinetics/state_kinetics.png` — state-occupancy plot
- `results/.../fit/fitted_coords*.npy` — per-frame fitted coordinates
- `results/.../fit/fitting_metadata.json` — tip, skip, correlations,
  per-frame matched conformer indices

## Reproducing the V2 αVβ3 v7 result

The validated v7 fit at `results/afm_pipeline/v7_smoothed_final/video2/`
was produced with the parameters above against
`data/runs/avb3/conformers/all_frames_bent_extended` and
`inbox/linz_avb3_video2.gif`. Mean per-frame correlation 0.94, BC
fraction 18.6%, EC 79%, mean fit time ~10 min on an A4500 (CPU
fitting takes ~4 h).

## Provenance of validated parameters

| Parameter | Validation feedback |
|---|---|
| Tip 2 nm | obj-013 (2026-04-11): identified pseudo-AFM tip was 6-12 nm — 3-6× too large vs Linz 1-2 nm |
| Skip 30 frames | obj-013: scan-stabilization observation in Linz data |
| `batch_size=1` | obj-016 (2026-04-16): `batch_size=16` truncated sampling to first 31 of 309 conformers |
| Head-tracked alignment | obj-014 (2026-04-14): fixes 10.4 nm cumulative drift from previous-frame anchoring |
| Tail-flip resolution | 2026-04-19: stronger settings (21-frame window, 0.05 hysteresis) caught 18-28% of frames missed by earlier 11-frame/0.1 settings |
| Rolling-median + head re-anchor | obj-021 (2026-04-22): 77% jitter reduction (V1 45→10Å, V2 38→9Å) vs per-frame fit |

## Dependencies

- Python: `mdtraj`, `torch`, `numpy`, `scipy`, `Pillow`
- afmfold at `$AFMFOLD_ROOT` (default `/home2/Documents/code/afmfold`)
- A GPU is helpful (4-12 min vs 4 h on CPU for SO(3) refinement)

## See Also

- `pipelines/sim-afm-video/` — forward-modelling arm; consumes the
  same `library.json` contract
- `pipelines/conformer-library/` — generates the libraries this
  pipeline consumes
- `pipelines.md` — top-level pipeline index
- `tasks/lessons.md` — gotchas: batch_size sampling bug, tip-size
  calibration, head-anchored vs previous-frame alignment
