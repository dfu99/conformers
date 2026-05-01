# sim-afm-video — PRIMARY PIPELINE

**Status: validated (V2 αVβ3, 2026-04-30)**

End-to-end pipeline that converts a conformer library (PDBs +
metadata) into a realistic stylized HS-AFM video. Used as the
forward-modelling arm of the AFMFold structure-to-image template
matching framework.

## Why this is a primary pipeline

The canonical sim-AFM render at `results/afm_pipeline/sim_afm/video2_stable_zblur/sim_afm_copper_v15.gif`
combines four years of feedback into a fixed contract: any future
generations of the project must reproduce these characteristics or
explicitly justify why they have been changed. This README + the
`scripts/` directory + `run_pipeline.sh` capture exactly what
generates that result.

## Scope

In: a conformer library (folder of PDBs + `library.json` metadata).
Out: a sim HS-AFM `.gif` plus intermediate npy artifacts for audit.

Not in scope: producing the conformer library itself. Use the
`pipelines/conformer-library` pipeline (steered MD via OpenMM with
distance CV bias) or any other conformer source — see Input Contract
below for the schema.

## Input Contract

A directory with:

```
my_library/
├── frame_0000.pdb
├── frame_0001.pdb
├── ...
├── frame_NNNN.pdb
└── library.json
```

`library.json` schema:

```json
{
  "topology_pdb": "frame_0000.pdb",
  "conformers": [
    {"file": "frame_0000.pdb", "cv0_A": 52.3, "cv1_A": 42.1, "cv2_A": 35.0, "state": "BC"},
    {"file": "frame_0001.pdb", "cv0_A": 53.1, "cv1_A": 42.5, "cv2_A": 34.9, "state": "BC"},
    ...
  ],
  "trajectory_order": [0, 1, 2, ..., N-1],
  "adjacency": [[0, 1], [1, 2], ...]
}
```

- `topology_pdb`: the PDB whose atom ordering defines the topology
  (all other PDBs must have identical atom counts). Optional — defaults
  to the first conformer.
- `conformers`: list of frames. Required field is `file`. CVs and
  state labels are optional but recommended for downstream analysis.
- `trajectory_order`: optional. Sequence of indices defining playback
  order. Each consecutive pair must be biophysically adjacent.
- `adjacency`: optional. Edge list of biophysical neighbors. Used
  if `trajectory_order` is absent — the loader walks the graph
  greedily from node 0.
- If neither is present, indices are played 0…N-1 in file order.

## Pipeline Stages

1. **Ingest** (`scripts/load_conformer_library.py`): load PDBs in
   playback order, stack into `(T, N, 3)` array.
2. **Stabilize orientation** (`scripts/stabilize_orientation.py`):
   - PCA-flatten so smallest principal axis = +z (surface normal)
   - Side-lock: keep αV chain on the same z-face frame to frame
   - Sign-align long axis between consecutive frames
   - Stepwise yaw with hysteresis (50° threshold, 20-frame dwell, 30° cap)
   - Gaussian-σ-8 head xy smoothing
3. **Forward-render** (`scripts/simulate_afm_video.py`): afmfold
   `generate_landscape` + `idilation` (rounded tip), then
   - z-axis Gaussian blur (cantilever-oscillation averaging)
   - z-noise (cantilever readout RMS)
   - per-frame normalize
4. **Stylize** (`scripts/stylize_sim_afm.py`): SINGLE-CANVAS render —
   - upscale to working canvas
   - substrate noise + anisotropic blur on same canvas
   - tail-direction-correlated surface slant
   - row-wise baseline jitter
   - partial-width flash streaks (rare)
   - localized horizontal distortions on molecule rows
   - copper colormap + soft-clip

The single-canvas constraint is critical: all stylization steps
operate on one array shape, so no inset crop, no paste between
mismatched sizes, no internal rectangular border possible.

## Validated Defaults (αVβ3, V2 trajectory)

| Parameter | Value | Reason |
|---|---|---|
| `WIDTH`, `HEIGHT` | 60 px | Fits 26-nm protein extent at 0.98 nm/px with margin |
| `RES_NM` | 0.98 | Matches Linz HS-AFM pixel size |
| `TIP_RADIUS_NM` | 3.0 | Matches Linz tip; from PI feedback iteration |
| `Z_BLUR_NM` | 0.8 | Cantilever oscillation averaging; PI-validated |
| `Z_NOISE_NM` | 0.35 | RMS noise floor; PI-validated |
| `HEAD_RES_A` / `HEAD_RES_B` | 1-440 / 1-350 | αVβ3 headpiece definition |
| `SIDE_CHAIN` | A | αV-up |
| `YAW_THRESHOLD_DEG` | 50 | Stepwise yaw hysteresis trigger |
| `YAW_DWELL_FRAMES` | 20 | Frames above threshold before commit |
| `YAW_STEP_CAP_DEG` | 30 | Max yaw change per commit |
| `HEAD_ANCHOR_SIGMA` | 8 | Frames; head xy Gaussian smoothing |
| `UPSCALE` | 5 | 60→300 px working canvas |

All overridable via env var when invoking `run_pipeline.sh`.

## Run

```bash
bash run_pipeline.sh data/my_library results/my_sim_afm
```

Output: `results/my_sim_afm/sim_afm.gif` plus intermediate files
under `ingest/`, `stable/`, `render/`.

## Reproducing the V2 αVβ3 result

The validated v15 render uses `data/runs/avb3/conformers/all_frames_bent_extended`
as the source library (615 frames). To reproduce:

```bash
# From repo root
python pipelines/sim-afm-video/scripts/load_conformer_library.py \
    --library data/runs/avb3/conformers/all_frames_bent_extended \
    --output /tmp/avb3_v15_input
# (requires a library.json — see examples/avb3_library_template.json)
bash pipelines/sim-afm-video/run_pipeline.sh \
    /tmp/avb3_v15_input \
    results/sim_afm_pipeline/avb3_v15
```

For real-AFM-fitted variants (where the conformer library is reordered
to match an experimental video frame-by-frame), use
`pipelines/avb3-conformers/scripts/fit_with_head_tracking.py` to
produce `fitted_coords*.npy` and feed that directly into stages 2-4.

## Provenance of validated parameters

| Parameter | Validation feedback iteration |
|---|---|
| Tip 3 nm | obj-022 v6 stiff-head test |
| Z-blur 0.8 nm + z-noise 0.35 nm | 2026-04-30, response to "z-axis too precise" |
| Side-lock + stepwise yaw | 2026-04-30, response to "molecule rolls to other side" |
| Single-canvas stylization | 2026-04-30, fix for "internal border cropping" |
| 60-px canvas | 2026-04-30, response to "molecule getting cut off at edge" |

## Dependencies

- Python: `mdtraj`, `torch`, `numpy`, `scipy`, `Pillow`
- afmfold at `$AFMFOLD_ROOT` (default `/home2/Documents/code/afmfold`)
- Active env: `venv_protenix` works; any env with the above will work

## See Also

- `pipelines/conformer-library/` — generates the input libraries
- `pipelines.md` — top-level pipeline index
- Earlier scripts at `pipelines/avb3-conformers/scripts/` (all
  copies of the canonical scripts in this directory; kept for
  historical context)
