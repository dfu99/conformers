# Pipelines Index

This file is an index of canonical, validated pipelines in this repo.
Pipelines marked **PRIMARY** are the ones that generate
publication-grade outputs and should be preserved unchanged across
future refactors. New work should plug into these contracts rather
than reinventing them.

## PRIMARY PIPELINES

The framework has *three* PRIMARY pipelines, forming a complete
forward + inverse workflow:

```
                   pipelines/conformer-library/
                       (steered MD)
                            │
                            ▼
                       library/
                  (PDBs + library.json)
                          ╱  ╲
                         ╱    ╲
                        ▼      ▼
        pipelines/         pipelines/
        sim-afm-video/     afm-overlay/
        (forward model)    (template-match real video)
              │                  │
              ▼                  ▼
         sim_afm.gif       pdb_projection_video.gif
                           + state kinetics
```

### `pipelines/sim-afm-video/` — sim HS-AFM video generator
**Status:** validated 2026-04-30 (V2 αVβ3, conformer library v7).
**Purpose:** Convert a conformer library (PDBs + metadata) into a
realistic stylized HS-AFM video for forward-modelling, template
matching, or qualitative comparison to experimental data.

**Input contract:** Folder of PDBs + `library.json` (per-frame CVs,
state, biophysical adjacency / trajectory order). See
`pipelines/sim-afm-video/README.md` § Input Contract.

**Pipeline stages:**
1. Ingest (`load_conformer_library.py`)
2. Stabilize orientation (PCA-flatten + side-lock + stepwise yaw +
   head xy anchor) (`stabilize_orientation.py`)
3. Forward-render through afmfold imaging model with z-blur
   (cantilever-oscillation averaging) and z-noise
   (`simulate_afm_video.py`)
4. Stylize on a single canvas — copper colormap, substrate noise,
   anisotropic blur, slant, row jitter, flash streaks
   (`stylize_sim_afm.py`)

**Output:** stylized GIF at the specified output directory plus
intermediate npy artifacts for audit.

**Validated reference render:**
`results/afm_pipeline/sim_afm/video2_stable_zblur/sim_afm_copper_v15.gif`

**Master entry:** `bash pipelines/sim-afm-video/run_pipeline.sh <library_dir> <output_dir>`

### `pipelines/afm-overlay/` — fit real HS-AFM video + render overlay
**Status:** validated 2026-04-22 (V1+V2 αVβ3 with v7 conformer library).
**Purpose:** Inference arm of AFMFold — match a real HS-AFM video
against the conformer library, fit per-frame rigid-body orientation,
render the fitted PDB overlay on top of the real video, and produce
state-kinetics analyses (BC/EC/EO/Intermediate occupancy).

**Input contract:** Same `library.json` schema as sim-afm-video,
plus a `.gif` of real HS-AFM data. See
`pipelines/afm-overlay/README.md` § Input Contract.

**Pipeline stages:**
1. Ingest conformer library
2. Build pseudo-AFM library for correlation matching (`process_frames_to_afm.py`)
3. Per-frame fit: correlation + SO(3) refinement + head tracking
   (`fit_with_head_tracking.py`)
4. Tail-flip resolution (`resolve_tail_flips.py`)
5. Temporal smoothing — rolling median + head re-anchor
   (`smooth_coords_temporal.py`, `reprocess_smooth_coords.py`)
6. State kinetics classification + plots (`analyze_state_kinetics.py`)
7. Overlay rendering (`render_projection_overlay.py`)

**Output:** PDB-overlay GIF + state-kinetics PNG + per-frame fitted
coords + correlation arrays.

**Validated reference fits:**
- `results/afm_pipeline/v7_smoothed_final/video1/` (379 frames, mean corr 0.97, BC 43.5%)
- `results/afm_pipeline/v7_smoothed_final/video2/` (1266 frames, mean corr 0.94, BC 18.6%)

**Master entry:** `bash pipelines/afm-overlay/run_pipeline.sh <afm_gif> <library_dir> <output_dir>`

### `pipelines/conformer-library/` — steered-MD conformer library
**Status:** validated 2026-04-22 (αVβ3 v7 library, 615 frames).
**Purpose:** Generate a CV-spanning conformer ensemble from a starting
structure via OpenMM steered MD with distance-CV bias. Output is
directly consumable by `pipelines/sim-afm-video`.

**Input contract:** Starting PDB + domain definitions (chain + resSeq
ranges) + steering preset (target CVs + force constant).

**Pipeline stages:**
1. Solvate, equilibrate, minimize *forcefield-only* (`run_domain_steering.py`)
2. Add CV bias as `CustomCentroidBondForce`, reinitialize Simulation,
   run production (`domain_steering.py`)
3. Extract protein-only PDB frames at fixed step intervals
   (`extract_frames_from_nc.py`)
4. Build `library.json` with per-frame CVs + state assignment
   + trajectory adjacency (`build_library_metadata.py`)

**Output:** `data/runs/<heterodimer>/conformers/<library_name>/` —
PDBs + `library.json` + `steering_metadata.json`.

**Validated reference library:**
`data/runs/avb3/conformers/all_frames_bent_extended/` (615 frames,
CV0 47.3-85.0 Å spanning bent → extended)

**Confirmed limitation:** classical SMD cannot open the αVβ3
headpiece (CV2) in 3-ns budget at any practical force constant.
See `tasks/objectives.yaml` obj-025.

## Supporting / experimental pipelines

These are useful but not primary; refactoring them is OK.

### `pipelines/avb3-conformers/scripts/` — AVB3 utilities
Mixed-purpose scripts: rigid-body fitting against real AFM, mechanical
sensitivity analysis (RMSF, cross-conformer std, bond-angle σ),
Matsumoto 2008 overlay, multi-integrin first-principles features.
Some scripts are duplicated under primary pipelines (originals
preserved here for git history).

### `pipelines/protenix-avb3-template/` — Protenix template inference
Template-guided structure prediction for αVβ3 conformers via
Protenix. Used as one of multiple starting structure sources.

### `pipelines/protenix-a5b1/` — Protenix A5B1 staged pipeline
α5β1 + tag attachment co-fold via Protenix on PACE A100. Sub-pipeline
unrelated to the AFMFold framework.

### `pipelines/boltz/`, `pipelines/afcluster/` — Boltz/AFCluster sweeps
Alternative structure-prediction backends. Used for comparison and
alternative conformer sources.

## Adding a new pipeline

If you build a new pipeline that produces a publication-grade output,
promote it to PRIMARY by:
1. Creating its own directory under `pipelines/`
2. Writing a README following the
   `pipelines/sim-afm-video/README.md` template (Status, Why-primary,
   Scope, Input contract, Stages, Validated defaults, Run, Reproducing,
   Provenance, Dependencies)
3. Adding an entry to this `pipelines.md`
4. Updating the `pipelines/sim-afm-video` and `pipelines/conformer-library`
   READMEs if they should consume your output

Non-primary pipelines just need a brief entry below the primary block.
