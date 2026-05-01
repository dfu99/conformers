# conformer-library — PRIMARY PIPELINE

**Status: validated (αVβ3 v7 library, 2026-04-22)**

Generate a steered-MD conformer library for a protein along a chosen
collective-variable (CV) axis (e.g., bent ↔ extended for integrins).
Outputs the `library.json` + per-frame PDBs that
`pipelines/sim-afm-video` consumes.

## Why this is a primary pipeline

The validated αVβ3 v7 library at
`data/runs/avb3/conformers/all_frames_bent_extended/` (615 frames
spanning CV0 47.3-85.0 Å) drove the v15 sim HS-AFM render. Future
heterodimers (αIIbβ3 next, then α5β1, αVβ6, αVβ8) will use the same
pipeline with adapted domain definitions. This README captures the
canonical workflow.

## Scope

In: a starting PDB (e.g., `1JV2.pdb` for αVβ3 bent), domain
definitions, target CV values + force constant.
Out: a directory of conformer PDBs + `library.json` consumable by
`pipelines/sim-afm-video`.

## Method

OpenMM steered MD with flat-bottom distance-CV bias on
domain-COM-pair distances. Two passes (one per direction):

1. **Bent direction** (`cv_distance_bent` preset): targets push the
   leg domain centroids closer together — generates BC conformers.
2. **Extended direction** (`cv_distance_extend` preset): targets pull
   them apart — generates EC conformers.

After each pass, frames are extracted from the production trajectory
at fixed intervals and converted to PDBs.

## Pipeline Stages

1. **Set up the system** (`scripts/run_domain_steering.py`):
   - Solvate the starting PDB in TIP3P water + Na+/Cl-
   - Equilibrate (NPT, 1 ns)
   - Minimize *forcefield-only* (lessons.md: minimization with steering forces hangs LBFGS for 30+ minutes)
2. **Apply steering** (`scripts/domain_steering.py`):
   - Add a `CustomCentroidBondForce` per CV pair (flat-bottom or harmonic)
   - Reinitialize the Simulation context
   - Run production with target schedule
3. **Extract frames** (`scripts/extract_frames_from_nc.py`):
   - Read `production.nc` with mdtraj
   - Strip waters, save protein-only PDBs at fixed step interval
4. **Build library metadata** (`scripts/build_library_metadata.py`):
   - For each frame, compute the CVs and assign state
   - Define `trajectory_order` (frame_0000 → frame_0001 → … is biophysically adjacent because consecutive frames in MD)
   - Write `library.json`

## Domain Definitions (αVβ3 example)

```python
AVB3_DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 956),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}
```

For other heterodimers, map the αVβ3 domain residue ranges via
sequence alignment. αIIbβ3 mapping in
`pipelines/avb3-conformers/scripts/map_aiib3_to_avb3_domains.py`.

## Validated Defaults

| Parameter | Value | Reason |
|---|---|---|
| Force constant `k` | 200 kJ/mol/nm² | Sufficient bias without unfolding (αVβ3 bent/extend) |
| Production length | 3 ns | Drives full BC → EC transition |
| Frame interval | 5000 steps | ~6 ps per frame, 500 frames per run |
| Force field | AMBER14SB + TIP3P | Standard for soluble proteins |
| Time step | 2 fs | Hydrogen bonds constrained |
| Temperature | 310 K | Body temperature |
| Pressure | 1 bar | NPT |
| Disk-pressure mitigations | `--report-interval 5000`, watchdog deletes equilibration files after production starts | RunPod 1.5 GB quota fix |

## Run on RunPod (A4500 or A40)

```bash
# (after rsync of project to pod, see tasks/lessons.md for setup gotchas)
ssh root@<pod_ip> -p <pod_port>
cd /workspace/conformers
pip install 'numpy<2'  # OpenMM ABI requirement
pip install --upgrade pyparsing  # mdtraj.reporters silent-failure fix
nohup bash pipelines/conformer-library/run_steering.sh \
    --start-pdb data/avb3/1JV2_clean.pdb \
    --preset cv_distance_bent \
    --output /workspace/runpod_bent \
    > steering.log 2>&1 & disown
```

## Reproducing the αVβ3 v7 library

```bash
# Bent steering
bash run_steering.sh \
    --start-pdb data/avb3/1JV2_clean.pdb \
    --preset cv_distance_bent \
    --target-cv0 4.0 --target-cv1 3.5 \
    --output data/runs/avb3/conformers/bent_run

# Extend steering (from a different seed)
bash run_steering.sh \
    --start-pdb data/avb3/1JV2_clean.pdb \
    --preset cv_distance_extend \
    --target-cv0 8.5 --target-cv1 8.0 \
    --output data/runs/avb3/conformers/extend_run

# Merge + build library metadata
python scripts/build_library_metadata.py \
    --frame-dirs data/runs/avb3/conformers/bent_run/frames data/runs/avb3/conformers/extend_run/frames \
    --domains avb3 \
    --output data/runs/avb3/conformers/all_frames_bent_extended
```

## Confirmed Negative Results

- **CV2 (head-head) opening with classical SMD**: even at k=1000
  kJ/mol/nm², SMD opens CV2 by only 0.07 Å/ps in αVβ3 (3 ns budget
  insufficient). See obj-025. To recover EO state coverage: αIIbβ3
  string-method templates (Ferg-Lab), metadynamics, or REMD.

## Output layout

```
data/runs/<heterodimer>/conformers/<library_name>/
├── frame_0000.pdb
├── frame_0001.pdb
├── ...
├── library.json              ← consumed by pipelines/sim-afm-video
└── steering_metadata.json    ← preset, k, target CVs, run timestamps
```

## Dependencies

- Python: `openmm` 8.x (with CUDA), `mdtraj`, `numpy<2`, `pyparsing` (latest)
- GPU: A4500 or larger (~10 h per direction for αVβ3)
- Storage: ~6 GB per 3-ns trajectory before frame extraction

## See Also

- `pipelines/sim-afm-video/` — consumes the libraries this pipeline produces
- `pipelines.md` — top-level pipeline index
- `tasks/lessons.md` — gotchas: numpy 2 ABI, pyparsing reporter bug,
  steering force minimization hang, RunPod disk quota
