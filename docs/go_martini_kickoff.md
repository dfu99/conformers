# Switching Gō-Martini parallel track kickoff

**Audit follow-up §11.5 P=2.** Cross-validation track for route A
(αIIbβ3 string method) per `docs/eo_coverage_strategy.md` route E.

Independent line of evidence for the EC→EO transition using
coarse-grained MD with switching Gō-contacts between two
endpoint conformations.

---

## Concept

Switching Gō-Martini (JCTC 2024, Pesce & Marrink) lets the protein
relax under Gō contacts of conformation A, then smoothly switches
to Gō contacts of conformation B. The trajectory traverses the
free-energy minimum-energy path between A and B at coarse-grained
resolution — typically 50× faster than atomistic MD.

Concept figure: `figures/go_martini_eo_trajectory.png` shows the
expected (CV0, CV2) trajectory the switching run should produce.

### Endpoints

- **Gō-A (bent)**: αVβ3 1JV2, CV0 ≈ 51, CV2 ≈ 36
- **Gō-B (EO target)**: literature EO cryo-EM reference. Eng et al.
  *Science* 2018 estimated αVβ3 EO has CV0 ≈ 92, CV2 ≈ 58. No
  high-resolution PDB available (per obj-041 null result), so use
  the cryo-EM-derived CV target as the Gō-B endpoint.

### Why CG and not atomistic for cross-validation

Independence is what we want. Atomistic string method (route A) has
a single source of error from one force-field family. CG Gō-Martini
has independent errors (different force-field, no Mg²⁺/Mn²⁺
explicit metal, simplified electrostatics). If the two methods
agree on the path topology and barrier locations, we have
evidence that the result is force-field-independent.

---

## Software stack

| component | tool | install needed |
|-----------|------|----------------|
| topology builder | Martinize3 | pip-installable |
| switching contacts | Switching Gō-contact generator | from JCTC 2024 supplementary |
| MD engine | GROMACS 2024+ with Martini patch | NEW dependency for project |
| CG-to-AA back-mapping | backward.py (Marrink lab) | pip-installable |
| CV analysis | reuse `pipelines/conformer-library/scripts/build_library_metadata.py` | already exists |

GROMACS install conflicts with the project's OpenMM-only stance.
Mitigation: install in a separate `venv_gomartini` env, isolated
from `venv_protenix` and `venv_boltz`.

---

## Day-by-day plan (2 weeks)

### Week 1

- **Day 1**: Install GROMACS 2024 + Martini patch in `venv_gomartini`.
- **Day 2**: Run Martinize3 on 1JV2 (αVβ3 bent) → CG topology.
- **Day 3**: Build the αVβ3 EO target structure. Two options:
  - (a) Homology model from αIIbβ3 extended via SWISS-MODEL (free)
        using the αV→αIIb residue mapping from
        `results/route_a/av_aiib_alignment.json`.
  - (b) Use the highest-CV0 v7 library frame as a near-EC endpoint
        and let the Gō-B contacts "open" the headpiece via CV2.
- **Day 4**: Generate switching Gō-contact tables (A + B).
- **Day 5**: Set up POPC bilayer + solvate + ionize + equilibrate
  for 50 ns. Sanity check: bilayer thickness ≈ 4.0 nm.
- **Day 6-7**: Production switching run, 500 ns CG (≈ 1 µs effective
  atomistic time). Estimated 6 hr on a single A40.

### Week 2

- **Day 8-9**: Back-map CG trajectory to atomistic via backward.py.
  Verify chain integrity + secondary structure with mdtraj.
- **Day 10**: Compute CV0/CV1/CV2 trajectory via existing
  `build_library_metadata.py`-style domain centroids.
- **Day 11**: Cross-compare to expected trajectory in
  `figures/go_martini_eo_trajectory.png`.
- **Day 12-13**: Compute MIDAS pocket SASA (reuse obj-039 script)
  on every 10th frame; verify EO endpoint *opens* the pocket
  (validates obj-039's negative result).
- **Day 14**: Write results doc + commit.

Total wall time: 14 days. Compute: 1× A40 for 6 hr (CG production)
+ ~12 hr post-processing on CPU. Cost: ~$200 GPU + free CPU.

---

## Expected deliverables

1. `results/afm_pipeline/go_martini/`:
   - `cg_trajectory.xtc` (raw CG trajectory)
   - `aa_backmapped.xtc` (atomistic back-map)
   - `cv_trajectory.npy` (CV0/CV1/CV2 per frame)
   - `midas_sasa_trajectory.npy`
2. `figures/go_martini_eo_trajectory.png` ← already exists as
   concept; will be **replaced** by actual trajectory.
3. `figures/go_martini_midas_sasa.png` — MIDAS opening through
   the trajectory.
4. obj-042 in `tasks/objectives.yaml` documenting the result.

---

## Acceptance criteria

The cross-validation **passes** if:

1. The switching trajectory traverses the expected (CV0, CV2)
   region monotonically (no large detours).
2. Final-frame CV0 reaches ≥ 80 Å AND CV2 reaches ≥ 50 Å
   (i.e., enters EO band).
3. MIDAS pocket SASA *increases* on extension (opposite of
   obj-039's negative result on the v7 EC frames) — confirming
   that the EO transition exposes the binding site.
4. CG trajectory + atomistic back-map remain stable (no chain
   breaks, secondary-structure preservation > 80 % in headpiece).

If 1 + 2 + 3 + 4 all pass, the cross-validation supports route A's
upcoming string-method result. If only some pass, the failure mode
informs route A's CV definition or membrane setup choice.

---

## Risk register (light, vs route A's full one)

| risk | P_fail | mitigation |
|------|--------|------------|
| GROMACS install conflicts | 10 % | isolated `venv_gomartini` env |
| CG → AA back-mapping artifacts | 20 % | secondary-structure check post-mapping |
| Switching too fast (kinetic trap) | 30 % | adjust switching schedule per JCTC 2024 §3.2 |
| Endpoint geometry too distant | 15 % | use intermediate library frame as switching anchor |
| Lipid composition wrong | 10 % | match route A's endothelial choice (POPC 50 % + cholesterol 25 %) |

Joint pass probability ≈ 50 % (independent risks).

---

## Cross-references

- `docs/eo_coverage_strategy.md` route E (this is the implementation)
- `docs/route_a_kickoff.md` (primary route, this is the parallel
  cross-validation)
- `docs/route_a_risk_register.md`
- `tasks/audit-2026-05-05.md` §11.5 P=2
- `figures/go_martini_eo_trajectory.png` (concept figure)
- obj-039 (RGD pocket buries on v7 extension — to be reversed by EO)
- obj-041 (no published αVβ3 EO PDB — drives need for CG modeling)

---

_Author: AFK audit deepening pass v3, 2026-05-05 late evening._
_Status: pre-kickoff. Activate after route-A starts to overlap
weeks; can optionally run independently if route-A is delayed._
