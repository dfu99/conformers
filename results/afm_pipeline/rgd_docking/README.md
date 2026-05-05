# RGD-binding-pocket accessibility (obj-039, 2026-05-05)

Reviewer-E PoC (audit-2026-05-05 §8 F4). Compares the 5 most-bent and 5
most-extended αVβ3 v7 library conformers on a geometric proxy for
RGD-binding feasibility.

## Why a proxy and not Vina/smina

The originally requested deliverable was Vina or smina docking. Two
practical blockers:

- αVβ3 RGD binding requires Mn²⁺ at the MIDAS site. AutoDock Vina has
  no native metal parameters, so Mn²⁺ would have to be added as a
  custom atom type or replaced with a bonded model — not a one-day
  PoC.
- Receptor PDBQT prep on a 1654-residue heterodimer needs rdkit/meeko
  with gemmi. The project's strict venv split (`venv_protenix` vs
  `venv_boltz`, CLAUDE.md §4) treats gemmi conflicts as a hard
  guardrail. Adding rdkit/meeko/gemmi to the working env risked
  breaking the protenix pipeline mid-week.

The geometric proxy (SASA + occlusion + steric-blocker distance)
answers the same first-pass question — *does extension expose the
RGD-binding pocket?* — without the dependency hazard.

## Method

`pipelines/conformer-library/scripts/rgd_pocket_accessibility.py`

1. Compute CV0 (αV head ↔ αV calf centroid) for all 615 v7 library
   PDBs. Cache to `library_cv0_cache.npz`.
2. Pick the 5 lowest-CV0 (most-bent) and 5 highest-CV0 (most-extended).
3. For each, compute three metrics:
   - **MIDAS pocket SASA** — Shrake-Rupley over the 10 canonical
     MIDAS-shell residues (αV D218; β3 D119, S121, S123, Y122, R214,
     N215, R216, E220, D251). Sum of per-atom SASA in Å².
   - **Pocket occlusion** — heavy-atom count inside a 12 Å sphere
     around the MIDAS centroid, normalized by 50 (empirical bent-state
     reference). Higher = more buried.
   - **MIDAS-leg distance** — distance from MIDAS centroid to the
     mean of (αV calf centroid, β3 tail centroid). Confirms the
     steering geometry.

## Result

| metric | bent (n=5) | extended (n=5) | Δ |
|--------|-----------|---------------|----|
| CV0 (Å) | 47.50 | 84.97 | +37.5 |
| MIDAS SASA (Å²) | 459.4 | 298.1 | **−35.1 %** |
| Occlusion (norm) | 3.78 | 4.19 | +0.41 |
| MIDAS-leg distance (Å) | 60.3 | 112.0 | +51.7 |

The legs retreat 52 Å from the MIDAS centroid (steering worked), but
the MIDAS pocket itself becomes *less* solvent-accessible — a 35 %
SASA drop and a 11 % rise in steric occlusion.

## Interpretation

The result is the textbook EC-vs-EO distinction:

- "Extended" in our v7 library means CV0 head-leg distance is large.
  This is structurally EC (extended-closed): legs splayed, head
  conformation unchanged.
- The high-affinity RGD-binding state is EO (extended-open), which
  requires headpiece opening — αV β-propeller / β3 βA-domain
  hybrid-domain swing. None of our SMD trajectories produced that
  (cf. obj-025: k=1000 in 3 ns gives only 0.9 Å CV2 motion).
- Without EO frames, real Vina docking on these conformers would
  return similar affinities for bent and extended (both MIDAS-closed)
  and the reviewer-E narrative ("extension exposes RGD") would not be
  supportable from this dataset.

## Implication for the project

This PoC sharpens the audit's #1 blocker (EO coverage gap) rather than
dissolving it. The right next step is not more docking work but
generating EO frames — via αIIbβ3 string method (compute-bound),
metadynamics with CV2 (tooling-bound), or replica exchange. Logged in
audit §3 and `tasks/planning.md` Confirmed Negative Results.

## Artifacts

- `pocket_accessibility.json` — per-PDB metrics, MIDAS-residue list,
  method note.
- `library_cv0_cache.npz` — CV0 for all 615 library frames.
- `rgd_pocket_accessibility.png` — 4-panel comparison figure.
- Top-level audit copy: `figures/rgd_docking_v1.png`.

## Cross-references

- audit-2026-05-05 §8 F4 (queued as obj-039 here)
- obj-025 EC→EO SMD negative result
- obj-038 ΔG(CV0) free-energy profile (companion reviewer-D answer)
