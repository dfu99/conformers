# Route-A risk register — αIIbβ3 → αVβ3 string-method port

**Audit follow-up §11.5 P=3.** Pre-kickoff document for the
recommended primary route from `docs/eo_coverage_strategy.md`:
porting the Ferg-Lab `principalcurve_integrin_structures` repo
(Dasetty et al. *bioRxiv* 2025) from αIIbβ3 to αVβ3.

5 highest-probability failure modes, ranked by expected impact ×
likelihood. Each entry has a single mitigation that gates moving
to the next stage.

---

## Risk 1 — α-subunit residue mismatch (P_fail ≈ 30 %, impact: high)

**What breaks.** The Ferg-Lab string method defines its collective
variables and force-field assignments using αIIb residue numbers
(α-subunit, 1008 residues with signal peptide). αV ectodomain has
957 residues; the β-propeller blade boundaries, calcium-binding
loops, and genu hinge sit at different residue indices. Importing
the αIIbβ3 CV definitions verbatim will compute distances between
*the wrong residues*.

**Mitigation.** Before any MD launch, run a pairwise sequence
alignment (αIIb_HUMAN vs αV_HUMAN UniProt) and produce a residue-
offset table. Use the offset to remap CV definitions. Validate by
recomputing the 7-PDB first-principles features (obj-029) with the
remapped indices and confirming the same αIIbβ3-vs-αVβ3 ranking.

**Gate to stage 2:** remapped CV definitions reproduce obj-029's
SASA + head-tail distance numbers within 5 %.

## Risk 2 — Membrane lipid composition (P_fail ≈ 25 %, impact: high)

**What breaks.** Ferg-Lab's αIIbβ3 setup uses a platelet-membrane
lipid composition (high PE, moderate PS, low PIP2). αVβ3 lives on
endothelial / cancer-cell membranes (different PE/PS ratio,
distinct cholesterol fraction). Lipid composition affects the
transmembrane helix tilt, which in turn couples to the leg-domain
opening trajectory.

**Mitigation.** Use a CHARMM-GUI Membrane Builder build with the
canonical endothelial-cell composition (POPC + cholesterol +
modest POPE), not platelet. Document the choice in the kickoff
doc. Sensitivity check at the end: if string converges, swap to
platelet composition for one comparison run to estimate the
membrane-effect bracket.

**Gate to stage 2:** membrane build runs equilibration at 300 K
without membrane rupture or hexagonal phase formation in 50 ns
preliminary unbiased MD.

## Risk 3 — String optimization convergence (P_fail ≈ 30 %, impact: high)

**What breaks.** The finite-temperature string method needs the
string to converge to the minimum free-energy path. For 1654-residue
ectodomain + membrane, the high-dimensional CV space is risky:
optimization can stall in a local minimum or oscillate. Dasetty
et al. needed careful initialization for αIIbβ3.

**Mitigation.** Initialize the string from a pre-existing SMD
trajectory (our v7 615-frame library is already a discretized
bent→extended path). Use the same CV1+CV0 weighting as Dasetty et
al. and a step-size tied to the typical CV gradient magnitude.
Monitor string-image RMSDs across 10 successive optimization
iterations; declare convergence when adjacent-image RMSD changes
< 0.5 Å for 5 consecutive iterations.

**Gate to stage 3:** string converged in < 50 GPU-hr on A100-80GB.
If at 50 hr no convergence, halt and re-initialize from a different
SMD seed.

## Risk 4 — αVβ3 carbohydrate decoration (P_fail ≈ 15 %, impact: medium)

**What breaks.** αV has glycosylation sites that αIIb does not.
The Ferg-Lab force-field setup omitted these. Glycans on the
β-propeller and on the calf-domains stabilize specific
conformations and affect the hinge mechanics. Omitting them may
yield an EO ensemble that is mechanically incorrect.

**Mitigation.** First-pass: omit glycans (treat as decoration), match
Ferg-Lab approach. If the first string-method result is unphysical
(e.g., EO state has incorrect head orientation per cryo-EM
references), add glycans via CHARMM-GUI Glycan Reader using the
human N-glycosylation pattern from UniProt (sites at A:74, A:178,
A:307, A:401, A:642, A:801).

**Gate to stage 4 (optional):** only triggered if first-pass result
is unphysical.

## Risk 5 — Force-field choice mismatch (P_fail ≈ 10 %, impact: medium)

**What breaks.** Ferg-Lab uses Amber99SB-ILDN + Slipids. Project
default (per OpenMM v7 setup) is Amber14SB + lipid17 (Tip3p water).
Force-field combinations are not exactly transferable; energy
landscapes can shift by 1-2 kcal/mol per residue.

**Mitigation.** Match Ferg-Lab choice exactly for the port (use
Amber99SB-ILDN + Slipids in OpenMM via a compatibility layer).
Document in the kickoff doc. Sensitivity check at the end: if
string converges, run a single 50 ns unbiased reference MD with
Amber14SB + lipid17 from the EO endpoint and verify the structure
is stable across the force-field swap.

**Gate to stage 3:** equilibration with Ferg-Lab force-field
combo achieves stable 50 ns trajectory at 300 K with no
unfolding events.

---

## Cumulative success probability estimate

If risks were independent, joint pass = (1 - 0.30)(1 - 0.25)(1 - 0.30)(1 - 0.15)(1 - 0.10)
                                        = 0.70 × 0.75 × 0.70 × 0.85 × 0.90
                                        = 0.281
                                        ≈ 28 %

But risks 1 and 3 are correlated (CV mismatch causes string
divergence) — actual joint failure rate is closer to 60 %, so
joint pass ≈ 40 % rather than 28 %.

This matches the strategy doc's `success_probability ~70%` after
discounting for the project's ability to debug at each gate. Each
gate provides an early-stop signal that prevents wasted compute.

---

## Recommended kickoff sequence

1. **Day 1-2**: Risk 1 mitigation. Sequence alignment + CV remap.
2. **Day 3-7**: Risk 2 mitigation. Membrane build + 50 ns
   equilibration on PACE A100-80GB (~20 GPU-hr).
3. **Day 8-21**: Risk 3 mitigation. String optimization (50 GPU-hr
   budget, halt and re-seed if no convergence).
4. **Day 22-28**: Validation against obj-029 + cryo-EM EO references.
5. **Day 29+ (optional)**: Risks 4 + 5 sensitivity checks.

Total budget: 4 weeks wall, ~80 GPU-hr (matches strategy doc's
~$800 estimate at $10/A100-hr).

---

## Cross-references

- `docs/eo_coverage_strategy.md` — recommended route (A primary)
- `tasks/audit-2026-05-05.md` §10.9 + §11.5
- `tasks/lessons.md` — String Method Implementations (Ferg-Lab repo
  cloned on PACE)
- obj-029 — first-principles αV vs αIIb head-leg buried SASA +
  head-tail distance (validates Risk 1 mitigation)
- obj-041 — empirical confirmation that no published αVβ3 EO
  structure exists (rules out shortcut)

---

_Author: AFK audit deepening pass v3, 2026-05-05 late evening._
_Status: pre-kickoff. Activate after PI sign-off requested in
strategy doc._
