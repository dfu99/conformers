# Multi-Integrin Pipeline Plan — 2026-04-23

## Scope

Reproduce the αVβ3 workflow (conformer library via MD steering + pseudo-AFM
imaging model + per-frame HS-AFM fitting) across additional integrin
heterodimers. The goal is to compare predicted bent/extended conformational
state distributions from *first principles* against published empirical
measurements *without* letting publications bias the conformer generation.

## Integrins to target and PDB starting structures

The 24 mammalian heterodimers fall into 4 ligand-specificity families
([Luo & Springer, Annu Rev Immunol 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC3039929/)).
The RGD-binders are the most biologically studied and most similar to
our αVβ3 work:

| Heterodimer  | Class           | Key PDBs (bent / extended)  | Size (res) | Notes |
|--------------|-----------------|-----------------------------|-----------:|-------|
| αVβ3         | RGD             | 1JV2 bent, 1L5G extended    | ~1600      | done (v7) |
| αIIbβ3       | RGD (platelet)  | 3FCS ectodomain, 2K9J TM    | ~1600      | extended available |
| α5β1         | RGD             | 3VI4 headpiece + peptide, 3VI3 apo headpiece | ~1700 | headpiece only |
| αVβ5         | RGD             | 4MMX, 4MMZ                  | ~1600      | |
| αVβ6         | RGD             | 4UM9 head + peptide, 5FFG   | ~1600      | |
| αVβ8         | RGD (atypical)  | 6OM2, 6UJC                  | ~1600      | constitutively extended |
| α4β1         | LDV (non-RGD)   | e.g., 1QC5 (β1 headpiece)   | ~1700      | |
| α5β1 (full)  | RGD             | 7NXD  (cryo-EM half-bent)   | ~1700      | corrected 2026-05-13: prior 6WOV ID was wrong (RyR2, not an integrin) |
| αLβ2         | leukocyte       | 6TR1 (full extended), 3HI6  | ~1700      | |
| αMβ2         | leukocyte       | 7USM, 7USL                  | ~1700      | |

Compute budget: ~10 h of A4500 steering per direction per heterodimer = ~20 h/integrin.
Including the pseudo-AFM library build and correlation matching, ~25 h/integrin.
For 5 heterodimers = ~125 h ≈ 5 days on the shared A4500.

## First-principles estimation of bent/extended distribution

The asks: estimate the expected HS-AFM-observable bent/extended
distribution from structural features alone, *without* reading the
published activation kinetics. Factors we can use from first principles:

1. **Knee geometry intrinsic energy**: domain-domain steric contact area
   in the bent state. Larger contact ⇒ more stabilized bent ⇒ more bent
   in the population.
2. **Hinge residue count at the genu**: charged/hydrophobic residues at
   the α-thigh / α-calf and β-head / β-tail interfaces. More contacts ⇒
   stabilized bent.
3. **β-tail length and mass distribution**: longer legs shift the
   free-energy balance toward the extended state (more entropy).
4. **Sequence conservation**: high conservation at the genu across all
   RGD integrins suggests a shared hinge mechanism; divergence may
   reflect different default populations.

Concrete prediction strategy (to execute):
1. Download each heterodimer's bent (or bent-proxy) PDB.
2. Compute genu domain-domain contact area (SASA difference).
3. Compute analogous CVs (α-head-thigh ↔ α-calf distance, etc.).
4. Rank predicted-bent-fraction from most-bent (most contact, smallest
   legs) to most-extended.
5. Run steering-MD pipeline (bent + extended) for each.
6. Build pseudo-AFM library and sim AFM.
7. *Only then* compare to published empirical measurements (Chen et al.
   2023 for α5β1/αVβ3, Li et al. 2024 for αIIbβ3 cascade, etc.).

## Proposed pipeline template (per integrin)

1. Download bent PDB, clean (strip water, select ectodomain).
2. Run `cv_distance_bent` steering (force constant 200, 3 ns, A4500).
3. Run `cv_distance_extend` steering from a different seed.
4. Extract ~300 bent + ~300 extended protein frames.
5. Compute CVs using domain definitions specific to that integrin.
6. Build pseudo-AFM library (50 epochs × 500, batch=1).
7. First-principles predicted distribution (ranked).
8. Generate sim-HS-AFM for each heterodimer.
9. Post-hoc literature comparison.

## Not-yet-decided

- Domain definitions for each integrin: need to map resSeq ranges.
  αVβ3 has `alpha_head_thigh`, `alpha_calf`, `beta_head_hybrid_egf1`,
  `beta_tail_egf2_3_4_btail`. These are αVβ3-specific. For α5β1, β6, etc.,
  we need to either (a) reuse by sequence alignment or (b) define new
  domains per integrin.
- Real HS-AFM data for non-αVβ3 integrins. We only have Linz data for
  αVβ3. For the sim-AFM-only comparison, this is fine — the
  first-principles prediction is what we're validating.
