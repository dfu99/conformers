# Literature Sweep — Integrin Mechanical Sensitivity + Dynamics

Background search done for AFK-LIT-1, focusing on whether our residue-level
flexibility findings converge with (or diverge from) known integrin biology.

## 1. Canonical hinge biology (established since Xiong et al. 2001)

Two well-known inter-domain flexibility points in integrins ([Xiong et al., Science 2001](https://pubmed.ncbi.nlm.nih.gov/11546839/), [Luo & Springer, Annu Rev Immunol 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC3039929/)):

1. **α-genu**: hinge between α-thigh and α-calf-1 — this is the textbook
   "knee" at which the α-subunit bends in the bent-closed state.
2. **β-knee**: analogous hinge between β-head+hybrid+EGF1 and the
   lower β-tail (EGF2–4 + β-tail). Extension requires simultaneous hinging
   at both.

Our pre-registered expectation for hinge identification matches: α-genu
around resSeq ~435 on chain A (boundary between α-head-thigh and α-calf)
and β-knee around resSeq ~352 on chain B (β-head / β-tail boundary).

## 2. Prior NMA on αVβ3 ectodomain

Matsumoto et al., ["Key interactions in integrin ectodomain responsible
for global conformational change detected by elastic network normal-mode
analysis"](https://pmc.ncbi.nlm.nih.gov/articles/PMC2527288/) (Biophys J
2008) applied ENM-NMA to the bent αVβ3 ectodomain and identified
"switch" residues that stabilize the bent state. Their low-frequency
modes reproduce the bent→extended transition and identify non-bonded
interactions (charged residues at the genu) that function as a "snap."

*Relevance to our work:* We can compare our per-residue angular variance
hinge candidates against their switch residue list. Agreement strengthens
the interpretation; disagreement flags where HS-AFM-derived flexibility
differs from the vacuum NMA picture.

## 3. Force-regulated bistability (recent HS-AFM / single-molecule)

Chen et al., ["Force-Regulated Spontaneous Conformational Changes of
Integrins α5β1 and αVβ3"](https://pubs.acs.org/doi/10.1021/acsnano.3c06253)
(ACS Nano 2023) show that αVβ3 is *bistable even without force* and
spontaneously transitions between bent and extended. This directly
supports our v7 result: video1 shows ~50/50 BC/EC population, i.e., a
spontaneously transitioning molecule.

Li et al., ["Ligand binding initiates single-molecule integrin
conformational activation"](https://www.sciencedirect.com/science/article/pii/S0092867424004756)
(Cell 2024) directly observe the conformational activation cascade from
bent→extended→open using single-molecule methods. Our data is compatible
with the bent→extended axis but does not resolve the closed→open
transition (CV2 constant at ~35 Å across our library — noted as the key
gap in `intuition.md`).

## 4. Scheuring lab HS-AFM capabilities

The [Scheuring Lab](https://scheuringlab.com/research/) established the
instrumental baseline: ~1 nm lateral, 0.1 nm vertical, 100 ms temporal
resolution. Our Linz data is at comparable resolution. They have
developed super-resolution LAFM ([Nature 2021](https://www.nature.com/articles/s41586-021-03551-x))
which could in principle improve our per-pixel fitting if the raw AFM
data is available at higher sampling.

## 5. DCC / flexibility analysis methodology

[MD-TASK documentation](https://md-task.readthedocs.io/en/latest/corr.html)
and [mDCC analysis (PLOS One 2014)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112419)
document the standard dynamic cross-correlation matrix approach. The
multi-modal extension handles transient correlations.

*Relevance:* our cross-conformer 615×615 RMSD matrix (AFK-MECH-2) is the
coarser cousin of DCC — it captures structural-state covariance rather
than per-frame dynamics. Per-residue variance across the library is
analogous to a one-dimensional flexibility fingerprint.

## 6. β3-integrin headpiece hinge angle MD

Literature search surfaced [Integrin Conformational Dynamics and
Mechanotransduction (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9688440/)
which notes MD simulations showed the Ångstrom-level structural pathway
of ligand-induced hinge-angle opening. This is another reference for
comparing our hinge candidates to published hinge residues.

## Synthesis and gaps

- Our approach differs from prior work in that it reads flexibility
  *out of HS-AFM imaging data* via the conformer library, rather than
  computing it from NMA or MD vacuum simulations. The 615-frame library
  captures the *sampled* conformational manifold (biased by steering
  targets and starting state), while NMA captures *linear-response*
  modes around a reference.
- Expected agreement: canonical hinges (α-genu, β-knee) should be top
  candidates in our angular-variance list (AFK-MECH-3 output).
- Expected *informative* disagreement: the C-terminal coil regions
  dominate our RMSF, which classical NMA usually assigns moderate
  flexibility. This may reflect the HS-AFM surface interaction damping
  long-range modes but not local coil motion.
- Next concrete experiment: overlay our top-10 hinge candidates with
  Matsumoto's switch residues. Whichever residues appear in both are
  robust, multi-method-validated hinges.

## References

- [Xiong et al., Science 2001 (αVβ3 crystal structure)](https://pubmed.ncbi.nlm.nih.gov/11546839/)
- [Luo & Springer, Annu Rev Immunol 2007 (Integrin structure review)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3039929/)
- [Matsumoto et al., Biophys J 2008 (ENM NMA of αVβ3)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2527288/)
- [Chen et al., ACS Nano 2023 (Force-regulated bistability)](https://pubs.acs.org/doi/10.1021/acsnano.3c06253)
- [Li et al., Cell 2024 (Single-molecule activation cascade)](https://www.sciencedirect.com/science/article/pii/S0092867424004756)
- [Scheuring Lab HS-AFM](https://scheuringlab.com/research/)
- [MD-TASK DCC docs](https://md-task.readthedocs.io/en/latest/corr.html)
- [Koike et al., PLOS One 2014 (mDCC analysis)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112419)
- [Integrin Conformational Dynamics and Mechanotransduction (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9688440/)
