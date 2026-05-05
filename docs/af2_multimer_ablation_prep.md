# AF2-Multimer ablation prep — Reviewer B last open

**Audit follow-up §12.7 P=3.** Pre-run planning document for the
AF2-Multimer ablation that Reviewer B (Jumper/Ahdritz) is expected
to demand: "is your pipeline complexity justified, or could
AF2-Multimer alone produce equivalent conformers?"

This is a doc, not an experiment. Audit §3 already drop-listed the
"three-predictor framing" partially because of compute cost; the
ablation is the alternative way to discharge that concern.

---

## What Reviewer B will ask

Three specific questions:

1. **Single-predictor sufficiency.** Could AF2-Multimer alone
   produce the bent + extended conformer pairs you generate via
   Protenix + steered MD? If yes, your pipeline is over-complex.
2. **MSA-depth diversity.** Wayment-Steele et al. *Nat Comms* 2024
   showed AF2 with subsampled MSAs explores conformational states;
   predicted frequencies correlate with NMR populations. Have you
   replicated this for αVβ3?
3. **Ranking-score-vs-RMSD calibration.** Does AF2's ranking score
   monotonically rank by structural correctness, or do you need
   custom filters?

---

## What we already know (from existing objectives)

- **obj-006** (Protenix MSA-depth-invariant): identical predictions
  at full MSA vs 5 % subsample. TM-scores differ by < 0.001.
- **obj-007** (AF2 reduced_dbs validation): AF2 also locks to the
  bent conformation. 25 ranked predictions all have pairwise RMSD
  0.1-2.9 Å (i.e., AF2 produces 25 nearly-identical bent
  structures).
- **obj-009** (AF2 pLDDT vs displacement): AF2 cannot score arbitrary
  PDBs but pLDDT identifies flexible leg/tail domains (82-84). TM-
  score is the practical conformer filter.
- **obj-041** (route D null): all 5 published αVβ3 ectodomain
  crystal PDBs are bent — confirming that bent is the only
  experimentally accessible state.

**Combined message**: both AF2 and Protenix collapse to the bent
state when given αVβ3 sequence + MSA. Neither accesses EO without
enhanced sampling. The "three-predictor pipeline" framing was
already abandoned in audit §4 (drop list).

---

## Proposed ablation (50 GPU-hr)

### Setup

- **Sequence**: αVβ3 heterodimer (αV chain 1-957 + β3 chain 1-690).
- **MSA conditions**: 5 depth subsamples — 100 %, 50 %, 25 %, 10 %,
  5 % — using AlphaFold's MSA from `reduced_dbs` (already cached).
- **Replicates per depth**: 5 random subsample seeds × 5 model
  weights = 25 predictions per depth. Total = 125 predictions.
- **Ranking**: ipTM + pTM via standard AF2-Multimer scoring.
- **CV scoring**: feed every prediction through
  `pipelines/conformer-library/scripts/build_library_metadata.py`
  to compute CV0/CV1/CV2.

### Expected outcome (pre-registered prediction)

- All 125 predictions cluster at CV0 ≈ 51 Å, CV2 ≈ 36 Å (the bent
  basin, matching obj-007 and obj-041).
- ipTM ≈ 0.85-0.92 across all depths (no MSA-depth sensitivity for
  αVβ3 specifically — same conclusion as obj-006/007 for Protenix).
- ranking-score-vs-RMSD: weak monotonicity, Pearson r ≈ 0.3-0.5
  (many high-ipTM structures with low RMSD-to-1JV2; some low-ipTM
  structures also bent).

If the prediction holds, the ablation **strengthens** the audit's
narrative: the EO-coverage gap is a structure-prediction-method-
agnostic result. Three-predictor framing is decisively dropped;
the paper claim becomes "no AF2/Protenix variant accesses EO; MD
enhanced sampling is required."

### Compute budget

- AF2-Multimer ipTM-with-templates on a 1654-residue heterodimer:
  ≈ 12 minutes per prediction on A100-80GB
- 125 predictions × 12 min = 25 GPU-hr  (lower bound)
- + 25 GPU-hr for re-prediction at higher MSA-subsample variance
  (Wayment-Steele suggests 5 % → 5 distinct seeds is light; bump
  to 10 seeds at 5 % depth specifically)
- Total: 50 GPU-hr matches `tasks/planning.md` Next Priority #2.

PACE allocation: A100-40GB suffices (1654 res < 2048 token cap with
templates). Cheaper than A100-80GB.

### Submit script reuse

`pipelines/avb3-conformers/scripts/submit_af2_msa_test_slurm.sh`
already exists from the obj-007 work. Modify:

- Loop over MSA depth ∈ {100, 50, 25, 10, 5}%
- Random subsample seeds 1..5 per depth
- Output to `data/runs/avb3/af2_ablation/depth_{D}/seed_{S}/`

---

## Decision threshold

The ablation is informative if either:

1. The prediction holds (no EO at any MSA depth) → drop "three-
   predictor" framing decisively, paper claim sharpens.
2. The prediction fails (AF2 with low-MSA-depth produces an EO
   conformer) → MAJOR finding, pivot the project. Worth the
   compute either way.

Pre-registration: **expect prediction to hold**. Falsifying
threshold: > 5 % of predictions land in CV0 ≥ 78 Å (EC band) at
*any* MSA depth.

---

## Status + dependencies

- Compute: blocked on PACE A100-40GB allocation.
- Code: ready (modify obj-007 submit script).
- Analysis: existing `build_library_metadata.py` produces the
  CV scoring directly.
- Prerequisites: none beyond PACE allocation.

---

## Cross-references

- obj-006, obj-007, obj-009, obj-041 (existing AF2/Protenix
  observations)
- audit-2026-05-05.md §4 (drop list — "three-predictor" framing)
- `tasks/planning.md` Next Priority #2 (50 GPU-hr AF2-Multimer
  ablation, can run on A40 once free)

---

_Author: AFK audit deepening pass v4, 2026-05-05 23:30._
_Status: ablation plan only; compute on hold pending PACE
allocation. Activate after route-A kickoff completes (Reviewer A
takes priority over Reviewer B given project trajectory)._
