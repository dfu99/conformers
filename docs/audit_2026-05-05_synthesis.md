# Audit 2026-05-05 — synthesis (one-page distillation, v20 end of day)

For PI / Slack / paper-prep reference. Captures the full day's
audit deepening (19 passes, obj-038 through obj-059) in a
single page. Composite figure: `figures/audit-2026-05-05.png`.

---

## Headline numbers

| metric | morning | end of day | net |
|--------|---------|-----------|-----|
| objectives complete                    | 35  | **57** | **+22** |
| reviewer concerns addressed (of 23)    | 5   | **13** | **+8** (160 % more) |
| reviewer concerns partial              | 6   | 3      | −3 |
| reviewer concerns open                 | 12  | 7      | −5 |
| audit blockers ranked                  | 4   | 5      | (added Reviewer-E #5) |
| audit blockers closed                  | 0   | **3 of 5** | #2 dynamics, #4 contact, #5 cryptic |
| audit blockers external-gated remaining| 4   | 2      | #1 EO, #3 PACE allocation |
| pre-kickoff docs ready                 | 0   | 7      | strategy + risk + 4 kickoffs + 2 prep |
| route-A code assets                    | 0   | 2 + alignment JSON | remap_cvs.py + topology stub |

Trajectory of reviewer panel:
**5/6/12 morning → 8/4/10 v3 → 10/5/8 v5 → 11/5/7 v11 →
12/4/7 v16 → 13/3/7 v17 (stable v18, v19, v20)**

---

## Four analytical threads — closure status

### Thread 1 — Dynamics + populations + kinetics (CLOSED)
*obj-038, 043, 044, 045, 046, 047, 053, 054*

ΔG(CV0) from 1645 unbiased frames; populations 25/46/24/5 % for
BC/Inter/EC/EO\* (Boltzmann + empirical agree within 1 %). FES drift
σ ≈ 8 Å is **intrinsic** — no metadata signal explains it (obj-046,
0/10 hypothesis tests reach p < 0.05). Stepwise re-organization
confirmed: per-frame CUSUM rejects stationarity at p < 1e-3 in both
videos (obj-047). 3-state HMM populations match Boltzmann; mean
dwell times BC 14.6 s, Inter 6.9 s, EC 11.4 s. ACF e-folding time
τ_e = 5.1 s (V1) / 7.8 s (V2) triangulates obj-049 dwell estimates.
Synthesis main-text figure: `figures/dynamics_synthesis_v1.png`.

### Thread 2 — HMM dimensional triangulation (CLOSED)
*obj-048 (1-D), 050, 051, 052, 055 (2-D), 059 (3-D)*

From-scratch implementation across all dimensionalities (no
`hmmlearn`). All three HMMs converge:
- 1-D obj-048 (CV0 only): BC/Inter/EC partition along extension axis,
  3-state log-L = -5286, mean dwells confirm 10-second timescale.
- 2-D obj-055 (CV0+CV2): 4-state preferred by ΔBIC = -393. Primary
  clustering axis is **CV2**, not CV0; same partition recovered as 3-D.
- 3-D obj-059 (CV0+CV1+CV2): same partition; CV1 confirmed redundant
  (within-state corr +0.69 to +0.82, matches obj-054). 3-state preferred
  by ΔBIC = +78; ΔAIC = -8 marginal. Adding CV1 just adds 28 free params
  and noise.

Markovianity validated (obj-049 KS p > 0.13 all three states). V1=V2
in dwell-time distributions (obj-050, 4/4 KS perm tests p > 0.25;
4/4 mean-ratio CIs overlap 1). Model selection well-characterized
(obj-051 multi-init + 4-state ablation).

### Thread 3 — Cryptic-binding triangulation (CLOSED)
*obj-039 (RGD-MIDAS), 056 (SASA scan), 057 (LIGSITE), 058 (Vina-proxy)*

Reviewer E's "does extension expose new ligand-binding sites?"
answered with a clean negative result via four independent
geometric metrics:

| obj | metric on K417-K422 candidate | direction on extension |
|-----|------------------------------|------------------------|
| 056 | per-residue ΔSASA            | **+237 Å²** (gains exposure) |
| 057 | LIGSITE ΔV pocket volume     | **−3153 Å³** (loses pocket-like) |
| 058 | Vina-proxy ΔS                | **+0.12 NS** (druggability stays flat) |

| obj | metric on canonical MIDAS pocket | direction on extension |
|-----|----------------------------------|------------------------|
| 039 | MIDAS SASA                       | -35 % (more occluded) |
| 057 | MIDAS pocket volume              | +471 Å³ (slightly deeper) |
| 058 | MIDAS Vina-proxy score           | -1.61 (clashes 1 → 5) |

**Combined picture**: αVβ3 BC→EC extension does NOT create a new
ligand-binding pocket. The K417-K422 region (β3 hybrid/EGF1 hinge)
is conformationally responsive (consistent with Springer's
documented hybrid swing-out site) but the response is "open to
bulk solvent" not "form discrete druggable pocket." The canonical
RGD-MIDAS pocket gets MORE buried (SASA-wise) while keeping its
cavity intact (volume) but with a more occluded mouth — EC = the
"extended-closed" state.

### Thread 4 — EO state coverage (BLOCKED, 4× empirically confirmed)
*obj-025 (SMD), 041 (no published EO crystals), 055 (2-D HMM), 059 (3-D HMM)*

Four independent empirical negatives now converge:
- **obj-025 (SMD)**: k=1000 force constant only opens CV2 by 0.9 Å in 3 ns.
- **obj-041 (crystallography)**: 5/5 published full-ectodomain αVβ3
  crystal structures (1JV2, 1L5G, 4G1E, 4G1M, 4MMX) sit in the BC
  band (CV0 ≈ 51-52, CV2 ≈ 36); even cilengitide-bound 1L5G is bent.
- **obj-055 (2-D HMM)**: 4-state seed at CV2=55 collapsed to 41.6 in EM.
- **obj-059 (3-D HMM)**: 4-state seed at CV2=55 collapsed to 42.2,
  π=0.002 (degenerate, near-singular).

Bayesian floor on EO penalty: ΔG_EO ≥ 2.02 kcal/mol (obj-043,
Jeffreys 95 % upper bound). The EO blocker is **multiply-confirmed
empirical**. Recovery requires Route A (αIIbβ3 string method on
PACE A100-80GB) or Gō-Martini parallel track. Route-A starter
scripts are ready; PI sign-off Monday is the unlock.

---

## What unblocks Monday

PI sign-off on a 4-week PACE A100-80GB block (~80 GPU-hr × $10 ≈
$800) starting 2026-05-12.

This single decision unblocks:
- **Route A** (αIIbβ3 string-method port) — primary EO-coverage path.
  Day-1 starter code ready: `pipelines/route_a/scripts/{remap_cvs.py,
  build_av_topology.py}`. Risk register has 5 named failure modes
  (~40 % joint pass after gating).
- **Route E** (Switching Gō-Martini) — parallel cross-validation track.
  14-day plan, ~$200 GPU on A40. Acceptance criteria: monotonic
  CV0+CV2 trajectory, MIDAS SASA increases on extension (would
  reverse obj-039 negative).
- **AF2-Multimer ablation** — Reviewer B last open of 2 (50 GPU-hr).
- **Committor analysis** — Reviewer D last open of 1 (~80 GPU-hr).

All 4 share the same PACE allocation cycle.

---

## What's still open after today

| concern | status | why open |
|---------|--------|----------|
| #1 EO coverage              | external-compute-gated | needs Route A or Gō-Martini |
| #3 PACE A100 allocation     | external-compute-gated | needs PI sign-off |
| Reviewer A 2 of 5 items     | open                   | gated on EO sampling |
| Reviewer B AF2 ablation     | open                   | 50 GPU-hr |
| Reviewer D committor        | partial                | ~80 GPU-hr unbiased MD |
| Reviewer E 1 of 5 items     | partial                | BIOVIA / docking-venv item |

All open items are gated on either (a) external compute or (b) PI
authorization. None are gated on local analytical work.

---

## Confidence statement

The αVβ3 → HS-AFM fitting pipeline is **quantitatively defended
end-to-end** as of today, across all four analytical axes:

- **Calibration**: Tip dilation matches analytic √3·(r+R) (F2);
  V1/V2 real-vs-random baseline ΔSNR > 0.30 in both (F2).
- **Free energy + populations**: ΔG(CV0) has tight bootstrap CI
  (F3, < 0.5 kcal/mol across [50,85] Å); per-block ΔG drift is
  intrinsic, not metadata-driven (obj-046).
- **Kinetics**: 3-state HMM is Markovian (obj-049 KS p > 0.13);
  V1 = V2 in rates (obj-050); two independent timescale estimators
  agree on 5–15 s (obj-049 dwell + obj-054 ACF τ_e).
- **Mechanical sensitivity**: Per-residue RMSF reproducible
  V1 ↔ V2 r = 0.998; Hertzian δ < 1 nm at HS-AFM forces (obj-040).
- **Cryptic binding**: Triangulated via SASA + LIGSITE + Vina-proxy.
  No new druggable site created by extension; canonical RGD pocket
  gets MORE buried (EC = "extended-closed").

The remaining argument the manuscript has to make is "we cannot
reach EO without enhanced sampling" — obj-025 + obj-041 + obj-055
+ obj-059 form a **four-line empirical confirmation**. The
two open external-compute items (Route A + Gō-Martini) are
the two parallel paths that would close it.

---

## Cross-references

- `tasks/audit-2026-05-05.md` — full audit document (§1–§28, 3275 lines)
- `figures/audit-2026-05-05.png` — composite v20 status board (20 tiles, 2-row clustered)
- `tasks/objectives.yaml` — obj-038 through obj-059 (22 entries)
- `tasks/planning.md` — Recently Completed + Next Priority + Audit status
- `tasks/lessons.md` — 3+ new lessons today (route D null, V1 ↔ V2
  reproducibility, BLOSUM62 per-position lookup mandatory, LIGSITE
  pure-numpy implementation, Vina-proxy weighting)
- 25+ commits today on the `main` branch

---

_Generated 2026-05-05 21:00 (v20 end-of-day). The next deepening
pass priorities (when the autochain advances) are (a) PI sign-off
capture for Monday's PACE allocation request, (b) day-1 execution
of Route A once allocation lands, (c) optional: 4-state per-state
rate-matrix consolidation across 1-D / 2-D / 3-D HMMs._
