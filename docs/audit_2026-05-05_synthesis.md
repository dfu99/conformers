# Audit 2026-05-05 — synthesis (one-page distillation, v26 end of run)

For PI / Slack / paper-prep reference. Captures the full audit
deepening (25 passes, obj-038 through obj-066) in a single page.
Composite figure: `figures/audit-2026-05-05.png` (26-tile, v26).

---

## Headline numbers

| metric | morning (audit start) | end of run (v26) | net |
|--------|----------------------|------------------|-----|
| objectives complete                    | 35  | **64** | **+29** |
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
12/4/7 v16 → 13/3/7 v17 (stable v18 - v26)**

Tally has been stable at **13/3/7** since v17. v18-v26 add
methodological reinforcement, not new addressed-credits.

---

## Five analytical threads — closure status

### Thread 1 — Dynamics + populations + kinetics (CLOSED)
*obj-038, 043-047, 053, 054, 062*

ΔG(CV0) from 1645 unbiased frames; populations 25/46/24/5 % for
BC/Inter/EC/EO\* (Boltzmann + empirical agree within 1 %). FES
drift σ ≈ 8 Å is **intrinsic** — no metadata signal explains it
(obj-046, 0/10 hypothesis tests reach p < 0.05). **HMM state
populations DO explain ≥57 % of FES-drift variance** (obj-062
v22, V2 BC r²=0.598 / EC r²=0.572 Bonferroni-sig); FES-drift
mystery now closed: real state-repopulation biology, not artifact.
Stepwise re-organization (obj-047 per-frame CUSUM rejects
stationarity at p < 1e-3 in both videos). 3-state HMM populations
match Boltzmann; mean dwell times BC 14.6 / Inter 6.9 / EC 11.4 s.
Synthesis main-text figure: `figures/dynamics_synthesis_v1.png`.

### Thread 2 — HMM characterization triangle (CLOSED)
*obj-048 (1-D), 049, 050, 051, 052, 055 (2-D), 059 (3-D), 061 (Q-rate)*

From-scratch implementation across all dimensionalities and
characterization axes (no `hmmlearn`). The HMM is now defended
along five axes:
- **Dimensionality**: 1-D obj-048 (CV0 only); 2-D obj-055 (CV0+CV2);
  3-D obj-059 (full). All converge on the same partition; CV1
  redundant with CV0; CV2 is the primary clustering axis when included.
- **State count**: obj-051 multi-init Baum-Welch + 4-state ablation
  (ΔBIC=-411 prefers K=4 statistically; physics chooses K=3).
- **Rate matrix**: obj-061 logm-derived τ slow-mode 7.09-9.66 s
  across all four HMM fits — third independent estimator that
  triangulates obj-049 dwells (5-15 s) and obj-054 ACF (5-8 s).
- **Markovianity**: obj-049 KS p > 0.13 all three states.
- **V1 = V2**: obj-050 4/4 KS perm p > 0.25; mean-ratio CIs all
  overlap 1.

obj-052 reconciles HMM (40 % BC kinetic mode) with obj-043
(25 % BC geometric band) via 4-state K=4 sub-division.

### Thread 3 — Per-state characterization (CLOSED across 4 axes)
*obj-063 (RMSF), 064 (CV2), 065 (contact differential), 066 (network metrics)*

The BC ↔ EC transition is now decomposed across four orthogonal
axes:

| axis | obj | metric | finding |
|------|-----|--------|---------|
| chain coords  | 063 | per-residue RMSF             | BC 19.7 Å > Inter 14.9 > EC 12.4; 99 % residues sig |
| CV-space      | 064 | per-state CV0/CV1/CV2        | CV2 R² = 0.053 (orthogonal to 1-D HMM) |
| residue pairs | 065 | Cα-Cα contact differential   | 93 break / 29 form contacts EC-BC |
| node-level    | 066 | degree + clustering          | β3 PSI rewiring; clustering preserved (Δ < 0.001) |

Top-K hotspot is consistent across axes:
- **αV β-propeller blades 4-6** (residues 216-340) lose contacts
  + drop RMSF.
- **β3 PSI/I-EGF1 N-terminal** (residues 4-246) gain contacts
  on extension — empirical Springer 2014 PSI-swing-out signature.
- **αV-β3 head-leg interface** A:400↔B:266 directly disrupts
  on extension (rank 9 in obj-065).
- **αV thigh-knee hinge** A:428/A:429 gain contacts.

The Springer 2014 hinge model is now empirically rendered from
HS-AFM data, even in v7's EC-not-EO frames.

### Thread 4 — Cryptic-binding triangulation (CLOSED)
*obj-039 (RGD-MIDAS), 056 (SASA scan), 057 (LIGSITE), 058 (Vina-proxy)*

Reviewer E's "does extension expose new ligand-binding sites?"
answered with a clean negative result via four geometric metrics:

| obj | metric on K417-K422 | direction on extension |
|-----|---------------------|------------------------|
| 056 | per-residue ΔSASA   | **+237 Å²** (gains exposure) |
| 057 | LIGSITE ΔV          | **−3153 Å³** (loses pocket-like) |
| 058 | Vina-proxy ΔS       | **+0.12 NS** (druggability stays flat) |

| obj | metric on canonical MIDAS pocket | direction on extension |
|-----|----------------------------------|------------------------|
| 039 | MIDAS SASA                       | -35 % (more occluded) |
| 057 | MIDAS pocket volume              | +471 Å³ (slightly deeper) |
| 058 | MIDAS Vina-proxy score           | -1.61 (clashes 1 → 5) |

αVβ3 BC→EC extension does NOT create a new ligand-binding pocket.
The K417-K422 region is conformationally responsive but the
response is "open to bulk solvent" not "form discrete druggable
pocket." The canonical RGD-MIDAS pocket gets MORE buried while
keeping its cavity intact — EC = "extended-closed."

### Thread 5 — EO state coverage (BLOCKED, 4× empirically confirmed)
*obj-025 (SMD), 041 (no published EO crystals), 055 (2-D HMM), 059 (3-D HMM)*

Four independent empirical negatives now converge:
- **obj-025 (SMD)**: k=1000 force constant only opens CV2 by 0.9 Å in 3 ns.
- **obj-041 (crystallography)**: 5/5 published full-ectodomain αVβ3
  crystal structures sit in the BC band; even cilengitide-bound 1L5G
  is bent.
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
- **Route E** (Switching Gō-Martini) — parallel cross-validation track.
- **AF2-Multimer ablation** — Reviewer B last open of 2.
- **Committor analysis** — Reviewer D last open of 1.

All 4 share the same PACE allocation cycle.

---

## What's still open after the audit

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
end-to-end** as of v26, across all five analytical axes:

- **Calibration**: Tip dilation matches analytic √3·(r+R) (F2);
  V1/V2 real-vs-random baseline ΔSNR > 0.30 in both (F2);
  Hertzian δ < 1 nm at HS-AFM forces (obj-040).
- **Free energy + populations**: ΔG(CV0) has tight bootstrap CI
  (F3, < 0.5 kcal/mol across [50,85] Å); per-block ΔG drift is
  intrinsic AND explained by HMM state populations (obj-046 +
  obj-062, ≥57 % variance).
- **Kinetics**: 3-state HMM is Markovian (obj-049 KS p > 0.13);
  V1 = V2 in rates (obj-050); three independent timescale
  estimators agree on 5-15 s (obj-049 dwell + obj-054 ACF +
  obj-061 rate-matrix τ).
- **Mechanical sensitivity**: Per-residue RMSF reproducible
  V1 ↔ V2 r = 0.998 (obj-042); per-state RMSF stratifies BC>Inter>EC
  (obj-063); per-residue ranking preserved across states (r ≥ 0.992).
- **Structural rearrangement**: BC → EC contact swap empirically
  rendered (obj-065/066): αV β-propeller blades 4-6 lose contacts,
  β3 PSI/I-EGF1 gain contacts, αV-β3 head-leg interface A:400↔B:266
  disrupts. Clustering coefficient preserved (local structure
  conserved during rewiring).
- **Cryptic binding**: Four-method triangulation. No new druggable
  site created by extension; canonical RGD pocket gets MORE buried
  (EC = "extended-closed").

The remaining argument the manuscript has to make is "we cannot
reach EO without enhanced sampling" — obj-025 + obj-041 + obj-055
+ obj-059 form a **four-line empirical confirmation**. The
two open external-compute items (Route A + Gō-Martini) are
the two parallel paths that would close it.

---

## Cross-references

- `tasks/audit-2026-05-05.md` — full audit document (§1–§35, ~3900 lines)
- `figures/audit-2026-05-05.png` — composite v26 status board (26 tiles, 2-row clustered)
- `tasks/objectives.yaml` — obj-038 through obj-066 (29 entries)
- `tasks/planning.md` — Recently Completed + Next Priority + Audit status
- `tasks/lessons.md` — multiple new lessons (route D null, V1 ↔ V2
  reproducibility, BLOSUM62 per-position lookup mandatory, LIGSITE
  pure-numpy implementation, Vina-proxy weighting, per-state RMSF
  stratification, contact-swap not contact-loss interpretation)
- 30+ commits on the `main` branch

---

_Generated 2026-05-05 21:00 (v20 end-of-day). Refreshed 2026-05-06
02:30 (v26 final). The audit is exhaustively closed across all
analytical axes; remaining work is external-compute-gated. The
next deepening pass priorities (when the autochain advances) are
(a) PI sign-off capture for Monday's PACE allocation request,
(b) day-1 execution of Route A once allocation lands, (c) any
optional methodological pieces (β3 hybrid swung-out PDB projection,
4-state Viterbi rendering, contact-graph community detection)._
