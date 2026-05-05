# EO-coverage strategy decision (audit-2026-05-05 §10.9 P=1)

**Purpose.** Pick a single route to generate αVβ3 extended-open (EO)
conformers so the v7 conformer library covers the full bent→EO
manifold and the RGD-affinity question (obj-039 negative result) can
be answered with real docking. Needs PI sign-off.

**Status of the EO gap.** The 615-frame v7 library reaches CV0 ∈
[47.3, 85.0] Å but never crosses CV2 ≥ 50 Å. obj-025 confirmed SMD
with CV2 head-piece bias cannot drive EC→EO in 3 ns even at k=1000
(rate ≈0.07 Å/ps; ≈320 ns to reach the 60 Å EO target). obj-039
confirmed downstream consequence: the v7 "extended" frames are EC
(extended-closed), MIDAS pocket SASA *drops* 35 % on extension, real
RGD docking would not surface a bent-vs-extended affinity contrast.

The EO gap is therefore the single largest determinant of whether
the project is publishable at the target venue.

---

## Routes considered (5)

### A. αIIbβ3 string method via Ferg-Lab/`principalcurve_integrin_structures`

The Dasetty et al. 2025 bioRxiv paper solved bent→extended for αIIbβ3
using the finite-temperature string method on top of OpenMM + SEEKR2.
Their code repo is already cloned on PACE (per `tasks/lessons.md`
String Method Implementations).

| dimension | value |
|-----------|-------|
| compute time | 2-3 weeks production (string optimization) + 1 week analysis |
| GPU $$ | ~$800 PACE A100-80GB at 8 nodes × 100 hr |
| feasibility | medium — full code path runs αIIbβ3, must port α-subunit definitions to αV |
| pre-existing assets | `principalcurve_integrin_structures` cloned; SBMOpenMM also available |
| success probability | **~70 %** — method validated in published αIIbβ3 work on the same problem class |
| novel risk | force-field re-parameterization for αV β-propeller (αIIb-specific contacts in their setup); membrane embedding rebuild for αVβ3 lipid composition |
| reviewer A weight | strong — same lab as Dasetty 2025 |

### B. Funnel Metadynamics with CV2 (Plumed) on αVβ3

Funnel-MD on αIIbβ3 has been demonstrated (PMC 2024) and would be a
direct adaptation. Funnel restraints reduce sampling volume, making
EC→EO transition cheaper than vanilla well-tempered MD.

| dimension | value |
|-----------|-------|
| compute time | 4 weeks: install GROMACS+Plumed (1w), parameter tuning on small system (1w), production (2w) |
| GPU $$ | ~$1200 PACE A100-80GB (50-100 ns × 6 funnel replicas) |
| feasibility | low-medium — currently use OpenMM exclusively. New GROMACS install collides with project venv discipline (CLAUDE.md §4); needs separate funnel-md env |
| pre-existing assets | none — fresh build |
| success probability | **~50 %** — convergence of well-tempered MD on a 1654-residue ectodomain in funnel geometry is non-trivial |
| novel risk | Plumed bias deposit may not converge; well-tempered Δ-T tuning is system-specific |
| reviewer D weight | strongest — direct ΔG path along CV2 |

### C. Replica exchange (T-REMD) on αVβ3

8-replica REMD spanning 300-450 K on the head-piece subdomain,
2-3 ns/replica/exchange. Sampling EC↔EO via thermal driving.

| dimension | value |
|-----------|-------|
| compute time | 1 week setup + 2 weeks production |
| GPU $$ | ~$1500-2000 (8 replicas × 50 ns × A100) |
| feasibility | medium — OpenMM-REMDTorch or similar exists; tempering range 300-450 K may unfold the β-propeller |
| pre-existing assets | OpenMM env already on RunPod; tempering scripts in `pipelines/avb3-conformers/scripts/` not REMD-aware |
| success probability | **~40 %** — REMD on a 1654-residue ectodomain is not a clean fit; usually applied to small domains |
| novel risk | tempering-induced unfolding; head-domain hybrid swing may need orthogonal CV |
| reviewer D weight | medium |

### D. Use published EO-state PDBs as templates (1L5G + variants)

**EMPIRICALLY INVALIDATED 2026-05-05 evening (obj-041).** Originally
proposed: αVβ3-RGD complex 1L5G (Xiong 2002, *Science*) is "extended-
open with cilengitide bound." Implemented `library_coverage_v3.py` to
test this: downloaded all 5 published full-ectodomain αVβ3 PDBs
(1JV2, 1L5G, 4G1E, 4G1M, 4MMX), computed CV0/CV1/CV2 with our domain
definitions. Result: **all 5 sit in the BC band (CV0 ≈ 51-52 Å,
CV2 ≈ 36 Å)**. Cilengitide-bound 1L5G has an *open headpiece*
internally but the legs are still folded over — full-ectodomain
αVβ3 is only ever crystallized bent.

| dimension | value (post-test) |
|-----------|-------|
| compute time | 1 day  ✓ already executed |
| GPU $$ | $0  ✓ |
| feasibility | high  ✓ pipeline runs cleanly |
| pre-existing assets | scoring script lands as `pipelines/conformer-library/scripts/library_coverage_v3.py` |
| success probability | **revised: 0 % for full-ectodomain EO** (5/5 PDBs are bent) |
| novel risk | none — empirical null result |
| reviewer A weight | weak — confirms the EO crystal gap rather than filling it |
| limitation | route D **does NOT produce EO endpoints** for αVβ3. Confirms enhanced sampling (A/B/C/E) is necessary, not optional |

This is a publishable null result and strengthens routes A/B/C/E.
Recommendation revised below.

### E. Switching Gō-Martini coarse-grained (JCTC 2024)

CG model that switches between bent and EO Gō contacts. Substantially
faster than atomistic, captures large-scale conformational transition.

| dimension | value |
|-----------|-------|
| compute time | 2 weeks: install GōMartini + parameterize switching contacts + production |
| GPU $$ | ~$200 (CG runs ~50× faster than atomistic) |
| feasibility | medium — needs Martinize3 + GōMartini build; we have OpenMM only |
| pre-existing assets | none |
| success probability | **~60 %** — published method, but integrin-scale switching not yet done |
| novel risk | back-mapping CG → atomistic for downstream pseudo-AFM |
| reviewer-set weight | medium — CG is generally accepted but Reviewer A may flag atomistic-only |

---

## Recommendation (revised after obj-041 null result)

Original recommendation was D → A. After empirically invalidating D,
the choice narrows. Revised:

**A as the primary route, with E as a parallel backup.** Concretely:

1. **Skipped (route D was tested and FAILED)**: 5/5 published αVβ3
   ectodomain crystal structures sit in the BC band. There is no
   `download an EO PDB` shortcut for αVβ3 at the full-ectodomain
   resolution. obj-041 documents this; route D is removed from the
   active path.

2. **Primary, immediate ramp (4 weeks setup, ~$800)**: port Ferg-Lab's
   αIIbβ3 string method to αVβ3. The repo is cloned on PACE; the
   primary remaining work is residue-mapping the α-subunit
   collective-variable definitions and confirming the membrane
   embedding works for αV β-propeller (which has different
   carbohydrate decoration than αIIb). Validates obj-029's
   first-principles αIIbβ3-vs-αVβ3 prediction by going in the
   reverse direction: structures-derived for αIIbβ3, MD-derived
   for αVβ3.

3. **Parallel backup (2 weeks, ~$200)**: build a Switching
   Gō-Martini coarse-grained model on αVβ3 using existing
   ectodomain bent state (1JV2) as Gō-A and a homology-mapped
   αIIbβ3-extended structure as Gō-B. Faster sampling but lower
   atomistic detail; serves as an independent line of evidence
   if route A hits a force-field showstopper.

**Why not B (metadynamics) or C (REMD) first.** Both require
significant infrastructure builds (GROMACS+Plumed for B; tempered
REMD wrappers for C) that compete with project velocity. The
existing Ferg-Lab path has working code, validated on the
biologically-related αIIbβ3, and PI sign-off can fund a single
PACE A100-80GB block to reproduce + adapt their result. Routes B
and C remain tertiary backups if both A and E hit a showstopper.

**Why D was eliminated.** obj-041's library_coverage_v3 figure
shows that all 5 published αVβ3 ectodomain crystal structures
sit in the BC band — there is no shortcut from RCSB to an EO
endpoint at full-ectodomain resolution. Even cilengitide-bound
1L5G (open headpiece, 2002 Xiong) crystallizes bent overall.

---

## Decision request to PI (revised after obj-041)

| | Route | Time | $$ | Risk | Closes |
|---|-------|------|----|----|------|
| ✅ recommended | **A as primary + E as parallel backup** | 6-8 wk + 2 wk parallel | ~$1000 total | medium | Reviewer A geometric, D ΔG, E RGD docking |
| if A stalls | A → B (metadynamics) | +4 wk | +$1200 | medium-high | same |
| budget-constrained | E (Gō-Martini) alone | 2 wk | $200 | medium | partial — atomistic still needed downstream |
| budget-unlimited | parallel A + B + E | 6-8 wk | $2200 | low | as above + cross-validation |
| ❌ tested + failed | D (download EO PDBs) | 1 day | $0 | n/a | n/a — empirical null (obj-041) |

PI sign-off needed on: (1) approve A: 4-week PACE A100-80GB block
(~80 GPU-hr × $10 = $800) starting next week; (2) optionally approve
parallel E (Gō-Martini, ~$200 in CG MD compute) as cross-validation
if R&D-budget allows.

---

## Cross-references

- Audit §3 (EO blocker statement)
- Audit §10.5 (obj-039 RGD PoC negative result)
- obj-025 (SMD cannot reach EO at k=1000)
- obj-029 (1JV2 + 3FCS + 5 other PDBs first-principles analysis;
  D is a direct extension)
- `tasks/lessons.md` String Method Implementations note about
  Ferg-Lab repo on PACE
- `tasks/lessons.md` Targeted MD / Metadynamics for Integrins
  section (PMC 2024 funnel-MD on αIIbβ3, ACS Nano 2024 SMD
  intermediate states on αVβ3)
- Dasetty et al. *bioRxiv* 2025 "finite-temperature string method
  bent→extended for αIIbβ3"

---

_Author: AFK audit deepening pass, 2026-05-05 evening._
_Sign-off requested from PI by 2026-05-08._
