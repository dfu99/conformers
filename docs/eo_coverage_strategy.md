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

αVβ3-RGD complex 1L5G (Xiong 2002, *Science*) is in an extended-open
configuration with cilengitide bound. Combined with 1JV2 (bent) we
have an experimental EO endpoint without any new MD.

| dimension | value |
|-----------|-------|
| compute time | 1 day: download, strip ligand, compute CVs, fit into pipeline |
| GPU $$ | $0 |
| feasibility | high — uses existing pipeline (`build_library_metadata.py`) |
| pre-existing assets | `multi_integrin_first_principles.py` already does this for 7 PDBs |
| success probability | **>95 %** — structures are crystallographically determined |
| novel risk | only one EO state, not a trajectory — cannot give ΔG along EC→EO |
| reviewer A weight | strong — uses canonical experimental reference |
| limitation | does not produce a *trajectory*; cannot answer reviewer D's "ΔG along EC→EO" |

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

## Recommendation (with rationale)

**D + A in sequence.** Concretely:

1. **Immediate (next week, $0)**: extend the 7-PDB
   `multi_integrin_first_principles.py` analysis to include 1L5G
   (αVβ3-cilengitide, EO endpoint), 2VDR (αVβ3-fibronectin), and any
   other RCSB αVβ3-RGD complex. Compute CV0/CV1/CV2 for each.
   Fold these as **endpoint EO templates** into the v7 library;
   regenerate library_coverage_v3.png with the EO panel populated.
   This closes Reviewer A's geometric coverage concern within a day,
   $0 compute. Caveat — these are static endpoints with bound
   cilengitide; we lose the trajectory needed for ΔG.

2. **Medium term (next 6-8 weeks, ~$1000)**: port Ferg-Lab's
   αIIbβ3 string method to αVβ3. The repo is cloned on PACE; the
   primary remaining work is residue-mapping the α-subunit
   collective-variable definitions and confirming the membrane
   embedding works for αV β-propeller (which has different
   carbohydrate decoration than αIIb). Validates obj-029's
   first-principles αIIbβ3-vs-αVβ3 prediction by going in the
   reverse direction: structures-derived for αIIbβ3, MD-derived
   for αVβ3, and comparing.

**Why not B (metadynamics) or C (REMD) first.** Both require
significant infrastructure builds (GROMACS+Plumed for B; tempered
REMD wrappers for C) that compete with project velocity. The
existing Ferg-Lab path has working code, validated on the
biologically-related αIIbβ3, and PI sign-off can fund a single
PACE A100-80GB block to reproduce + adapt their result. Routes B
and C remain backups if A hits a force-field showstopper.

**Why D alone is not enough.** 1L5G gives one EO endpoint structure,
not the EC→EO trajectory. The ΔG bootstrap (F3) confirmed that the
1645 unbiased fitted frames give a tight free-energy estimate for
the *populated* range — but ΔG along EC→EO requires sampling that
range, which only enhanced sampling (A/B/C/E) can provide.

---

## Decision request to PI

| | Route | Time | $$ | Risk | Closes |
|---|-------|------|----|----|------|
| ✅ recommended | **D + A** in sequence | 1d + 6-8 wk | ~$1000 total | medium | Reviewer A geometric, Reviewer D ΔG, Reviewer E RGD docking |
| backup if A stalls | A → B | +4 wk | +$1200 | medium-high | same |
| fastest feasible-only | D alone | 1d | $0 | low | Reviewer A geometric only |
| if budget unconstrained | parallel A + B | 6-8 wk | $2000 | low | as above + ΔG cross-validation |

PI sign-off needed on: (1) approve D immediately for the next
working day; (2) approve A start once D ships, with a 4-week PACE
A100-80GB block (~80 GPU-hr × $10 = $800).

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
