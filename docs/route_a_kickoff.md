# Route-A kickoff — αIIbβ3 string-method port to αVβ3

**Audit follow-up §11.5 P=1.** Pre-kickoff issues list for porting the
Ferg-Lab `principalcurve_integrin_structures` repo (Dasetty et al.
*bioRxiv* 2025) from αIIbβ3 to αVβ3.

Companion documents:
- `docs/eo_coverage_strategy.md` — recommendation rationale
- `docs/route_a_risk_register.md` — 5 named failure modes with gates

This doc is the actionable issues list to start day 1 once PI approves.

---

## 1. αV vs αIIb α-subunit residue mapping

Pairwise BLOSUM62 alignment of αV chain A from **1JV2** (927 CAs,
resSeq 1-956) against αIIb chain A from **3FCS** (913 CAs,
resSeq 1-959). Global alignment, gap_open=−10, gap_extend=−1.

- **Alignment length**: 949 columns
- **Sequence identity**: **38.6 %** — moderate; well within homology-
  port range
- **Alignment score**: 1507 (BLOSUM62, gap_open=−10, gap_extend=−1)

Saved as `results/route_a/av_aiib_alignment.json` for downstream
CV-remapping scripts.

### Landmark residue table (αV → αIIb)

| αV resSeq | αV aa | αIIb resSeq | αIIb aa | Δ | landmark |
|-----------|-------|-------------|---------|----|----------|
| 1 | F | 1 | L | +0 | β-propeller blade 1 (W1) |
| 62 | S | 60 | A | −2 | β-propeller blade 2 (W2) |
| 124 | P | 126 | P | +2 | β-propeller blade 3 (W3) |
| 180 | Q | 192 | L | +12 | β-propeller blade 4 (W4) |
| 218 | D | 231 | F | +13 | RGD-Arg pocket residue |
| 243 | V | 257 | A | +14 | β-propeller blade 5 (W5) |
| 303 | R | 316 | S | +13 | β-propeller blade 6 (W6) |
| 367 | E | 383 | P | +16 | β-propeller blade 7 (W7) |
| 435 | Y | 448 | Y | +13 | thigh end |
| 442 | T | 455 | K | +13 | genu hinge |
| 505 | K | 516 | R | +11 | calf-1 mid |
| 605 | K | 611 | Q | +6 | calf-2 start |
| 702 | H | 708 | G | +6 | calf-2 mid |
| 740 | A | 746 | A | +6 | calf-2 end |
| 820 | M | 824 | L | +4 | membrane-proximal mid |
| 940 | D | 944 | L | +4 | C-terminal |

### Offset structure summary

The offset varies through the structure:

- **N-terminal (residues 1-65)**: Δ ≈ −2 to 0  (propeller blades W1-W2)
- **Mid-propeller (residues 65-435)**: Δ stabilizes at **+12 to +16**
  — αIIb has a 13-residue insertion in W4 / between W4-W5
- **Genu + calf-1 (residues 442-605)**: Δ ≈ +11 to +13
- **Calf-2 (residues 605-740)**: Δ contracts to +6
- **Membrane-proximal (740-940)**: Δ further contracts to +4

**Practical implication**: a single global Δ offset will NOT work.
The CV-remap script must use the per-position lookup from
`av_aiib_alignment.json`.

### Critical residue substitutions to flag

The RGD-Arg pocket residue **αV D218 → αIIb F231 (D vs F)** is a
non-conservative substitution. αIIb uses a different residue
(literature: αIIb D224) for RGD-Arg coordination than αV. This
means the αIIbβ3 string method's CV definitions for ligand-binding
contacts cannot be ported verbatim — they must be replaced by the
αV equivalents.

The MIDAS β-subunit residues (β3 D119, S121, S123, etc.) are *the
same* in both heterodimers because both share the β3 chain. β-chain
CV definitions can be ported as-is.

---

## 2. Force-field gotchas

### 2.1 Force-field family choice

| | Ferg-Lab default | Conformers project default |
|---|------------------|-----------------------------|
| protein | Amber99SB-ILDN | Amber14SB |
| lipid | Slipids | lipid17 (CHARMM-GUI default) |
| water | TIP3P | TIP3P |
| ions | Joung-Cheatham | Joung-Cheatham |

**Decision**: match Ferg-Lab exactly for the port (Amber99SB-ILDN +
Slipids). Sensitivity check after string convergence by running a
50 ns reference trajectory under Amber14SB + lipid17 from the EO
endpoint and verifying structural stability.

### 2.2 Membrane composition

| | Ferg-Lab (αIIbβ3) | recommended (αVβ3) |
|---|-------------------|--------------------|
| host cell | platelet | endothelial |
| POPC | 30 % | 50 % |
| POPE | 25 % | 15 % |
| POPS | 15 % | 5 % |
| cholesterol | 25 % | 25 % |
| PIP2 | 5 % | 0 % |

Endothelial-cell composition matches the cells αVβ3 actually lives
on. PIP2 is platelet-specific.

### 2.3 Glycosylation

αV has 6 N-glycosylation sites (UniProt: A:74, A:178, A:307, A:401,
A:642, A:801). αIIb has 7 (different positions). Ferg-Lab omits
glycans. Day-1 plan: **omit glycans**, match Ferg-Lab. If string
result is unphysical, add glycans via CHARMM-GUI Glycan Reader
(documented in risk register §4).

### 2.4 Calcium / Mg²⁺ handling

αV β-propeller has 4 calcium-binding loops. αIIb has 4 + 1 (W3
extra). Ferg-Lab parameterizes Ca²⁺ as bonded to coordinating
residues (Aqvist parameters). Port must:

1. Identify the 4 αV Ca²⁺-binding loops (residues 220-225, 281-286,
   348-353, 412-417 per Xiong 2001).
2. Drop the αIIb extra W3 site.
3. Use Aqvist parameters as Ferg-Lab did.

### 2.5 β3 MIDAS metal

Both heterodimers share β3, so the MIDAS Mn²⁺ / Mg²⁺ parameterization
is **identical**. No port work needed. Use Ferg-Lab's β3 setup
verbatim.

### 2.6 N-terminus of αIIb is heavy / light chain split

αIIb is post-translationally cleaved at αIIb K877 into a heavy
chain (1-859) + light chain (878-1008). The disulfide between
C826-C880 holds them together. Ferg-Lab models this as a single
peptide with a special-case bond. αV is **not** cleaved — the
α-subunit is one continuous polypeptide. Port must:

1. Remove the αIIb heavy/light bond from the topology.
2. Make αV chain A fully continuous.

### 2.7 PSI-domain disulfide in β3 (shared)

β3 has the αIIb-shared β3 PSI-domain disulfide network. No port
work for the β-chain.

---

## 3. Day 1-2 deliverables to verify port readiness

Before any MD launch:

1. ✓ Done: `results/route_a/av_aiib_alignment.json` — saved.
2. **TODO**: write `pipelines/route_a/scripts/remap_cvs.py` — read
   Ferg-Lab CV definitions in αIIb numbering and emit αV-numbered
   versions using the alignment table.
3. **TODO**: write `pipelines/route_a/scripts/build_av_topology.py`
   — apply force-field gotchas §2.1, §2.4, §2.6 to produce a
   ready-to-equilibrate αVβ3 + membrane topology.
4. **TODO**: 50 ns equilibration on PACE A100-80GB to verify the
   topology + membrane build are stable (Risk-register Risk 2 gate).

---

## 4. Compute budget and PACE allocation

**What this is — and is not.** This block is *not* another structure
predictor (Protenix/Boltz/AFCluster/AF2 all collapse to bent; that is
settled — see obj-041, obj-055). It is **enhanced-sampling molecular
dynamics** — the finite-temperature string method — which physically
drives αVβ3 along the bent→extended-open (EO) path and yields a
free-energy profile. The EO state is absent from every source we have:
the HS-AFM movies top out short of it (obj-055), no EO crystal exists
(obj-041), and plain steered MD is ~100× too slow within budget
(obj-025). The string method is the one remaining route, and it is the
only one that produces the ΔG curve reviewers keep demanding.

**GPU sizing — A100-80GB is NOT required (corrected 2026-06-27).** The
80 GB figure was inherited from our structure-prediction jobs
(Protenix/Boltz/AF2), which need 40–80 GB to hold attention over the
full sequence. This is **MD**, a much lighter memory regime: the system
is αVβ3 + endothelial bilayer + explicit solvent ≈ ~1 M atoms, and MD
holds only positions, forces, and neighbor lists — **comfortably under
16 GB VRAM**. Cost here is driven by *throughput* (ns/day), not VRAM.
The finite-temperature string method also parallelizes as ~8
independent images (1 image/GPU), so this is "several modest cards,"
never one large one. Per-GPU need ≈ 8 GB; 16 GB gives headroom.
- **Stage 1 card**: `V100-16GB` (PACE cheapest; `-C V100-16GB`) — the
  setup-and-check gate is a single short equilibration, not production.
- **Production card (Stages 2–3)**: `A100-40GB` or `L40S` — full MD
  throughput at lower cost than A100-80GB. Lock the final choice from a
  real ns/day benchmark measured during Stage 1, not a guess.
- The $ figures below are rough proxies at the old A100-80GB rate;
  dropping the card tier either lowers the bill or buys more sampling
  for the same spend.

**Staged budget with early-stop gates (cap the downside).** Do NOT
commit the full ~$800 up front. The run is gated (see
`docs/route_a_risk_register.md`); each gate is a cheap decision point:

| Stage | Work | GPU-hr | ~$ | Card | Decision gate |
|------|------|-------:|---:|------|---------------|
| 1 | Topology build + remapped-CV check + membrane equilibration | ~20 | ~$200 | V100-16GB | Does the ported setup reproduce obj-029 geometry at 300 K? |
| 2 | String optimization (production) | ~50 | ~$500 | A100-40GB / L40S | String converged in < 50 GPU-hr? |
| 3 | Validation + ΔG / committor analysis | ~10 | ~$100 | A100-40GB / L40S | EO endpoint reaches CV0 ≥ 80, CV2 ≥ 50? |

- **Account**: `gts-yke8`. Cheaper cards reduce per-hour cost; the
  staged $ caps above are upper bounds at the old rate.
- **Recommended approval**: Stage 1 only (~$200, ~1 week). If it clears,
  continue to Stages 2–3. If Stage 1 fails, stop — $200 spent, not
  $800 — and publish around the bent→EC transition already fully
  characterized (reviewer tally 13/3/7 of 23), documenting EO as a hard
  limit. Real downside is capped at ~$200.
- **Success odds**: ~70 % the method works in principle (validated in
  published αIIbβ3 work); ~40 % joint pass after all porting risks
  (`docs/route_a_risk_register.md`). The gates exist precisely because
  the joint odds are a coin-flip — they convert a $800 gamble into a
  $200 test.
- **Submit script**: copy `pipelines/protenix-avb3-template/scripts/submit_slurm.sh`
  as starting template, modify for the string-method run loop.

**Compute target: RunPod is the cheaper option (added 2026-06-28).** PI
flagged an already-running RunPod pod. Read-only inspection
(`ssh -p 11363 root@213.173.109.6`): **RTX 2000 Ada 16 GB**, 48 CPU,
251 GB RAM, persistent **2.3 PB `/workspace`** network volume (ideal for
multi-week runs that must survive pod restarts). 16 GB fits the MD
comfortably. *Constraint*: the GPU is currently at 99 % utilization
(memory free, compute saturated) by the caDNAgentic oxDNA job — do NOT
contend with it; wait for a free GPU window or spin a second pod.
- **Stage 1 on RunPod ≈ $0**: the pod is already paid for; the
  build-and-check gate is short. Replaces the ~$200 PACE Stage-1.
- **Production on RunPod**: rent on the same `/workspace` volume. RunPod
  $/hr ≪ the ~$10/hr PACE A100-80GB proxy → real production likely a few
  hundred $, not $800. (A40 already in-budget per
  `docs/go_martini_kickoff.md`.) **Card choice — for MD, bandwidth >
  VRAM > FLOPS** (all candidates have ample VRAM; we need only a few GB).
  Budget-card MD ranking (PI-suggested cheaper options, 2026-06-28):
  - `RTX A5000` (24 GB, 768 GB/s) — fastest of the budget set.
  - `RTX A4500` (20 GB, 640 GB/s) — close second.
  - `RTX A4000` (16 GB, 448 GB/s) — value pick; cheapest, exactly fits.
  - `RTX 4000 Ada` (20 GB, ~360 GB/s) — Ada-efficient but bandwidth-limited.
  - `L4` (24 GB, ~300 GB/s, 72 W) — **avoid for MD**; inference card,
    bandwidth-starved despite newest gen + most VRAM.
  - Premium tier: A40 (48 GB, 696 GB/s) / L40S (48 GB, 864 GB/s).
    **Key insight — the A5000 (768 GB/s) out-bandwidths the A40 and ≈
    ties it for MD**; the A40's only edge is 48 GB VRAM we don't use.
    Rough MD throughput (A4000 = 1.0): A4000 1.0 / A4500 ~1.3 / A5000
    ~1.6 / A40 ~1.7 / L40S ~2.3. Premium gap is only ~1.5–2.3×, not
    5–10×; budget cards usually match/beat premium on $/ns. A40/L40S buy
    wall-clock turnaround, not value. **A5000 = sweet spot** (premium-
    class speed at budget price); pay up for L40S only if turnaround
    matters.
  - **Decide by $/ns, not spec sheet**: Stage 1 benchmarks ns/day across
    tiers (A4000 + A4500 + A5000 + L40S) on the same `/workspace` volume
    (RunPod card-swap = minutes); production tier locked from real
    nanoseconds-per-dollar before committing to the weeks-long Stage 2.
- **Setup**: isolated env on `/workspace` (OpenMM + Ferg-Lab string
  stack + force fields), ~1 hr one-time. `nvcc`/conda not in base PATH —
  use a conda env that ships the CUDA runtime (OpenMM does).
- PACE (`gts-yke8`, A100-40GB/L40S) remains the fallback if RunPod
  capacity/stability is insufficient for the weeks-long Stage-2 run.

---

## 5. PI sign-off requested

1. Approve **Stage 1 only** (~$0 on the existing RunPod pod, or ~$200 on
   PACE): topology build + remapped-CV reproduction check + membrane
   equilibration. Stages 2–3 gated on Stage 1 passing its decision check.
   *Preferred: RunPod* (existing pod for Stage 1; rent A40/L40S for
   production on the same persistent volume — cheaper than PACE).
2. Approve the use of the cloned Ferg-Lab string-method repo.
3. Optionally approve the parallel Switching Gō-Martini track
   (~$200) per audit §11.5 P=2 — independent cross-validation.

---

## 6. Cross-references

- `docs/eo_coverage_strategy.md` — strategy decision
- `docs/route_a_risk_register.md` — 5 risks + gates
- `tasks/audit-2026-05-05.md` §11.1-§11.5
- obj-029 — first-principles αV vs αIIb (validates Risk 1 mitigation)
- obj-041 — empirical confirmation that route D is invalidated
- `results/route_a/av_aiib_alignment.json` — alignment landmarks

---

_Author: AFK audit deepening pass v3, 2026-05-05 late evening._
_Revised 2026-06-27: budget restructured into staged early-stop gates
(Stage 1 ~$200 decision point) after PI questioned the flat $800 ask and
the "is it another model?" framing. This is enhanced-sampling MD, not a
predictor._
_Revised 2026-06-27 (2): GPU spec corrected after PI asked "do you really
need A100-80GB?" — no. MD needs < 16 GB; Stage 1 on V100-16GB, production
on A100-40GB / L40S. 80 GB was inherited from the prediction jobs._
_Revised 2026-06-28: RunPod added as preferred (cheaper) compute target.
Existing pod = RTX 2000 Ada 16 GB + 2.3 PB persistent volume; Stage 1 ≈ $0
once a GPU slot frees (currently 99 % busy with the oxDNA job). PACE is
the fallback._
_Status: pre-kickoff. Awaiting PI go-ahead to set up Stage 1 on RunPod._
