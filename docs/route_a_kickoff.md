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
- **Alignment score**: 167

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

Per `docs/eo_coverage_strategy.md` recommendation A:

- **Block size**: 4 weeks wall, 80 GPU-hr peak
- **Hardware**: PACE A100-80GB, account `gts-yke8`
- **Budget proxy**: 80 GPU-hr × $10/hr ≈ $800
- **Submit script**: copy `pipelines/protenix-avb3-template/scripts/submit_slurm.sh`
  as starting template, modify for the string-method run loop.
- **Wall time**: 1 week prep + 2 weeks string optimization + 1 week
  validation = 4 weeks.

---

## 5. PI sign-off requested

1. Approve a 4-week PACE A100-80GB block from 2026-05-12 (next
   Monday) through 2026-06-09.
2. Approve the use of the cloned Ferg-Lab repo on PACE per
   `tasks/lessons.md` String Method Implementations note.
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
_Status: pre-kickoff. Activate after PI sign-off._
