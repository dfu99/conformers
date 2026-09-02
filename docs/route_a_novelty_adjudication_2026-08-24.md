## (1) VERDICT

**Narrowly novel in claim-shape, not currently defensible** — nothing published scores named αV genu ion pairs against clamped load or mutates them, but the pairs themselves form only in a model missing the ion that the literature says holds that site, so the novel part is the part most likely to be an artifact.

---

## (2) WHAT THE VERIFICATION CHANGED

**The prior sweep's two load-bearing threats both collapsed as novelty threats — and both got worse as validity threats.**

*Li 2024 (doi:10.1007/s10237-023-01805-3, PMID 38308770) — MIS-STATED, and unverifiable in the published body.* The published text is paywalled; the "lock" quote is **PREPRINT_ONLY** (biorxiv 10.1101/2023.04.10.536312). Worse, the prior sweep truncated it mid-clause. Actual: "...which could lock the α subunit **in the opposite bent conformation**." Not the extended state. It is the **isolated α subunit** (Methods: "The α and β subunits alone were obtained by removing the other subunit from integrin αVβ3 (3IJE)"), unloaded, in a direction the intact heterodimer never takes. The sweep also welded two unrelated passages together: the D457–K688 salt bridge (Results, α-alone) and a Discussion proposal that is a **Ca²⁺ coordination** credited to Chen 2011, not a salt bridge — and cited it as preprint Fig 7A when the published figure is 6c. K459, E598, D595 appear **zero times** in the full preprint text and zero times in the published SI (retrieved from static-content.springer.com). **This paper does not threaten obj-083 at all.** The sweep's verdict rested substantially on a citation that, read in full, says something close to the opposite.

*Cormier 2018 (doi:10.1038/s41594-018-0093-x, PMC6214843) — WRONG SYSTEM, and the D595/E598 attribution was ours.* It is an **αvβ8** cryo-EM paper at **6.4 Å global / 4.8 Å headpiece**; Fig 4b calls interacting residues "putative." αvβ3 is only the mutant-validation system. It names **zero** salt bridges, ion pairs, or charged contacts anywhere. The claim that "the knee loop containing D595/E598 does not participate" is a residue mapping *we* supplied — the paper writes the generic knee-region exclusion and never writes D595 or E598.

**Three genuine scoops the prior sweep missed entirely:**

- **Kolasangiani 2025** (doi:10.1016/j.str.2024.11.016, PMC12335895) already publishes **constant 66.4 pN on full-length αVβ3, n=3 replicas**, plus force-release re-bending (αVβ3 bends 4.3 nm in 100 ns). The sweep claimed "no replicate statistics" — false; they ran our n, in explicit solvent, at 100× our system size, with public trajectories and scripts.
- **Driscoll 2021** (doi:10.1016/j.bpj.2021.09.010, PMC8553792) publishes **αVβ3 genu angle and genu separation vs applied force, WT vs point mutants** (122°±7.8 WT vs 166°±18.3 S243E). That is our F-half readout family, already owned.
- **Chen Y 2017** (doi:10.1016/j.matbio.2016.07.002, PMC5237428) already runs "single point mutation abolishes force-regulated bending in αVβ3": L138I and D723R bend "≥35-fold lower than that of WT."

**Two structural findings the sweep never checked, both damaging:**

- **One of the four pairs is not cross-knee.** By CATH/SCOP/Pfam/InterPro on 1JV2 chain A, E598 is interior to the genu disulfide loop (1JV2 SSBOND: CYS A596–CYS A602) and K459 is thigh — same side of the bend. Only D457–K688, K459–E636, D595–K688 span thigh↔calf-1.
- **None of the four pairs exist in either starting structure.** In 1JV2: 22.95 / 19.36 / 25.42 / 7.39 Å. In `results/route_a/extended_seed.pdb`: 21.82 / 6.39 / 15.67 / 7.39 Å. In `extended_state_b.pdb` (production seed): 2.64 / 2.83 / 2.92 / 2.67 Å. All four form during equilibration, with K459 NZ travelling 13–19 Å.

---

## (3) WHAT IS ACTUALLY UNCLAIMED

One sentence, and it must carry its own caveat:

> In a calcium-free implicit-solvent model of extended αVβ3, three thigh–calf-1 ion pairs at the αV genu remain clasped throughout a 60 pN clamp and 60→0 pN ramp but rupture in every unloaded replicate (Fisher p = 0.029), and single alanine substitutions at K459 or E598 abolish this load-dependent protection.

No retrieved paper scores a **named αV genu ion pair against applied force**, and none mutates a genu lock residue. Chen 2011's D457 result is explicitly force-off ("Free MD began after turning off the pulling force on the βA domain"). Chen 2024's named H-bonds are scored vs **extension**, not force, and are bent-state headpiece–tailpiece clasps — opposite polarity. Kolasangiani 2025 and Driscoll 2021 have no residue-resolved interaction analysis (Driscoll's CG model is ~8 Cα/bead — it has no side chains to score).

---

## (4) WHAT CANNOT BE CLAIMED

| Claim | Foreclosed by |
|---|---|
| Constant-force MD of αVβ3 holding the extended state; force-release causing re-bending | **Kolasangiani 2025**, doi:10.1016/j.str.2024.11.016 (66.4 pN, n=3, 4.3 nm re-bend) |
| "Force stabilizes the extended αVβ3 state" | **Chen Y 2017**, doi:10.1016/j.matbio.2016.07.002 — "force accelerates unbending, decelerates bending, stabilizes the extended state" |
| Genu angle / genu separation as a force-dependent readout, WT vs mutant | **Driscoll 2021**, doi:10.1016/j.bpj.2021.09.010 |
| "A single point mutation abolishes force-regulated conformational change in αVβ3" as a novel form | **Chen Y 2017** (L138I, D723R, ≥35-fold) |
| A cross-knee αv thigh↔calf-1 contact holding the extended leg open | **Cormier 2018**, doi:10.1038/s41594-018-0093-x (S546C/N685C disulfide → "only extended but no bent integrins," validated *in αvβ3*) |
| Something at the αV genu involving D457 stabilizes the extended state | **Chen 2011**, doi:10.1371/journal.pcbi.1001086 (D457↔genu Ca²⁺, <4 Å) |
| Contact-resolved WT-vs-mutant SMD in integrins as a method | **Montes 2024**, doi:10.1016/j.bpj.2024.06.008; **Jiang 2021**, doi:10.3389/fmolb.2021.638396 (per-named-pair occupancy at 0/25/50/75 pN, n=3) |
| Conserved ion pairs constraining a loaded integrin element | **Zhang 2018**, doi:10.1111/febs.14335 (β2 αI α7 helix; note: Gly not Ala, and **single** mutations there failed) |
| That 60 pN is inside the experimentally measured window | **Chen Y 2017**: αVβ3 bending/unbending measured "over forces ranging from 0–45 pN." 50 pN clamps exist for lifetimes only. Report the clamp as *above* the measured window. |

---

## (5) THE KILLER OBJECTION

**Verification makes it substantially worse, and it is now residue-specific rather than generic.**

- 1JV2 genu Ca²⁺ (A4008) first shell, measured from `files.rcsb.org/download/1JV2.pdb`: **E636 OE1 at 2.09 Å** — the shortest contact in the entire shell — plus D599 2.11, V601 O 2.54, G597 O 2.55. **D595 at 4.5 Å and E598 at 4.7 Å are second shell.** So one of our scored partners was a direct ligand of the ion we deleted, and two more sit in its second shell.
- **Xiong 2001** (doi:10.1126/science.1064535, PMC2885948) anticipated exactly our conformation: "The base of the thigh domain and the top of the calf-1 domain both contain acidic patches, and these would be expected to face each other in an extended integrin. A well-coordinated calcium ion lies at the genu; in an extended structure it may help to neutralize the negative charge at this interface." We removed the neutralizer and set κ = 0. Lysine recruitment into that interface is the *predicted* consequence of our own deletion.
- **Xie 2004** (doi:10.1073/pnas.0406680101, PMC524430): the E636–Ca²⁺ coordination is likely **maintained** through extension, and "there are no hydrogen bonds or Ca²⁺ coordinations from the thigh domain to the genu." Our K459(thigh)–E636 pair is precisely the contact this paper says does not exist.
- Geometry of the migration: K459 NZ sits **24.36 Å** from the genu-Ca ligand centroid in 1JV2, **5.43 Å** in `extended_state_b.pdb`, **6.14 Å** in `wt_ramp/final`. It moved into the outer shell of the vacated pocket.
- And none of the four pairs exist in 1JV2 *or* in `extended_seed.pdb` — they are entirely products of equilibration in the ion-free, unscreened model.

**One real mitigation:** Chen 2011's partially-extended αVβ3 collapsed back toward bent within 21 ns at zero force **despite 6 bound Ca²⁺** ("the Cα RMSD dropped to ~4 Å within 21 ns"). So our unloaded collapse at 1.2–9.6 ns is not per se evidence of a broken model — cite it.

**Survivable? Only by one computation.** Rebuild the extended seed with the genu Ca²⁺ transferred from 1JV2 A4008 (superpose the G597–C602 loop), set GB-OBC2 κ for 150 mM ionic strength, re-equilibrate, and re-run the loaded/unloaded WT contrast. Three outcomes: (a) all three cross-knee pairs still form and still show load protection → the claim is real and much stronger; (b) K459 no longer reaches E636 but D457–K688 and D595–K688 survive → drop that pair, claim narrows but holds; (c) nothing forms → the result is an ion-vacancy artifact and should not be written up.

---

## (6) CHEAPEST PATH TO A DEFENSIBLE RESULT

Baseline: 26 GPU-h / 12 runs × 10 ns ≈ **2.2 GPU-h per run**.

**0 — Zero GPU, do first (hours of CPU).** Highest value per unit cost.
- Measure the four pairs in **PDB 6DJP** (Cormier's deposited extended αv-leg model). If they are formed there, the morph is validated for free; if not, that is decisive and costs nothing.
- Score the four pairs on the **public explicit-solvent αVβ3 trajectories**: `github.com/tamarabidone/alphaV_vs_alphaIIB` (Kolasangiani 2025, 66.4 pN clamp, 150 mM NaCl, n=3 extension + n=3 bending) and the DynMoCo force-clamp set (doi:10.1016/j.bpj.2026.03.034, PMC13075543). This is an independent extended αVβ3 under load with proper screening. If our pairs appear there, the killer objection collapses immediately — and if we do not run this, a referee will, in an afternoon.
- Re-label: report **three** thigh↔calf-1 pairs, not four. E598–K459 is genu-loop↔thigh (CATH/SCOP/Pfam/InterPro on 1JV2 A). Rewrite "E598A abolishes cross-knee protection" accordingly.

**1 — Calcium + screening control. ~13 GPU-h.** Genu Ca²⁺ restored, κ set for 150 mM, n=3 loaded + n=3 unloaded WT, 10 ns. This is the experiment that decides publication. Do not run anything below until this returns.

**2 — Mutant arms in the corrected system. ~26 GPU-h.** Only if step 1 survives. Consider swapping E598A (inside the C596–C602 genu loop, 4.7 Å from the ion — it perturbs a calcium site the model doesn't contain) for **E636A or K688A**, which are unambiguously cross-knee.

**3 — Optional, ~9 GPU-h.** A 20 pN clamp arm (n=3), landing inside Chen Y 2017's measured 0–45 pN window and near the ~20 pN bending catch peak. Cheap insurance against "your clamp is off-scale."

**Total: ~13 GPU-h to know whether there is a paper; ~39 GPU-h to have one.** Spend nothing beyond step 0 until step 1 reports.