# Head Scientist — Per-Project Critique for Top-Venue Submission
## 2026-03-22

---

## 1. WorldNN — C_i Coordination Quality (Target: ICLR 2027)

**VERDICT: Do not submit now. Revise first. Strong workshop paper; main-track with second task.**

**Novelty:** C_i (cosine alignment between degraded-perception policy and optimal-state action) fills a genuine gap. No prior work defines this specific metric across a designed perception-degradation ladder. Closest: asymmetric actor-critic (privileged training), DreamerV3 prediction probes, quasimetric RL. C_i differs by operating in real state space, not latent.

**Rigor Gaps:**
- Only one task (rock-push, 4D/2D) — threshold values likely task-specific
- 3 seeds/condition, C_i≥0.8 bucket has only N=5
- C_i is post-hoc (requires optimal action) — no proxy estimator proposed
- No ablation of metric definition (why cosine vs L2?)
- Threshold claim is correlational, not causal (needs factorial ANOVA)

**Framing Fix:** Rename "coordination quality" → "sensorimotor alignment" to avoid MARL confusion. Lead with "blind cat" hook more aggressively.

**Experiments Needed:**
1. Second task (8D+ locomotion or manipulation) — non-negotiable
2. C_i dynamics during training — makes it actionable as diagnostic
3. Self-supervised proxy for C_i without oracle access
4. Increase seeds to 5 for threshold-critical conditions
5. Formal interaction test (log-linear model with p-value)

**Target:** ICLR 2027 (Oct 2026 deadline) if second task replicates.

---

## 2. CorticalNN — Bio Topology Negative Result (Target: Workshop/Journal)

**VERDICT: Do not target NeurIPS/ICLR main. Target Neuro-AI workshop or PLOS CompBio.**

**Novelty:** Partially known (Goulas 2021 Bio2Art found similar null result for cortical topology RNNs). Your result differs: generative growth model, clean topology-vs-weight-init ablation, larger scale. A 2026 biorxiv preprint claims evolutionary topology wins — making your counter-result timely.

**Rigor Gaps:**
- Density confound: bio networks sparser, Kaiming init favors sparse random
- Single task (MNIST) — may be task-specific null
- "Spatial constraint hurts" asserted not demonstrated — need structural feature mediation analysis
- 50-sample survey OOM'd, weight init ablation failed

**Framing Fix:** Don't frame as "our method failed." Frame as: "Random growth rules cannot produce useful frozen reservoirs — spatial constraint is a bottleneck." The density threshold (bio wins >1000 connections) is publishable.

**Experiments Needed:**
1. Structural feature mediation (partial correlation of Gini/clustering/LCC vs advantage)
2. CIFAR-10 or sequential MNIST (one non-trivial task)
3. Density-stratified figure with confidence bands
4. Complete the 50-sample survey (fix OOM)

**Target:** NeurIPS Neuro-AI workshop, PLOS CompBio, or Neural Networks.

---

## 3. cadnano2 — AI-Assisted DNA Origami Design (Target: DNA 32)

**VERDICT: Submit Track B now; pursue Track A if RLVR yields numbers by April 3.**

**Novelty:** Partially novel. scadnano (2020) and inSēquio (2024) provide scriptable DNA design. MENDEL does automated sequential commands. Your differentiator: NL intent → constrained tool calls → verification reward → local model training. That pipeline is not in the literature.

**Rigor Gaps:**
- No wet-lab validation (oxDNA simulations exist but no folding confirmation)
- Agent capability benchmarks thin — no success rate on held-out tasks
- No expert-vs-agent time comparison baseline

**Framing Fix:** "Constraint-aware agentic interface enabling non-expert DNA origami design with verifiable correctness" — better than "parametric design tool." The RLVR loop is the most differentiated piece.

**Experiments Needed:**
1. Agent success rate on held-out scaffold routing tasks (before/after RLVR)
2. At least one 6-helix bundle oxDNA equilibration
3. User study or timing comparison: NL agent vs GUI

**Target:** DNA 32 Track B (April 3 deadline). Track A if RLVR produces numbers.

---

## 4. FIND-SNP — GPU FFS for SNP Discrimination (Target: NAR or J Phys Chem B)

**VERDICT: 2-3 months from publishable. NUPACK simulator story is closer than oxDNA FFS.**

**Novelty:** oxDNA FFS for strand displacement is established (Machinek 2014, Zhang 2022 already scanned mismatches with 124x discrimination). Your 1.4x approach flux ratio is NOT the discrimination signal — mismatches slow branch migration, not approach. GPU-batch FFS at scale across 94 sequences would be novel.

**Rigor Gaps:**
- Only approach flux measured — full FFS with all interfaces needed for rate constants
- No temperature sweep
- No comparison of oxDNA rates to experimental dataset
- Topology bugs (obj-023/024) raise validity concerns

**Framing Fix:** Lead with the branching-probability mechanism (obj-017: P(forward) drops 0.50→0.06, 161x slowdown). That's a genuine mechanistic insight. FFS is validation, not the headline.

**Experiments Needed:**
1. Complete full FFS on ≥10 PM + 10 SNP sequences
2. Validate against Machinek/Broadwater experimental data
3. Temperature sweep (40-60°C)
4. GPU vs CPU wall-clock comparison

**Target:** NAR (methods+biophysics hybrid) or J Phys Chem B letter.

---

## 5. RL-Arm — Muscle Dynamics in RL (Target: Workshop)

**VERDICT: Not ready. 2/5 readiness. Workshop-level at best.**

**Novelty:** Not novel at CoRL/RSS level. OpenSim-RL (NeurIPS 2018) established field. Action Space Design (CoRL 2025) and musculoskeletal arm papers (Springer 2025) directly overlap.

**Rigor Gaps (fatal):**
- Tau ablation is confounded (trajectory shifted simultaneously, obj-013)
- No baseline against torque-controlled or PD-controlled arm
- No comparison to MyoSuite/OpenSim-RL — reviewers will ask why
- 2-DOF is kinematically degenerate for a robotics claim
- No real-robot transfer or physical plausibility argument

**Experiments Needed:**
1. Clean tau ablation (fix trajectory, vary tau in {0, 0.02, 0.05, 0.1})
2. PPO vs SAC baseline
3. Mass generalization curve
4. Extend to 3-DOF minimum

**Target:** IROS or RA-L (letters format). More realistically: RSS/CoRL workshop.

---

## 6. conformers — MSA-Subsampled Conformer Validation (Target: AI4Science Workshop)

**VERDICT: Not ready. 4-6 weeks with AF2 test + one baseline comparison.**

**Novelty:** MSA-subsampling for conformers is saturated (AFCluster Nature 2023, SPEACH_AF 2022, AFsample2 2025). What's potentially novel: steered MD pull trajectories + pseudo-AFM image generation + TM-score validity filtering as an integrated pipeline for integrin conformers.

**Rigor Gaps:**
- Only 19 frames scored — TM=0.06 at frame 300 is likely a destroyed structure, not a valid conformer
- No AFCluster comparison baseline
- Pseudo-AFM images not quantitatively evaluated
- AF2 MSA-depth test still pending
- Single protein system

**Experiments Needed:**
1. Complete AF2 MSA-depth sensitivity test
2. Validate against PDB structures (1JV2, 3IJE, 8IJ5) for AVB3 bent/extended
3. Run AFCluster as direct comparison
4. Second protein system for generality

**Target:** NeurIPS/ICML AI4Science workshop or ML4Molecules workshop.

---

## Summary Ranking (publication readiness)

| Project | Readiness | Target | Timeline |
|---------|-----------|--------|----------|
| cadnano2 | 3.5/5 | DNA 32 Track B | Now (April 3) |
| WorldNN | 3/5 | ICLR 2027 | 3-4 months |
| FIND-SNP | 2.5/5 | NAR / J Phys Chem B | 2-3 months |
| CorticalNN | 2/5 | Neuro-AI workshop | 1-2 months |
| conformers | 1.5/5 | AI4Science workshop | 4-6 weeks |
| RL-Arm | 1/5 | IROS/workshop | 3+ months |
