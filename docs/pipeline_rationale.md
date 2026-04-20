# Pipeline Rationale — Why Protenix + Boltz + AFCluster + Steered MD?

This document anticipates the question: "Why this four-part pipeline instead
of just AlphaFold-Multimer with MSA subsampling?"

## The short answer

Because single-model structure predictors collapse onto one preferred
conformation — for αVβ3 that's the extended state. We needed a
conformational *ensemble* spanning bent ↔ extended, and no single predictor
gives that. Each of the four components contributes something the others
lack:

| Component | Role | What it gives |
|-----------|------|---------------|
| Protenix  | Primary co-fold + template-conditioned inference | MSA-depth-invariant bent seed (with template biasing) |
| Boltz     | Alternative inference + ranking | Cross-validation of Protenix predictions; handles partner-swap more flexibly |
| AFCluster | MSA clustering → coevolution-diverse sweep | Broader conformational sampling than single-MSA approaches |
| Steered MD (CV-distance extend / bent) | Transitional conformers between end states | Fills CV0 = 47.3–85.0 Å continuum |

## Component-by-component justification

### Why Protenix (and not just AlphaFold-Multimer)?

The publicly released AlphaFold-Multimer collapses to the bent αVβ3 state
and is **depth-invariant** — MSA subsampling does not produce extended
decoys for this protein. (Observation from obj-006, 2026-03-20: TM-score
spread across 25 predictions was 0.1–2.9Å, essentially monomodal.)

Protenix with template conditioning gave us a workable bent seed. We
still rely on Protenix for initial co-fold (not conformational diversity).

### Why Boltz?

Two complementary roles:
1. **Template-free alternative** when we want to check whether a Protenix
   output is model-dependent.
2. **Partner-swap flexibility** — Boltz handles multi-chain rearrangement
   (A5B1 + SpyTag + Streptavidin) more gracefully in our tests.

Boltz alone does not solve the ensemble problem — it also prefers
extended αVβ3 as its default state.

### Why AFCluster?

AFCluster clusters the MSA by coevolutionary signal, then runs folding on
each cluster. For proteins where MSA subsampling *does* produce
conformational diversity (flhAC, various enzymes), this is the cleanest
way to generate alternatives. For αVβ3 it was a **negative result** —
AFCluster confirmed that the MSA signal is not the bottleneck for αVβ3
conformational diversity. Worth keeping as a documented negative control,
not as active pipeline infrastructure.

### Why steered MD?

This is the component that actually generates the bent↔extended continuum.
CV-distance-based flat-bottom biases on inter-domain distances move the
structure along a reaction coordinate (α-leg extension, β-leg extension,
headpiece separation). We run it twice:

- `cv_distance_extend`: target 20/18/16 nm → produces conformers at
  CV0 52.9–85.0 Å (generally extended half)
- `cv_distance_bent`: target 4.0/3.5/2.0 nm → produces conformers at
  CV0 47.3–63.7 Å (generally bent half)

The 615-frame union covers the full observed AFM range. Without steered
MD, our conformer library would be stuck at the Protenix/Boltz mode and
we could not do library-based template matching.

## What an AlphaFold-Multimer ablation would cost

Reviewer B's implicit request: run AF2-Multimer with several MSA depths
(default, 10%, 5%, 1%), record ranking_score and CV distribution per
sample, plot against Protenix+Boltz+AFCluster+MD output. Estimated
compute:

- AF2-Multimer on αVβ3 (~1600 residues) with full MSA: ~2 h / sample on
  A100 80GB.
- 5 MSA depths × 5 samples each = 25 samples × 2 h = 50 GPU-hours.
- Cost on RunPod A4500 (~3x slower): 150 h ≈ $150 at $1/h.

This is tractable but not on the main path while we still lack PACE
access. Listed as follow-up experiment in `tasks/planning.md`.

## What this pipeline enables that single-predictor approaches don't

1. **Library-based template matching** of HS-AFM frames (the central
   contribution of this project) — needs >300 conformers spanning the CV
   manifold. Single-predictor approaches give ~10 structural modes.
2. **Frame-by-frame CV recovery** — correlation matching assigns each AFM
   frame a CV triplet by finding the nearest pseudo-AFM image. This is
   only meaningful with a library that actually contains the right CVs.
3. **State-occupancy recovery** — v6 (no bent frames) vs v7 (+bent frames)
   diverges for video1 (12% → 43% BC) but not for video2 (12% → 19% BC).
   The pipeline is sensitive to which end states are available in the
   library — a property not testable with a single-mode predictor.
