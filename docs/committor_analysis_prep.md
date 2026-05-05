# Committor analysis prep — Reviewer D last open

**Audit follow-up §12.7 P=3.** Pre-run planning document for the
committor analysis Reviewer D (Elber/Chipot) is expected to demand
after the ΔG(CV0) profile (obj-038) lands.

This is a planning doc, not an experiment. The ΔG profile (obj-038)
shows a saddle near CV0 ≈ 67 Å between BC (≈55 Å minimum) and EC
(≈80 Å minimum). Reviewer D will ask: "from a saddle seed, do
unbiased trajectories commit to BC and EC with the symmetry your
ΔG predicts?"

---

## What Reviewer D will ask

Two specific committor concerns:

1. **Saddle seed validity.** Where is the BC↔EC barrier in
   *unbiased* MD? Does a frame at ΔG-saddle CV0 ≈ 67 Å actually
   sit on the dividing surface, or is it kinetically biased toward
   one basin?
2. **Committor symmetry.** From identical seeds, what fraction of
   short trajectories commit to EC vs BC? Symmetric committor (50/50)
   confirms saddle. Asymmetric committor (e.g., 80/20) means the
   saddle estimate is off.

---

## What we already have

- **obj-038**: ΔG(CV0) profile from 1645 fitted-trajectory frames.
  Free-energy minimum at CV0 ≈ 70-75 Å (Intermediate); saddle
  inferred at ≈ 67 Å between BC and Intermediate, ≈ 78 Å between
  Intermediate and EC.
- **F3 bootstrap**: 95 % CI < 0.5 kcal/mol across [50, 85] Å.
- **v7 library**: 615 frames spanning CV0 [47.3, 85.0] Å — has
  frames in the CV0 ≈ 67 Å window (12 frames in 65-69 Å bin).
- **Already committed**: `pipelines/conformer-library/scripts/`
  contains `domain_steering.py` for biased MD; can be modified to
  remove bias for unbiased seeded runs.

---

## Proposed analysis (~30 GPU-hr)

### Setup

- **Saddle seeds**: pick 10 frames from v7 library with CV0 ∈
  [66.5, 67.5] Å (the 1-Å window around the ΔG saddle). Spread
  across orientations to avoid local seed bias.
- **Per-seed replicates**: 5 distinct random velocity seeds per
  saddle frame → 50 trajectories total.
- **Length**: 5 ns per trajectory (committor experiments need
  short, many; not long, few).
- **Termination**: stop early if CV0 enters BC basin (CV0 < 60 Å)
  or EC basin (CV0 > 75 Å) for 100 ps continuous.
- **MD details**: same OpenMM setup as `domain_steering.py` but
  with `apply_steering_preset(...)` skipped.

### Expected outcome (pre-registered prediction)

If obj-038's ΔG estimate is correct:

- Committor fraction toward EC: ≈ 0.5-0.6 (slightly biased toward
  Intermediate/EC because the ΔG profile shows EC basin minimum
  closer than BC minimum).
- 50 trajectories should split ≈ 25 EC / 25 BC, with a few
  trajectories that don't decide in 5 ns and stay in Intermediate.

If obj-038's ΔG is *miscalibrated* (e.g., HS-AFM populations are
biased by surface-binding effects):

- Committor will be strongly asymmetric (e.g., 40/10/0 or 5/45/0).
- Falsifies obj-038's interpretation of the fitted-trajectory
  populations as a true thermodynamic ΔG.

### Compute budget

- 50 trajectories × 5 ns × ~3 ns/day on A100-80GB = ~80 GPU-hr
- Achievable with 1 week of PACE A100-80GB time
- Can be parallelized across 8 jobs simultaneously → 1 day wall

### Submit script outline

```bash
# pipelines/committor/scripts/submit_committor_slurm.sh
#SBATCH -A gts-yke8
#SBATCH -p gpu-a100
#SBATCH -C A100-80GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.fu@emory.edu
#SBATCH -t 24:00:00
#SBATCH --array=0-49

source venv_protenix/bin/activate
SEED_FRAME=$((SLURM_ARRAY_TASK_ID / 5))
VEL_SEED=$((SLURM_ARRAY_TASK_ID % 5))

python pipelines/committor/scripts/run_committor.py \
    --seed-pdb data/runs/avb3/conformers/all_frames_bent_extended/saddle_${SEED_FRAME}.pdb \
    --vel-seed $VEL_SEED \
    --max-ns 5.0 \
    --bc-threshold 60 \
    --ec-threshold 75 \
    --output data/runs/avb3/committor/seed_${SEED_FRAME}/vel_${VEL_SEED}/
```

### Analysis script outline

```python
# pipelines/committor/scripts/analyze_committor.py
# Read all 50 trajectories' final commit state.
# Build histogram: count(EC) / count(EC + BC) = committor fraction.
# Also: time-to-commit distribution (mean, std, max).
# Plot: committor histogram + time-to-commit histogram.
```

---

## Decision threshold

Pre-registered:

- **Validating** (committor in [0.40, 0.65]): obj-038's ΔG
  interpretation is supported. Reviewer D's concern resolves
  with a single short experiment.
- **Falsifying** (committor outside [0.30, 0.75]): obj-038's
  ΔG is biased by HS-AFM surface effects. The paper claim
  weakens — but this is itself a publishable finding (HS-AFM
  surface biases the equilibrium population).

---

## Risk register (light)

| risk | P_fail | mitigation |
|------|--------|------------|
| Saddle seed kinetically biased | 25 % | spread 10 seeds across orientations |
| 5 ns too short to commit | 30 % | extend to 10 ns for non-committers |
| Wall-time exceeded on PACE | 10 % | array job; each task independent |
| OpenMM unbiased setup divergence | 15 % | reuse domain_steering.py minus bias |
| MD instabilities at saddle | 5 % | small Δt + Langevin γ=1 |

---

## Status + dependencies

- Compute: blocked on PACE A100-80GB allocation (same queue as
  route-A; coordinate timing).
- Code: needs `pipelines/committor/` directory + 2 scripts (run +
  analyze). Both small (~150 lines each).
- Prerequisites: route-A kickoff (which uses the same allocation).

---

## Cross-references

- obj-038 (ΔG(CV0) profile)
- F3 bootstrap (CI on ΔG)
- audit-2026-05-05.md §6 Reviewer D analysis (`Open (3): ΔG
  profile, committor at midpoint, why CV0 over leg-knee
  dihedral`) — committor was the second open item.
- `pipelines/avb3-conformers/scripts/domain_steering.py` (template
  for the unbiased setup)

---

_Author: AFK audit deepening pass v4, 2026-05-05 23:45._
_Status: planning only; compute on hold. Activate after route-A
kickoff completes (would share PACE allocation)._
