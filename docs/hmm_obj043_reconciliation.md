# HMM ↔ obj-043 state-fraction reconciliation

Audit-2026-05-05 §20.5 v12 P=4 follow-up to obj-051. Resolves the
16-percentage-point BC discrepancy between obj-048's 3-state HMM
stationary distribution (BC = 40 %) and obj-043's hard-threshold
partition (BC = 25 %).

## The two pictures

| state | obj-043 hard threshold (CV0) | obj-048 HMM stationary π |
|-------|------------------------------|---------------------------|
| BC    | ≤ 65 Å → 26.4 %              | 40.5 %                    |
| Intermediate | 65 – 78 Å → 45.1 %    | 35.2 %                    |
| EC    | 78 – 85 Å → 26.0 %           | 24.4 %                    |
| EO*   | > 85 Å → 2.6 %               | (folded into HMM EC)      |

## Where the discrepancy lives — Viterbi cross-tab

|              | obj-043 BC ≤65 | obj-043 Inter 65-78 | obj-043 EC 78-85 | obj-043 EO* >85 |
|--------------|---------------:|---------------------:|------------------:|----------------:|
| HMM BC       | 433            | **215**              | 8                 | 2               |
| HMM Inter    | 1              | 517                  | 69                | 0               |
| HMM EC       | 0              | 10                   | 350               | **40**          |

## The 215-frame BC↔Intermediate transition zone

The single largest off-diagonal cell is **215 frames assigned HMM-BC
but obj-043-Intermediate** (13.1 % of all frames). Their CV0 range is
**[65.0, 77.7] Å, mean 67.8 Å** — sitting just above the obj-043
65-Å BC/Intermediate boundary.

This is exactly the **"high-BC" sub-state at mean 66.0 Å** revealed by
the obj-051 4-state HMM. The 4-state model's split:

  - low-BC: 54.5 Å (very-bent, "true" BC)
  - high-BC: 66.0 Å (BC↔Intermediate transition zone)
  - Intermediate: 75.0 Å
  - EC: 82.5 Å

shows the transition zone explicitly. The 3-state HMM lumps low-BC +
high-BC into a single "BC" mode (μ = 61, σ = 8.8 Å). The obj-043 hard
threshold instead splits the same frames at exactly 65 Å and assigns
them to Intermediate.

## EO* frames

obj-043's 42 EO* frames (CV0 > 85 Å) split: 40 → HMM EC, 2 → HMM BC.
The 40 EO* frames absorbed into HMM EC explain why HMM EC and obj-043
EC fractions agree at 24-26 % despite the HMM having no EO state.

## Resolution

**Both pictures are internally consistent**; they answer different
questions.

  - **obj-043** answers "what fraction of time is CV0 in the BC band?"
    → use the hard threshold.
  - **obj-048** answers "what fraction of time is the protein in the
    BC kinetic mode?" → use HMM stationary π.

The 215-frame difference is the BC↔Intermediate transition zone at
CV0 ≈ 65–78 Å. The 4-state obj-051 HMM names it explicitly as
"high-BC". For the paper, **report obj-043 fractions for the
geometric / FES interpretation** (matches Springer 2014 BC band) and
**report obj-048 stationary π for the kinetic interpretation**
(matches the Markov-chain dwell-time framework). Both are correct.
