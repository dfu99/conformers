# Matsumoto 2008 switch residues vs our v7 hinge-σ

σ = std of ∠(CA_{i-1}, CA_i, CA_{i+1}) across 615 v7 conformers.
Total CA residues with σ measured: 1652.

| Matsumoto residue | role | σ (°) | rank | percentile | nearest top-20 (±10) | classification |
|---|---|---|---|---|---|---|
| B:633 (Arg633) | primary_snap | 7.6 | 286 | 82.7% | — | near_hit |
| B:375 (Leu375) | snap_sandwich | 12.2 | 47 | 97.2% | B:374 (Δ=1, σ=18.0°) | direct_hit |
| B:389 (Leu389) | snap_sandwich | 7.0 | 383 | 76.9% | — | miss |
| B:374 (Cys374) | interaction_A | 18.0 | 11 | 99.4% | B:374 (Δ=0, σ=18.0°) | direct_hit |
| B:388 (Gly388) | interaction_A | 4.8 | 1013 | 38.7% | — | miss |
| B:404 (Arg404) | hybrid_egf | 8.0 | 240 | 85.5% | — | near_hit |
| B:364 (Glu364) | hybrid_egf | 4.5 | 1141 | 31.0% | B:374 (Δ=10, σ=18.0°) | neighborhood_hit |
| B:543 (Ser543) | interaction_B | 9.5 | 127 | 92.4% | — | near_hit |
| B:548 (Leu548) | interaction_B | 5.7 | 706 | 57.3% | — | miss |
| B:549 (Cys549) | interaction_B | 5.7 | 714 | 56.8% | — | miss |
| B:550 (Asp550) | interaction_B | 4.9 | 979 | 40.8% | — | miss |
| A:305 (Ser305) | alpha_constraint | 4.8 | 1022 | 38.2% | — | miss |

## Summary

- direct_hit (≥95th pct): 2
- near_hit (80–95th pct): 3
- neighborhood_hit (top-20 within ±10): 1
- miss: 6
- absent (residue not in our PDB): 0
