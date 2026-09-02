# Extended conformations across integrin variants

Genu-hinge incremental morph (vacuum ff14SB, CPU-only), generalised from the
alphaVbeta3 route-A recipe by alignment-transferred domain boundaries.

| variant | start PDB | conf. | genu bent→ext (°) | extent bent→ext (Å) | Rg bent→ext (Å) | steps | acceptance |
|---|---|---|---|---|---|---|---|
| alphaVbeta3 | 1JV2 | high | 41.2 → 176.2 | 130.0 → 209.9 | 39.1 → 66.6 | 16 | PASS |
| alphaIIbbeta3 | 3FCS | high | 36.3 → 179.4 | 118.5 → 207.3 | 39.0 → 63.5 | 17 | PASS |
| alpha5beta1 | 7NXD | high | 78.0 → 178.2 | 160.2 → 214.7 | 48.7 → 63.3 | 12 | PASS |
| alphaMbeta2 | 7USM | low | 48.2 → 172.3 | 109.0 → 207.9 | 39.8 → 61.4 | 15 | PASS |
| alphaXbeta2 | 5ES4 | low | 12.9 → 175.5 | 137.3 → 238.6 | 54.8 → 80.3 | 20 | PASS |

## alphaVbeta3 positive control

Published route-A endpoint (obj-075): Rg 67.32 A, extent 211.36 A.
This run: Rg 66.57 A, extent 209.9 A (deviation 1.1% / 0.7%).
Reproduces published endpoint: **True**

## Comparison with experimental extended structures

| variant | reference | exp genu (°) | morph genu (°) | exp extent (Å) | morph extent (Å) | |Δextent| (Å) |
|---|---|---|---|---|---|---|
| alphaVbeta3 | 8XEN | 145.3 | 176.2 | 214.5 | 209.9 | 4.6 |
