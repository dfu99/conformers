# Multi-integrin first-principles bent-state features

Bent-state geometry from RCSB-deposited heterodimer ectodomain PDBs. Features are computed without reference to experimental activation kinetics or HS-AFM data.

| PDB | name | full ectodomain | α / β residues | Rg (Å) | long axis (Å) | head-tail (Å) | α head↔tail (Å) | β head↔tail (Å) | head↔head (Å) | head-leg buried SASA (Å²) |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1JV2 | αVβ3 (bent crystal) | yes | 934/542 | 39.23 | 129.22 | 42.61 | 55.55 | 31.07 | 40.63 | 13442 |
| 3FCS | αIIbβ3 (bent crystal) | yes | 1154/756 | 39.14 | 119.06 | 37.2 | 49.95 | 26.88 | 43.26 | 15298 |
| 3VI4 | α5β1 headpiece (bound peptide) | no | 609/435 | 34.94 | 116.7 | 20.05 | 50.76 | 12.02 | 57.3 | — |
| 4UM9 | αVβ6 head + peptide | no | 686/520 | 35.03 | 115.36 | 28.37 | 51.3 | 6.96 | 52.73 | — |
| 5FFG | αVβ6 headpiece (apo) | no | 767/325 | 32.84 | 111.89 | 32.26 | 51.24 | 14.28 | 39.27 | — |
| 6OM2 | αVβ8 headpiece + LAP | no | 680/323 | 33.34 | 113.34 | 36.0 | 52.45 | 21.11 | 41.61 | — |
| 6UJC | αVβ8 headpiece + Fab | no | 598/363 | 33.66 | 115.74 | 32.75 | 51.14 | 15.72 | 43.53 | — |

## Predicted bent-fraction ranking

Bent-fraction prediction uses head-leg buried SASA as the primary stabilizer (more contact ⇒ stable bent ⇒ higher experimental bent population), tied with smaller head-tail distance (more compact ⇒ bent). Headpiece-only structures are excluded since they cannot constrain bent-fraction.

| rank | PDB | name | head-leg buried SASA (Å²) | head-tail (Å) | predicted bent fraction (relative) |
|---:|---|---|---:|---:|---|
| 1 | 3FCS | αIIbβ3 (bent crystal) | 15298 | 37.2 | highest |
| 2 | 1JV2 | αVβ3 (bent crystal) | 13442 | 42.61 | #2 |
