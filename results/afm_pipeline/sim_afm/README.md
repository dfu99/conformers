# Simulated HS-AFM from v7 Fitted Trajectories

Forward-rendered each v7 fitted PDB through the afmfold imaging pipeline
(tip R=1.5nm, angle=20°, noise 0.1nm, resolution 0.98nm/px).

## Self-consistency metric: sim vs real AFM correlation
- video1: mean 0.824, median 0.846, std 0.089
- video2: mean 0.722, median 0.724, std 0.130

Context: random-conformer baseline is corr 0.65 (V1), 0.43 (V2).

## Interpretation
Above-random agreement shows the library + imaging model jointly reproduce
the real HS-AFM images frame-by-frame. Not perfect — residual gap suggests
(a) our tip model is simpler than the real one (no contact mechanics),
(b) scanning artifacts in real AFM (drift, tip convolution asymmetry),
(c) library completeness still imperfect (no EO states, CV0 gap to fully
bent).

## Files
- `video{1,2}/sim_afm.gif` — simulated AFM video
- `video{1,2}/sim_vs_real.gif` — side-by-side (real | sim) comparison
- `video{1,2}/sim_vs_real_correlation.npy` — per-frame cosine similarity
- `video{1,2}/sim_vs_real_corr.png` — time-course plot
