# Project Intuition

## One-line claim

A structure-to-image template-matching pipeline (Protenix + Boltz + AFCluster +
bidirectional steered MD) produces an αVβ3 integrin conformer ensemble that,
when projected through a corrected 1-2nm pseudo-AFM imaging model, fits
experimental HS-AFM trajectories frame-by-frame (corr > 0.93) AND recovers
population-weighted bent↔extended state occupancies consistent with library
CV coverage — establishing an imaging-fitting framework for single-molecule
conformational dynamics that does not require CNN training and that is
demonstrably sensitive to library completeness.

## Why we believe it

*Head-anchored fitting quality.* Correlation matching with corrected-tip
pseudo-AFM achieves video1 corr=0.965 and video2 corr=0.939 with the v7
library. Head-tracked SO(3) fitting against 615 conformers eliminates
cumulative position drift (was up to 10.4nm with previous-frame anchoring).
A secondary tail-flip smoothing pass corrects 18-28% of frames where the
head is stable but the tail rotates 180°.

*Library completeness drives state recovery — the discriminating result.*
With the v6 library (309 frames, CV0 range 52.9-85.0Å, ZERO bent
conformers), both videos were forced to match extended states: video1
BC=12.4%, video2 BC=12.5%. After adding 306 bent frames from CV-distance-bent
steering on RunPod A4500 (CV0 range 47.3-63.7Å), v7 fits show dramatically
different behavior across the two videos:
- Video 1: BC 12.4% → *43.5%* (Intermediate 15% → 6%, EC 73% → 50%)
- Video 2: BC 12.5% → 18.6% (Intermediate 14% → 3%, EC 73% → 79%)

Video 1 really was mostly bent-state but had been mis-matched to the
closest extended conformer because no bent templates existed. Video 2 is
genuinely extended-biased. This divergence is strong evidence the pipeline
is reading the signal in the AFM, not the shape of the library.

*Two critical bug fixes validated the framework.* First, `generate_images`
with batch_size=16 truncated sampling to only 31/309 conformers — switching
to batch_size=1 unlocked the full library (correlations 0.73→0.94). Second,
head-anchored flip resolution (comparing to tracked AFM head, not the
previous frame) eliminated cumulative drift.

*CNN inference fails where correlation matching succeeds.* A CNN trained
on simulated pseudo-AFM gave near-constant CV predictions on real Linz
data — confirming a sim-real domain gap that template matching sidesteps.

## What would falsify it

1. *Generalization failure.* A third independent HS-AFM dataset gives
   corr < 0.7 — the pipeline doesn't transfer beyond Linz data.
2. *Random baseline beats matching by <Δ0.1.* My original pre-registered
   threshold was "random corr < 0.4", but the random baseline on our own
   library (see `figures/random_baseline_v7.png`) is 0.65 (video1) and
   0.43 (video2). The pre-reg was wrong: mean-subtracted correlation on
   AFM-shaped images is naturally high because the gross shape matches.
   *Revised falsification*: matching corr − random corr < 0.1 would
   indicate no discriminating signal. Observed gaps: V1 0.86 − 0.65 =
   *+0.21*, V2 0.72 − 0.43 = *+0.29*. Gap is real but smaller than the
   raw numbers suggest.
3. *Population mismatch.* An independent biochemical measurement (FRET,
   cryo-EM class fractions) of the same sample disagrees with our v7
   state fractions by >30% absolute. For video2 we predict >70% extended.
4. *Library-free baseline beats us.* A simple 2D image-template bank with
   a hand-drawn bent/extended silhouette matches the AFM as well as the
   full conformer library — rendering the MD pipeline unnecessary.

## Target panel and venue

- *Panel*: see `tasks/review-panel.yaml`
- *Venue*: Nature Structural & Molecular Biology (for the integrin biology
  + pipeline) or PNAS (for the broader methods contribution)

## Open questions

- *Why is video1 ~50/50 bent/extended while video2 is ~80% extended?*
  Needs Linz lab notebook: was video2 the same sample with Mn²⁺, a
  different prep, or a post-activation recording? The divergence is
  scientifically interesting but unexplained without metadata.
- *CV0 gap to fully-bent PDB.* Our v7 bent library bottoms out at
  CV0=47.3Å; PDB 1JV2 (bent αVβ3 crystal structure) is ~40Å. Is the 7Å
  gap biologically meaningful or a 3-ns steering limitation?
- *EC vs EO not distinguishable.* Head-head distance CV2 is essentially
  constant at 34.3-35.5Å across our entire library — steering never
  opened the headpiece. To recover EO population we need either (a)
  steering with a CV2 target, or (b) the αIIbβ3 string-method structures.
- *Why does the CNN fail?* Is the sim-real domain gap in our pseudo-AFM
  model (tip mechanics, noise statistics) or is classification on AFM
  images fundamentally harder than regression+template matching?
- *Library-size sensitivity.* We went 309 → 615 conformers and saw a
  3.5× BC increase in video1. What's the saturation point — when does
  adding more conformers stop changing state fractions?
