#!/usr/bin/env bash
# Independent replicates of the snap-back suite.
#
# The first pass (obj-082) ran ONE 20 ns trajectory per system. It produced a clean
# separation -- E598A fell to ~117 deg while WT, K459A and the double all settled at
# ~144-147 deg -- but with n=1 the 28-31 deg spread among WT/K459A/double is well
# inside plausible run-to-run variance, so nothing in it is claimable yet. That
# matters more than usual here because the result CONTRADICTS obj-081: the lock
# energetics predicted K459A would snap back ~1.8x harder than E598A, and the
# single run showed the inverse (ratio 0.53).
#
# snapback_md.py calls setVelocitiesToTemperature() with no explicit seed, so each
# invocation draws fresh velocities -- a different --tag is a genuinely independent
# replicate with no code change. run() skips any tag already carrying summary.json,
# so this is resumable.
#
# Usage (from /workspace/route_a):  bash scripts/run_snapback_replicates.sh [ns] [n_extra]
set -euo pipefail
NS="${1:-20}"
EXTRA="${2:-2}"          # extra replicates per system, on top of the existing run
ROOT=/workspace/route_a
PY=/workspace/envs/route_a/bin/python
FINAL="$ROOT/results/route_a/snapback"     # MooseFS network mount
LOCAL=/root/snapback_runs                  # container-local overlay disk
cd "$ROOT"
export CUDA_VISIBLE_DEVICES=0
mkdir -p "$LOCAL"

# Write the live trajectory to LOCAL disk, not to /workspace.
#
# The first replicate attempt died 400 ps in with OSError errno 5 (EIO) inside
# fh.flush(). /workspace is a MooseFS network mount shared across the RunPod
# cluster, and under the load the oxDNA FFS campaign generates (load average
# ~175-200) it intermittently drops a write. A 4 h trajectory that flushes every
# 100 ps has ~200 chances to hit that, and one hit kills the whole run. The
# container overlay disk has no such failure mode, so each system runs to
# completion locally and is copied to /workspace only once, at the end.
run() { # tag  mutations
  if [[ -f "$FINAL/$1/summary.json" || -f "$LOCAL/$1/summary.json" ]]; then
    echo "=== $1 already done — skipping ==="; return
  fi
  echo "=== $1 ($2) — ${NS} ns on A5000 (local scratch) ==="
  rm -rf "${LOCAL:?}/$1"
  "$PY" scripts/snapback_md.py \
        --pdb results/route_a/extended_state_b.pdb \
        --mutations "$2" --tag "$1" --ns "$NS" \
        --platform CUDA --report-ps 100 \
        --out-dir "$LOCAL"
  # publish; keep the local copy if the mount is unhappy right now. Copy the CONTENTS: `cp -r
  # src dst` nests into dst/$(basename src) when dst already exists, which a half-finished
  # earlier publish leaves behind.
  mkdir -p "$FINAL/$1"
  if cp -r "$LOCAL/$1/." "$FINAL/$1/" 2>/dev/null; then
    echo "=== $1 published to $FINAL/$1 ==="
  else
    echo "=== $1 WARNING: copy to $FINAL failed; result retained at $LOCAL/$1 ==="
  fi
}

# Replicate the three systems that carry the falsification: the WT baseline and the
# two singles whose rank order is in dispute. The double tracked WT exactly in run 1;
# it is deferred until the singles are resolved.
for r in $(seq 2 $((EXTRA + 1))); do
  run "wt_r${r}"    WT
  run "e598a_r${r}" A:598:ALA
  run "k459a_r${r}" A:459:ALA
done

echo "ALL REPLICATES DONE"
