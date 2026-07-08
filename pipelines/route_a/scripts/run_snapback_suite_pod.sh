#!/usr/bin/env bash
# Snap-back suite on the RunPod A5000 (pod paths + pod env). Runs highest-contrast
# systems first; skips any already finished (resumable). Called by the watcher once the
# A5000 slot is free, or manually after confirming no oxDNA is running.
set -euo pipefail
NS="${1:-20}"
ROOT=/workspace/route_a
PY=/workspace/envs/route_a/bin/python
cd "$ROOT"
export CUDA_VISIBLE_DEVICES=0

run() { # tag  mutations
  if [[ -f "results/route_a/snapback/$1/summary.json" ]]; then
    echo "=== $1 already done — skipping ==="; return
  fi
  echo "=== $1 ($2) — ${NS} ns on A5000 ==="
  "$PY" scripts/snapback_md.py \
        --pdb results/route_a/extended_state_b.pdb \
        --mutations "$2" --tag "$1" --ns "$NS" \
        --platform CUDA --report-ps 100 \
        --out-dir results/route_a/snapback
}

run wt      WT                    # control (should hold extended)
run double  A:459:ALA,A:598:ALA   # break the core cross-knee lock (biggest expected effect)
run k459a   A:459:ALA             # single
run e598a   A:598:ALA             # single
"$PY" scripts/plot_snapback.py --dir results/route_a/snapback \
      --out results/route_a/snapback/route_a_snapback.png || true
echo "ALL DONE"
