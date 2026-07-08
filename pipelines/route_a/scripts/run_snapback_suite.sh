#!/usr/bin/env bash
# ⚠️ DEPRECATED / DO NOT USE — PI directive 2026-07-08: "Do not use local GPU. Just wait
# for a GPU slot on the A5000." The snap-back runs on the RunPod A5000 via
# scripts/run_snapback_suite_pod.sh, launched by scripts/wait_and_run_snapback.sh once the
# oxDNA jobs clear. This local-GPU version is kept only as a record; do not run it.
#
# Snap-back suite: WT control + linchpin mutants, implicit-solvent MD on the local GPU.
# RUN ONLY after receiving "GPU ACCESS GRANTED" (sets CUDA_VISIBLE_DEVICES=0).
# Uses the isolated CUDA-enabled OpenMM env (conda env: snapback).
#
#   bash pipelines/route_a/scripts/run_snapback_suite.sh [ns]
#
# WT should stay extended (knee ~175°); mutants that release the genu salt-bridge lock
# should show the knee angle fall (re-bending). The WT-vs-mutant differential is the signal.
set -euo pipefail
cd "$(dirname "$0")/../../.."          # repo root
NS="${1:-20}"
PY=/home/dan/anaconda3/envs/snapback/bin/python
export CUDA_VISIBLE_DEVICES=0

run() { # tag  mutations
  if [[ -f "results/route_a/snapback/$1/summary.json" ]]; then
    echo "=== $1 already done (summary.json present) — skipping ==="
    return
  fi
  echo "=== $1 ($2) — ${NS} ns ==="
  "$PY" pipelines/route_a/scripts/snapback_md.py \
        --mutations "$2" --tag "$1" --ns "$NS" \
        --platform CUDA --report-ps 100
}

# Ordered for maximum information first: WT (control) then the double mutant give the
# strongest expected contrast, so even a short GPU window yields the key comparison.
# Each system is skipped if already finished (summary.json present) -> resumable suite.
run wt      WT                    # control (should hold extended)
run double  A:459:ALA,A:598:ALA   # break the core cross-knee lock hard (biggest effect)
run k459a   A:459:ALA             # single: break E598-K459 + K459-E636 (K459 in two bridges)
run e598a   A:598:ALA             # single: break the top-ranked lock (E598-K459 only)

$PY pipelines/route_a/scripts/plot_snapback.py
echo "ALL DONE"
