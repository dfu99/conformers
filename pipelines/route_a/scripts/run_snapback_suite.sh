#!/usr/bin/env bash
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
  echo "=== $1 ($2) — ${NS} ns ==="
  "$PY" pipelines/route_a/scripts/snapback_md.py \
        --mutations "$2" --tag "$1" --ns "$NS" \
        --platform CUDA --report-ps 100
}

run wt      WT                    # control
run k459a   A:459:ALA             # break E598-K459 + K459-E636 (K459 is in two bridges)
run e598a   A:598:ALA             # break the top-ranked lock partner
run double  A:459:ALA,A:598:ALA   # break the core cross-knee lock hard

$PY pipelines/route_a/scripts/plot_snapback.py
echo "ALL DONE"
