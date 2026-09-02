#!/usr/bin/env bash
# Force-ramp suite (obj-083): hold the knee at a load, then release it and watch where it lets go.
#
# obj-082 ran these same three systems at F = 0 and they were indistinguishable, because bent is
# the ground state at zero load -- there was no barrier to cross. Under load the readout becomes
# F½ (the force at which the knee stops holding) rather than a 9° smear against ±10-17° scatter.
#
# Two design points that are NOT free choices, both from the 2026-08-17 audit:
#
#  1. RAMP DOWN, not up. The structure starts extended, so an up-ramp from 0 spends its low-force
#     half in exactly the obj-082 regime: the knee collapses unaided (half the drop by 3.8 ± 2.4
#     ns across the 9 zero-force replicates) before the load arrives, and the resulting "F½"
#     is the collapse clock times the ramp rate -- 12 ± 7 pN with no force applied, which is
#     indistinguishable from the literature answer. Equilibrating AT the holding load and ramping
#     down measures the release point instead, which is what F½ actually means.
#  2. REPLICATES. n = 1 is what made obj-082 uninterpretable: pooled within-genotype SD is 12.8°
#     of knee across its replicates, so a single trace per genotype cannot rank genotypes. Three
#     shorter runs beat one long one here -- the ordering, not the absolute pN, is the result.
#
# Each run aborts on its own if the requested pN does not reach the atoms (verify_pull) or if the
# knee ranges do not resolve to the expected residue counts (make_observables), so neither a
# silent null nor a wrong-hinge readout can get as far as producing a trajectory.
#
# Usage (from /workspace/route_a):  bash scripts/run_force_ramp.sh [ramp_ns] [hold_pN] [reps] [equil_ns]
set -euo pipefail
NS="${1:-8}"        # ramp length per run
MAXPN="${2:-60}"    # holding load: equilibrate here, then release to 0
REPS="${3:-3}"      # replicates per genotype
EQ="${4:-2}"        # equilibration under the holding load, before recording
RETRIES="${6:-3}"   # attempts per run before giving up on it
RETRY_WAIT="${7:-300}"
SUF="${5:-}"        # tag suffix. Set MAXPN to ~0 and SUF to _f0 for the matched no-load control
                    # arm: same code, same 2 ns pre-equilibration, same observables, no force.
ROOT=/workspace/route_a
PY=/workspace/envs/route_a/bin/python
FINAL="$ROOT/results/route_a/force_ramp"    # MooseFS network mount
LOCAL=/root/force_ramp_runs                 # container-local overlay disk
cd "$ROOT"
export CUDA_VISIBLE_DEVICES=0
mkdir -p "$LOCAL" "$FINAL"

# GPU preflight. The 2026-08-21 campaign died on its first run with
# CUDA_ERROR_MPS_SERVER_NOT_READY and `set -e` discarded all 20 GPU-h of queued work. Cause: the
# pod's shared MPS server had entered FAULT (`get_server_status` said so) and refused every new
# client, while the tenant job that attached before the fault kept running. Compute mode is
# Default, so MPS is optional -- pointing CUDA at a pipe directory with no daemon opens a normal
# context instead. Restarting the shared daemon is NOT an option: it would kill the other
# tenant's job. Prefer MPS when it works, fall back when it does not, never start blind.
gpu_ok() {
  CUDA_VISIBLE_DEVICES=0 "$PY" -c 'import openmm as mm
s = mm.System(); s.addParticle(1.0)
mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("CUDA"),
           {"Precision": "mixed"})' >/dev/null 2>&1
}
if gpu_ok; then
  echo "=== GPU preflight OK (MPS healthy) ==="
else
  export CUDA_MPS_PIPE_DIRECTORY=/tmp/no-mps-route-a
  if gpu_ok; then
    echo "=== GPU preflight: MPS unavailable, bypassing it (compute mode is Default) ==="
  else
    echo "=== GPU preflight FAILED: no usable CUDA context with or without MPS. Not starting. ==="
    exit 1
  fi
fi

TOTAL_NS=$(python3 -c "print(3*$REPS*($NS+$EQ))")
echo "=== ${REPS} reps x 3 genotypes x (${EQ} ns hold @ ${MAXPN} pN + ${NS} ns ramp -> 0) "
echo "=== ${TOTAL_NS} ns total, ~$(python3 -c "print(round(0.22*$TOTAL_NS,1))") GPU-h at the obj-082 rate, ~55 MB of frames per run"

# Same local-disk-then-publish pattern as run_snapback_replicates.sh: /workspace is a MooseFS
# mount that intermittently drops a write under cluster load, and a multi-hour run that flushes
# every report chunk gets hundreds of chances to hit it. One EIO kills the whole trajectory.
FAILED=()
run() { # tag  mutations
  if [[ -f "$FINAL/$1/summary.json" || -f "$LOCAL/$1/summary.json" ]]; then
    echo "=== $1 already done — skipping ==="; return 0
  fi
  echo "=== $1 ($2) — ${EQ} ns @ ${MAXPN} pN then ${NS} ns ramp ${MAXPN}->0 pN on A5000 ==="
  # Retry, and never let one run abort the campaign. The 2026-08-21 attempt died on its FIRST
  # run with CUDA_ERROR_MPS_SERVER_NOT_READY -- the pod's MPS daemon was mid-restart -- and
  # `set -e` threw away all 20 GPU-h of queued work over a hiccup that cleared by itself.
  local attempt=1
  while :; do
    rm -rf "${LOCAL:?}/$1"
    if "$PY" scripts/snapback_md.py \
          --pdb results/route_a/extended_state_b.pdb \
          --mutations "$2" --tag "$1" --ns "$NS" \
          --equil-ns "$EQ" --ramp-from-pn "$MAXPN" --ramp-pn 0 \
          --platform CUDA --report-ps 100 \
          --out-dir "$LOCAL"; then
      break
    fi
    if [[ $attempt -ge $RETRIES ]]; then
      echo "=== $1 FAILED after ${attempt} attempts — moving on to the next system ==="
      FAILED+=("$1"); return 0
    fi
    echo "=== $1 attempt ${attempt} failed; retrying in ${RETRY_WAIT}s ==="
    attempt=$((attempt + 1)); sleep "$RETRY_WAIT"
  done
  # copy the CONTENTS into the target dir: `cp -r src dst` nests into dst/$(basename src) when
  # dst already exists, which a half-finished earlier publish leaves behind.
  mkdir -p "$FINAL/$1"
  if cp -r "$LOCAL/$1/." "$FINAL/$1/" 2>/dev/null; then
    echo "=== $1 published to $FINAL/$1 ==="
  else
    echo "=== $1 WARNING: copy to $FINAL failed; result retained at $LOCAL/$1 ==="
  fi
}

# All three genotypes at replicate 1, then replicate 2, then 3 -- an interrupted campaign leaves
# a complete comparison at lower n rather than three runs of WT and nothing to compare them to.
for r in $(seq 1 "$REPS"); do
  sfx=""; [[ $r -gt 1 ]] && sfx="_r${r}"
  run "wt_ramp${SUF}${sfx}"    WT
  run "k459a_ramp${SUF}${sfx}" A:459:ALA
  run "e598a_ramp${SUF}${sfx}" A:598:ALA
done

if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "CAMPAIGN INCOMPLETE — these systems never produced a trajectory: ${FAILED[*]}"
  exit 1
fi
echo "ALL RAMPS DONE — analyse with: python scripts/analyze_force_ramp.py --dir $FINAL"
