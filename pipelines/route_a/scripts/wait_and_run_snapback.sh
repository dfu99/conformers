#!/usr/bin/env bash
# Wait for the A5000 slot to open (oxDNA FFS jobs fully gone), THEN run the snap-back
# suite. Conservative by design: requires oxDNA to be absent for NEED consecutive minutes
# before taking the GPU, so it can never contend with a running or between-runs campaign.
# Launch under setsid so it survives ssh disconnect:
#   setsid bash /workspace/route_a/scripts/wait_and_run_snapback.sh </dev/null >/dev/null 2>&1 &
set -u
ROOT=/workspace/route_a
LOG="$ROOT/results/route_a/snapback/watcher.log"
NEED="${1:-30}"          # consecutive oxDNA-free minutes required before launching
mkdir -p "$ROOT/results/route_a/snapback"
idle=0
echo "$(date -u) watcher up; will launch after ${NEED} min with no oxDNA" >> "$LOG"
while true; do
  if pgrep -f "[o]xDNA" >/dev/null 2>&1; then
    if [ "$idle" -ne 0 ]; then echo "$(date -u) oxDNA active again — reset" >> "$LOG"; fi
    idle=0
  else
    idle=$((idle + 1))
    if [ $((idle % 10)) -eq 0 ] || [ "$idle" -ge "$NEED" ]; then
      echo "$(date -u) oxDNA-free ${idle}/${NEED} min" >> "$LOG"
    fi
  fi
  if [ "$idle" -ge "$NEED" ]; then
    # final confirmation the GPU has no oxDNA compute app
    if nvidia-smi --query-compute-apps=process_name --format=csv,noheader 2>/dev/null | grep -qi oxDNA; then
      echo "$(date -u) oxDNA back on GPU at launch — reset" >> "$LOG"; idle=0; sleep 60; continue
    fi
    echo "$(date -u) === A5000 slot free: launching snap-back suite ===" >> "$LOG"
    cd "$ROOT"
    CUDA_VISIBLE_DEVICES=0 bash "$ROOT/scripts/run_snapback_suite_pod.sh" 20 >> "$LOG" 2>&1
    echo "$(date -u) === snap-back suite finished ===" >> "$LOG"
    break
  fi
  sleep 60
done
