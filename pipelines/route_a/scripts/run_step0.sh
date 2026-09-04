#!/usr/bin/env bash
# Step-0 external pair check, wrapped so it restores its own deps on the pod.
#
# WHY THIS WRAPPER EXISTS: the gpuq runner's default environment restore (derived
# from requirements.txt, bin/runpod-runner.sh:419) has never actually executed —
# the daemon has been up since 2026-08-24 and the feature landed 2026-08-26, so
# 0 of 52 job logs since then carry a '--- setup:' line. Until that daemon is
# restarted, a job on this queue gets whatever the pod happens to have. numpy is
# in the image; matplotlib and MDAnalysis are not.
#
# Audit first, install second, exit 0 when already correct — the same contract
# tasks/gpu-setup.sh has to meet. Once the runner is restarted this logic belongs
# in tasks/gpu-setup.sh and this wrapper can go away.
set -euo pipefail

need=()
python3 -c 'import matplotlib' 2>/dev/null || need+=(matplotlib)
python3 -c 'import MDAnalysis'  2>/dev/null || need+=(MDAnalysis)
python3 -c 'import numpy'       2>/dev/null || need+=(numpy)

if (( ${#need[@]} )); then
    echo "restoring: ${need[*]}"
    pip install -q --no-input "${need[@]}"
else
    echo "environment already correct"
fi

exec python3 pipelines/route_a/scripts/check_pairs_external.py "$@"
