#!/usr/bin/env bash
# Find the extended conformation of every integrin variant, on the RunPod A5000 pod.
#
# Generalises the alphaVbeta3 route-A recipe (obj-074..077) to all full-ectodomain
# integrin heterodimers in the PDB. Every stage is CPU-only, so this never touches
# the GPU and cannot contend with the oxDNA FFS campaign that holds the A5000.
#
# Stages:
#   1 discover   query RCSB for integrin beta-subunit entities, verify each hit
#                against its own coordinates, keep only full ectodomains (legs
#                present) -- never trust a recalled PDB ID (tasks/lessons.md, obj-070)
#   2 survey     emit every accepted structure, transfer the validated alphaVbeta3
#                domain split onto each by per-module local alignment, and measure
#                the genu angle -- the validated primary reaction coordinate
#   3 select     per variant, pick the most-bent SUFFICIENTLY COMPLETE structure as
#                the morph start; flag any whose leg-module alignment is too weak to
#                trust; register experimentally extended structures as validation targets
#   4 morph      incremental genu-hinge swing + vacuum minimisation to the extended
#                endpoint, emitting every intermediate frame (string-method initial path)
#   5 summarise  acceptance checks, alphaVbeta3 positive control, comparison against
#                experimental extended cryo-EM structures, figure + table
#
# Usage (from /workspace/route_a on the pod):
#   bash scripts/run_all_variants.sh [threads_per_variant]
set -euo pipefail

PY=/workspace/envs/route_a/bin/python
THREADS="${1:-16}"
RES=results/variants
mkdir -p "$RES" logs/variants

echo "=== 1/5 discover full-ectodomain integrin structures ==="
$PY -u scripts/discover_integrin_ectodomains.py \
    --pdb-dir inputs/variants \
    --out-json "$RES/discovered_ectodomains.json"

echo "=== 2/5 survey: transfer boundaries onto every candidate ==="
$PY -u scripts/select_variant_endpoints.py \
    --discovered "$RES/discovered_ectodomains.json" \
    --all --out-json "$RES/endpoints_all.json"
$PY -u scripts/map_domain_boundaries.py \
    --discovered "$RES/endpoints_all.json" \
    --out-json "$RES/boundaries_all.json" \
    --out-dir inputs/variants/prepared

echo "=== 3/5 select bent starts by measured genu angle ==="
$PY -u scripts/select_variant_endpoints.py \
    --discovered "$RES/discovered_ectodomains.json" \
    --from-boundaries "$RES/boundaries_all.json" \
    --out-json "$RES/variant_endpoints.json"

echo "=== 4/5 morph each variant bent -> extended (parallel, CPU-only) ==="
keys=$($PY -c "
import json
ep=json.load(open('$RES/variant_endpoints.json'))
for v,r in sorted(ep['representatives'].items()):
    print(r['survey_key'], v)
")
while read -r key label; do
    [ -z "$key" ] && continue
    setsid nohup $PY -u scripts/morph_extend_generic.py \
        --variant "$key" --label "$label" \
        --boundaries "$RES/boundaries_all.json" \
        --out-dir "$RES" --threads "$THREADS" \
        > "logs/variants/${label}.log" 2>&1 < /dev/null &
    echo "  launched $label ($key)"
done <<< "$keys"

while [ "$(pgrep -fc morph_extend_generic || true)" != "0" ]; do sleep 30; done
echo "  all morphs finished"

echo "=== 5/5 summarise + validate ==="
$PY -u scripts/summarize_variant_extension.py \
    --endpoints "$RES/variant_endpoints.json" \
    --boundaries "$RES/boundaries_all.json" \
    --results-dir "$RES" \
    --out-json "$RES/extension_summary.json" \
    --out-md "$RES/extension_summary.md" \
    --out-png figures/variant_extension_summary.png

echo "done -- see $RES/extension_summary.md"
