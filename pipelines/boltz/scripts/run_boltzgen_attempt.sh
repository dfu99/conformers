#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF2'
Usage:
  run_boltzgen_attempt.sh \
    --design-spec <path.yaml> \
    --outdir <dir> \
    [--protocol protein-anything] \
    [--num-designs 10000] \
    [--budget 200] \
    [--extra "<extra args>"]

Notes:
- Uses official BoltzGen style:
  boltzgen run <spec.yaml> --output <dir> --protocol <name> --num_designs <N> --budget <K>
- If `boltzgen` executable is unavailable, falls back to `python3 -m boltzgen`.
EOF2
}

DESIGN_SPEC=""
OUTDIR=""
PROTOCOL="protein-anything"
NUM_DESIGNS=10000
BUDGET=200
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --design-spec) DESIGN_SPEC="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --protocol) PROTOCOL="$2"; shift 2 ;;
    --num-designs) NUM_DESIGNS="$2"; shift 2 ;;
    --budget) BUDGET="$2"; shift 2 ;;
    --extra) EXTRA_ARGS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$DESIGN_SPEC" || -z "$OUTDIR" ]]; then
  echo "Missing required args." >&2
  usage
  exit 1
fi
if [[ ! -f "$DESIGN_SPEC" ]]; then
  echo "Design spec not found: $DESIGN_SPEC" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

if command -v boltzgen >/dev/null 2>&1; then
  echo "Running: boltzgen run"
  # shellcheck disable=SC2086
  boltzgen run "$DESIGN_SPEC" \
    --output "$OUTDIR" \
    --protocol "$PROTOCOL" \
    --num_designs "$NUM_DESIGNS" \
    --budget "$BUDGET" \
    $EXTRA_ARGS
  exit $?
fi

if command -v python3 >/dev/null 2>&1; then
  echo "Running: python3 -m boltzgen run"
  # shellcheck disable=SC2086
  python3 -m boltzgen run "$DESIGN_SPEC" \
    --output "$OUTDIR" \
    --protocol "$PROTOCOL" \
    --num_designs "$NUM_DESIGNS" \
    --budget "$BUDGET" \
    $EXTRA_ARGS
  exit $?
fi

echo "ERROR: could not find boltzgen executable or python3." >&2
exit 1
