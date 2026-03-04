#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_boltzgen_attempt.sh \
    --design-spec <path.yaml> \
    --outdir <dir> \
    [--extra "<extra args>"]

Notes:
- This is a thin launcher because BoltzGen CLI variants differ by install.
- It tries, in order:
  1) boltzgen run <spec> --out_dir <outdir>
  2) python -m boltzgen run <spec> --out_dir <outdir>
EOF
}

DESIGN_SPEC=""
OUTDIR=""
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --design-spec) DESIGN_SPEC="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
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
  echo "Running: boltzgen run ..."
  # shellcheck disable=SC2086
  boltzgen run "$DESIGN_SPEC" --out_dir "$OUTDIR" $EXTRA_ARGS
  exit $?
fi

if command -v python3 >/dev/null 2>&1; then
  echo "Running: python -m boltzgen run ..."
  # shellcheck disable=SC2086
  python3 -m boltzgen run "$DESIGN_SPEC" --out_dir "$OUTDIR" $EXTRA_ARGS
  exit $?
fi

echo "ERROR: could not find boltzgen executable or python3." >&2
exit 1
