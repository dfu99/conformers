#!/usr/bin/env bash
set -euo pipefail

# Helper for provisioning Protenix template mmcif database on cluster.
# Preferred for large datasets: symlink an existing mmcif directory.
#
# Usage examples:
#   ./download_mmcif.sh --root "$HOME/scratch/Protenix" --link-from "/path/to/existing/mmcif"
#   ./download_mmcif.sh --root "$HOME/scratch/Protenix" --download

ROOT_DIR="${HOME}/scratch/Protenix"
MMCIF_DIR=""
MODE=""

usage() {
  cat <<'EOF'
Usage:
  download_mmcif.sh --root <protenix_root> --link-from <existing_mmcif_dir>
  download_mmcif.sh --root <protenix_root> --download

Options:
  --root <path>       Protenix runtime root (will host mmcif at <root>/mmcif)
  --link-from <path>  Existing mmcif directory to symlink into <root>/mmcif
  --download          Download mmcif snapshot via rsync into <root>/mmcif
  -h, --help          Show this help text
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)
      ROOT_DIR="$2"
      shift 2
      ;;
    --link-from)
      MMCIF_DIR="$2"
      MODE="link"
      shift 2
      ;;
    --download)
      MODE="download"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$MODE" ]]; then
  echo "ERROR: Must choose either --link-from or --download" >&2
  usage
  exit 1
fi

TARGET="${ROOT_DIR}/mmcif"
mkdir -p "${ROOT_DIR}"

if [[ "$MODE" == "link" ]]; then
  if [[ -z "$MMCIF_DIR" || ! -d "$MMCIF_DIR" ]]; then
    echo "ERROR: --link-from directory does not exist: ${MMCIF_DIR}" >&2
    exit 1
  fi
  ln -sfn "$MMCIF_DIR" "$TARGET"
  echo "Linked mmcif:"
  echo "  ${TARGET} -> ${MMCIF_DIR}"
  exit 0
fi

# Download mode
if ! command -v rsync >/dev/null 2>&1; then
  echo "ERROR: rsync is required for --download mode." >&2
  exit 1
fi

mkdir -p "$TARGET"
echo "Downloading mmcif into ${TARGET} (this is very large)..."
rsync -avz --delete \
  --info=progress2 \
  rsync://rsync.rcsb.org/pub/pdb/data/structures/divided/mmCIF/ \
  "${TARGET}/"

echo "Download complete: ${TARGET}"
