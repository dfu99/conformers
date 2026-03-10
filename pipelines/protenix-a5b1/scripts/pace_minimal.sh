#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../../.." && pwd)"

PACE_HOST="${PACE_HOST:-pace}"
PACE_REMOTE_ROOT="${PACE_REMOTE_ROOT:-}"
PACE_CONTROL_PERSIST="${PACE_CONTROL_PERSIST:-15m}"
PACE_CONTROL_PATH="${PACE_CONTROL_PATH:-$HOME/.ssh/cm-%r@%h:%p}"
PACE_GIT_PULL="${PACE_GIT_PULL:-1}"
PACE_ACCOUNT="${PACE_ACCOUNT:-gts-yke8}"

DEFAULT_SUBMIT_SCRIPT="pipelines/protenix-a5b1/scripts/submit_complete_tagged_pipeline_slurm.sh"

DEFAULT_REMOTE_LOG_DIR_REL="logs/protenix-a5b1"
DEFAULT_LOG_PREFIX="a5b1_tagged_full"
DEFAULT_RESULTS_REL="data/runs/a5b1/staged_attachment/outputs/final"

DEFAULT_LOCAL_LOG_DIR="$REPO_ROOT/logs/protenix-a5b1"
DEFAULT_LOCAL_RESULTS_DIR="$REPO_ROOT/data/runs/a5b1/staged_attachment/outputs/final"

SMOKE_LOG_DIR_REL="logs/pace-smoke"
SMOKE_RESULTS_REL="data/runs/pace-smoke"
SMOKE_LOG_PREFIX="pace_smoke"

usage() {
  cat <<'EOF'
Usage:
  pace_minimal.sh check
  pace_minimal.sh submit [submit_script_rel]
  pace_minimal.sh watch <job_id> [poll_seconds]
  pace_minimal.sh fetch <job_id> [remote_log_dir_or_rel] [log_prefix] [remote_results_rel] [local_log_dir] [local_results_dir]
  pace_minimal.sh smoke [poll_seconds] [sleep_seconds]

Environment:
  PACE_HOST             SSH host alias for PACE (default: pace)
  PACE_REMOTE_ROOT      Remote conformers root (default on remote: $HOME/scratch/conformers)
  PACE_ACCOUNT          SLURM account for smoke submits (default: gts-yke8)
  PACE_GIT_PULL         If 1, run git pull --ff-only before submit (default: 1)
  PACE_CONTROL_PERSIST  SSH ControlPersist value (default: 15m)
  PACE_CONTROL_PATH     SSH ControlPath (default: ~/.ssh/cm-%r@%h:%p)
EOF
}

q() {
  printf "%q" "$1"
}

ssh_opts=(
  -o BatchMode=yes
  -o ControlMaster=auto
  -o "ControlPersist=$PACE_CONTROL_PERSIST"
  -o "ControlPath=$PACE_CONTROL_PATH"
)

run_ssh() {
  ssh "${ssh_opts[@]}" "$PACE_HOST" "$@"
}

rsync_ssh_cmd() {
  printf "ssh -o BatchMode=yes -o ControlMaster=auto -o ControlPersist=%s -o ControlPath=%s" \
    "$PACE_CONTROL_PERSIST" "$PACE_CONTROL_PATH"
}

ensure_control_socket_dir() {
  local control_dir
  control_dir="$(dirname -- "$PACE_CONTROL_PATH")"
  mkdir -p "$control_dir"
}

extract_job_id() {
  local raw="$1"
  local line
  line="$(printf "%s\n" "$raw" | tail -n 1)"
  line="${line%%;*}"
  if [[ "$line" =~ ^[0-9]+$ ]]; then
    printf "%s" "$line"
    return 0
  fi
  return 1
}

resolve_remote_root() {
  run_ssh \
    "PACE_REMOTE_ROOT=$(q "$PACE_REMOTE_ROOT") bash -lc 'set -euo pipefail; printf \"%s\" \"\${PACE_REMOTE_ROOT:-\$HOME/scratch/conformers}\"'"
}

check_remote() {
  ensure_control_socket_dir
  run_ssh \
    "PACE_REMOTE_ROOT=$(q "$PACE_REMOTE_ROOT") bash -lc 'set -euo pipefail
remote_root=\"\${PACE_REMOTE_ROOT:-\$HOME/scratch/conformers}\"
printf \"host=%s\\n\" \"\$(hostname)\"
printf \"remote_root=%s\\n\" \"\$remote_root\"
test -d \"\$remote_root\"
command -v sbatch >/dev/null
command -v squeue >/dev/null
command -v sacct >/dev/null
echo \"status=ok\"'"
}

submit_job() {
  local submit_script_rel="${1:-$DEFAULT_SUBMIT_SCRIPT}"
  local raw
  local job_id

  ensure_control_socket_dir

  raw="$(
    run_ssh \
      "PACE_REMOTE_ROOT=$(q "$PACE_REMOTE_ROOT") PACE_GIT_PULL=$(q "$PACE_GIT_PULL") SUBMIT_SCRIPT_REL=$(q "$submit_script_rel") bash -lc 'set -euo pipefail
remote_root=\"\${PACE_REMOTE_ROOT:-\$HOME/scratch/conformers}\"
cd \"$remote_root\"
if [[ \"$PACE_GIT_PULL\" == \"1\" ]]; then
  git pull --ff-only
fi
sbatch --parsable \"$SUBMIT_SCRIPT_REL\"'"
  )"

  if ! job_id="$(extract_job_id "$raw")"; then
    printf "ERROR: could not parse job id from submit output.\n" >&2
    printf "%s\n" "$raw" >&2
    exit 1
  fi

  mkdir -p "$DEFAULT_LOCAL_LOG_DIR"
  printf "%s\n" "$job_id" > "$DEFAULT_LOCAL_LOG_DIR/.last_job_id"
  printf "%s\n" "$raw"
  printf "JOB_ID=%s\n" "$job_id"
}

watch_job() {
  local job_id="${1:-}"
  local poll_seconds="${2:-20}"

  if [[ -z "$job_id" ]]; then
    printf "ERROR: watch requires <job_id>\n" >&2
    usage
    exit 1
  fi

  ensure_control_socket_dir

  run_ssh \
    "JOB_ID=$(q "$job_id") POLL_SECONDS=$(q "$poll_seconds") bash -lc 'set -euo pipefail
while true; do
  ts=\$(date \"+%Y-%m-%dT%H:%M:%S%z\")
  line=\$(squeue -h -j \"\$JOB_ID\" -o \"%i|%t|%M|%D|%R\" || true)
  if [[ -n \"\$line\" ]]; then
    printf \"%s|SQUEUE|%s\\n\" \"\$ts\" \"\$line\"
    sleep \"\$POLL_SECONDS\"
    continue
  fi

  final=\$(sacct -n -P -j \"\$JOB_ID\" --format=JobIDRaw,State,ExitCode,Elapsed,End | awk -F \"|\" -v id=\"\$JOB_ID\" \"\\\$1==id {print; exit}\")
  if [[ -n \"\$final\" ]]; then
    printf \"%s|SACCT|%s\\n\" \"\$ts\" \"\$final\"
    break
  fi

  printf \"%s|WAIT|No squeue/sacct record yet\\n\" \"\$ts\"
  sleep \"\$POLL_SECONDS\"
done'"
}

fetch_artifacts() {
  local job_id="${1:-}"
  local remote_log_dir_input="${2:-$DEFAULT_REMOTE_LOG_DIR_REL}"
  local log_prefix="${3:-$DEFAULT_LOG_PREFIX}"
  local remote_results_rel="${4:-$DEFAULT_RESULTS_REL}"
  local local_log_dir="${5:-$DEFAULT_LOCAL_LOG_DIR}"
  local local_results_dir="${6:-$DEFAULT_LOCAL_RESULTS_DIR}"

  local remote_root
  local remote_log_dir
  local remote_results_path
  local rsync_ssh
  local remote_file
  local ext
  local pulled_any=0

  if [[ -z "$job_id" ]]; then
    printf "ERROR: fetch requires <job_id>\n" >&2
    usage
    exit 1
  fi

  ensure_control_socket_dir
  mkdir -p "$local_log_dir" "$local_results_dir"
  rsync_ssh="$(rsync_ssh_cmd)"
  remote_root="$(resolve_remote_root)"

  if [[ "$remote_log_dir_input" = /* ]]; then
    remote_log_dir="$remote_log_dir_input"
  else
    remote_log_dir="$remote_root/$remote_log_dir_input"
  fi

  for ext in log err out; do
    remote_file="$remote_log_dir/${log_prefix}_${job_id}.${ext}"
    if run_ssh "test -f $(q "$remote_file")"; then
      rsync -av -e "$rsync_ssh" "$PACE_HOST:$remote_file" "$local_log_dir/"
      pulled_any=1
    else
      printf "Missing remote log: %s\n" "$remote_file"
    fi
  done

  if [[ "$remote_results_rel" = /* ]]; then
    remote_results_path="$remote_results_rel"
  else
    remote_results_path="$remote_root/$remote_results_rel"
  fi

  if run_ssh "test -d $(q "$remote_results_path")"; then
    rsync -av -e "$rsync_ssh" "$PACE_HOST:$remote_results_path/" "$local_results_dir/"
    pulled_any=1
  else
    printf "Missing remote results directory: %s\n" "$remote_results_path"
  fi

  if [[ "$pulled_any" -eq 0 ]]; then
    printf "No artifacts were copied for job %s.\n" "$job_id"
    exit 1
  fi

  printf "Fetched artifacts for job %s.\n" "$job_id"
  printf "  Local logs:    %s\n" "$local_log_dir"
  printf "  Local results: %s\n" "$local_results_dir"
}

smoke_check() {
  local poll_seconds="${1:-10}"
  local sleep_seconds="${2:-60}"
  local raw
  local job_id
  local remote_root
  local remote_log_dir
  local local_log_dir
  local local_results_dir

  ensure_control_socket_dir

  raw="$(
    run_ssh \
      "PACE_REMOTE_ROOT=$(q "$PACE_REMOTE_ROOT") PACE_ACCOUNT=$(q "$PACE_ACCOUNT") bash -lc 'set -euo pipefail
remote_root=\"\${PACE_REMOTE_ROOT:-\$HOME/scratch/conformers}\"
log_dir=\"\$remote_root/$SMOKE_LOG_DIR_REL\"
result_dir=\"\$remote_root/$SMOKE_RESULTS_REL\"
mkdir -p \"\$log_dir\" \"\$result_dir\"
account=\"\${PACE_ACCOUNT:-gts-yke8}\"

wrap_cmd=\"set -euo pipefail; mkdir -p \\\"\$result_dir\\\"; result=\\\"\$result_dir/smoke_\\\${SLURM_JOB_ID}.txt\\\"; echo \\\"hello world from pace smoke\\\"; echo \\\"hello=world\\\" > \\\"\\\$result\\\"; echo \\\"job_id=\\\${SLURM_JOB_ID}\\\" >> \\\"\\\$result\\\"; echo \\\"host=\\\$(hostname)\\\" >> \\\"\\\$result\\\"; date -Is >> \\\"\\\$result\\\"; sleep $sleep_seconds\"
sbatch --parsable -A \"\$account\" --job-name=pace_smoke --time=00:03:00 --mem=1G --output \"\$log_dir/${SMOKE_LOG_PREFIX}_%j.log\" --error \"\$log_dir/${SMOKE_LOG_PREFIX}_%j.err\" --wrap \"\$wrap_cmd\"'"
  )"

  if ! job_id="$(extract_job_id "$raw")"; then
    printf "ERROR: could not parse job id from smoke submit output.\n" >&2
    printf "%s\n" "$raw" >&2
    exit 1
  fi

  printf "Smoke job submitted: %s\n" "$job_id"
  watch_job "$job_id" "$poll_seconds"

  remote_root="$(resolve_remote_root)"
  remote_log_dir="$remote_root/$SMOKE_LOG_DIR_REL"
  local_log_dir="$REPO_ROOT/logs/pace-smoke"
  local_results_dir="$REPO_ROOT/data/runs/pace-smoke"

  fetch_artifacts \
    "$job_id" \
    "$remote_log_dir" \
    "$SMOKE_LOG_PREFIX" \
    "$SMOKE_RESULTS_REL" \
    "$local_log_dir" \
    "$local_results_dir"
}

main() {
  local cmd="${1:-}"
  shift || true

  case "$cmd" in
    check)
      check_remote
      ;;
    submit)
      submit_job "${1:-$DEFAULT_SUBMIT_SCRIPT}"
      ;;
    watch)
      watch_job "${1:-}" "${2:-20}"
      ;;
    fetch)
      fetch_artifacts \
        "${1:-}" \
        "${2:-$DEFAULT_REMOTE_LOG_DIR_REL}" \
        "${3:-$DEFAULT_LOG_PREFIX}" \
        "${4:-$DEFAULT_RESULTS_REL}" \
        "${5:-$DEFAULT_LOCAL_LOG_DIR}" \
        "${6:-$DEFAULT_LOCAL_RESULTS_DIR}"
      ;;
    smoke)
      smoke_check "${1:-10}" "${2:-60}"
      ;;
    -h|--help|help)
      usage
      ;;
    *)
      usage >&2
      exit 1
      ;;
  esac
}

main "$@"
