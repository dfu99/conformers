#!/bin/bash
set -euo pipefail
#
# AVB3 Conformer Generation + Pseudo-AFM Image Pipeline
#
# Stages:
#   1. PULL  — RoyalMD steered MD on AVB3 (PACE GPU job)
#   2. SPLIT — Convert .nc trajectory → individual PDB frames
#   3. RELAX — RoyalMD equilibrium MD on selected frames (PACE GPU jobs)
#   4. COLLECT — Gather all relaxed PDB frames into one directory
#   5. AFM   — Compute CVs and generate pseudo-AFM images (CPU OK)
#
# Usage:
#   bash run_conformer_pipeline.sh <stage> [options]
#
# Examples:
#   bash run_conformer_pipeline.sh pull          # submit pull job to PACE
#   bash run_conformer_pipeline.sh split         # convert + split pull trajectory
#   bash run_conformer_pipeline.sh relax         # submit relax jobs for selected frames
#   bash run_conformer_pipeline.sh collect       # gather all PDB frames
#   bash run_conformer_pipeline.sh afm           # generate pseudo-AFM images
#   bash run_conformer_pipeline.sh all-local     # run split + collect + afm locally
#

# ── Paths ──────────────────────────────────────────────────────────────────────
CONFORMERS_ROOT="${CONFORMERS_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}"
ROYALMD_ROOT="${ROYALMD_ROOT:-$(cd "$CONFORMERS_ROOT/../RoyalMD" && pwd)}"
AFMFOLD_ROOT="${AFMFOLD_ROOT:-$(cd "$CONFORMERS_ROOT/../afmfold" && pwd)}"

PIPELINE_DIR="$CONFORMERS_ROOT/pipelines/avb3-conformers"
WORK_DIR="${WORK_DIR:-$CONFORMERS_ROOT/data/runs/avb3/conformers}"

# Input PDB for pull phase
INPUT_PDB="${INPUT_PDB:-$ROYALMD_ROOT/test_systems/AVB3_clean.pdb}"

# ── Pull parameters ────────────────────────────────────────────────────────────
PULL_FORCE_PN="${PULL_FORCE_PN:-2.0}"
PULL_PRODUCTION_TIME="${PULL_PRODUCTION_TIME:-1000}"  # ps
PULL_REPORT_INTERVAL="${PULL_REPORT_INTERVAL:-1500}"  # steps

# ── Split parameters ──────────────────────────────────────────────────────────
STRIDE="${STRIDE:-1}"
MAX_PULL_FRAME="${MAX_PULL_FRAME:-50}"  # don't keep frames beyond this (unrealistic)

# ── Relax parameters ──────────────────────────────────────────────────────────
# Comma-separated list of pull-trajectory frame indices to seed relaxation runs
RELAX_FRAMES="${RELAX_FRAMES:-10,20,30,40,50}"
RELAX_PRODUCTION_TIME="${RELAX_PRODUCTION_TIME:-1000}"  # ps
RELAX_REPORT_INTERVAL="${RELAX_REPORT_INTERVAL:-1500}"

# ── AFM parameters ────────────────────────────────────────────────────────────
PROTEIN_NAME="${PROTEIN_NAME:-avb3}"
AFM_WIDTH="${AFM_WIDTH:-35}"
AFM_HEIGHT="${AFM_HEIGHT:-35}"
AFM_RESOLUTION="${AFM_RESOLUTION:-0.98}"
AFM_EPOCHS="${AFM_EPOCHS:-1}"
AFM_DATASET_SIZE="${AFM_DATASET_SIZE:-1000}"
AFM_DEVICE="${AFM_DEVICE:-cpu}"

# ── Derived paths ──────────────────────────────────────────────────────────────
PULL_DIR="$WORK_DIR/pull"
PULL_FRAMES_DIR="$WORK_DIR/pull_frames"
RELAX_DIR="$WORK_DIR/relax"
COLLECT_DIR="$WORK_DIR/all_frames"
AFM_DIR="$WORK_DIR/pseudo_afm"

# ═══════════════════════════════════════════════════════════════════════════════
stage_pull() {
    echo "=== Stage 1: PULL (steered MD) ==="
    mkdir -p "$PULL_DIR"

    if [[ -f "$PULL_DIR/production.nc" ]]; then
        echo "Pull trajectory already exists at $PULL_DIR/production.nc"
        echo "Set WORK_DIR to a new directory or delete to rerun."
        return 0
    fi

    echo "Input PDB: $INPUT_PDB"
    echo "Force: ${PULL_FORCE_PN} pN, Production: ${PULL_PRODUCTION_TIME} ps"
    echo ""
    echo "Run on PACE with:"
    echo "  sbatch $PIPELINE_DIR/scripts/submit_pull_slurm.sh"
    echo ""
    echo "Or run directly (GPU required):"
    echo "  python $ROYALMD_ROOT/RoyalMD.py $INPUT_PDB \\"
    echo "    --out-dir $PULL_DIR \\"
    echo "    --pulling-enabled --pulling-force-pn $PULL_FORCE_PN \\"
    echo "    --production-time $PULL_PRODUCTION_TIME \\"
    echo "    --report-interval $PULL_REPORT_INTERVAL"
}

# ═══════════════════════════════════════════════════════════════════════════════
stage_split() {
    echo "=== Stage 2: SPLIT (trajectory → PDB frames) ==="

    if [[ ! -f "$PULL_DIR/production.nc" ]]; then
        echo "ERROR: Pull trajectory not found at $PULL_DIR/production.nc" >&2
        echo "Run the pull stage first." >&2
        return 1
    fi

    # Find topology file
    local topology=""
    for candidate in "$PULL_DIR/equilibrated.cif" "$PULL_DIR/solvated.pdb" "$PULL_DIR/minimized.cif"; do
        if [[ -f "$candidate" ]]; then
            topology="$candidate"
            break
        fi
    done
    if [[ -z "$topology" ]]; then
        echo "ERROR: No topology file found in $PULL_DIR" >&2
        return 1
    fi

    mkdir -p "$PULL_FRAMES_DIR"

    echo "Converting .nc → PDB trajectory (stride=$STRIDE)..."
    python3 "$ROYALMD_ROOT/scripts/visualize_nc.py" \
        "$PULL_DIR/production.nc" \
        --top "$topology" \
        --stride "$STRIDE"

    # visualize_nc.py writes trajectory_all.pdb in the .nc directory
    local traj_pdb="$PULL_DIR/trajectory_all.pdb"
    if [[ ! -f "$traj_pdb" ]]; then
        # fallback: might be named differently
        traj_pdb="$(ls "$PULL_DIR"/production_trajectory*.pdb 2>/dev/null | head -1 || true)"
    fi
    if [[ -z "$traj_pdb" || ! -f "$traj_pdb" ]]; then
        echo "ERROR: Could not find converted PDB trajectory in $PULL_DIR" >&2
        return 1
    fi

    echo "Splitting into individual frames (max frame: $MAX_PULL_FRAME)..."
    python3 "$ROYALMD_ROOT/scripts/split_pdb_trajectory.py" \
        "$traj_pdb" \
        --start 0 \
        --end "$MAX_PULL_FRAME" \
        --prefix pull

    # Move splits to our frames directory
    local split_src="$PULL_DIR/splits"
    if [[ -d "$split_src" ]]; then
        mv "$split_src"/pull_*.pdb "$PULL_FRAMES_DIR/" 2>/dev/null || true
        echo "Frames saved to $PULL_FRAMES_DIR/"
        echo "  $(ls "$PULL_FRAMES_DIR"/*.pdb 2>/dev/null | wc -l) frames extracted"
    else
        echo "ERROR: splits directory not found after split_pdb_trajectory.py" >&2
        return 1
    fi
}

# ═══════════════════════════════════════════════════════════════════════════════
stage_relax() {
    echo "=== Stage 3: RELAX (equilibrium MD on selected frames) ==="

    IFS=',' read -ra FRAMES <<< "$RELAX_FRAMES"
    echo "Relaxing frames: ${FRAMES[*]}"
    echo ""

    for frame_idx in "${FRAMES[@]}"; do
        local frame_pdb
        # Try to find the frame PDB (handle different naming conventions)
        frame_pdb=$(printf "%s/pull_pull_%03d.pdb" "$PULL_FRAMES_DIR" "$frame_idx")
        if [[ ! -f "$frame_pdb" ]]; then
            frame_pdb=$(printf "%s/pull_frame_%03d.pdb" "$PULL_FRAMES_DIR" "$frame_idx")
        fi
        if [[ ! -f "$frame_pdb" ]]; then
            # glob fallback
            frame_pdb=$(ls "$PULL_FRAMES_DIR"/*"$(printf '%03d' "$frame_idx")"*.pdb 2>/dev/null | head -1 || true)
        fi

        if [[ -z "$frame_pdb" || ! -f "$frame_pdb" ]]; then
            echo "WARNING: Frame $frame_idx PDB not found in $PULL_FRAMES_DIR, skipping."
            continue
        fi

        local relax_out="$RELAX_DIR/seed_$(printf '%03d' "$frame_idx")"
        if [[ -f "$relax_out/production.nc" ]]; then
            echo "  Frame $frame_idx: relaxation already exists, skipping."
            continue
        fi

        mkdir -p "$relax_out"
        echo "  Frame $frame_idx: $frame_pdb → $relax_out"
        echo "    Submit to PACE or run:"
        echo "    python $ROYALMD_ROOT/RoyalMD.py $frame_pdb \\"
        echo "      --out-dir $relax_out --no-pulling \\"
        echo "      --production-time $RELAX_PRODUCTION_TIME \\"
        echo "      --report-interval $RELAX_REPORT_INTERVAL"
    done

    echo ""
    echo "To submit all relax jobs to PACE at once:"
    echo "  RELAX_FRAMES=$RELAX_FRAMES bash $PIPELINE_DIR/scripts/submit_relax_slurm.sh"
}

# ═══════════════════════════════════════════════════════════════════════════════
stage_collect() {
    echo "=== Stage 4: COLLECT (gather all PDB frames) ==="
    mkdir -p "$COLLECT_DIR"

    local count=0

    # Collect pull frames (up to MAX_PULL_FRAME)
    if [[ -d "$PULL_FRAMES_DIR" ]]; then
        for f in "$PULL_FRAMES_DIR"/*.pdb; do
            [[ -f "$f" ]] || continue
            local base
            base=$(basename "$f")
            cp "$f" "$COLLECT_DIR/$base"
            count=$((count + 1))
        done
        echo "  Collected $count pull frames"
    fi

    # Collect relaxation frames
    local relax_count=0
    if [[ -d "$RELAX_DIR" ]]; then
        for seed_dir in "$RELAX_DIR"/seed_*; do
            [[ -d "$seed_dir" ]] || continue
            local seed_name
            seed_name=$(basename "$seed_dir")

            # Convert relaxation .nc → frames if not already done
            if [[ -f "$seed_dir/production.nc" && ! -d "$seed_dir/splits" ]]; then
                local topology=""
                for candidate in "$seed_dir/equilibrated.cif" "$seed_dir/solvated.pdb" "$seed_dir/minimized.cif"; do
                    if [[ -f "$candidate" ]]; then
                        topology="$candidate"
                        break
                    fi
                done
                if [[ -n "$topology" ]]; then
                    echo "  Converting $seed_name trajectory..."
                    python3 "$ROYALMD_ROOT/scripts/visualize_nc.py" \
                        "$seed_dir/production.nc" \
                        --top "$topology" \
                        --stride "$STRIDE"

                    local traj_pdb
                    traj_pdb=$(ls "$seed_dir"/trajectory_all.pdb "$seed_dir"/production_trajectory*.pdb 2>/dev/null | head -1 || true)
                    if [[ -n "$traj_pdb" && -f "$traj_pdb" ]]; then
                        python3 "$ROYALMD_ROOT/scripts/split_pdb_trajectory.py" \
                            "$traj_pdb" \
                            --prefix "${seed_name}"
                    fi
                fi
            fi

            # Copy split frames
            if [[ -d "$seed_dir/splits" ]]; then
                for f in "$seed_dir/splits"/*.pdb; do
                    [[ -f "$f" ]] || continue
                    local base
                    base=$(basename "$f")
                    cp "$f" "$COLLECT_DIR/${seed_name}_${base}"
                    relax_count=$((relax_count + 1))
                done
            fi
        done
        echo "  Collected $relax_count relaxation frames"
    fi

    local total
    total=$(ls "$COLLECT_DIR"/*.pdb 2>/dev/null | wc -l)
    echo "  Total frames in $COLLECT_DIR: $total"
}

# ═══════════════════════════════════════════════════════════════════════════════
stage_afm() {
    echo "=== Stage 5: AFM (compute CVs + generate pseudo-AFM images) ==="

    local frame_count
    frame_count=$(ls "$COLLECT_DIR"/*.pdb 2>/dev/null | wc -l)
    if [[ "$frame_count" -eq 0 ]]; then
        echo "ERROR: No PDB frames found in $COLLECT_DIR" >&2
        echo "Run the collect stage first." >&2
        return 1
    fi
    echo "  Found $frame_count PDB frames"

    mkdir -p "$AFM_DIR"

    # Run the processing script
    python3 "$PIPELINE_DIR/scripts/process_frames_to_afm.py" \
        --frames-dir "$COLLECT_DIR" \
        --output-dir "$AFM_DIR" \
        --afmfold-root "$AFMFOLD_ROOT" \
        --protein-name "$PROTEIN_NAME" \
        --width "$AFM_WIDTH" \
        --height "$AFM_HEIGHT" \
        --resolution-nm "$AFM_RESOLUTION" \
        --epochs "$AFM_EPOCHS" \
        --dataset-size "$AFM_DATASET_SIZE" \
        --device "$AFM_DEVICE"

    echo "  Pseudo-AFM images saved to $AFM_DIR/"
}

# ═══════════════════════════════════════════════════════════════════════════════
usage() {
    echo "Usage: $0 <stage>"
    echo ""
    echo "Stages:"
    echo "  pull       Submit steered MD job (prints PACE command)"
    echo "  split      Convert pull .nc → individual PDB frames"
    echo "  relax      Submit relaxation jobs for selected frames"
    echo "  collect    Gather all frames into one directory"
    echo "  afm        Compute CVs and generate pseudo-AFM images"
    echo "  all-local  Run split + collect + afm (post-PACE processing)"
    echo "  status     Show pipeline state"
    echo ""
    echo "Environment variables:"
    echo "  WORK_DIR             Working directory (default: data/runs/avb3/conformers)"
    echo "  INPUT_PDB            Input PDB for pull (default: RoyalMD/test_systems/AVB3_clean.pdb)"
    echo "  PULL_FORCE_PN        Pulling force in pN (default: 2.0)"
    echo "  RELAX_FRAMES         Comma-separated frame indices to relax (default: 10,20,30,40,50)"
    echo "  MAX_PULL_FRAME       Max pull frame to keep (default: 50)"
    echo "  PROTEIN_NAME         Protein name for domain CVs (default: avb3)"
    echo "  AFM_DEVICE           torch device for AFM gen (default: cpu)"
}

stage_status() {
    echo "=== Pipeline Status ==="
    echo "Work directory: $WORK_DIR"
    echo ""

    echo "Pull:"
    if [[ -f "$PULL_DIR/production.nc" ]]; then
        echo "  ✓ Trajectory exists"
    else
        echo "  ✗ Not run"
    fi

    echo "Split:"
    local pull_frame_count
    pull_frame_count=$(ls "$PULL_FRAMES_DIR"/*.pdb 2>/dev/null | wc -l)
    echo "  $pull_frame_count frames extracted"

    echo "Relax:"
    if [[ -d "$RELAX_DIR" ]]; then
        for seed_dir in "$RELAX_DIR"/seed_*; do
            [[ -d "$seed_dir" ]] || continue
            local name
            name=$(basename "$seed_dir")
            if [[ -f "$seed_dir/production.nc" ]]; then
                echo "  ✓ $name complete"
            else
                echo "  ✗ $name incomplete"
            fi
        done
    else
        echo "  Not started"
    fi

    echo "Collect:"
    local total
    total=$(ls "$COLLECT_DIR"/*.pdb 2>/dev/null | wc -l)
    echo "  $total frames collected"

    echo "AFM:"
    local img_count
    img_count=$(ls "$AFM_DIR"/image_*.npy 2>/dev/null | wc -l)
    echo "  $img_count image batches generated"
}

# ═══════════════════════════════════════════════════════════════════════════════
STAGE="${1:-}"
case "$STAGE" in
    pull)      stage_pull ;;
    split)     stage_split ;;
    relax)     stage_relax ;;
    collect)   stage_collect ;;
    afm)       stage_afm ;;
    all-local) stage_split && stage_collect && stage_afm ;;
    status)    stage_status ;;
    *)         usage; exit 1 ;;
esac
