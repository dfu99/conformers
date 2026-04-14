#!/usr/bin/env python3
"""Reprocess fitted coordinates with head-anchored flip resolution.

Reads raw fitted coords and head positions from a previous run,
applies the improved flip resolution that anchors to tracked AFM head
positions (instead of previous-frame comparison which drifts), and
saves corrected smooth coords. No GPU or re-fitting needed.

Usage:
  python reprocess_smooth_coords.py \
    --fitted-dir results/afm_pipeline/video2_v4 \
    --head-residues-a 1-440 --head-residues-b 1-350
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fitted-dir", type=Path, required=True)
    p.add_argument("--head-residues-a", type=str, default="1-440")
    p.add_argument("--head-residues-b", type=str, default="1-350")
    p.add_argument("--diagnostic-frames", type=str, default=None,
                   help="Comma-separated frame numbers for diagnostic plot")
    return p.parse_args()


def main():
    args = parse_args()
    d = args.fitted_dir

    # Load existing data
    coords_raw = np.load(str(d / "fitted_coords.npy"))
    rotations_raw = np.load(str(d / "fitted_rotations.npy"))
    head_px = np.load(str(d / "head_positions_px.npy"))
    correlations = np.load(str(d / "fitted_correlations.npy"))

    with open(d / "fitting_metadata.json") as f:
        meta = json.load(f)
    skip = meta["skip_frames"]
    resolution_nm = meta.get("tip_radius_nm", 2.0)  # fallback
    # Resolution is 0.98 nm/px for this data
    resolution_nm = 0.98
    grid_size = 35

    n_frames = len(coords_raw)
    print(f"Loaded {n_frames} frames (skip={skip})")
    print(f"Coords shape: {coords_raw.shape}")

    # Load old smooth for comparison
    old_smooth_path = d / "fitted_coords_smooth.npy"
    old_smooth = np.load(str(old_smooth_path)) if old_smooth_path.exists() else None

    # Load topology and get head indices
    ref = md.load(str(d / "topology.pdb"))
    a_start, a_end = map(int, args.head_residues_a.split("-"))
    b_start, b_end = map(int, args.head_residues_b.split("-"))
    head_indices = ref.topology.select(
        f"(chainid 0 and resSeq {a_start} to {a_end}) or "
        f"(chainid 1 and resSeq {b_start} to {b_end})"
    )
    print(f"Head domain: {len(head_indices)} atoms")

    # Convert head pixel positions to nm (same as fit_with_head_tracking.py)
    half_w = grid_size / 2.0
    half_h = grid_size / 2.0
    head_positions_nm = np.zeros((n_frames, 3), dtype=np.float32)
    head_positions_nm[:, 0] = (head_px[:, 1] - half_w) * resolution_nm
    head_positions_nm[:, 1] = (head_px[:, 0] - half_h) * resolution_nm
    head_positions_nm[:, 2] = 0.0

    # Import the improved flip resolution
    from fit_with_head_tracking import resolve_flips_head_anchored

    print("\nRunning head-anchored flip resolution...")
    smooth_rot, smooth_coords = resolve_flips_head_anchored(
        rotations_raw, coords_raw, head_indices, head_positions_nm
    )

    # Compare old vs new
    if old_smooth is not None:
        old_head_xy = old_smooth[:, head_indices, :2].mean(axis=1)
        new_head_xy = smooth_coords[:, head_indices, :2].mean(axis=1)
        raw_head_xy = coords_raw[:, head_indices, :2].mean(axis=1)

        old_drift = np.linalg.norm(old_head_xy - raw_head_xy, axis=1)
        new_drift = np.linalg.norm(new_head_xy - raw_head_xy, axis=1)

        print(f"\nHead drift from raw fitted positions:")
        print(f"  Old smooth: mean={old_drift.mean():.2f}nm, "
              f"max={old_drift.max():.2f}nm, >1nm: {(old_drift > 1).sum()} frames")
        print(f"  New smooth: mean={new_drift.mean():.2f}nm, "
              f"max={new_drift.max():.2f}nm, >1nm: {(new_drift > 1).sum()} frames")

    # Backup old and save new
    if old_smooth_path.exists():
        backup = d / "fitted_coords_smooth_old.npy"
        np.save(str(backup), np.load(str(old_smooth_path)))
        print(f"\nBacked up old smooth coords to {backup}")

    old_rot_path = d / "fitted_rotations_smooth.npy"
    if old_rot_path.exists():
        np.save(str(d / "fitted_rotations_smooth_old.npy"),
                np.load(str(old_rot_path)))

    np.save(str(d / "fitted_coords_smooth.npy"), smooth_coords)
    np.save(str(d / "fitted_rotations_smooth.npy"), smooth_rot)
    print(f"Saved corrected smooth coords to {d}")

    # Diagnostic plot
    if args.diagnostic_frames:
        diag_frames = [int(x) for x in args.diagnostic_frames.split(",")]
    else:
        diag_frames = None

    if diag_frames and old_smooth is not None:
        _plot_diagnostic(d, coords_raw, old_smooth, smooth_coords,
                         head_indices, head_positions_nm, correlations,
                         diag_frames, skip, grid_size, resolution_nm)


def _plot_diagnostic(output_dir, raw, old_smooth, new_smooth,
                     head_indices, head_nm, correlations,
                     frames, skip, grid_size, resolution_nm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(frames)
    fig, axes = plt.subplots(n, 3, figsize=(12, 3.0 * n))
    if n == 1:
        axes = axes[np.newaxis, :]

    half = grid_size / 2.0

    for row, f in enumerate(frames):
        idx = f - skip
        if idx < 0 or idx >= len(raw):
            continue

        for col, (label, coords) in enumerate([
            ("Raw (fitted)", raw[idx]),
            ("Old smooth", old_smooth[idx]),
            ("New smooth\n(head-anchored)", new_smooth[idx]),
        ]):
            ax = axes[row, col]
            # Plot all atoms XY
            xy = coords[:, :2]
            cols_px = xy[:, 0] / resolution_nm + half
            rows_px = xy[:, 1] / resolution_nm + half

            # Head vs tail
            all_idx = np.arange(len(coords))
            tail_idx = np.setdiff1d(all_idx, head_indices)

            ax.scatter(cols_px[tail_idx], rows_px[tail_idx], s=0.3, alpha=0.3,
                       c="gray", label="Tail")
            ax.scatter(cols_px[head_indices], rows_px[head_indices], s=0.3,
                       alpha=0.3, c="steelblue", label="Head")

            # Head centroid
            hc = coords[head_indices].mean(axis=0)
            hc_px = hc[:2] / resolution_nm + half
            ax.plot(hc_px[0], hc_px[1], "r+", markersize=10, markeredgewidth=2)

            # Tracked AFM head position
            target = head_nm[idx, :2] / resolution_nm + half
            ax.plot(target[0], target[1], "gx", markersize=10, markeredgewidth=2)

            ax.set_xlim(0, grid_size)
            ax.set_ylim(grid_size, 0)
            ax.set_aspect("equal")
            ax.grid(True, alpha=0.2)

            if row == 0:
                ax.set_title(label, fontsize=10, fontweight="bold")
            if col == 0:
                cc = correlations[idx]
                ax.set_ylabel(f"Frame {f}\n(r={cc:.3f})", fontsize=9)

    fig.suptitle("Flip Resolution Diagnostic\n"
                 "Red +: PDB head COM  |  Green x: tracked AFM head",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = output_dir / "flip_resolution_diagnostic.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved diagnostic: {out}")

    # Also plot head drift time series
    fig, ax = plt.subplots(figsize=(12, 4))
    raw_head = raw[:, head_indices, :2].mean(axis=1)
    old_head = old_smooth[:, head_indices, :2].mean(axis=1)
    new_head = new_smooth[:, head_indices, :2].mean(axis=1)
    target_xy = head_nm[:, :2]

    old_drift = np.linalg.norm(old_head - target_xy, axis=1)
    new_drift = np.linalg.norm(new_head - target_xy, axis=1)
    raw_drift = np.linalg.norm(raw_head - target_xy, axis=1)

    frame_nums = np.arange(len(raw)) + skip
    ax.plot(frame_nums, old_drift, "r-", alpha=0.6, linewidth=0.8,
            label=f"Old smooth (mean={old_drift.mean():.2f}nm)")
    ax.plot(frame_nums, new_drift, "g-", alpha=0.6, linewidth=0.8,
            label=f"New head-anchored (mean={new_drift.mean():.2f}nm)")
    ax.plot(frame_nums, raw_drift, "b-", alpha=0.3, linewidth=0.5,
            label=f"Raw fitted (mean={raw_drift.mean():.2f}nm)")

    for f in frames:
        ax.axvline(f, color="k", alpha=0.3, linestyle="--", linewidth=0.8)

    ax.set_xlabel("Frame")
    ax.set_ylabel("Head drift from AFM target (nm)")
    ax.set_title("Head Position Drift: Old vs New Smooth Coordinates")
    ax.legend(fontsize=9)
    ax.set_ylim(0, None)
    fig.tight_layout()
    out2 = output_dir / "head_drift_comparison.png"
    fig.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved drift comparison: {out2}")


if __name__ == "__main__":
    main()
