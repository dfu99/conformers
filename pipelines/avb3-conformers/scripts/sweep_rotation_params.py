#!/usr/bin/env python3
"""Sweep rotation fitting parameters to minimize worst-case frame mismatch.

Pre-computes per-frame rotation scores at fine angular resolution, then
evaluates many (max_delta, smooth_sigma) combinations. Reports the combo
that maximizes the minimum per-frame alignment score (minimax).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from scipy.ndimage import gaussian_filter1d

# Import functions from render script
spec = importlib.util.spec_from_file_location(
    "rao", Path(__file__).parent / "render_atomistic_overlay.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def score_trajectory(all_scores: np.ndarray, angles_deg: np.ndarray,
                     max_delta: float, smooth_sigma: float,
                     initial_angle: float) -> tuple[np.ndarray, np.ndarray]:
    """Run forward propagation + smoothing and return (smoothed_angles, per_frame_scores)."""
    n = all_scores.shape[0]
    n_angles = len(angles_deg)
    angle_step = angles_deg[1] - angles_deg[0]

    # Forward propagation
    raw_angles = np.zeros(n, dtype=np.float32)
    prev = initial_angle

    for i in range(n):
        scores = all_scores[i]
        if scores.max() <= 0:
            raw_angles[i] = prev
            continue
        diffs = np.abs((angles_deg - prev + 180) % 360 - 180)
        allowed = diffs <= max_delta
        if not allowed.any():
            raw_angles[i] = prev
            continue
        masked = np.where(allowed, scores, -1.0)
        raw_angles[i] = angles_deg[int(np.argmax(masked))]
        prev = raw_angles[i]

    # Smoothing
    if smooth_sigma > 0:
        sin_t = np.sin(np.deg2rad(raw_angles))
        cos_t = np.cos(np.deg2rad(raw_angles))
        sin_s = gaussian_filter1d(sin_t, sigma=smooth_sigma)
        cos_s = gaussian_filter1d(cos_t, sigma=smooth_sigma)
        smoothed = np.rad2deg(np.arctan2(sin_s, cos_s)) % 360
    else:
        smoothed = raw_angles.copy()

    # Evaluate: interpolate score at smoothed angle for each frame
    per_frame_scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        if all_scores[i].max() <= 0:
            per_frame_scores[i] = 0.0
            continue
        # Find nearest scored angle
        diffs = np.abs((angles_deg - smoothed[i] + 180) % 360 - 180)
        nearest_idx = int(np.argmin(diffs))
        per_frame_scores[i] = all_scores[i, nearest_idx]

    return smoothed, per_frame_scores


def main():
    base = Path("/home2/Documents/code/conformers")

    # Load video2 data
    matched_indices = np.load(str(base / "results/afm_pipeline/video2_v2/matched_indices.npy"))
    real_frames = mod.extract_gif_frames(base / "inbox/linz_avb3_video2.gif")
    n_frames = len(matched_indices)
    n_conformers = 264

    # GIF frame indices (same subsampling as main script)
    max_gif = 150
    step = max(1, n_frames // max_gif)
    gif_indices = list(range(0, n_frames, step))[:max_gif]
    n = len(gif_indices)

    # Load renders
    render_dir = base / "results/afm_pipeline/video2_v2/chimerax_renders"
    rendered = {}
    for p in render_dir.glob("render_*.png"):
        idx = int(p.stem.split("_")[1])
        rendered[idx] = p

    # Fine angular resolution: 5° steps
    n_angles = 72
    angles_deg = np.linspace(0, 360, n_angles, endpoint=False)
    render_size = 400

    # --- Step 1: Pre-compute all scores (the expensive part) ---
    print(f"Pre-computing scores: {n} frames × {n_angles} angles...")
    all_scores = np.zeros((n, n_angles), dtype=np.float32)

    for i, frame_i in enumerate(gif_indices):
        conf_idx = int(matched_indices[frame_i]) % n_conformers
        if conf_idx not in rendered:
            continue
        render_rgba = mod.load_render(rendered[conf_idx], target_size=render_size)
        afm_gray, _ = mod._prepare_afm(real_frames[frame_i], render_size)
        threshold = np.median(afm_gray) + 0.15 * (afm_gray.max() - afm_gray.min())
        afm_mask = (afm_gray > threshold).astype(np.float32)
        if afm_mask.sum() < 100:
            continue
        pil_rgba = Image.fromarray((render_rgba * 255).astype(np.uint8), mode="RGBA")
        for j, a in enumerate(angles_deg):
            all_scores[i, j] = mod._score_rotation(afm_mask, pil_rgba, a)
        if (i + 1) % 30 == 0:
            print(f"  {i+1}/{n}")

    # Save scores for re-use
    np.save(str(base / "results/afm_pipeline/video2_v2/rotation_scores.npy"), all_scores)
    print("Scores saved.")

    # Calibration: initial angle
    n_cal = min(10, n)
    cal_pos = np.linspace(0, n - 1, n_cal, dtype=int)
    votes = np.zeros(n_angles, dtype=np.float32)
    for pos in cal_pos:
        s = all_scores[pos]
        if s.max() > 0:
            votes += s / s.max()
    initial_angle = angles_deg[int(np.argmax(votes))]
    print(f"Calibrated initial: {initial_angle:.0f}°")

    # --- Step 2: Parameter sweep ---
    max_deltas = [30, 45, 60, 75, 90, 120, 150]
    smooth_sigmas = [0, 1, 2, 3, 5, 8, 12, 20]

    results = []
    print(f"\nSweeping {len(max_deltas)} × {len(smooth_sigmas)} = "
          f"{len(max_deltas)*len(smooth_sigmas)} combinations...")

    for md in max_deltas:
        for ss in smooth_sigmas:
            angles, scores = score_trajectory(all_scores, angles_deg, md, ss, initial_angle)
            # Only count frames with signal (scores > 0 means the frame had data)
            valid = scores[all_scores.max(axis=1) > 0]
            if len(valid) == 0:
                continue
            worst = float(valid.min())
            p5 = float(np.percentile(valid, 5))
            mean = float(valid.mean())
            # Count "misaligned" frames: score below 50% of frame's own best
            frame_bests = all_scores.max(axis=1)
            ratios = np.where(frame_bests > 0, scores / (frame_bests + 1e-8), 1.0)
            n_misaligned = int((ratios[frame_bests > 0] < 0.7).sum())

            results.append({
                "max_delta": md, "smooth_sigma": ss,
                "worst": worst, "p5": p5, "mean": mean,
                "n_misaligned": n_misaligned,
            })

    # Sort by: fewest misaligned, then best worst-case
    results.sort(key=lambda r: (r["n_misaligned"], -r["worst"]))

    print(f"\n{'max_delta':>10} {'sigma':>6} {'worst':>7} {'p5':>7} "
          f"{'mean':>7} {'misaligned':>10}")
    print("-" * 55)
    for r in results[:20]:
        print(f"{r['max_delta']:>10} {r['smooth_sigma']:>6} {r['worst']:>7.3f} "
              f"{r['p5']:>7.3f} {r['mean']:>7.3f} {r['n_misaligned']:>10}")

    # --- Step 3: Visualize sweep ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Heatmap: worst-case score
    worst_grid = np.full((len(max_deltas), len(smooth_sigmas)), np.nan)
    mis_grid = np.full((len(max_deltas), len(smooth_sigmas)), np.nan)
    p5_grid = np.full((len(max_deltas), len(smooth_sigmas)), np.nan)
    for r in results:
        i = max_deltas.index(r["max_delta"])
        j = smooth_sigmas.index(r["smooth_sigma"])
        worst_grid[i, j] = r["worst"]
        mis_grid[i, j] = r["n_misaligned"]
        p5_grid[i, j] = r["p5"]

    for ax, grid, title, cmap in [
        (axes[0], worst_grid, "Worst-case score\n(higher=better)", "RdYlGn"),
        (axes[1], p5_grid, "5th percentile score\n(higher=better)", "RdYlGn"),
        (axes[2], mis_grid, "Misaligned frames (<70% of best)\n(lower=better)", "RdYlGn_r"),
    ]:
        im = ax.imshow(grid, aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(smooth_sigmas)))
        ax.set_xticklabels(smooth_sigmas)
        ax.set_yticks(range(len(max_deltas)))
        ax.set_yticklabels([f"±{d}°" for d in max_deltas])
        ax.set_xlabel("Smooth sigma")
        ax.set_ylabel("Max delta")
        ax.set_title(title, fontsize=10)
        plt.colorbar(im, ax=ax)
        # Annotate cells
        for ii in range(grid.shape[0]):
            for jj in range(grid.shape[1]):
                if not np.isnan(grid[ii, jj]):
                    val = grid[ii, jj]
                    fmt = f"{val:.0f}" if title.startswith("Mis") else f"{val:.2f}"
                    ax.text(jj, ii, fmt, ha="center", va="center", fontsize=6)

    # Mark best
    best = results[0]
    bi = max_deltas.index(best["max_delta"])
    bj = smooth_sigmas.index(best["smooth_sigma"])
    for ax in axes:
        ax.plot(bj, bi, "s", ms=18, mec="black", mew=2, mfc="none")

    fig.suptitle(f"Rotation Parameter Sweep — Video 2\n"
                 f"Best: max_delta=±{best['max_delta']}°, sigma={best['smooth_sigma']} "
                 f"(worst={best['worst']:.3f}, misaligned={best['n_misaligned']})",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out_path = base / "results/afm_pipeline/video2_v2/rotation_param_sweep.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSweep plot: {out_path}")

    # --- Step 4: Generate trajectory comparison for best vs current ---
    best_md, best_ss = best["max_delta"], best["smooth_sigma"]
    current_md, current_ss = 150, 3  # current settings

    best_angles, best_scores = score_trajectory(
        all_scores, angles_deg, best_md, best_ss, initial_angle)
    curr_angles, curr_scores = score_trajectory(
        all_scores, angles_deg, current_md, current_ss, initial_angle)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    x = np.arange(n)
    ax1.plot(x, curr_angles, "r-", alpha=0.5, lw=1, label=f"Current (±{current_md}°, σ={current_ss})")
    ax1.plot(x, best_angles, "b-", lw=1.5, label=f"Best (±{best_md}°, σ={best_ss})")
    ax1.set_ylabel("Rotation angle (°)")
    ax1.legend(fontsize=8)
    ax1.set_title("Rotation Trajectory Comparison")

    ax2.plot(x, curr_scores, "r-", alpha=0.5, lw=1, label=f"Current per-frame score")
    ax2.plot(x, best_scores, "b-", lw=1.5, label=f"Best per-frame score")
    ax2.axhline(0.7 * all_scores.max(axis=1).mean(), color="gray", ls="--", lw=0.8,
                label="~70% threshold")
    ax2.set_ylabel("Correlation score")
    ax2.set_xlabel("GIF frame index")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    traj_path = base / "results/afm_pipeline/video2_v2/rotation_trajectory_comparison.png"
    fig.savefig(traj_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Trajectory plot: {traj_path}")

    print(f"\n=== RECOMMENDATION ===")
    print(f"  max_delta = {best_md}°")
    print(f"  smooth_sigma = {best_ss}")
    print(f"  worst-case score = {best['worst']:.3f}")
    print(f"  misaligned frames = {best['n_misaligned']}")


if __name__ == "__main__":
    main()
