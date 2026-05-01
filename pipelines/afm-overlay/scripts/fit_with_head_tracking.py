#!/usr/bin/env python3
"""Improved rigid-body fitting with head tracking and tip radius comparison.

Improvements over fit_rigid_body.py:
  1. Skip initial frames (scan window stabilization)
  2. Track head position in AFM images for alignment (instead of COM)
  3. Align PDB head domain to tracked AFM head position
  4. Compare multiple tip radii to find optimal fit

Usage:
  # Compare tip radii on a subset first:
  python fit_with_head_tracking.py --gif video.gif --output-dir out/ \
    --frames-dir frames/ --matched-indices matched.npy \
    --afmfold-root /path/to/afmfold \
    --skip-frames 30 --compare-tips

  # Full fitting with specific tip:
  python fit_with_head_tracking.py --gif video.gif --output-dir out/ \
    --frames-dir frames/ --matched-indices matched.npy \
    --afmfold-root /path/to/afmfold \
    --skip-frames 30 --tip-radius 2.04
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import mdtraj as md
from PIL import Image
from scipy.ndimage import gaussian_filter, center_of_mass, label


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gif", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--frames-dir", type=Path, required=True,
                   help="Directory of protein-only PDB conformer frames")
    p.add_argument("--matched-indices", type=Path, required=True,
                   help="matched_indices.npy from correlation matching")
    p.add_argument("--afmfold-root", type=Path, required=True)
    p.add_argument("--n-conformers", type=int, default=264)
    p.add_argument("--width", type=int, default=35)
    p.add_argument("--height", type=int, default=35)
    p.add_argument("--resolution-nm", type=float, default=0.98)
    p.add_argument("--tip-radius", type=float, default=None,
                   help="Tip radius in pixels. If omitted, uses --compare-tips result")
    p.add_argument("--tip-angle", type=float, default=20.0)
    p.add_argument("--n-rotations", type=int, default=2048)
    p.add_argument("--rot-batch", type=int, default=32)
    p.add_argument("--frame-batch", type=int, default=20)
    p.add_argument("--device", default="cuda")
    # New arguments
    p.add_argument("--skip-frames", type=int, default=30,
                   help="Skip first N frames (scan window adjustment)")
    p.add_argument("--compare-tips", action="store_true",
                   help="Compare tip radii (1,2,3,5 nm) on a subset before full fitting")
    p.add_argument("--tip-compare-subsample", type=int, default=20,
                   help="Use every Nth frame for tip comparison")
    p.add_argument("--head-residues-a", type=str, default="1-440",
                   help="Alpha chain head domain residue range")
    p.add_argument("--head-residues-b", type=str, default="1-350",
                   help="Beta chain head domain residue range")
    p.add_argument("--head-sigma", type=float, default=2.0,
                   help="Gaussian sigma for AFM head detection smoothing")
    p.add_argument("--head-threshold-pct", type=float, default=85,
                   help="Percentile threshold for head region detection")
    return p.parse_args()


# ---------------------------------------------------------------------------
# GIF / image utilities
# ---------------------------------------------------------------------------

def extract_gif_frames(gif_path: Path) -> list[np.ndarray]:
    gif = Image.open(gif_path)
    frames = []
    try:
        while True:
            frames.append(np.array(gif.convert("L"), dtype=np.float32))
            gif.seek(gif.tell() + 1)
    except EOFError:
        pass
    return frames


def resize_to_grid(frames: list[np.ndarray], h: int, w: int) -> np.ndarray:
    from skimage.transform import resize as sk_resize
    out = np.zeros((len(frames), h, w), dtype=np.float32)
    for i, f in enumerate(frames):
        fh, fw = f.shape
        if fh != fw:
            s = min(fh, fw)
            y0, x0 = (fh - s) // 2, (fw - s) // 2
            f = f[y0:y0 + s, x0:x0 + s]
        out[i] = sk_resize(f, (h, w), anti_aliasing=True, preserve_range=True)
    return out


# ---------------------------------------------------------------------------
# Head tracking in AFM images
# ---------------------------------------------------------------------------

def track_head_positions(afm_frames: np.ndarray, sigma: float = 2.0,
                         threshold_pct: float = 85) -> np.ndarray:
    """Track the head (brightest/tallest feature) in each AFM frame.

    Returns (N, 2) array of (row, col) pixel positions of the head center.
    Uses Gaussian smoothing + thresholding + largest-component centroid.
    """
    N, H, W = afm_frames.shape
    positions = np.zeros((N, 2), dtype=np.float32)

    for i in range(N):
        frame = afm_frames[i]
        # Smooth to reduce noise
        smoothed = gaussian_filter(frame, sigma=sigma)
        # Threshold: keep pixels above percentile
        thresh = np.percentile(smoothed, threshold_pct)
        binary = smoothed >= thresh
        # Find largest connected component (the head)
        labeled, n_features = label(binary)
        if n_features == 0:
            # Fallback: use overall centroid
            positions[i] = center_of_mass(smoothed)
            continue
        # Find largest component
        component_sizes = np.bincount(labeled.ravel())[1:]  # skip background
        largest = np.argmax(component_sizes) + 1
        # Centroid of largest component, weighted by intensity
        mask = labeled == largest
        weighted = smoothed * mask
        cy, cx = center_of_mass(weighted)
        positions[i] = [cy, cx]

    return positions


def smooth_head_trajectory(positions: np.ndarray, sigma: float = 5.0) -> np.ndarray:
    """Temporal smoothing of head positions to remove jitter."""
    smoothed = np.zeros_like(positions)
    smoothed[:, 0] = gaussian_filter(positions[:, 0], sigma=sigma)
    smoothed[:, 1] = gaussian_filter(positions[:, 1], sigma=sigma)
    return smoothed


# ---------------------------------------------------------------------------
# PDB head domain selection
# ---------------------------------------------------------------------------

def get_head_atom_indices(topology, res_range_a: str, res_range_b: str) -> np.ndarray:
    """Get atom indices for the integrin head domain.

    Args:
        topology: mdtraj Topology
        res_range_a: "start-end" residue range for chain A head
        res_range_b: "start-end" residue range for chain B head
    """
    a_start, a_end = map(int, res_range_a.split("-"))
    b_start, b_end = map(int, res_range_b.split("-"))
    sel = topology.select(
        f"(chainid 0 and resSeq {a_start} to {a_end}) or "
        f"(chainid 1 and resSeq {b_start} to {b_end})"
    )
    return sel


# ---------------------------------------------------------------------------
# SO(3) fitting with head-based centering
# ---------------------------------------------------------------------------

def fit_orientations_head(
    target_images: np.ndarray,
    xyz: "torch.Tensor",
    head_positions_nm: "torch.Tensor",
    head_atom_indices: np.ndarray,
    xedges: "torch.Tensor",
    yedges: "torch.Tensor",
    tip: "torch.Tensor",
    n_rotations: int,
    rot_batch: int,
    height: int,
    width: int,
    device: "torch.device",
):
    """SO(3) fitting with head-domain-based centering.

    Instead of centering PDB COM at image COM, centers the PDB head domain
    at the tracked head position in the AFM image.
    """
    import torch
    from afmfold.images import sample_uniform_so3, apply_rotations, generate_landscape, idilation
    from afmfold.rigid_body_fitting import compute_correlation_coefficient

    N = xyz.shape[0]
    target = torch.from_numpy(target_images).to(device)
    head_idx = torch.from_numpy(head_atom_indices).long().to(device)

    tmed = torch.median(target)
    tmax = torch.max(target)

    best_cc = torch.full((N,), -1.0, device=device)
    best_rot = torch.zeros((N, 3, 3), device=device)
    best_coords = torch.zeros_like(xyz)
    best_pseudo = torch.zeros((N, height, width), device=device)

    n_steps = (n_rotations + rot_batch - 1) // rot_batch

    for step in range(n_steps):
        rots = sample_uniform_so3(rot_batch, device=device)
        rotated = apply_rotations(xyz, rots)  # (N, rot_batch, natoms, 3)

        # Head-based centering: align PDB head centroid to AFM head position
        head_com = rotated[:, :, head_idx, :].mean(dim=2, keepdim=True)  # (N, rot_batch, 1, 3)
        # Translate so head_com -> head_positions_nm
        target_pos = head_positions_nm.unsqueeze(1).unsqueeze(2)  # (N, 1, 1, 3)
        centered = rotated - head_com + target_pos

        # Z-shift: minimum z to 0
        z_min = centered[..., 2].min(dim=-1, keepdim=True).values
        centered = centered.clone()
        centered[..., 2] = centered[..., 2] - z_min

        flat = centered.reshape(-1, centered.shape[-2], 3)
        landscape, _, _ = generate_landscape(flat, xedges, yedges)
        landscape = landscape.reshape(-1, height, width)
        pseudo = idilation(landscape, tip)

        pmin = pseudo.reshape(len(pseudo), -1).min(dim=1).values.view(-1, 1, 1)
        pmax = pseudo.reshape(len(pseudo), -1).max(dim=1).values.view(-1, 1, 1)
        pseudo_norm = (pseudo - pmin) / (pmax - pmin + 1e-6) * (tmax - tmed) + tmed
        pseudo_norm = pseudo_norm.reshape(N, rot_batch, height, width)

        cc = compute_correlation_coefficient(
            target.unsqueeze(1), pseudo_norm
        ).squeeze(1)

        batch_best_cc, batch_best_idx = cc.max(dim=-1)
        improved = batch_best_cc > best_cc
        if improved.any():
            imp_idx = torch.where(improved)[0]
            best_cc[imp_idx] = batch_best_cc[imp_idx]
            for f in imp_idx:
                bi = batch_best_idx[f]
                best_rot[f] = rots[bi]
                best_coords[f] = centered[f, bi]
                best_pseudo[f] = pseudo_norm[f, bi]

        if (step + 1) % 16 == 0 or step == n_steps - 1:
            print(f"    step {step+1}/{n_steps}: "
                  f"mean_cc={best_cc.mean():.4f} "
                  f"min={best_cc.min():.4f} max={best_cc.max():.4f}")

    return (best_rot.cpu().numpy(), best_cc.cpu().numpy(),
            best_coords.cpu().numpy(), best_pseudo.cpu().numpy())


# ---------------------------------------------------------------------------
# Tip radius comparison
# ---------------------------------------------------------------------------

def compare_tip_radii(
    afm_grid: np.ndarray,
    xyz_tensor,
    head_positions_nm,
    head_atom_indices: np.ndarray,
    xedges, yedges,
    tip_radii_nm: list[float],
    tip_angle: float,
    resolution_nm: float,
    n_rotations: int,
    rot_batch: int,
    height: int,
    width: int,
    device,
    subsample: int = 20,
):
    """Compare fitting quality across multiple tip radii on a frame subset."""
    import torch
    from afmfold.images import generate_tip_shape

    N = afm_grid.shape[0]
    subset_indices = list(range(0, N, subsample))
    n_sub = len(subset_indices)
    print(f"\nTip radius comparison on {n_sub} frames (every {subsample}th)")

    sub_images = afm_grid[subset_indices]
    sub_xyz = xyz_tensor[subset_indices]
    sub_head = head_positions_nm[subset_indices]

    results = {}
    batch_size = 10  # Process comparison frames in small batches to avoid OOM
    for tip_nm in tip_radii_nm:
        tip_px = tip_nm / resolution_nm
        tip = generate_tip_shape(tip_px, tip_angle, device=device)
        print(f"\n  Tip {tip_nm:.1f} nm ({tip_px:.2f} px), kernel size {tip.shape}:")

        all_cc = []
        for b_start in range(0, n_sub, batch_size):
            b_end = min(b_start + batch_size, n_sub)
            _, batch_cc, _, _ = fit_orientations_head(
                sub_images[b_start:b_end],
                sub_xyz[b_start:b_end],
                sub_head[b_start:b_end],
                head_atom_indices,
                xedges, yedges, tip,
                n_rotations, rot_batch, height, width, device,
            )
            all_cc.append(batch_cc)
            torch.cuda.empty_cache()
        cc = np.concatenate(all_cc)
        results[tip_nm] = {
            "mean_cc": float(cc.mean()),
            "median_cc": float(np.median(cc)),
            "min_cc": float(cc.min()),
            "max_cc": float(cc.max()),
            "std_cc": float(cc.std()),
            "tip_px": tip_px,
            "correlations": cc,
        }
        print(f"    → mean={cc.mean():.4f} median={np.median(cc):.4f} "
              f"min={cc.min():.4f} max={cc.max():.4f}")

    # Rank by mean correlation
    ranked = sorted(results.items(), key=lambda x: x[1]["mean_cc"], reverse=True)
    print(f"\n{'='*60}")
    print("TIP RADIUS COMPARISON RESULTS (ranked by mean correlation):")
    for i, (tip_nm, r) in enumerate(ranked):
        marker = " ★ BEST" if i == 0 else ""
        print(f"  {tip_nm:.1f} nm: mean={r['mean_cc']:.4f} "
              f"median={r['median_cc']:.4f} std={r['std_cc']:.4f}{marker}")
    print(f"{'='*60}")

    best_tip_nm = ranked[0][0]
    return results, best_tip_nm


# ---------------------------------------------------------------------------
# Temporal smoothing (flip resolution)
# ---------------------------------------------------------------------------

def resolve_flips(rotations: np.ndarray, coords: np.ndarray,
                  head_indices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Resolve 180-degree orientation flips between consecutive frames.

    Tries 180-degree rotations around x, y, z axes and picks the orientation
    closest to the previous frame.
    """
    from scipy.spatial.transform import Rotation as R

    N = len(rotations)
    smooth_rot = rotations.copy()
    smooth_coords = coords.copy()

    flip_axes = [
        R.from_rotvec([np.pi, 0, 0]).as_matrix(),
        R.from_rotvec([0, np.pi, 0]).as_matrix(),
        R.from_rotvec([0, 0, np.pi]).as_matrix(),
    ]

    n_flips = 0
    for i in range(1, N):
        # Current and previous head centroids
        prev_head = smooth_coords[i-1, head_indices].mean(axis=0)
        curr_head = smooth_coords[i, head_indices].mean(axis=0)
        best_dist = np.linalg.norm(curr_head - prev_head)
        best_rot = smooth_rot[i]
        best_coords = smooth_coords[i]

        for flip in flip_axes:
            # Apply 180-degree flip to current rotation
            flipped_rot = flip @ smooth_rot[i]
            # Recompute coordinates with flipped rotation
            # Use original xyz and re-apply rotation
            flipped_coords = smooth_coords[i].copy()
            com = flipped_coords.mean(axis=0)
            flipped_coords = (flipped_coords - com) @ flip.T + com
            # Z-shift
            flipped_coords[:, 2] -= flipped_coords[:, 2].min()

            flipped_head = flipped_coords[head_indices].mean(axis=0)
            dist = np.linalg.norm(flipped_head - prev_head)
            if dist < best_dist:
                best_dist = dist
                best_rot = flipped_rot
                best_coords = flipped_coords

        if not np.array_equal(best_rot, smooth_rot[i]):
            n_flips += 1
        smooth_rot[i] = best_rot
        smooth_coords[i] = best_coords

    pct = 100 * n_flips / max(1, N - 1)
    print(f"  Flip resolution: {n_flips}/{N-1} frames flipped ({pct:.1f}%)")
    return smooth_rot, smooth_coords


def resolve_flips_head_anchored(
    rotations: np.ndarray,
    coords: np.ndarray,
    head_indices: np.ndarray,
    head_positions_nm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Resolve flips using tracked AFM head positions as anchor.

    Unlike resolve_flips (which compares to the previous frame and can drift),
    this version compares each frame's head centroid to the tracked AFM head
    position, then re-centers to eliminate any residual XY offset.
    """
    from scipy.spatial.transform import Rotation as R

    N = len(rotations)
    smooth_rot = rotations.copy()
    smooth_coords = coords.copy()

    flip_axes = [
        R.from_rotvec([np.pi, 0, 0]).as_matrix(),
        R.from_rotvec([0, np.pi, 0]).as_matrix(),
        R.from_rotvec([0, 0, np.pi]).as_matrix(),
    ]

    n_flips = 0
    for i in range(N):
        target_xy = head_positions_nm[i, :2]
        curr_head = smooth_coords[i, head_indices].mean(axis=0)
        best_dist = np.linalg.norm(curr_head[:2] - target_xy)
        best_rot = smooth_rot[i]
        best_coords = smooth_coords[i]

        for flip in flip_axes:
            flipped_rot = flip @ smooth_rot[i]
            flipped_coords = smooth_coords[i].copy()
            com = flipped_coords.mean(axis=0)
            flipped_coords = (flipped_coords - com) @ flip.T + com
            flipped_coords[:, 2] -= flipped_coords[:, 2].min()

            flipped_head = flipped_coords[head_indices].mean(axis=0)
            dist = np.linalg.norm(flipped_head[:2] - target_xy)
            if dist < best_dist:
                best_dist = dist
                best_rot = flipped_rot
                best_coords = flipped_coords

        if not np.array_equal(best_rot, smooth_rot[i]):
            n_flips += 1
        smooth_rot[i] = best_rot
        smooth_coords[i] = best_coords

        # Re-center head to tracked position
        head_com = smooth_coords[i, head_indices].mean(axis=0)
        smooth_coords[i, :, 0] += target_xy[0] - head_com[0]
        smooth_coords[i, :, 1] += target_xy[1] - head_com[1]
        smooth_coords[i, :, 2] -= smooth_coords[i, :, 2].min()

    pct = 100 * n_flips / max(1, N)
    print(f"  Flip resolution (head-anchored): {n_flips}/{N} frames flipped ({pct:.1f}%)")
    return smooth_rot, smooth_coords


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_tip_comparison(results: dict, output_path: Path):
    """Plot tip radius comparison results."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tip_nms = sorted(results.keys())
    means = [results[t]["mean_cc"] for t in tip_nms]
    medians = [results[t]["median_cc"] for t in tip_nms]
    stds = [results[t]["std_cc"] for t in tip_nms]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Bar chart of mean/median correlation
    ax = axes[0]
    x = np.arange(len(tip_nms))
    w = 0.35
    ax.bar(x - w/2, means, w, label="Mean", color="steelblue", yerr=stds, capsize=4)
    ax.bar(x + w/2, medians, w, label="Median", color="salmon")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{t:.0f} nm" for t in tip_nms])
    ax.set_ylabel("Correlation Coefficient")
    ax.set_title("Fitting Quality vs Tip Radius")
    ax.legend()
    ax.set_ylim(0, 1)
    best_idx = np.argmax(means)
    ax.annotate("BEST", xy=(best_idx - w/2, means[best_idx] + stds[best_idx] + 0.02),
                ha="center", fontsize=10, fontweight="bold", color="green")

    # Box plot of per-frame correlations
    ax = axes[1]
    data = [results[t]["correlations"] for t in tip_nms]
    bp = ax.boxplot(data, labels=[f"{t:.0f} nm" for t in tip_nms], patch_artist=True)
    colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_ylabel("Correlation Coefficient")
    ax.set_title("Per-Frame Correlation Distribution")
    ax.set_xlabel("Tip Radius")

    fig.suptitle("Tip Radius Comparison for SO(3) Rigid-Body Fitting",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved tip comparison plot: {output_path}")


def plot_head_tracking(positions: np.ndarray, smoothed: np.ndarray,
                       output_path: Path, skip_frames: int = 0):
    """Plot head position trajectory."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    frames = np.arange(len(positions)) + skip_frames

    ax = axes[0]
    ax.plot(frames, positions[:, 1], ".", alpha=0.3, markersize=2, label="Raw X")
    ax.plot(frames, smoothed[:, 1], "-", linewidth=1.5, label="Smooth X")
    ax.plot(frames, positions[:, 0], ".", alpha=0.3, markersize=2, label="Raw Y")
    ax.plot(frames, smoothed[:, 0], "-", linewidth=1.5, label="Smooth Y")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Position (pixels)")
    ax.set_title("Head Position Trajectory")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.scatter(positions[:, 1], positions[:, 0], c=frames, cmap="viridis",
               s=3, alpha=0.5)
    ax.plot(smoothed[:, 1], smoothed[:, 0], "r-", linewidth=1.5, alpha=0.8)
    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")
    ax.set_title("Head Position (color = time)")
    ax.set_aspect("equal")
    cb = plt.colorbar(ax.collections[0], ax=ax, label="Frame")

    fig.suptitle("AFM Head Position Tracking", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved head tracking plot: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    afmfold_src = str(args.afmfold_root / "src")
    if afmfold_src not in sys.path:
        sys.path.insert(0, afmfold_src)

    import torch
    from afmfold.images import generate_tip_shape
    from afmfold.rigid_body_fitting import threshold_and_mask

    device = torch.device(args.device)

    # --- Load data ---
    matched_indices = np.load(str(args.matched_indices))
    n_total = len(matched_indices)
    print(f"Total frames from correlation matching: {n_total}")

    # Skip initial frames
    if args.skip_frames > 0:
        print(f"Skipping first {args.skip_frames} frames (scan window adjustment)")
        matched_indices = matched_indices[args.skip_frames:]
        n_frames = len(matched_indices)
        print(f"Fitting frames {args.skip_frames} to {n_total-1} ({n_frames} frames)")
    else:
        n_frames = n_total

    conformer_indices = (matched_indices % args.n_conformers).astype(int)

    # Extract and resize GIF frames
    print("Extracting GIF frames...")
    all_real_frames = extract_gif_frames(args.gif)
    # Only use frames after skip
    real_frames_subset = all_real_frames[args.skip_frames:]
    afm_grid = resize_to_grid(real_frames_subset, args.height, args.width)
    print(f"  Resized {len(afm_grid)} frames to {args.height}x{args.width}")

    # --- Head tracking in AFM images ---
    print("\nTracking head position in AFM images...")
    head_positions_px = track_head_positions(
        afm_grid, sigma=args.head_sigma, threshold_pct=args.head_threshold_pct
    )
    head_positions_smooth_px = smooth_head_trajectory(head_positions_px, sigma=5.0)
    print(f"  Head position range: "
          f"x=[{head_positions_smooth_px[:,1].min():.1f}, {head_positions_smooth_px[:,1].max():.1f}] "
          f"y=[{head_positions_smooth_px[:,0].min():.1f}, {head_positions_smooth_px[:,0].max():.1f}] px")
    drift = np.sqrt(np.sum(np.diff(head_positions_smooth_px, axis=0)**2, axis=1))
    print(f"  Mean frame-to-frame drift: {drift.mean():.2f} px ({drift.mean()*args.resolution_nm:.2f} nm)")

    # Convert head positions to nm coordinates (matching the pseudo-AFM grid)
    head_positions_nm = np.zeros((n_frames, 3), dtype=np.float32)
    half_w = args.width / 2
    half_h = args.height / 2
    head_positions_nm[:, 0] = (head_positions_smooth_px[:, 1] - half_w) * args.resolution_nm
    head_positions_nm[:, 1] = (head_positions_smooth_px[:, 0] - half_h) * args.resolution_nm
    head_positions_nm[:, 2] = 0.0

    # Plot head tracking
    plot_head_tracking(head_positions_px, head_positions_smooth_px,
                       args.output_dir / "head_tracking.png",
                       skip_frames=args.skip_frames)

    # --- Load conformer PDBs ---
    print("\nLoading conformer PDBs...")
    pdb_files = sorted(args.frames_dir.glob("*.pdb"))
    ref = md.load(str(pdb_files[0]))
    keep = ref.topology.select("protein")
    ref = ref.atom_slice(keep)
    n_atoms = ref.n_atoms
    topology = ref.topology

    conformer_xyz = {}
    for pf in pdb_files:
        idx = int(pf.stem.split("_")[-1])
        t = md.load(str(pf))
        t = t.atom_slice(t.topology.select("protein"))
        if t.n_atoms == n_atoms:
            conformer_xyz[idx] = t.xyz[0]
    print(f"  Loaded {len(conformer_xyz)} conformers ({n_atoms} atoms)")

    # Head domain atom indices
    head_atom_indices = get_head_atom_indices(
        topology, args.head_residues_a, args.head_residues_b
    )
    print(f"  Head domain: {len(head_atom_indices)} atoms "
          f"({100*len(head_atom_indices)/n_atoms:.1f}% of structure)")

    # Build per-frame xyz
    all_xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    for i in range(n_frames):
        cidx = conformer_indices[i]
        if cidx in conformer_xyz:
            all_xyz[i] = conformer_xyz[cidx]
        else:
            nearest = min(conformer_xyz.keys(), key=lambda k: abs(k - cidx))
            all_xyz[i] = conformer_xyz[nearest]

    # Setup grid
    xedges = args.resolution_nm * torch.arange(
        -int(args.width / 2) - 1, args.width - int(args.width / 2),
        device=device)
    yedges = args.resolution_nm * torch.arange(
        -int(args.height / 2) - 1, args.height - int(args.height / 2),
        device=device)

    xyz_tensor = torch.from_numpy(all_xyz).to(device)
    head_nm_tensor = torch.from_numpy(head_positions_nm).to(device)

    # --- Tip radius comparison ---
    if args.compare_tips:
        tip_radii_nm = [1.0, 2.0, 3.0, 5.0]
        tip_results, best_tip_nm = compare_tip_radii(
            afm_grid, xyz_tensor, head_nm_tensor, head_atom_indices,
            xedges, yedges,
            tip_radii_nm, args.tip_angle, args.resolution_nm,
            args.n_rotations, args.rot_batch,
            args.height, args.width, device,
            subsample=args.tip_compare_subsample,
        )
        plot_tip_comparison(tip_results, args.output_dir / "tip_comparison.png")

        # Save comparison results
        np.savez(str(args.output_dir / "tip_comparison.npz"),
                 tip_radii_nm=tip_radii_nm,
                 **{f"cc_{t:.0f}nm": tip_results[t]["correlations"] for t in tip_radii_nm})

        if args.tip_radius is None:
            args.tip_radius = best_tip_nm / args.resolution_nm
            print(f"\nUsing best tip radius: {best_tip_nm:.1f} nm ({args.tip_radius:.2f} px)")

    if args.tip_radius is None:
        args.tip_radius = 2.04  # Default 2nm
        print(f"No tip radius specified, using default: 2.0 nm ({args.tip_radius:.2f} px)")

    tip = generate_tip_shape(args.tip_radius, args.tip_angle, device=device)
    tip_nm = args.tip_radius * args.resolution_nm
    print(f"\nFitting with tip radius {tip_nm:.1f} nm ({args.tip_radius:.2f} px), "
          f"kernel size {tip.shape}")

    # --- Full fitting ---
    fitted_rotations = np.zeros((n_frames, 3, 3), dtype=np.float32)
    fitted_coords = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    fitted_correlations = np.zeros(n_frames, dtype=np.float32)
    fitted_pseudo = np.zeros((n_frames, args.height, args.width), dtype=np.float32)

    t0 = time.time()
    print(f"\nFull SO(3) fitting: {args.n_rotations} rotations/frame, "
          f"head-based centering, {n_frames} frames")

    for batch_start in range(0, n_frames, args.frame_batch):
        batch_end = min(batch_start + args.frame_batch, n_frames)
        n_batch = batch_end - batch_start
        batch_num = batch_start // args.frame_batch + 1
        total_batches = (n_frames + args.frame_batch - 1) // args.frame_batch

        print(f"\nBatch [{batch_start}:{batch_end}] ({batch_num}/{total_batches})")

        rot, cc, coords, pseudo = fit_orientations_head(
            afm_grid[batch_start:batch_end],
            xyz_tensor[batch_start:batch_end],
            head_nm_tensor[batch_start:batch_end],
            head_atom_indices,
            xedges, yedges, tip,
            args.n_rotations, args.rot_batch,
            args.height, args.width, device,
        )

        fitted_rotations[batch_start:batch_end] = rot
        fitted_coords[batch_start:batch_end] = coords
        fitted_correlations[batch_start:batch_end] = cc
        fitted_pseudo[batch_start:batch_end] = pseudo

        elapsed = time.time() - t0
        rate = batch_end / elapsed
        eta = (n_frames - batch_end) / rate if rate > 0 else 0
        print(f"  → corr: mean={cc.mean():.4f} | "
              f"elapsed={elapsed:.0f}s, ETA={eta:.0f}s")

    # --- Temporal smoothing ---
    print("\nApplying head-anchored flip resolution...")
    smooth_rot, smooth_coords = resolve_flips_head_anchored(
        fitted_rotations, fitted_coords, head_atom_indices, head_positions_nm
    )

    # --- Save results ---
    # Save with metadata about what was skipped
    np.save(str(args.output_dir / "fitted_rotations.npy"), fitted_rotations)
    np.save(str(args.output_dir / "fitted_coords.npy"), fitted_coords)
    np.save(str(args.output_dir / "fitted_correlations.npy"), fitted_correlations)
    np.save(str(args.output_dir / "fitted_pseudo_afm.npy"), fitted_pseudo)
    np.save(str(args.output_dir / "fitted_rotations_smooth.npy"), smooth_rot)
    np.save(str(args.output_dir / "fitted_coords_smooth.npy"), smooth_coords)
    np.save(str(args.output_dir / "head_positions_px.npy"), head_positions_smooth_px)

    # Save topology
    ref.save_pdb(str(args.output_dir / "topology.pdb"))

    # Save metadata
    import json
    meta = {
        "skip_frames": args.skip_frames,
        "n_frames_fitted": n_frames,
        "n_frames_total": n_total,
        "frame_range": [args.skip_frames, n_total],
        "tip_radius_nm": float(tip_nm),
        "tip_radius_px": float(args.tip_radius),
        "tip_angle": float(args.tip_angle),
        "n_rotations": args.n_rotations,
        "head_residues_a": args.head_residues_a,
        "head_residues_b": args.head_residues_b,
        "n_head_atoms": int(len(head_atom_indices)),
        "correlation_mean": float(fitted_correlations.mean()),
        "correlation_std": float(fitted_correlations.std()),
        "correlation_min": float(fitted_correlations.min()),
        "correlation_max": float(fitted_correlations.max()),
        "elapsed_seconds": float(time.time() - t0),
    }
    with open(args.output_dir / "fitting_metadata.json", "w") as f:
        json.dump(meta, f, indent=2)

    total_time = time.time() - t0
    print(f"\n{'='*60}")
    print(f"HEAD-TRACKING RIGID-BODY FITTING COMPLETE")
    print(f"  Frames: {n_frames} (skipped first {args.skip_frames})")
    print(f"  Tip radius: {tip_nm:.1f} nm ({args.tip_radius:.2f} px)")
    print(f"  Head atoms: {len(head_atom_indices)}")
    print(f"  Correlation: mean={fitted_correlations.mean():.4f} "
          f"std={fitted_correlations.std():.4f}")
    print(f"  Time: {total_time:.0f}s ({total_time/60:.1f} min)")
    print(f"  Output: {args.output_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
