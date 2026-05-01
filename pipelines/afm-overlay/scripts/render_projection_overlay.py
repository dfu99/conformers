#!/usr/bin/env python3
"""Render PDB overlays via direct 2D projection (no ChimeraX).

Projects fitted 3D coordinates to 2D pixel positions on the AFM image,
preserving correct spatial alignment. This replaces the ChimeraX-based
rendering which auto-centers every frame, breaking the positioning.

The projection is top-down (x,y → col,row), with z used for coloring
(height-based brightness). Atoms are drawn as circles colored by chain.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from PIL import Image
from scipy.ndimage import gaussian_filter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gif", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--fitted-dir", type=Path, required=True,
                   help="Directory from fit_with_head_tracking.py")
    p.add_argument("--matched-cvs", type=Path, required=True)
    p.add_argument("--resolution-nm", type=float, default=0.98)
    p.add_argument("--grid-size", type=int, default=35)
    p.add_argument("--render-size", type=int, default=400,
                   help="Output image size (square)")
    p.add_argument("--max-gif-frames", type=int, default=150)
    p.add_argument("--n-grid-frames", type=int, default=16)
    p.add_argument("--fps", type=int, default=8)
    p.add_argument("--atom-selection", default="name CA",
                   help="MDTraj atom selection for rendering")
    p.add_argument("--atom-radius", type=float, default=3.0,
                   help="Atom circle radius in pixels")
    p.add_argument("--alpha", type=float, default=0.7,
                   help="Overlay opacity")
    return p.parse_args()


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


def prepare_afm(frame: np.ndarray, size: int) -> np.ndarray:
    """Crop to square, resize, apply colormap. Returns RGB float [0,1]."""
    pil = Image.fromarray(frame.astype(np.uint8) if frame.max() > 1
                          else (frame * 255).astype(np.uint8))
    w, h = pil.size
    if w != h:
        s = min(w, h)
        left, top = (w - s) // 2, (h - s) // 2
        pil = pil.crop((left, top, left + s, top + s))
    pil = pil.resize((size, size), Image.LANCZOS)
    gray = np.array(pil, dtype=np.float32) / 255.0
    colored = plt.cm.afmhot(gray)[:, :, :3]
    return colored


def project_atoms_to_image(
    coords_nm: np.ndarray,
    chain_ids: np.ndarray,
    grid_size: int,
    resolution_nm: float,
    render_size: int,
    atom_radius: float = 3.0,
    alpha: float = 0.7,
) -> np.ndarray:
    """Project 3D coordinates to 2D overlay image (RGBA).

    Coordinates are in nm on the pseudo-AFM grid (centered at 0,0).
    Returns (render_size, render_size, 4) RGBA float array.
    """
    half = grid_size / 2.0

    # Convert nm to pixel positions on the grid
    cols = coords_nm[:, 0] / resolution_nm + half  # x → col
    rows = coords_nm[:, 1] / resolution_nm + half  # y → row
    z = coords_nm[:, 2]

    # Scale to render size
    scale = render_size / grid_size
    px_cols = cols * scale
    px_rows = rows * scale

    # Normalize z for brightness (higher = brighter)
    z_min, z_max = z.min(), z.max()
    if z_max > z_min:
        z_norm = (z - z_min) / (z_max - z_min)
    else:
        z_norm = np.ones_like(z)

    # Create RGBA overlay
    overlay = np.zeros((render_size, render_size, 4), dtype=np.float32)

    # Draw atoms as filled circles on the overlay
    # Sort by z so higher atoms are drawn on top
    order = np.argsort(z)

    # Chain colors
    chain_colors = {
        0: np.array([0.39, 0.58, 0.93]),   # cornflowerblue for chain A
        1: np.array([0.98, 0.50, 0.45]),    # salmon for chain B
    }

    fig, ax = plt.subplots(figsize=(render_size / 100, render_size / 100), dpi=100)
    ax.set_xlim(0, render_size)
    ax.set_ylim(render_size, 0)  # flip y for image coords
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)

    for i in order:
        color = chain_colors.get(chain_ids[i], np.array([0.5, 0.5, 0.5]))
        brightness = 0.4 + 0.6 * z_norm[i]
        c = color * brightness
        circle = Circle((px_cols[i], px_rows[i]), atom_radius,
                         facecolor=(*c, alpha), edgecolor="none")
        ax.add_patch(circle)

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, transparent=True,
                bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    buf.seek(0)
    img = Image.open(buf).convert("RGBA").resize((render_size, render_size), Image.LANCZOS)
    return np.array(img, dtype=np.float32) / 255.0


def composite_overlay(afm_colored: np.ndarray, overlay_rgba: np.ndarray) -> np.ndarray:
    """Alpha-blend overlay on AFM image."""
    rgb = overlay_rgba[:, :, :3]
    alpha = overlay_rgba[:, :, 3:4]
    return np.clip(afm_colored * (1 - alpha) + rgb * alpha, 0, 1)


def classify_state(cvs: np.ndarray) -> str:
    mean_cv = cvs.mean()
    if mean_cv < 60:
        return "Compact"
    elif mean_cv < 70:
        return "Intermediate"
    else:
        return "Extended"


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load fitting metadata
    meta_path = args.fitted_dir / "fitting_metadata.json"
    skip_frames = 0
    if meta_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)
        skip_frames = meta.get("skip_frames", 0)
        print(f"Skip frames: {skip_frames}")

    # Load fitted coordinates
    smooth_path = args.fitted_dir / "fitted_coords_smooth.npy"
    raw_path = args.fitted_dir / "fitted_coords.npy"
    if smooth_path.exists():
        fitted_coords = np.load(str(smooth_path))
        print("Using smoothed coordinates")
    else:
        fitted_coords = np.load(str(raw_path))
    correlations = np.load(str(args.fitted_dir / "fitted_correlations.npy"))
    topology_path = args.fitted_dir / "topology.pdb"
    n_frames = len(fitted_coords)
    print(f"Loaded {n_frames} fitted frames")

    # Load CVs (slice if needed)
    all_cvs = np.load(str(args.matched_cvs))
    if skip_frames > 0:
        matched_cvs = all_cvs[skip_frames:]
    else:
        matched_cvs = all_cvs

    # Load topology and get atom selection + chain IDs
    ref = md.load(str(topology_path))
    atom_sel = ref.topology.select(args.atom_selection)
    print(f"Rendering {len(atom_sel)} atoms ({args.atom_selection})")

    # Get chain ID for each selected atom
    chain_ids = np.array([ref.topology.atom(i).residue.chain.index for i in atom_sel])
    unique_chains = np.unique(chain_ids)
    print(f"Chains: {unique_chains}")

    # Extract GIF frames
    all_real_frames = extract_gif_frames(args.gif)
    if skip_frames > 0:
        real_frames = all_real_frames[skip_frames:]
    else:
        real_frames = all_real_frames
    print(f"Using {len(real_frames)} GIF frames (skipped {skip_frames})")

    # Select frames for output
    step = max(1, n_frames // args.max_gif_frames)
    gif_indices = list(range(0, n_frames, step))[:args.max_gif_frames]
    print(f"Selected {len(gif_indices)} frames for GIF")

    # --- Create grid ---
    n_show = min(args.n_grid_frames, len(gif_indices))
    grid_step = max(1, len(gif_indices) // n_show)
    grid_selected = gif_indices[::grid_step][:n_show]

    n_cols = min(4, len(grid_selected))
    n_rows = (len(grid_selected) + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols * 2, figsize=(n_cols * 4, n_rows * 2.2))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    for idx, fi in enumerate(grid_selected):
        row = idx // n_cols
        col_base = (idx % n_cols) * 2

        afm = prepare_afm(real_frames[fi], args.render_size)
        coords_sel = fitted_coords[fi, atom_sel]
        overlay = project_atoms_to_image(
            coords_sel, chain_ids,
            args.grid_size, args.resolution_nm, args.render_size,
            args.atom_radius, args.alpha
        )
        comp = composite_overlay(afm, overlay)

        # Real AFM
        ax_real = axes[row, col_base]
        ax_real.imshow(afm, interpolation="bilinear")
        ax_real.set_title(f"Frame {fi + skip_frames}", fontsize=8, fontweight="bold")
        ax_real.axis("off")

        # Overlay
        ax_ov = axes[row, col_base + 1]
        ax_ov.imshow(comp, interpolation="bilinear")
        cvs = matched_cvs[fi]
        state = classify_state(cvs)
        cc = correlations[fi]
        ax_ov.set_title(f"{state} (r={cc:.2f})\n"
                        f"CVs: {cvs[0]:.0f}/{cvs[1]:.0f}/{cvs[2]:.0f}\u00c5",
                        fontsize=6)
        ax_ov.axis("off")

    for idx in range(len(grid_selected), n_rows * n_cols):
        row = idx // n_cols
        for c in range(2):
            axes[row, (idx % n_cols) * 2 + c].axis("off")

    fig.suptitle("HS-AFM with 2D Projection Overlay (head-tracked)\n"
                 "(Real AFM | Fitted Conformer \u2014 CA atoms, colored by chain)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    grid_path = args.output_dir / "pdb_projection_grid.png"
    fig.savefig(grid_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved grid: {grid_path}")

    # --- Create animated GIF ---
    cv_names = ["\u03b1H\u2194\u03b1T", "\u03b2H\u2194\u03b1T", "\u03b1H\u2194\u03b2T"]
    state_colors = {"Compact": "#3366CC", "Intermediate": "#FF9900",
                    "Extended": "#CC3333"}

    pil_frames = []
    for i, fi in enumerate(gif_indices):
        afm = prepare_afm(real_frames[fi], args.render_size)
        coords_sel = fitted_coords[fi, atom_sel]
        overlay = project_atoms_to_image(
            coords_sel, chain_ids,
            args.grid_size, args.resolution_nm, args.render_size,
            args.atom_radius, args.alpha
        )
        comp = composite_overlay(afm, overlay)

        cvs = matched_cvs[fi]
        cc = correlations[fi]
        state = classify_state(cvs)

        fig = plt.figure(figsize=(10, 3.5))
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.6], wspace=0.15)

        # Panel 1: Real AFM
        ax1 = fig.add_subplot(gs[0])
        ax1.imshow(afm, interpolation="bilinear")
        ax1.set_title("HS-AFM", fontsize=9, fontweight="bold")
        ax1.axis("off")

        # Panel 2: Projection overlay
        ax2 = fig.add_subplot(gs[1])
        ax2.imshow(comp, interpolation="bilinear")
        ax2.set_title(f"PDB Overlay (r={cc:.2f})", fontsize=9, fontweight="bold")
        ax2.axis("off")

        # Panel 3: CV bar chart
        ax3 = fig.add_subplot(gs[2])
        bars = ax3.barh(cv_names, cvs, color=state_colors[state], height=0.6)
        ax3.set_xlim(45, 85)
        ax3.set_xlabel("Distance (\u00c5)", fontsize=8)
        ax3.set_title(state, fontsize=10, fontweight="bold",
                      color=state_colors[state])
        for bar, v in zip(bars, cvs):
            ax3.text(v + 0.3, bar.get_y() + bar.get_height() / 2,
                     f"{v:.1f}", va="center", fontsize=7)
        ax3.tick_params(labelsize=7)

        display_frame = fi + skip_frames
        n_total = len(real_frames) + skip_frames
        fig.suptitle(f"Frame {display_frame}/{n_total}  |  "
                     f"\u03b1v\u03b23 Integrin Conformer",
                     fontsize=10, fontweight="bold", y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120, bbox_inches="tight",
                    facecolor="white")
        plt.close(fig)
        buf.seek(0)
        pil_frames.append(Image.open(buf).convert("RGB"))

        if (i + 1) % 30 == 0:
            print(f"  GIF frame {i + 1}/{len(gif_indices)}")

    gif_path = args.output_dir / "pdb_projection_video.gif"
    duration_ms = int(1000 / args.fps)
    pil_frames[0].save(
        gif_path, save_all=True, append_images=pil_frames[1:],
        duration=duration_ms, loop=0,
    )
    print(f"Saved GIF: {gif_path} ({len(pil_frames)} frames)")

    # Optimize GIF
    opt_path = args.output_dir / "pdb_projection_video_opt.gif"
    gif_img = Image.open(gif_path)
    opt_frames = []
    try:
        while True:
            frame = gif_img.copy().convert("RGB")
            w, h = frame.size
            frame = frame.resize((int(w * 0.6), int(h * 0.6)), Image.LANCZOS)
            frame = frame.quantize(colors=64, method=Image.Quantize.MEDIANCUT)
            opt_frames.append(frame)
            gif_img.seek(gif_img.tell() + 1)
    except EOFError:
        pass
    opt_frames[0].save(opt_path, save_all=True, append_images=opt_frames[1:],
                       duration=duration_ms, loop=0, optimize=True)
    import os
    orig_mb = os.path.getsize(gif_path) / 1e6
    opt_mb = os.path.getsize(opt_path) / 1e6
    print(f"Optimized: {orig_mb:.1f}MB -> {opt_mb:.1f}MB")

    print(f"\n{'='*60}")
    print("PROJECTION OVERLAY COMPLETE")
    print(f"  Grid: {grid_path}")
    print(f"  Video: {gif_path}")
    print(f"  Optimized: {opt_path}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
