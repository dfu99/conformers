#!/usr/bin/env python3
"""Visualize AFMFold pipeline results: overlay fitted PDB on HS-AFM frames.

Two matching strategies:
  1. CNN-predicted CVs → nearest training conformer
  2. Direct image correlation between real AFM and pseudo-AFM library

Outputs:
  - Grid of AFM frames with PDB silhouette overlay
  - Animated GIF video showing frame-by-frame conformer evolution
  - CV prediction summary plots
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gif", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--frames-dir", type=Path, required=True,
                   help="Directory of protein-only PDB conformer frames.")
    p.add_argument("--training-data-dir", type=Path, required=True,
                   help="Directory with pseudo-AFM image_*.npy and label_*.npy.")
    p.add_argument("--cv-distances", type=Path, default=None,
                   help="Pre-computed CV distances .npy for the conformer frames.")
    p.add_argument("--predicted-cvs", type=Path, default=None,
                   help="CNN-predicted CVs .npy from predict_from_afm_gif.py.")
    p.add_argument("--image-size", type=int, default=35)
    p.add_argument("--n-grid-frames", type=int, default=20,
                   help="Number of frames to show in grid.")
    p.add_argument("--device", default="cpu")
    return p.parse_args()


# ── GIF frame extraction ──────────────────────────────────────────────────

def extract_gif_frames(gif_path: Path) -> list[np.ndarray]:
    from PIL import Image
    gif = Image.open(gif_path)
    frames = []
    try:
        while True:
            frames.append(np.array(gif.convert("L"), dtype=np.float32))
            gif.seek(gif.tell() + 1)
    except EOFError:
        pass
    print(f"Extracted {len(frames)} frames ({frames[0].shape}) from {gif_path.name}")
    return frames


def preprocess_frame(frame: np.ndarray, target_size: int) -> np.ndarray:
    from PIL import Image
    h, w = frame.shape
    if h != w:
        s = min(h, w)
        y0, x0 = (h - s) // 2, (w - s) // 2
        frame = frame[y0:y0+s, x0:x0+s]
    img = Image.fromarray(frame)
    img = img.resize((target_size, target_size), Image.BILINEAR)
    arr = np.array(img, dtype=np.float32)
    arr = arr - arr.min()
    if arr.max() > 0:
        arr = arr / arr.max()
    return arr


# ── Load pseudo-AFM library ───────────────────────────────────────────────

def load_pseudo_afm_library(data_dir: Path):
    """Load all pseudo-AFM images and their CV labels."""
    img_files = sorted(data_dir.glob("image_*.npy"))
    lbl_files = sorted(data_dir.glob("label_*.npy"))
    images = np.concatenate([np.load(str(f)) for f in img_files], axis=0)
    labels = np.concatenate([np.load(str(f)) for f in lbl_files], axis=0)
    print(f"Pseudo-AFM library: {images.shape[0]} images, {labels.shape[1]} CVs")
    return images, labels


# ── Image correlation matching ─────────────────────────────────────────────

def correlation_match(real_frames: list[np.ndarray], pseudo_images: np.ndarray,
                      pseudo_labels: np.ndarray, image_size: int,
                      device: str = "cpu") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Match each real AFM frame to best pseudo-AFM by correlation coefficient.

    Returns:
        matched_cvs: (N, D) CVs of best-matching pseudo-AFM
        matched_indices: (N,) index into pseudo library
        correlations: (N,) best correlation per frame
    """
    import torch

    # Preprocess real frames to match pseudo-AFM dimensions
    real = np.stack([preprocess_frame(f, image_size) for f in real_frames])
    n_real = len(real)
    n_pseudo = len(pseudo_images)

    # Normalize pseudo images to [0,1]
    pseudo_norm = pseudo_images.copy()
    for i in range(n_pseudo):
        mn, mx = pseudo_norm[i].min(), pseudo_norm[i].max()
        if mx > mn:
            pseudo_norm[i] = (pseudo_norm[i] - mn) / (mx - mn)

    # Compute correlation in batches on GPU/CPU
    dev = torch.device(device)
    real_t = torch.from_numpy(real).float().to(dev)          # (N_real, H, W)
    pseudo_t = torch.from_numpy(pseudo_norm).float().to(dev)  # (N_pseudo, H, W)

    # Flatten for batch correlation
    real_flat = real_t.reshape(n_real, -1)       # (N_real, H*W)
    pseudo_flat = pseudo_t.reshape(n_pseudo, -1)  # (N_pseudo, H*W)

    # Normalize for cosine similarity = correlation coefficient
    real_flat = real_flat - real_flat.mean(dim=1, keepdim=True)
    pseudo_flat = pseudo_flat - pseudo_flat.mean(dim=1, keepdim=True)
    real_norm = real_flat / (real_flat.norm(dim=1, keepdim=True) + 1e-8)
    pseudo_norm_t = pseudo_flat / (pseudo_flat.norm(dim=1, keepdim=True) + 1e-8)

    # Batch correlation: (N_real, N_pseudo)
    batch_size = 50
    matched_indices = np.zeros(n_real, dtype=np.int64)
    correlations = np.zeros(n_real, dtype=np.float32)

    for start in range(0, n_real, batch_size):
        end = min(start + batch_size, n_real)
        cc = torch.mm(real_norm[start:end], pseudo_norm_t.T)  # (batch, N_pseudo)
        best_idx = cc.argmax(dim=1).cpu().numpy()
        best_val = cc.max(dim=1).values.cpu().numpy()
        matched_indices[start:end] = best_idx
        correlations[start:end] = best_val

    matched_cvs = pseudo_labels[matched_indices]
    print(f"Correlation matching: {n_real} frames → {n_pseudo} pseudo-AFM images")
    print(f"  Correlation range: [{correlations.min():.3f}, {correlations.max():.3f}]")
    print(f"  Matched CV ranges:")
    for i in range(matched_cvs.shape[1]):
        print(f"    CV{i}: [{matched_cvs[:,i].min():.1f}, {matched_cvs[:,i].max():.1f}] Å")

    return matched_cvs, matched_indices, correlations


# ── PDB to 2D silhouette projection ───────────────────────────────────────

def pdb_to_silhouette(pdb_path: Path, image_size: int = 200,
                      sigma: float = 3.0) -> np.ndarray:
    """Project PDB structure to 2D density map (top-down z-projection)."""
    import mdtraj as md
    from scipy.ndimage import gaussian_filter

    t = md.load(str(pdb_path))
    # Use CA atoms for cleaner silhouette
    ca_idx = t.topology.select("name CA")
    if len(ca_idx) == 0:
        ca_idx = t.topology.select("protein")
    xyz = t.xyz[0, ca_idx]  # (N, 3) in nm

    # Project to x-y plane
    x, y = xyz[:, 0], xyz[:, 1]
    # Center and scale
    x = x - x.mean()
    y = y - y.mean()
    scale = max(x.max() - x.min(), y.max() - y.min())
    if scale > 0:
        x = (x / scale) * 0.7 * image_size + image_size / 2
        y = (y / scale) * 0.7 * image_size + image_size / 2

    # Rasterize as density
    density = np.zeros((image_size, image_size), dtype=np.float32)
    for xi, yi in zip(x, y):
        ix, iy = int(round(yi)), int(round(xi))
        if 0 <= ix < image_size and 0 <= iy < image_size:
            density[ix, iy] += 1.0

    density = gaussian_filter(density, sigma=sigma)
    if density.max() > 0:
        density /= density.max()
    return density


# ── Grid visualization ─────────────────────────────────────────────────────

def create_overlay_grid(real_frames, matched_indices, pseudo_images,
                        matched_cvs, correlations, pdb_files,
                        output_path: Path, n_show: int = 20,
                        image_size: int = 35):
    """Create a grid showing real AFM with pseudo-AFM overlay and PDB silhouette."""
    n_frames = len(real_frames)
    step = max(1, n_frames // n_show)
    selected = list(range(0, n_frames, step))[:n_show]

    n_cols = min(5, len(selected))
    n_rows = (len(selected) + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols * 3, figsize=(n_cols * 6, n_rows * 2.2))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    cv_names = ["αH↔αT", "βH↔αT", "αH↔βT"]

    for idx, frame_i in enumerate(selected):
        row = idx // n_cols
        col_base = (idx % n_cols) * 3

        real = preprocess_frame(real_frames[frame_i], image_size)
        pseudo_idx = matched_indices[frame_i]
        pseudo = pseudo_images[pseudo_idx]
        pseudo_n = (pseudo - pseudo.min()) / (pseudo.max() - pseudo.min() + 1e-8)

        # Real AFM frame
        ax_real = axes[row, col_base]
        ax_real.imshow(real, cmap="afmhot", interpolation="nearest")
        ax_real.set_title(f"Frame {frame_i}", fontsize=7)
        ax_real.axis("off")

        # Best-matching pseudo-AFM
        ax_pseudo = axes[row, col_base + 1]
        ax_pseudo.imshow(pseudo_n, cmap="afmhot", interpolation="nearest")
        cc = correlations[frame_i]
        ax_pseudo.set_title(f"Match (r={cc:.2f})", fontsize=7)
        ax_pseudo.axis("off")

        # Overlay
        ax_overlay = axes[row, col_base + 2]
        ax_overlay.imshow(real, cmap="gray", interpolation="nearest")
        ax_overlay.imshow(pseudo_n, cmap="hot", alpha=0.4, interpolation="nearest")
        cvs = matched_cvs[frame_i]
        cv_str = "/".join(f"{v:.0f}" for v in cvs)
        ax_overlay.set_title(f"CVs: {cv_str}Å", fontsize=6)
        ax_overlay.axis("off")

    # Turn off unused axes
    for idx in range(len(selected), n_rows * n_cols):
        row = idx // n_cols
        for c in range(3):
            axes[row, (idx % n_cols) * 3 + c].axis("off")

    fig.suptitle("HS-AFM → Pseudo-AFM Correlation Matching\n"
                 "(Real | Best Match | Overlay + CVs)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"Saved grid: {output_path}")


def create_pdb_overlay_grid(real_frames, matched_indices, pseudo_labels,
                            pdb_files, output_path: Path,
                            n_show: int = 16, image_size: int = 200):
    """Grid with real AFM frames and PDB silhouette contour overlaid."""
    import mdtraj as md
    from scipy.ndimage import gaussian_filter

    n_frames = len(real_frames)
    step = max(1, n_frames // n_show)
    selected = list(range(0, n_frames, step))[:n_show]

    # Precompute PDB silhouettes for unique matched conformers
    unique_indices = np.unique([matched_indices[i] for i in selected])
    # Map pseudo-AFM index back to conformer frame index
    # Each pseudo-AFM batch has 66 images (one per conformer frame)
    # So pseudo_idx % 66 gives the conformer frame index
    n_conformers = len(pdb_files)
    silhouette_cache = {}

    for pseudo_idx in unique_indices:
        frame_idx = int(pseudo_idx) % n_conformers
        frame_idx = min(frame_idx, n_conformers - 1)
        if frame_idx not in silhouette_cache:
            pdb = pdb_files[frame_idx]
            silhouette_cache[frame_idx] = pdb_to_silhouette(pdb, image_size=image_size)

    n_cols = min(4, len(selected))
    n_rows = (len(selected) + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes[np.newaxis, :]
    elif n_cols == 1:
        axes = axes[:, np.newaxis]

    for idx, frame_i in enumerate(selected):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]

        # Real AFM frame (full resolution)
        real = real_frames[frame_i]
        h, w = real.shape
        if h != w:
            s = min(h, w)
            y0, x0 = (h - s) // 2, (w - s) // 2
            real = real[y0:y0+s, x0:x0+s]
        real_norm = (real - real.min()) / (real.max() - real.min() + 1e-8)

        ax.imshow(real_norm, cmap="afmhot", interpolation="bilinear",
                  extent=[0, image_size, image_size, 0])

        # PDB silhouette overlay
        pseudo_idx = matched_indices[frame_i]
        frame_idx = int(pseudo_idx) % n_conformers
        frame_idx = min(frame_idx, n_conformers - 1)
        if frame_idx in silhouette_cache:
            sil = silhouette_cache[frame_idx]
            ax.contour(sil, levels=[0.2, 0.5, 0.8],
                       colors=["cyan"], linewidths=[0.5, 0.8, 1.2],
                       extent=[0, image_size, image_size, 0], alpha=0.7)

        cvs = pseudo_labels[pseudo_idx]
        ax.set_title(f"t={frame_i}\nCVs: {cvs[0]:.0f}/{cvs[1]:.0f}/{cvs[2]:.0f}Å",
                     fontsize=7)
        ax.axis("off")

    for idx in range(len(selected), n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].axis("off")

    fig.suptitle("HS-AFM with Fitted PDB Silhouette Overlay",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"Saved PDB overlay grid: {output_path}")


# ── CV trajectory and statistics ───────────────────────────────────────────

def plot_cv_trajectories(matched_cvs: np.ndarray, correlations: np.ndarray,
                         output_path: Path, cnn_cvs: np.ndarray = None):
    """Plot matched CVs over time with correlation quality."""
    n_frames = len(matched_cvs)
    cv_names = ["αHead↔αTail", "βHead↔αTail", "αHead↔βTail"]
    n_panels = matched_cvs.shape[1] + 1  # CVs + correlation

    fig, axes = plt.subplots(n_panels, 1, figsize=(12, 2.5 * n_panels), sharex=True)
    t = np.arange(n_frames)

    for j in range(matched_cvs.shape[1]):
        ax = axes[j]
        ax.plot(t, matched_cvs[:, j], "b-", linewidth=0.5, alpha=0.7,
                label="Correlation match")
        if cnn_cvs is not None:
            ax.plot(t, cnn_cvs[:, j], "r--", linewidth=0.5, alpha=0.5,
                    label="CNN prediction")
        ax.set_ylabel(f"{cv_names[j]} (Å)", fontsize=9)
        ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=7)

    ax = axes[-1]
    ax.plot(t, correlations, "g-", linewidth=0.5, alpha=0.7)
    ax.set_ylabel("Correlation r", fontsize=9)
    ax.set_xlabel("HS-AFM Frame")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)

    fig.suptitle("Conformational State Trajectory from HS-AFM",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved CV trajectory: {output_path}")


# ── Animated GIF ───────────────────────────────────────────────────────────

def create_animated_gif(real_frames, matched_indices, pseudo_images,
                        matched_cvs, correlations, output_path: Path,
                        image_size: int = 35, fps: int = 10,
                        max_frames: int = 200):
    """Create animated GIF showing real AFM with pseudo-AFM overlay per frame."""
    from PIL import Image
    import io

    n_frames = len(real_frames)
    step = max(1, n_frames // max_frames)
    selected = list(range(0, n_frames, step))[:max_frames]

    pil_frames = []
    for i, frame_i in enumerate(selected):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 3))

        real = preprocess_frame(real_frames[frame_i], image_size)
        pseudo_idx = matched_indices[frame_i]
        pseudo = pseudo_images[pseudo_idx]
        pseudo_n = (pseudo - pseudo.min()) / (pseudo.max() - pseudo.min() + 1e-8)

        ax1.imshow(real, cmap="afmhot", interpolation="nearest")
        ax1.set_title("HS-AFM", fontsize=9)
        ax1.axis("off")

        ax2.imshow(pseudo_n, cmap="afmhot", interpolation="nearest")
        ax2.set_title(f"Matched (r={correlations[frame_i]:.2f})", fontsize=9)
        ax2.axis("off")

        ax3.imshow(real, cmap="gray", interpolation="nearest")
        ax3.imshow(pseudo_n, cmap="hot", alpha=0.45, interpolation="nearest")
        cvs = matched_cvs[frame_i]
        ax3.set_title(f"CVs: {cvs[0]:.0f}/{cvs[1]:.0f}/{cvs[2]:.0f}Å", fontsize=8)
        ax3.axis("off")

        fig.suptitle(f"Frame {frame_i}/{n_frames}", fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.92])

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        pil_frames.append(Image.open(buf).convert("RGB"))

    # Save as GIF
    duration_ms = int(1000 / fps)
    pil_frames[0].save(
        output_path,
        save_all=True,
        append_images=pil_frames[1:],
        duration=duration_ms,
        loop=0,
    )
    print(f"Saved animated GIF: {output_path} ({len(pil_frames)} frames)")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Extract GIF frames
    print("=" * 60)
    print("Loading HS-AFM frames...")
    print("=" * 60)
    real_frames = extract_gif_frames(args.gif)

    # 2. Load pseudo-AFM library
    print("\n" + "=" * 60)
    print("Loading pseudo-AFM library...")
    print("=" * 60)
    pseudo_images, pseudo_labels = load_pseudo_afm_library(args.training_data_dir)

    # 3. Correlation matching
    print("\n" + "=" * 60)
    print("Running correlation matching...")
    print("=" * 60)
    matched_cvs, matched_indices, correlations = correlation_match(
        real_frames, pseudo_images, pseudo_labels,
        args.image_size, device=args.device)

    # Save matched results
    np.save(str(args.output_dir / "matched_cvs.npy"), matched_cvs)
    np.save(str(args.output_dir / "matched_indices.npy"), matched_indices)
    np.save(str(args.output_dir / "correlations.npy"), correlations)

    # 4. Load PDB files for silhouette overlay
    pdb_files = sorted(args.frames_dir.glob("*.pdb"))
    print(f"\nConformer PDB library: {len(pdb_files)} structures")

    # 5. Visualizations
    print("\n" + "=" * 60)
    print("Creating visualizations...")
    print("=" * 60)

    # Grid: real | pseudo | overlay
    create_overlay_grid(
        real_frames, matched_indices, pseudo_images,
        matched_cvs, correlations, pdb_files,
        args.output_dir / "correlation_grid.png",
        n_show=args.n_grid_frames, image_size=args.image_size)

    # PDB silhouette overlay grid
    create_pdb_overlay_grid(
        real_frames, matched_indices, pseudo_labels,
        pdb_files, args.output_dir / "pdb_overlay_grid.png",
        n_show=16, image_size=200)

    # CV trajectory
    cnn_cvs = None
    if args.predicted_cvs and args.predicted_cvs.exists():
        cnn_cvs = np.load(str(args.predicted_cvs))
    plot_cv_trajectories(
        matched_cvs, correlations,
        args.output_dir / "cv_trajectory.png",
        cnn_cvs=cnn_cvs)

    # Animated GIF
    create_animated_gif(
        real_frames, matched_indices, pseudo_images,
        matched_cvs, correlations,
        args.output_dir / "conformer_video.gif",
        image_size=args.image_size, fps=8)

    # Summary
    print("\n" + "=" * 60)
    print("VISUALIZATION COMPLETE")
    print("=" * 60)
    print(f"  Frames analyzed: {len(real_frames)}")
    print(f"  Correlation grid: {args.output_dir / 'correlation_grid.png'}")
    print(f"  PDB overlay grid: {args.output_dir / 'pdb_overlay_grid.png'}")
    print(f"  CV trajectory:    {args.output_dir / 'cv_trajectory.png'}")
    print(f"  Animated video:   {args.output_dir / 'conformer_video.gif'}")


if __name__ == "__main__":
    main()
