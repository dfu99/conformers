#!/usr/bin/env python3
"""Render atomistic PDB models via ChimeraX and overlay on HS-AFM frames.

Two modes:
  1. Fitted mode (--fitted-dir): Uses per-frame fitted coordinates from
     fit_rigid_body.py (exhaustive SO(3) search, 2048+ rotations/frame).
  2. Legacy mode (--rotations-npy): Uses epoch-based rotation lookup from
     the 50-rotation pseudo-AFM library.

Outputs:
  - pdb_atomistic_grid.png: Grid of selected frames with atomistic overlay
  - pdb_atomistic_video.gif: Animated GIF with atomistic PDB per frame
"""
from __future__ import annotations

import argparse
import io
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gif", type=Path, required=True,
                   help="Input HS-AFM GIF")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--frames-dir", type=Path, required=True,
                   help="Directory of protein-only PDB conformer frames")
    p.add_argument("--matched-indices", type=Path, required=True,
                   help="matched_indices.npy from correlation matching")
    p.add_argument("--matched-cvs", type=Path, required=True,
                   help="matched_cvs.npy from correlation matching")
    p.add_argument("--correlations", type=Path, default=None,
                   help="correlations.npy from correlation matching (legacy mode)")
    p.add_argument("--rotations-npy", type=Path, default=None,
                   help="rotations.npy: SO(3) matrices (n_epochs, 3, 3) (legacy mode)")
    p.add_argument("--fitted-dir", type=Path, default=None,
                   help="Directory from fit_rigid_body.py with fitted_coords.npy etc.")
    p.add_argument("--render-size", type=int, default=400,
                   help="ChimeraX render resolution (square)")
    p.add_argument("--max-gif-frames", type=int, default=150,
                   help="Max frames in animated GIF")
    p.add_argument("--n-grid-frames", type=int, default=16,
                   help="Frames to show in grid")
    p.add_argument("--fps", type=int, default=8)
    p.add_argument("--style", default="cartoon",
                   choices=["cartoon", "surface", "ribbon", "sphere"],
                   help="ChimeraX representation style")
    p.add_argument("--n-conformers", type=int, default=264,
                   help="Number of conformer frames in PDB library")
    return p.parse_args()


def extract_gif_frames(gif_path: Path) -> list[np.ndarray]:
    """Extract all frames from a GIF as grayscale arrays."""
    gif = Image.open(gif_path)
    frames = []
    try:
        while True:
            frames.append(np.array(gif.convert("L"), dtype=np.float32))
            gif.seek(gif.tell() + 1)
    except EOFError:
        pass
    print(f"Extracted {len(frames)} frames from {gif_path.name}")
    return frames


def pre_rotate_pdb(pdb_path: Path, rotation_matrix: np.ndarray,
                   output_path: Path) -> bool:
    """Apply SO(3) rotation to PDB coordinates and save.

    Matches the pseudo-AFM rendering: xyz_rotated = xyz @ R.T,
    then center and shift Z minimum to 0.
    """
    t = md.load(str(pdb_path))
    keep = t.topology.select("protein")
    t = t.atom_slice(keep)
    xyz = t.xyz[0]  # (N_atoms, 3)
    # Apply rotation (same as pseudo-AFM: xyz @ R.T)
    xyz_rotated = xyz @ rotation_matrix.T
    # Center
    xyz_rotated -= xyz_rotated.mean(axis=0)
    # Z-shift: minimum z to 0
    xyz_rotated[:, 2] -= xyz_rotated[:, 2].min()
    t.xyz[0] = xyz_rotated
    t.save_pdb(str(output_path))
    return True


def batch_render_rotated(pdb_dir: Path, full_indices: list[int],
                         rotations: np.ndarray, n_conformers: int,
                         render_dir: Path, render_size: int = 400,
                         style: str = "cartoon") -> dict[int, Path]:
    """Batch-render conformers at their matched SO(3) orientations.

    For each unique full_index, pre-rotates the PDB with the epoch's
    rotation matrix, then renders via ChimeraX.

    Returns mapping: full_index -> rendered PNG path
    """
    render_dir.mkdir(parents=True, exist_ok=True)
    tmp_pdb_dir = render_dir / "rotated_pdbs"
    tmp_pdb_dir.mkdir(exist_ok=True)

    # Determine what needs rendering
    to_render = []
    for fidx in full_indices:
        out_path = render_dir / f"render_f{fidx:05d}.png"
        if out_path.exists():
            continue
        conf_idx = fidx % n_conformers
        epoch = fidx // n_conformers
        pdb_path = pdb_dir / f"frame_{conf_idx:04d}.pdb"
        if not pdb_path.exists():
            continue
        if epoch >= len(rotations):
            continue
        to_render.append(fidx)

    if to_render:
        # Pre-rotate PDB files
        print(f"  Pre-rotating {len(to_render)} PDB files...")
        rotated_pdbs = {}
        for fidx in to_render:
            conf_idx = fidx % n_conformers
            epoch = fidx // n_conformers
            pdb_path = pdb_dir / f"frame_{conf_idx:04d}.pdb"
            rot_pdb = tmp_pdb_dir / f"rotated_f{fidx:05d}.pdb"
            if pre_rotate_pdb(pdb_path, rotations[epoch], rot_pdb):
                rotated_pdbs[fidx] = rot_pdb

        # Style commands
        style_cmds = {
            "cartoon": [
                "cartoon",
                "color /A cornflowerblue",
                "color /B salmon",
                "lighting soft",
                "set silhouetteWidth 2",
                "set silhouette true",
            ],
            "surface": [
                "surface",
                "color /A cornflowerblue",
                "color /B salmon",
                "transparency 20",
                "lighting soft",
            ],
        }
        cmds = style_cmds.get(style, style_cmds["cartoon"])

        # Write ChimeraX batch script
        py_script = render_dir / "batch_render.py"
        py_lines = ["from chimerax.core.commands import run", ""]
        for fidx, rot_pdb in rotated_pdbs.items():
            out_path = render_dir / f"render_f{fidx:05d}.png"
            py_lines.append(f'run(session, "open {rot_pdb}")')
            for cmd in cmds:
                py_lines.append(f'run(session, "{cmd}")')
            # No 'view orient' — use pre-rotated coordinates directly
            py_lines.append('run(session, "view")')
            py_lines.append(
                f'run(session, "save {out_path} width {render_size} '
                f'height {render_size} transparentBackground true supersample 3")'
            )
            py_lines.append('run(session, "close all")')

        py_script.write_text("\n".join(py_lines))

        print(f"Rendering {len(to_render)} orientation-specific conformers via ChimeraX "
              f"({len(full_indices) - len(to_render)} cached)...")
        result = subprocess.run(
            ["chimerax", "--offscreen", "--exit", "--script", str(py_script)],
            capture_output=True, text=True, timeout=1800
        )
        if result.returncode != 0:
            print(f"ChimeraX stderr (last 500 chars): {result.stderr[-500:]}")

    # Collect all rendered files
    rendered = {}
    for fidx in full_indices:
        out_path = render_dir / f"render_f{fidx:05d}.png"
        if out_path.exists():
            rendered[fidx] = out_path

    print(f"Successfully rendered {len(rendered)}/{len(full_indices)} conformers")
    return rendered


def batch_render_fitted(topology_path: Path, fitted_coords: np.ndarray,
                        frame_indices: list[int], render_dir: Path,
                        render_size: int = 400,
                        style: str = "cartoon") -> dict[int, Path]:
    """Render frames using per-frame fitted coordinates from rigid-body fitting.

    Each frame has its own pre-rotated, centered, z-shifted coordinates.
    No epoch/conformer lookup needed — coordinates are ready to write as PDBs.

    Returns mapping: frame_index -> rendered PNG path
    """
    render_dir.mkdir(parents=True, exist_ok=True)
    tmp_pdb_dir = render_dir / "fitted_pdbs"
    tmp_pdb_dir.mkdir(exist_ok=True)

    ref = md.load(str(topology_path))

    to_render = []
    for fi in frame_indices:
        out_path = render_dir / f"render_frame{fi:05d}.png"
        if not out_path.exists():
            to_render.append(fi)

    if to_render:
        print(f"  Writing {len(to_render)} fitted PDB files...")
        fitted_pdbs = {}
        for fi in to_render:
            pdb_path = tmp_pdb_dir / f"fitted_frame{fi:05d}.pdb"
            t = ref[0]  # single-frame copy
            t.xyz[0] = fitted_coords[fi]
            t.save_pdb(str(pdb_path))
            fitted_pdbs[fi] = pdb_path

        style_cmds = {
            "cartoon": [
                "cartoon",
                "color /A cornflowerblue",
                "color /B salmon",
                "lighting soft",
                "set silhouetteWidth 2",
                "set silhouette true",
            ],
            "surface": [
                "surface",
                "color /A cornflowerblue",
                "color /B salmon",
                "transparency 20",
                "lighting soft",
            ],
        }
        cmds = style_cmds.get(style, style_cmds["cartoon"])

        py_script = render_dir / "batch_render_fitted.py"
        py_lines = ["from chimerax.core.commands import run", ""]
        for fi, pdb_path in fitted_pdbs.items():
            out_path = render_dir / f"render_frame{fi:05d}.png"
            py_lines.append(f'run(session, "open {pdb_path}")')
            for cmd in cmds:
                py_lines.append(f'run(session, "{cmd}")')
            py_lines.append('run(session, "view")')
            py_lines.append(
                f'run(session, "save {out_path} width {render_size} '
                f'height {render_size} transparentBackground true supersample 3")'
            )
            py_lines.append('run(session, "close all")')

        py_script.write_text("\n".join(py_lines))

        print(f"Rendering {len(to_render)} fitted conformers via ChimeraX "
              f"({len(frame_indices) - len(to_render)} cached)...")
        result = subprocess.run(
            ["chimerax", "--offscreen", "--exit", "--script", str(py_script)],
            capture_output=True, text=True, timeout=1800
        )
        if result.returncode != 0:
            print(f"ChimeraX stderr (last 500 chars): {result.stderr[-500:]}")

    rendered = {}
    for fi in frame_indices:
        out_path = render_dir / f"render_frame{fi:05d}.png"
        if out_path.exists():
            rendered[fi] = out_path

    print(f"Successfully rendered {len(rendered)}/{len(frame_indices)} conformers")
    return rendered


def load_render(path: Path, target_size: int = None) -> np.ndarray:
    """Load a ChimeraX render as RGBA numpy array."""
    img = Image.open(path).convert("RGBA")
    if target_size:
        img = img.resize((target_size, target_size), Image.LANCZOS)
    return np.array(img, dtype=np.float32) / 255.0


def _prepare_afm(afm_frame: np.ndarray, size: int) -> tuple[np.ndarray, np.ndarray]:
    """Crop to square, resize, colormap. Returns (grayscale_0_1, colored_rgb)."""
    afm_pil = Image.fromarray(afm_frame.astype(np.uint8) if afm_frame.max() > 1
                               else (afm_frame * 255).astype(np.uint8))
    aw, ah = afm_pil.size
    if aw != ah:
        s = min(aw, ah)
        left = (aw - s) // 2
        top = (ah - s) // 2
        afm_pil = afm_pil.crop((left, top, left + s, top + s))
    afm_pil = afm_pil.resize((size, size), Image.LANCZOS)
    gray = np.array(afm_pil, dtype=np.float32) / 255.0
    colored = plt.cm.afmhot(gray)[:, :, :3]
    return gray, colored


def composite_overlay(afm_frame: np.ndarray, render_rgba: np.ndarray,
                      alpha_boost: float = 1.2) -> np.ndarray:
    """Composite transparent PDB render over AFM frame.

    Returns RGB image as float [0,1] array of shape (H, W, 3).
    """
    h, w = render_rgba.shape[:2]
    _, afm_colored = _prepare_afm(afm_frame, h)

    rgb = render_rgba[:, :, :3]
    alpha = np.clip(render_rgba[:, :, 3:4] * alpha_boost, 0, 1)

    composite = afm_colored * (1 - alpha) + rgb * alpha
    return np.clip(composite, 0, 1)


def classify_state(cvs: np.ndarray) -> str:
    """Classify conformational state from mean CV distance."""
    mean_cv = cvs.mean()
    if mean_cv < 60:
        return "Compact"
    elif mean_cv < 70:
        return "Intermediate"
    else:
        return "Extended"


def create_atomistic_grid(real_frames: list[np.ndarray],
                          matched_cvs: np.ndarray,
                          correlations: np.ndarray,
                          rendered: dict[int, Path],
                          frame_indices: list[int],
                          render_key_fn,
                          output_path: Path,
                          n_show: int = 16,
                          render_size: int = 400):
    """Grid: real AFM | atomistic overlay per frame."""
    n_total = len(frame_indices)
    step = max(1, n_total // n_show)
    selected_pos = list(range(0, n_total, step))[:n_show]

    n_cols = min(4, len(selected_pos))
    n_rows = (len(selected_pos) + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols * 2, figsize=(n_cols * 4, n_rows * 2.2))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    for idx, pos in enumerate(selected_pos):
        frame_i = frame_indices[pos]
        row = idx // n_cols
        col_base = (idx % n_cols) * 2

        real = real_frames[frame_i]
        rkey = render_key_fn(frame_i)

        # Real AFM
        ax_real = axes[row, col_base]
        h, w = real.shape
        if h != w:
            s = min(h, w)
            y0, x0 = (h - s) // 2, (w - s) // 2
            real_sq = real[y0:y0+s, x0:x0+s]
        else:
            real_sq = real
        ax_real.imshow(real_sq, cmap="afmhot", interpolation="bilinear")
        ax_real.set_title(f"Frame {frame_i}", fontsize=8, fontweight="bold")
        ax_real.axis("off")

        # Atomistic overlay — orientation baked into render
        ax_overlay = axes[row, col_base + 1]
        if rkey in rendered:
            render_rgba = load_render(rendered[rkey], target_size=render_size)
            comp = composite_overlay(real, render_rgba)
            ax_overlay.imshow(comp, interpolation="bilinear")
        else:
            ax_overlay.imshow(real_sq, cmap="afmhot", interpolation="bilinear")

        cvs = matched_cvs[frame_i]
        state = classify_state(cvs)
        cc = correlations[frame_i]
        ax_overlay.set_title(f"{state} (r={cc:.2f})\n"
                             f"CVs: {cvs[0]:.0f}/{cvs[1]:.0f}/{cvs[2]:.0f}Å",
                             fontsize=6)
        ax_overlay.axis("off")

    for idx in range(len(selected_pos), n_rows * n_cols):
        row = idx // n_cols
        for c in range(2):
            axes[row, (idx % n_cols) * 2 + c].axis("off")

    fig.suptitle("HS-AFM with Atomistic PDB Overlay (ChimeraX)\n"
                 "(Real AFM | Fitted Conformer — SO(3) oriented)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved atomistic grid: {output_path}")


def create_atomistic_gif(real_frames: list[np.ndarray],
                         matched_cvs: np.ndarray,
                         correlations: np.ndarray,
                         rendered: dict[int, Path],
                         frame_indices: list[int],
                         render_key_fn,
                         output_path: Path,
                         render_size: int = 400,
                         fps: int = 8):
    """Animated GIF: real AFM | atomistic overlay | CV bar chart per frame."""
    n_all = len(real_frames)
    cv_names = ["\u03b1H\u2194\u03b1T", "\u03b2H\u2194\u03b1T", "\u03b1H\u2194\u03b2T"]
    state_colors = {"Compact": "#3366CC", "Intermediate": "#FF9900",
                    "Extended": "#CC3333"}

    pil_frames = []
    for i, frame_i in enumerate(frame_indices):
        rkey = render_key_fn(frame_i)
        cvs = matched_cvs[frame_i]
        cc = correlations[frame_i]
        state = classify_state(cvs)

        fig = plt.figure(figsize=(10, 3.5))
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.6], wspace=0.15)

        # Panel 1: Real AFM
        ax1 = fig.add_subplot(gs[0])
        real = real_frames[frame_i]
        h, w = real.shape
        if h != w:
            s = min(h, w)
            y0, x0 = (h - s) // 2, (w - s) // 2
            real_sq = real[y0:y0+s, x0:x0+s]
        else:
            real_sq = real
        ax1.imshow(real_sq, cmap="afmhot", interpolation="bilinear")
        ax1.set_title("HS-AFM", fontsize=9, fontweight="bold")
        ax1.axis("off")

        # Panel 2: Atomistic overlay — orientation baked into render
        ax2 = fig.add_subplot(gs[1])
        if rkey in rendered:
            render_rgba = load_render(rendered[rkey], target_size=render_size)
            comp = composite_overlay(real, render_rgba)
            ax2.imshow(comp, interpolation="bilinear")
        else:
            ax2.imshow(real_sq, cmap="afmhot", interpolation="bilinear")
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
            ax3.text(v + 0.3, bar.get_y() + bar.get_height()/2,
                     f"{v:.1f}", va="center", fontsize=7)
        ax3.tick_params(labelsize=7)

        fig.suptitle(f"Frame {frame_i}/{n_all}  |  \u03b1v\u03b23 Integrin Conformer",
                     fontsize=10, fontweight="bold", y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120, bbox_inches="tight",
                    facecolor="white")
        plt.close(fig)
        buf.seek(0)
        pil_frames.append(Image.open(buf).convert("RGB"))

        if (i + 1) % 30 == 0:
            print(f"  GIF frame {i+1}/{len(frame_indices)}")

    duration_ms = int(1000 / fps)
    pil_frames[0].save(
        output_path,
        save_all=True,
        append_images=pil_frames[1:],
        duration=duration_ms,
        loop=0,
    )
    print(f"Saved atomistic GIF: {output_path} ({len(pil_frames)} frames)")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    use_fitted = args.fitted_dir is not None

    if use_fitted:
        # Fitted mode: per-frame coordinates from exhaustive SO(3) fitting
        fitted_coords = np.load(str(args.fitted_dir / "fitted_coords.npy"))
        correlations = np.load(str(args.fitted_dir / "fitted_correlations.npy"))
        topology_path = args.fitted_dir / "topology.pdb"
        n_frames = len(fitted_coords)
        print(f"Fitted mode: {n_frames} frames with per-frame SO(3) orientations")
        print(f"  Correlation: mean={correlations.mean():.4f} "
              f"min={correlations.min():.4f} max={correlations.max():.4f}")
    else:
        # Legacy mode: epoch-based rotation lookup
        if args.rotations_npy is None or args.correlations is None:
            raise ValueError("Legacy mode requires --rotations-npy and --correlations")
        matched_indices = np.load(str(args.matched_indices))
        correlations = np.load(str(args.correlations))
        rotations = np.load(str(args.rotations_npy))
        n_frames = len(matched_indices)
        print(f"Legacy mode: {n_frames} frame matches, "
              f"{len(rotations)} rotation matrices")

    matched_cvs = np.load(str(args.matched_cvs))

    # Extract GIF frames
    real_frames = extract_gif_frames(args.gif)

    # Select frames for GIF
    step = max(1, n_frames // args.max_gif_frames)
    gif_frame_indices = list(range(0, n_frames, step))[:args.max_gif_frames]

    if use_fitted:
        # Render using fitted coordinates — keyed by frame index
        render_dir = args.output_dir / "chimerax_renders_fitted"
        rendered = batch_render_fitted(
            topology_path, fitted_coords, gif_frame_indices,
            render_dir, render_size=args.render_size, style=args.style
        )
        # For grid/GIF: lookup by frame index directly
        render_key_fn = lambda frame_i: frame_i
    else:
        # Render using epoch rotation lookup — keyed by full_index
        needed_full = sorted(set(
            int(matched_indices[i]) for i in gif_frame_indices
        ))
        print(f"Need {len(needed_full)} unique (conformer, orientation) renders")
        render_dir = args.output_dir / "chimerax_renders_rotated"
        rendered = batch_render_rotated(
            args.frames_dir, needed_full, rotations, args.n_conformers,
            render_dir, render_size=args.render_size, style=args.style
        )
        render_key_fn = lambda frame_i: int(matched_indices[frame_i])

    # Create grid
    create_atomistic_grid(
        real_frames, matched_cvs, correlations,
        rendered, gif_frame_indices, render_key_fn,
        args.output_dir / "pdb_atomistic_grid.png",
        n_show=args.n_grid_frames, render_size=args.render_size
    )

    # Create animated GIF
    create_atomistic_gif(
        real_frames, matched_cvs, correlations,
        rendered, gif_frame_indices, render_key_fn,
        args.output_dir / "pdb_atomistic_video.gif",
        render_size=args.render_size, fps=args.fps
    )

    print("\n" + "=" * 60)
    print(f"ATOMISTIC OVERLAY COMPLETE ({'fitted SO(3)' if use_fitted else 'legacy'})")
    print(f"  Grid:  {args.output_dir / 'pdb_atomistic_grid.png'}")
    print(f"  Video: {args.output_dir / 'pdb_atomistic_video.gif'}")
    print("=" * 60)


if __name__ == "__main__":
    main()
