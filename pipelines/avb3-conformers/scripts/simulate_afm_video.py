#!/usr/bin/env python3
"""Forward-render fitted PDB trajectory → simulated HS-AFM video.

Takes the v7 fitted coordinates (already head-anchored to tracked AFM
position) and runs each frame through the pseudo-AFM pipeline with the
tip/resolution parameters we calibrated against the real data. Output is
a simulated AFM GIF alongside per-frame sim-vs-real correlation.

Validates pipeline self-consistency: if the sim AFM matches the real,
the conformer library + imaging model are jointly adequate.
"""
from __future__ import annotations

import argparse
import io
from pathlib import Path

import numpy as np
import torch
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fitted-dir", type=Path, required=True,
                   help="Directory from fit_with_head_tracking.py containing "
                        "fitted_coords_smooth.npy and topology.pdb")
    p.add_argument("--gif", type=Path, required=True, help="Real AFM GIF for comparison")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--afmfold-root", type=Path, default=Path("/home2/Documents/code/afmfold"))
    p.add_argument("--width", type=int, default=35)
    p.add_argument("--height", type=int, default=35)
    p.add_argument("--resolution-nm", type=float, default=0.98)
    p.add_argument("--tip-radius", type=float, default=1.5, help="nm")
    p.add_argument("--tip-angle", type=float, default=20.0, help="deg")
    p.add_argument("--noise-nm", type=float, default=0.1)
    p.add_argument("--max-frames", type=int, default=0, help="0 = all")
    p.add_argument("--device", default="cpu")
    p.add_argument("--n-grid-frames", type=int, default=16)
    return p.parse_args()


def extract_gif_frames(gif_path, target_size):
    gif = Image.open(gif_path)
    frames = []
    try:
        while True:
            arr = np.array(gif.convert("L"), dtype=np.float32)
            h, w = arr.shape
            if h != w:
                s = min(h, w)
                y0, x0 = (h - s) // 2, (w - s) // 2
                arr = arr[y0:y0+s, x0:x0+s]
            img = Image.fromarray(arr)
            img = img.resize((target_size, target_size), Image.BILINEAR)
            arr = np.array(img, dtype=np.float32)
            arr -= arr.min()
            if arr.max() > 0:
                arr /= arr.max()
            frames.append(arr)
            gif.seek(gif.tell() + 1)
    except EOFError:
        pass
    return np.stack(frames)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    import sys
    sys.path.insert(0, str(args.afmfold_root / "src"))
    from afmfold.images import generate_landscape, idilation, generate_tip_shape

    # Load fitted trajectory
    fitted = np.load(str(args.fitted_dir / "fitted_coords_smooth.npy"))
    import json
    with open(args.fitted_dir / "fitting_metadata.json") as f:
        meta = json.load(f)
    skip = meta.get("skip_frames", 0)
    n = fitted.shape[0]
    print(f"Loaded {n} fitted frames (skip_frames={skip})")

    # Load real AFM for comparison (skip the same initial frames)
    real = extract_gif_frames(args.gif, args.width)
    real_fitted = real[skip:skip + n]
    print(f"Real AFM fitted-aligned: {real_fitted.shape}")

    # Determine image center to compare (head tracked to this region)
    # fitted_coords are in nm with head at tracked AFM head position.
    # Need to shift so the head pixel maps to image center.
    head_positions_px = np.load(str(args.fitted_dir / "head_positions_px.npy"))
    if len(head_positions_px) != n:
        head_positions_px = head_positions_px[:n]

    # Sim rendering setup
    device = torch.device(args.device)
    # Grid edges: must match afmfold convention exactly
    # For W=35: xedges = -17..+17 giving 35 edges (wait — W+1=36, actually -18..+17 = 36 entries)
    # Match afmfold process_frames_to_afm exactly:
    w_half = int(args.width / 2)
    h_half = int(args.height / 2)
    xedges = args.resolution_nm * torch.arange(
        -w_half - 1, args.width - w_half, dtype=torch.float32, device=device)
    yedges = args.resolution_nm * torch.arange(
        -h_half - 1, args.height - h_half, dtype=torch.float32, device=device)
    tip = generate_tip_shape(
        radius=torch.tensor(args.tip_radius, dtype=torch.float32),
        angle=torch.tensor(args.tip_angle, dtype=torch.float32),
        device=device)

    # Render each frame
    print(f"Rendering {n} sim-AFM frames at {args.width}x{args.height} "
          f"(tip R={args.tip_radius}nm, angle={args.tip_angle}°)...")
    N = n if args.max_frames == 0 else min(n, args.max_frames)
    sim_images = np.zeros((N, args.height, args.width), dtype=np.float32)
    # fitted_coords are already centered — head in range [-7, +7] nm from
    # origin. Sim grid is also centered on (0, 0). To reproduce the real
    # AFM's head position (head at px (head_positions_px[i])), we need to
    # shift so that the head ends up at the right grid location.
    #
    # In the real AFM pre-processing, head_positions_px[i] is in the
    # 35x35 image coord. Pixel (px, py) maps to grid nm
    # ((px - width/2) * res_nm, (py - height/2) * res_nm).
    # So the head should be at nm ((head_px - width/2) * res_nm, ...).
    for i in range(N):
        coords = fitted[i]
        # Find head centroid in current coords
        # (head is already in the expected coord system; we just need to
        # translate it from current position to target grid position)
        shifted = coords.copy()
        shifted[:, 2] -= shifted[:, 2].min()

        xyz = torch.from_numpy(shifted).float().unsqueeze(0).to(device)
        pure, _, _ = generate_landscape(xyz, xedges, yedges)
        pure = pure.reshape((-1, args.height, args.width))
        afm_img = idilation(pure, tip)
        if args.noise_nm > 0:
            noise = torch.from_numpy(
                np.random.normal(0, args.noise_nm, afm_img.shape).astype(np.float32)
            ).to(device)
            afm_img = afm_img + noise
        arr = afm_img[0].detach().cpu().numpy()
        arr -= arr.min()
        if arr.max() > 0:
            arr /= arr.max()
        sim_images[i] = arr

        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{N}")

    np.save(str(args.output_dir / "sim_images.npy"), sim_images)
    print(f"Saved sim_images.npy")

    # Per-frame correlation between sim and real
    real_flat = real_fitted[:N].reshape(N, -1)
    sim_flat = sim_images.reshape(N, -1)
    real_flat = real_flat - real_flat.mean(axis=1, keepdims=True)
    sim_flat = sim_flat - sim_flat.mean(axis=1, keepdims=True)
    real_norm = np.linalg.norm(real_flat, axis=1, keepdims=True)
    sim_norm = np.linalg.norm(sim_flat, axis=1, keepdims=True)
    real_norm[real_norm < 1e-9] = 1.0
    sim_norm[sim_norm < 1e-9] = 1.0
    correlations = np.sum(
        (real_flat / real_norm) * (sim_flat / sim_norm), axis=1
    )
    np.save(str(args.output_dir / "sim_vs_real_correlation.npy"), correlations)
    print(f"Sim-vs-real corr: mean={correlations.mean():.3f}, "
          f"std={correlations.std():.3f}, median={np.median(correlations):.3f}")

    # Save sim GIF + side-by-side GIF
    def to_pil(arr, scale=8):
        arr = (arr * 255).clip(0, 255).astype(np.uint8)
        img = Image.fromarray(arr, mode="L")
        return img.resize((arr.shape[1] * scale, arr.shape[0] * scale), Image.NEAREST)

    pil_sim = [to_pil(sim_images[i]) for i in range(N)]
    pil_sim[0].save(str(args.output_dir / "sim_afm.gif"), save_all=True,
                     append_images=pil_sim[1:], duration=80, loop=0)
    print(f"Saved sim_afm.gif")

    # Side-by-side (real | sim)
    pil_side = []
    for i in range(N):
        r = to_pil(real_fitted[i])
        s = to_pil(sim_images[i])
        combined = Image.new("L", (r.width + s.width, r.height))
        combined.paste(r, (0, 0))
        combined.paste(s, (r.width, 0))
        pil_side.append(combined)
    pil_side[0].save(str(args.output_dir / "sim_vs_real.gif"), save_all=True,
                      append_images=pil_side[1:], duration=80, loop=0)
    print(f"Saved sim_vs_real.gif")

    # Correlation time-course
    fig, ax = plt.subplots(figsize=(10, 4))
    t = np.arange(N) * 0.1
    ax.plot(t, correlations, color="#2b8cbe", linewidth=0.8)
    ax.axhline(correlations.mean(), color="red", linestyle="--",
                label=f"mean={correlations.mean():.3f}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Sim-vs-Real cosine correlation")
    ax.set_title(f"Forward-rendered AFM vs Linz real — {args.fitted_dir.name}")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(str(args.output_dir / "sim_vs_real_corr.png"), dpi=150)
    print(f"Saved sim_vs_real_corr.png")


if __name__ == "__main__":
    main()
