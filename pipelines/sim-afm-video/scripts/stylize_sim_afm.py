#!/usr/bin/env python3
"""Single-canvas HS-AFM stylization for forward-rendered sim images.

Loads sim_images.npy (the output of simulate_afm_video.py) and applies
realism layers — substrate noise, anisotropic blur, surface slant,
row-wise baseline jitter, partial-width flash streaks, localized
distortions, copper colormap — all on ONE canvas at one resolution.

This deliberately avoids the v11 bug: that script generated noise on
a larger canvas and pasted the molecule render into a smaller inset,
producing a visible rectangular border whenever the molecule extended
past the inset edge. Here every operation happens on the same array.

Usage:
    python stylize_sim_afm.py \\
        --sim-images results/afm_pipeline/sim_afm/video2/sim_images.npy \\
        --fitted-dir results/afm_pipeline/v7_smoothed_final/video2 \\
        --output results/afm_pipeline/sim_afm/video2/sim_afm_copper_v12.gif \\
        --upscale 8
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, zoom


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sim-images", type=Path, required=True,
                   help="(N, H, W) np array of base sim AFM images in [0,1]")
    p.add_argument("--fitted-dir", type=Path, default=None,
                   help="Optional: provides fitted_coords for tail-direction slant")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--upscale", type=int, default=8,
                   help="Bilinear upscale factor for the working canvas")
    p.add_argument("--substrate-noise-std", type=float, default=0.022)
    p.add_argument("--substrate-baseline", type=float, default=0.15)
    p.add_argument("--blur-sigma-x", type=float, default=1.0,
                   help="Blur σ_x (working pixels). Default isotropic "
                        "(σx=σy=1.0). The earlier σx=1.2/σy=0.7 squashed "
                        "the molecule in x once orientation was locked.")
    p.add_argument("--blur-sigma-y", type=float, default=1.0,
                   help="Blur σ_y (working pixels). Set isotropic with "
                        "blur-sigma-x by default.")
    p.add_argument("--slant-amplitude", type=float, default=0.045)
    p.add_argument("--row-jitter-std", type=float, default=0.008)
    p.add_argument("--flash-streak-prob", type=float, default=0.02)
    p.add_argument("--soft-clip", type=float, default=0.92)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--max-frames", type=int, default=0,
                   help="0 = all frames")
    p.add_argument("--frame-duration-ms", type=int, default=80)
    return p.parse_args()


COPPER = np.array(
    [
        [0, 0, 0],
        [50, 32, 16],
        [110, 70, 38],
        [170, 105, 60],
        [222, 145, 86],
        [255, 188, 130],
        [255, 226, 178],
    ],
    dtype=np.float32,
)


def copper_colormap(arr: np.ndarray) -> np.ndarray:
    """Map a 2D float32 array in [0,1] to a (H, W, 3) uint8 RGB image
    using a piecewise-linear copper colormap."""
    arr_clip = np.clip(arr, 0.0, 1.0)
    n = COPPER.shape[0]
    pos = arr_clip * (n - 1)
    lo = np.floor(pos).astype(int).clip(0, n - 1)
    hi = np.minimum(lo + 1, n - 1)
    frac = (pos - lo)[..., None]
    rgb = COPPER[lo] * (1 - frac) + COPPER[hi] * frac
    return rgb.clip(0, 255).astype(np.uint8)


def soft_clip(arr: np.ndarray, threshold: float) -> np.ndarray:
    """Smoothly compress values above `threshold` toward 1.0."""
    out = arr.copy()
    mask = arr > threshold
    if mask.any():
        excess = arr[mask] - threshold
        compressed = threshold + (1.0 - threshold) * np.tanh(
            excess / (1.0 - threshold)
        )
        out[mask] = compressed
    return out


def estimate_tail_direction(fitted_coords_path: Path):
    """Approximate per-frame tail-direction angle (rad) for slant tilt.

    Returns None if not available.
    """
    if fitted_coords_path is None or not fitted_coords_path.exists():
        return None
    fc = np.load(str(fitted_coords_path))  # [T, N_atoms, 3]
    # Use the mean Z of the first 1/3 of atoms (head) vs last 1/3 (tail)
    # to get the head→tail vector projected onto XY.
    n_atoms = fc.shape[1]
    head = fc[:, : n_atoms // 3].mean(axis=1)
    tail = fc[:, -n_atoms // 3 :].mean(axis=1)
    direction = tail[:, :2] - head[:, :2]
    angle = np.arctan2(direction[:, 1], direction[:, 0])
    return angle


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    sim = np.load(args.sim_images)
    if sim.ndim != 3:
        raise ValueError(f"Expected (N, H, W), got {sim.shape}")
    N_full, H_low, W_low = sim.shape
    N = N_full if args.max_frames <= 0 else min(N_full, args.max_frames)
    print(f"Stylizing {N} frames, base resolution {H_low}x{W_low}")

    H = H_low * args.upscale
    W = W_low * args.upscale
    print(f"Working canvas: {H}x{W} (upscale {args.upscale})")

    # Optional tail-direction slant.
    tail_angle = None
    if args.fitted_dir is not None:
        fpath = args.fitted_dir / "fitted_coords_median_reanchored.npy"
        if not fpath.exists():
            fpath = args.fitted_dir / "fitted_coords_smooth.npy"
        tail_angle = estimate_tail_direction(fpath)
        if tail_angle is not None:
            print(f"Got tail-direction angles (T={len(tail_angle)})")

    # Pre-build one large noise field per frame at full canvas size, so
    # the molecule signal and the noise share the same canvas — no inset
    # crops, no pastes, no border possible.
    frames_rgb = []
    yy, xx = np.indices((H, W))
    for i in range(N):
        # Upscale base sim to the working canvas (bilinear via scipy.zoom).
        base_low = sim[i].astype(np.float32)
        upscaled = zoom(base_low, (args.upscale, args.upscale), order=1)
        upscaled = upscaled[:H, :W]  # safety crop if zoom rounded up

        # Substrate noise on the SAME canvas.
        noise = rng.normal(0, args.substrate_noise_std, (H, W)).astype(np.float32)
        # Smooth the noise to match the molecule's anisotropic point spread.
        noise = gaussian_filter(noise, sigma=(args.blur_sigma_y * 0.6,
                                              args.blur_sigma_x * 0.6))
        # Combine: molecule ON TOP of substrate baseline + noise.
        canvas = args.substrate_baseline + noise + upscaled

        # Anisotropic blur to imitate scan-direction PSF — applied to
        # the whole canvas (substrate + molecule together), so noise and
        # signal share the same response.
        canvas = gaussian_filter(canvas, sigma=(args.blur_sigma_y,
                                                args.blur_sigma_x))

        # Tail-direction-correlated surface slant.
        if tail_angle is not None and i < len(tail_angle):
            ang = tail_angle[i]
            grad_x = np.cos(ang)
            grad_y = np.sin(ang)
            # Slant amplitude in [-amp, +amp] across the canvas.
            slant = (
                args.slant_amplitude
                * (grad_x * (xx / W - 0.5) + grad_y * (yy / H - 0.5))
            )
            canvas = canvas + slant.astype(np.float32)

        # Row-wise baseline jitter (all rows, full width — single canvas).
        row_jitter = rng.normal(0, args.row_jitter_std, H).astype(np.float32)
        canvas = canvas + row_jitter[:, None]

        # Partial-width flash streaks: rare bright horizontal stripes that
        # span only part of the row width. Applied across the SAME canvas.
        if rng.random() < args.flash_streak_prob:
            streak_row = rng.integers(0, H)
            streak_height = rng.integers(1, max(2, H // 30))
            streak_width = rng.integers(W // 5, W)
            streak_x0 = rng.integers(0, max(1, W - streak_width))
            streak_amp = rng.uniform(0.04, 0.10)
            canvas[streak_row : streak_row + streak_height,
                   streak_x0 : streak_x0 + streak_width] += streak_amp

        # Localized horizontal distortion on rows that touch the molecule.
        # Apply as a small per-row x-shift drawn from a smoothed random
        # field — only where there's actual molecule signal.
        molecule_mask = upscaled > (0.05 * upscaled.max() + 1e-6)
        if molecule_mask.any():
            rows_with_molecule = np.where(molecule_mask.any(axis=1))[0]
            if len(rows_with_molecule):
                shifts = rng.normal(0, 1.5, len(rows_with_molecule))
                shifts = gaussian_filter(shifts, sigma=2.0)
                for rr_idx, r in enumerate(rows_with_molecule):
                    s = int(round(shifts[rr_idx]))
                    if s != 0:
                        canvas[r] = np.roll(canvas[r], s)

        canvas = soft_clip(canvas, args.soft_clip)
        canvas = np.clip(canvas, 0.0, 1.0)

        rgb = copper_colormap(canvas)
        frames_rgb.append(Image.fromarray(rgb))

        if (i + 1) % 200 == 0:
            print(f"  styled {i+1}/{N}")

    print(f"Saving {len(frames_rgb)} frames to {args.output}")
    frames_rgb[0].save(
        str(args.output),
        save_all=True,
        append_images=frames_rgb[1:],
        duration=args.frame_duration_ms,
        loop=0,
        optimize=False,
    )
    # Also dump four sample frames as PNGs so the inset-border can be
    # visually audited.
    audit_dir = args.output.with_suffix("")
    audit_dir.mkdir(exist_ok=True)
    for label, idx in zip(("f100", "f400", "f700", "f1100"), (100, 400, 700, 1100)):
        if idx < len(frames_rgb):
            frames_rgb[idx].save(str(audit_dir / f"{label}.png"))
    print(f"Wrote audit PNGs to {audit_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
