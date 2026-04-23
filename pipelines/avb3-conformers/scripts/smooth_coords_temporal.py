#!/usr/bin/env python3
"""Temporal coordinate smoothing for HS-AFM fitted overlays.

Applies a rolling-median filter on the fitted coordinates directly, then
re-anchors each frame's head centroid back to the pre-temporal head
position (which was already locked to the tracked AFM head by
head-anchored SO(3) fitting).

This is a post-processing step on the output of fit_with_head_tracking.py.
It replaces `fitted_coords_smooth.npy` with a temporally smoothed version
that reduces frame-to-frame jitter by ~77% on our data without drifting
the overlay off the image.

Reasoning:
  1. Per-frame SO(3) fit over-rotates to match AFM noise, producing
     visible overlay twist-ups between adjacent frames.
  2. A rolling median on fitted_coords directly (per atom, per dim)
     removes these high-frequency artifacts. Median is preferred over
     mean because it preserves sharp transitions (BC→EC) while smoothing
     one-off jitter.
  3. The median smears head XY across its window (~1Å mean on our data)
     which is invisible on a 379-frame video but builds visible drift on
     a 1266-frame one. Re-anchor step eliminates it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md
from scipy.ndimage import median_filter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fitted-dir", type=Path, required=True,
                   help="Output dir from fit_with_head_tracking.py")
    p.add_argument("--window", type=int, default=7,
                   help="Rolling median window size (frames, odd). "
                        "Larger = smoother but blurs real transitions.")
    p.add_argument("--head-residues-a", default="1-440")
    p.add_argument("--head-residues-b", default="1-350")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite fitted_coords_smooth.npy in place "
                        "(preserves the pre-temporal version as "
                        "fitted_coords_smooth_pretemporal.npy)")
    return p.parse_args()


def parse_range(s):
    lo, hi = s.split("-")
    return int(lo), int(hi)


def main():
    args = parse_args()
    fd = args.fitted_dir

    pre_path = fd / "fitted_coords_smooth.npy"
    pre_backup = fd / "fitted_coords_smooth_pretemporal.npy"
    coords = np.load(str(pre_path)).astype(np.float32)
    print(f"Loaded {coords.shape[0]} frames from {pre_path}")

    if not pre_backup.exists():
        np.save(str(pre_backup), coords)
        print(f"Backed up original to {pre_backup}")

    # Median filter along time axis (axis 0), per atom, per dim
    print(f"Applying rolling-median (window={args.window})...")
    smoothed = median_filter(coords, size=(args.window, 1, 1), mode="nearest")

    # Re-anchor head to pre-temporal per-frame positions
    topo = md.load(str(fd / "topology.pdb")).topology
    lo_a, hi_a = parse_range(args.head_residues_a)
    lo_b, hi_b = parse_range(args.head_residues_b)
    head_a = topo.select(f"chainid 0 and resSeq {lo_a} to {hi_a}")
    head_b = topo.select(f"chainid 1 and resSeq {lo_b} to {hi_b}")
    head_idx = np.concatenate([head_a, head_b])
    print(f"Head atoms: {len(head_idx)}")

    target = coords[:, head_idx].mean(axis=1)[:, :2]  # pre-temporal head XY
    out = smoothed.copy()
    for i in range(len(out)):
        cur = out[i, head_idx].mean(axis=0)[:2]
        out[i, :, 0] += target[i, 0] - cur[0]
        out[i, :, 1] += target[i, 1] - cur[1]
        out[i, :, 2] -= out[i, :, 2].min()

    # Jitter diagnostics
    def jit(c):
        return np.linalg.norm(c[1:] - c[:-1], axis=-1).mean(axis=1) * 10
    print(f"Jitter: per-frame {jit(coords).mean():.2f}Å → "
          f"smoothed {jit(out).mean():.2f}Å "
          f"(-{100*(jit(coords).mean()-jit(out).mean())/jit(coords).mean():.1f}%)")

    if args.overwrite:
        np.save(str(pre_path), out)
        print(f"Overwrote {pre_path}")
    else:
        out_path = fd / "fitted_coords_median_reanchored.npy"
        np.save(str(out_path), out)
        print(f"Saved {out_path} (use --overwrite to replace fitted_coords_smooth.npy)")


if __name__ == "__main__":
    main()
