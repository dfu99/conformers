#!/usr/bin/env python3
"""Resolve tail flips in existing head-anchored fitted coordinates.

The head-anchored flip resolution in fit_with_head_tracking.py only checks
head position. After the head is locked to the AFM tracked position, the
remaining degree of freedom is rotation around the vertical (z) axis
through the head, which can flip the tail 180° without moving the head.

This script adds a second-stage smoothing that:
  1. Identifies the head centroid (already aligned to AFM)
  2. Computes tail position relative to head for each frame
  3. For each frame, tries flipping 180° around z-axis through head
  4. Picks the orientation whose tail position is closest to the local
     temporal median (rolling window of neighboring frames)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import mdtraj as md


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fitted-dir", type=Path, required=True,
                   help="Directory from fit_with_head_tracking.py")
    p.add_argument("--output-dir", type=Path, default=None,
                   help="Output dir (defaults to fitted-dir)")
    p.add_argument("--head-residues-a", default="1-440")
    p.add_argument("--head-residues-b", default="1-350")
    p.add_argument("--window", type=int, default=11,
                   help="Temporal window (odd) for median reference")
    return p.parse_args()


def parse_residue_range(s: str) -> tuple[int, int]:
    lo, hi = s.split("-")
    return int(lo), int(hi)


def get_tail_atoms(topology: md.Topology, head_a_range, head_b_range,
                   atom_sel: str = "name CA") -> tuple[np.ndarray, np.ndarray]:
    """Return (head_atom_indices, tail_atom_indices) for CA atoms."""
    ca = topology.select(atom_sel)
    head_idx = []
    tail_idx = []
    for i in ca:
        atom = topology.atom(i)
        chain_id = atom.residue.chain.index  # 0=alpha, 1=beta
        resid = atom.residue.resSeq
        if chain_id == 0:
            lo, hi = head_a_range
        else:
            lo, hi = head_b_range
        if lo <= resid <= hi:
            head_idx.append(i)
        else:
            tail_idx.append(i)
    return np.array(head_idx), np.array(tail_idx)


def flip_around_head_z(coords: np.ndarray, head_idx: np.ndarray) -> np.ndarray:
    """Rotate coords 180° around z-axis through head centroid."""
    out = coords.copy()
    head_xy = out[head_idx, :2].mean(axis=0)
    # 180° z-rotation: (x,y,z) → (-x+2cx, -y+2cy, z)
    out[:, 0] = 2 * head_xy[0] - out[:, 0]
    out[:, 1] = 2 * head_xy[1] - out[:, 1]
    return out


def resolve_tail_flips(coords: np.ndarray, head_idx: np.ndarray,
                        tail_idx: np.ndarray, window: int = 11) -> tuple:
    """Pass over the trajectory and flip frames whose tail disagrees with
    the temporal-median tail direction.

    Strategy:
      1. Compute each frame's tail centroid minus head centroid (vector).
      2. For each frame, compute the median tail-vector direction in a
         local window of neighbors (both flipped and unflipped candidates).
      3. Pick the orientation whose tail direction is closer to that median.
    """
    n = len(coords)
    flipped = np.zeros(n, dtype=bool)
    out_coords = coords.copy()
    half = window // 2

    # Iteratively resolve: flip frames whose tail direction conflicts with
    # the median of their neighbors. Repeat until stable.
    # Lower hysteresis + more iterations to catch subtle flips.
    for iteration in range(20):
        any_change = False
        # Compute current tail vectors
        tail_vecs = np.zeros((n, 2))
        for i in range(n):
            head_xy = out_coords[i, head_idx, :2].mean(axis=0)
            tail_xy = out_coords[i, tail_idx, :2].mean(axis=0)
            tail_vecs[i] = tail_xy - head_xy

        # Normalize to unit vectors
        norms = np.linalg.norm(tail_vecs, axis=1, keepdims=True)
        norms[norms < 1e-9] = 1.0
        unit_vecs = tail_vecs / norms

        # Use circular median (via vector mean of unit vectors)
        for i in range(n):
            lo = max(0, i - half)
            hi = min(n, i + half + 1)
            # Neighbors excluding current frame
            neighbor_idx = [j for j in range(lo, hi) if j != i]
            if not neighbor_idx:
                continue
            neighbor_mean = unit_vecs[neighbor_idx].mean(axis=0)
            nm_norm = np.linalg.norm(neighbor_mean)
            if nm_norm < 1e-6:
                continue
            neighbor_mean /= nm_norm

            curr_dot = np.dot(unit_vecs[i], neighbor_mean)
            flipped_dot = -curr_dot

            # Lower hysteresis for more aggressive flipping
            if flipped_dot > curr_dot + 0.05:
                out_coords[i] = flip_around_head_z(out_coords[i], head_idx)
                flipped[i] = not flipped[i]
                any_change = True

        if not any_change:
            print(f"  Converged after {iteration + 1} iteration(s)")
            break
    else:
        print(f"  Did not fully converge in 20 iterations")

    n_flipped = int(flipped.sum())
    pct = 100 * n_flipped / max(1, n)
    print(f"  Tail flip resolution: {n_flipped}/{n} frames flipped ({pct:.1f}%)")
    return out_coords, flipped


def main():
    args = parse_args()
    out_dir = args.output_dir or args.fitted_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    coords = np.load(str(args.fitted_dir / "fitted_coords_smooth.npy"))
    topology = md.load(str(args.fitted_dir / "topology.pdb")).topology

    head_a = parse_residue_range(args.head_residues_a)
    head_b = parse_residue_range(args.head_residues_b)
    head_idx, tail_idx = get_tail_atoms(topology, head_a, head_b)
    print(f"Head atoms: {len(head_idx)}, Tail atoms: {len(tail_idx)}")
    print(f"Loaded {len(coords)} frames from {args.fitted_dir}")
    print(f"Resolving tail flips with window={args.window}...")

    # Note: coords here are full atom arrays, head_idx/tail_idx are over full atoms
    # but we got them from CA selection. Need to convert to be over all atoms
    # indexed the same as coords. The fitted_coords are already full-atom.
    # Actually the head/tail idx from CA selection is fine because coords has
    # the full atom set from topology.pdb.

    new_coords, flipped = resolve_tail_flips(coords, head_idx, tail_idx,
                                              window=args.window)

    # Save updated coords with tail-flip-resolved smoothing
    out_path = out_dir / "fitted_coords_smooth_tail.npy"
    np.save(str(out_path), new_coords)
    print(f"Saved: {out_path}")

    flipped_path = out_dir / "tail_flipped_mask.npy"
    np.save(str(flipped_path), flipped)
    print(f"Saved: {flipped_path}")


if __name__ == "__main__":
    main()
