#!/usr/bin/env python3
"""Stabilize per-frame conformer orientation before pseudo-AFM rendering.

The base fitted_coords are head-anchored in xy but free to roll about
the long axis between frames (causing the molecule to "flip to its
other side") and to undergo continuous yaw drift. Real HS-AFM shows:
- One face of the molecule rests flat on the substrate; the other
  face stays up.
- In-plane rotation (yaw) happens in stepwise jumps, not smooth drift.
- Head xy position drifts much less than this rendering implies.

Pipeline:
1. *PCA-flatten*: rotate each frame so the smallest principal axis
   (perpendicular to the molecular plane) points to +z. The molecule
   lies flat on the substrate.
2. *Side-lock*: pick a reference signed direction in frame 0 (e.g., a
   user-specified atom relative to the COM). For every subsequent
   frame, if the reference direction's sign flipped, rotate 180° about
   the long axis so the same face stays up.
3. *Sign-align long axis*: PCA eigenvectors are sign-ambiguous; align
   each frame's long axis to be acute relative to the previous frame.
4. *Stepwise yaw lock*: track the in-plane long-axis angle. Hysteresis:
   only update committed yaw if the running angle differs by more than
   threshold for >dwell frames, then cap step at max_step.
5. *Head xy-anchor*: smooth the head xy trajectory with a Gaussian
   kernel and re-anchor each frame to the smoothed reference.

Output: fitted_coords_stable.npy (same shape as input).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import mdtraj as md
from scipy.ndimage import gaussian_filter1d


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fitted-coords", type=Path, required=True)
    p.add_argument("--topology", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--head-residues-a", default="1-440",
                   help="α-chain residue range used as the head reference")
    p.add_argument("--head-residues-b", default="1-350",
                   help="β-chain residue range used as the head reference")
    p.add_argument("--side-lock-chain", default="A",
                   help="Use this chain's COM as the up-side reference")
    p.add_argument("--yaw-step-threshold", type=float, default=50.0,
                   help="Degrees the running mean must deviate before committing")
    p.add_argument("--yaw-step-cap", type=float, default=30.0,
                   help="Maximum committed yaw change per step (deg)")
    p.add_argument("--yaw-dwell-frames", type=int, default=20,
                   help="Min consecutive frames above threshold before committing")
    p.add_argument("--head-anchor-sigma", type=float, default=8.0,
                   help="Gaussian σ (frames) for head xy smoothing")
    return p.parse_args()


def parse_range(spec: str):
    lo, hi = (int(x) for x in spec.split("-"))
    return lo, hi


def rotation_matrix(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rodrigues rotation."""
    a = axis / (np.linalg.norm(axis) + 1e-12)
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    K = np.array(
        [[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]], dtype=np.float64
    )
    return np.eye(3) + s * K + (1 - c) * (K @ K)


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    coords = np.load(str(args.fitted_coords)).astype(np.float64)
    if coords.ndim != 3 or coords.shape[2] != 3:
        raise ValueError(f"fitted_coords shape {coords.shape}")
    T, N, _ = coords.shape
    print(f"Loaded {T} frames, {N} atoms")

    topo = md.load(str(args.topology)).topology
    ca_idx = np.array([a.index for a in topo.atoms if a.name == "CA"], dtype=int)
    ca_chain = np.array([topo.atom(int(i)).residue.chain.index for i in ca_idx])
    ca_resseq = np.array([topo.atom(int(i)).residue.resSeq for i in ca_idx])
    ca_chain_id = np.array([
        topo.atom(int(i)).residue.chain.chain_id or
        chr(ord('A') + topo.atom(int(i)).residue.chain.index)
        for i in ca_idx
    ])

    a_lo, a_hi = parse_range(args.head_residues_a)
    b_lo, b_hi = parse_range(args.head_residues_b)
    head_mask = (
        ((ca_chain == 0) & (ca_resseq >= a_lo) & (ca_resseq <= a_hi))
        | ((ca_chain == 1) & (ca_resseq >= b_lo) & (ca_resseq <= b_hi))
    )
    head_ca = ca_idx[head_mask]
    print(f"Head reference: {len(head_ca)} CAs")

    side_chain_mask = ca_chain_id == args.side_lock_chain
    side_chain_ca = ca_idx[side_chain_mask]
    print(f"Side-lock reference: {len(side_chain_ca)} CAs (chain {args.side_lock_chain})")

    out = np.zeros_like(coords)
    long_axis_unit_prev = None
    yaw_committed_deg = 0.0
    yaw_running_running = 0.0
    yaw_dwell_count = 0
    yaw_history_deg = np.zeros(T)
    side_sign_ref = None
    head_xy_trajectory = np.zeros((T, 2))

    for t in range(T):
        frame = coords[t].copy()
        com = frame.mean(axis=0)
        rel = frame - com

        # Step 1: PCA-flatten
        cov = rel.T @ rel
        eig_w, eig_v = np.linalg.eigh(cov)  # ascending eigenvalues
        normal = eig_v[:, 0]                # smallest variance = surface normal
        long_axis = eig_v[:, -1]            # largest variance = long axis
        # Build a rotation that takes (long, mid, normal) to (x, y, z).
        mid_axis = eig_v[:, 1]
        # Right-handed frame: ensure long × mid = normal (otherwise flip mid).
        if np.dot(np.cross(long_axis, mid_axis), normal) < 0:
            mid_axis = -mid_axis
        R_pca = np.column_stack([long_axis, mid_axis, normal]).T  # rows are basis vectors

        # Sign-align long axis with previous frame's long axis (acute angle).
        if long_axis_unit_prev is not None:
            new_long = R_pca[0]
            if np.dot(new_long, long_axis_unit_prev) < 0:
                # flip x and y rows together to keep right-handed-ness;
                # equivalently rotate 180° about the normal (z).
                R_pca[:2] = -R_pca[:2]
        long_axis_unit_prev = R_pca[0].copy()

        flat = (R_pca @ rel.T).T  # PCA-aligned coords centered at origin

        # Step 2: Side-lock — keep the same chain on the +z (up) side.
        if len(side_chain_ca) > 0:
            side_com = flat[side_chain_ca].mean(axis=0)
            side_sign_now = np.sign(side_com[2]) if side_com[2] != 0 else 1.0
            if side_sign_ref is None:
                side_sign_ref = side_sign_now
            elif side_sign_now != side_sign_ref:
                # Rotate 180° about x (long axis) so z flips while x stays.
                R_flip = rotation_matrix(np.array([1.0, 0.0, 0.0]), np.pi)
                flat = (R_flip @ flat.T).T

        # Step 3: head-up — also enforce head z > tail z so the head rests
        # on the substrate face (head is the readout).
        if len(head_ca) > 0:
            head_z = flat[head_ca, 2].mean()
            tail_mask = ~np.isin(np.arange(N), head_ca)
            tail_z = flat[tail_mask, 2].mean()
            if head_z < tail_z:
                R_flip = rotation_matrix(np.array([1.0, 0.0, 0.0]), np.pi)
                flat = (R_flip @ flat.T).T

        # Translate head centroid back to (0, 0) in xy; we'll re-anchor
        # to a smoothed trajectory in step 5.
        if len(head_ca) > 0:
            head_xy = flat[head_ca, :2].mean(axis=0)
            flat[:, 0] -= head_xy[0]
            flat[:, 1] -= head_xy[1]
            head_xy_trajectory[t] = head_xy

        # Step 4: stepwise yaw lock. The instantaneous yaw is the angle
        # of the long axis in xy. Use hysteresis: maintain a "committed"
        # yaw and only update it when the running angle has been outside
        # the threshold for `dwell` consecutive frames, then jump by at
        # most `cap` degrees toward it.
        long_xy = R_pca[0, :2]  # long-axis projection in original frame, before flips
        # After all the flips, the canonical long axis is +x in `flat`,
        # so the running yaw is best estimated by looking at where the
        # long axis points BEFORE flatten relative to the world frame.
        # That's equivalent to atan2 of the original long_axis_unit_prev's
        # xy projection.
        yaw_inst_rad = np.arctan2(long_axis_unit_prev[1], long_axis_unit_prev[0])
        yaw_inst_deg = np.rad2deg(yaw_inst_rad)
        # Bring yaw_inst into [-180, 180] continuous with previous committed.
        delta = yaw_inst_deg - yaw_committed_deg
        delta = ((delta + 180.0) % 360.0) - 180.0
        yaw_running_running = 0.7 * yaw_running_running + 0.3 * delta
        if abs(yaw_running_running) > args.yaw_step_threshold:
            yaw_dwell_count += 1
            if yaw_dwell_count >= args.yaw_dwell_frames:
                step = max(-args.yaw_step_cap, min(args.yaw_step_cap, yaw_running_running))
                yaw_committed_deg += step
                yaw_running_running = 0.0
                yaw_dwell_count = 0
        else:
            yaw_dwell_count = 0
        yaw_history_deg[t] = yaw_committed_deg

        # Apply committed yaw about z. (In `flat`, long axis is +x;
        # rotating by yaw_committed gives the locked in-plane orientation.)
        yaw_rad = np.deg2rad(yaw_committed_deg)
        cz, sz = np.cos(yaw_rad), np.sin(yaw_rad)
        Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
        flat = (Rz @ flat.T).T

        out[t] = flat

    # Step 5: head xy anchoring with Gaussian smoothing.
    if args.head_anchor_sigma > 0:
        # Smooth head_xy_trajectory in time, then translate each frame so
        # the head sits at the smoothed xy at that frame (relative to its
        # current head xy, which is (0, 0) after step 4).
        smooth_x = gaussian_filter1d(head_xy_trajectory[:, 0],
                                     sigma=args.head_anchor_sigma)
        smooth_y = gaussian_filter1d(head_xy_trajectory[:, 1],
                                     sigma=args.head_anchor_sigma)
        for t in range(T):
            out[t, :, 0] += smooth_x[t]
            out[t, :, 1] += smooth_y[t]

    # Lift the molecule so its lowest z is 0 (substrate plane).
    for t in range(T):
        out[t, :, 2] -= out[t, :, 2].min()

    np.save(str(args.output), out.astype(np.float32))
    print(f"Wrote {args.output}")
    yaw_path = args.output.with_name(args.output.stem + "_yaw.npy")
    np.save(str(yaw_path), yaw_history_deg.astype(np.float32))
    print(f"Wrote yaw history to {yaw_path}")
    print(f"Yaw stats: range [{yaw_history_deg.min():.1f}, {yaw_history_deg.max():.1f}] deg, "
          f"std {yaw_history_deg.std():.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
