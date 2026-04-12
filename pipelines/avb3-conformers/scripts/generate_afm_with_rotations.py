#!/usr/bin/env python3
"""Generate pseudo-AFM library with saved SO(3) rotation matrices.

Wraps the standard afmfold generate_images but intercepts and saves
the rotation matrices used for each epoch. These rotations can then
be applied in ChimeraX to render conformers at the exact same
3D orientation as their matched pseudo-AFM image.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frames-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--afmfold-root", type=Path, required=True)
    p.add_argument("--protein-name", default="avb3")
    p.add_argument("--width", type=int, default=35)
    p.add_argument("--height", type=int, default=35)
    p.add_argument("--resolution-nm", type=float, default=0.98)
    p.add_argument("--noise-nm", type=float, default=0.1)
    p.add_argument("--tip-radius-px", type=float, default=3.57,
                   help="Fixed tip radius in pixels (3.5nm at 0.98nm/px = 3.57)")
    p.add_argument("--tip-angle", type=float, default=20.0)
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--device", default="cuda")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Add afmfold to path
    afmfold_src = str(args.afmfold_root / "src")
    if afmfold_src not in sys.path:
        sys.path.insert(0, afmfold_src)

    # Import and load frames
    import mdtraj as md

    pdb_files = sorted(args.frames_dir.glob("*.pdb"))
    ref = md.load(str(pdb_files[0]))
    keep = ref.topology.select("protein")
    ref = ref.atom_slice(keep)
    frames = [ref]
    for pf in pdb_files[1:]:
        t = md.load(str(pf))
        t = t.atom_slice(t.topology.select("protein"))
        if t.n_atoms != ref.n_atoms:
            continue
        t.unitcell_lengths = None
        t.unitcell_angles = None
        frames.append(t)
    frames[0].unitcell_lengths = None
    frames[0].unitcell_angles = None
    traj = md.join(frames)
    print(f"Loaded {traj.n_frames} frames")

    # Compute CVs
    from afmfold.domain import get_domain_pairs, compute_domain_distance
    domain_pairs = get_domain_pairs(args.protein_name)
    distances = []
    for d1, d2 in domain_pairs:
        distances.append(compute_domain_distance(traj, d1, d2))
    cv_array = np.concatenate(distances, axis=1) * 10.0  # nm → Å
    np.save(str(args.output_dir / "domain_distances.npy"), cv_array)

    # Monkey-patch sample_uniform_so3 to save rotation matrices
    import afmfold.images as img_mod
    saved_rotations = []
    original_sample = img_mod.sample_uniform_so3

    def patched_sample(*a, **kw):
        rot = original_sample(*a, **kw)
        saved_rotations.append(rot.detach().cpu().numpy().copy())
        return rot

    img_mod.sample_uniform_so3 = patched_sample

    # Generate images (rotations are intercepted)
    from afmfold.images import generate_images
    print(f"Generating {args.epochs} epochs × {traj.n_frames} conformers "
          f"with tip={args.tip_radius_px:.1f}px...")

    generate_images(
        traj, args.resolution_nm, args.width, args.height,
        args.epochs, traj.n_frames,
        distance=cv_array, batch_size=1, min_z=0.0,
        noise_nm=args.noise_nm,
        max_tip_radius=args.tip_radius_px,
        min_tip_radius=args.tip_radius_px,
        max_tip_angle=args.tip_angle,
        min_tip_angle=args.tip_angle,
        save_dir=str(args.output_dir),
        device=args.device,
    )

    # Restore original function
    img_mod.sample_uniform_so3 = original_sample

    # Save rotation matrices: one 3x3 matrix per epoch
    rotations = np.array([r[0] for r in saved_rotations])  # (n_epochs, 3, 3)
    rot_path = args.output_dir / "rotations.npy"
    np.save(str(rot_path), rotations)
    print(f"Saved {len(rotations)} rotation matrices to {rot_path}")
    print(f"  Shape: {rotations.shape}")

    # Verify
    n_img = len(list(args.output_dir.glob("image_*.npy")))
    print(f"Generated {n_img} image epochs, {len(rotations)} rotation matrices")


if __name__ == "__main__":
    main()
