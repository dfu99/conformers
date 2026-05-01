#!/usr/bin/env python3
"""Process collected PDB frames into pseudo-AFM images via afmfold.

Steps:
1. Load all PDB frames as an MDTraj trajectory (requires identical topology).
2. Compute inter-domain collective variables (CVs) for the given protein.
3. Save combined trajectory (PDB + DCD) and CV distances.
4. Call afmfold's generate_images to produce pseudo-AFM images + labels.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import mdtraj as md
import numpy as np


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frames-dir", type=Path, required=True,
                   help="Directory of PDB frames (all must share topology).")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Directory for output images/labels.")
    p.add_argument("--afmfold-root", type=Path, required=True,
                   help="Root of the afmfold repository.")
    p.add_argument("--protein-name", default="avb3",
                   help="Protein name for domain CV lookup (default: avb3).")
    p.add_argument("--width", type=int, default=35)
    p.add_argument("--height", type=int, default=35)
    p.add_argument("--resolution-nm", type=float, default=0.98)
    p.add_argument("--noise-nm", type=float, default=0.1)
    p.add_argument("--min-tip-radius", type=float, default=6.0)
    p.add_argument("--max-tip-radius", type=float, default=12.0)
    p.add_argument("--min-tip-angle", type=float, default=10.0)
    p.add_argument("--max-tip-angle", type=float, default=30.0)
    p.add_argument("--epochs", type=int, default=1)
    p.add_argument("--dataset-size", type=int, default=1000)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--device", default="cpu")
    p.add_argument("--strip-water", action="store_true", default=True,
                   help="Strip water/ions before processing (default: True).")
    return p.parse_args()


def load_frames(frames_dir: Path, strip_water: bool = True) -> md.Trajectory:
    """Load all PDBs in a directory as a single trajectory."""
    pdb_files = sorted(frames_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in {frames_dir}")

    print(f"Loading {len(pdb_files)} PDB frames...")

    # Load first frame as reference topology
    ref = md.load(str(pdb_files[0]))
    if strip_water:
        keep = ref.topology.select("protein")
        ref = ref.atom_slice(keep)

    frames = [ref]
    skipped = 0
    for pf in pdb_files[1:]:
        try:
            t = md.load(str(pf))
            if strip_water:
                keep = t.topology.select("protein")
                t = t.atom_slice(keep)
            # Verify atom count matches reference
            if t.n_atoms != ref.n_atoms:
                print(f"  WARNING: {pf.name} has {t.n_atoms} atoms "
                      f"(expected {ref.n_atoms}), skipping.")
                skipped += 1
                continue
            # Clear unitcell to avoid join errors from mixed CRYST1 records
            t.unitcell_lengths = None
            t.unitcell_angles = None
            frames.append(t)
        except Exception as e:
            print(f"  WARNING: Failed to load {pf.name}: {e}")
            skipped += 1

    # Clear unitcell on reference too
    frames[0].unitcell_lengths = None
    frames[0].unitcell_angles = None
    traj = md.join(frames)
    print(f"  Loaded {traj.n_frames} frames ({skipped} skipped)")
    return traj


def compute_cvs(traj: md.Trajectory, protein_name: str,
                afmfold_root: Path) -> np.ndarray:
    """Compute inter-domain distance CVs using afmfold domain definitions."""
    # Add afmfold to path so we can import it
    src_path = str(afmfold_root / "src")
    if src_path not in sys.path:
        sys.path.insert(0, src_path)

    from afmfold.domain import get_domain_pairs, compute_domain_distance

    domain_pairs = get_domain_pairs(protein_name)
    print(f"Computing {len(domain_pairs)} inter-domain distances...")

    distances = []
    for d1, d2 in domain_pairs:
        dist = compute_domain_distance(traj, d1, d2)  # shape: [B, 1], nm
        distances.append(dist)

    # Stack and convert nm → Angstroms
    cv_array = np.concatenate(distances, axis=1) * 10.0  # [B, num_pairs]
    print(f"  CV shape: {cv_array.shape} (frames x domain pairs)")
    return cv_array


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load frames
    traj = load_frames(args.frames_dir, strip_water=args.strip_water)

    # 2. Save combined trajectory
    traj_pdb = args.output_dir / "combined_trajectory.pdb"
    traj_dcd = args.output_dir / "combined_trajectory.dcd"
    traj[0].save_pdb(str(traj_pdb))  # topology reference
    traj.save_dcd(str(traj_dcd))
    print(f"Saved trajectory: {traj_pdb}, {traj_dcd}")

    # 3. Compute CVs
    cv_array = compute_cvs(traj, args.protein_name, args.afmfold_root)
    distance_path = args.output_dir / "domain_distances.npy"
    np.save(str(distance_path), cv_array)
    print(f"Saved CVs: {distance_path}")

    # 4. Generate pseudo-AFM images
    afmfold_src = str(args.afmfold_root / "src")
    if afmfold_src not in sys.path:
        sys.path.insert(0, afmfold_src)

    from afmfold.images import generate_images

    print(f"Generating pseudo-AFM images ({args.dataset_size} samples, "
          f"{args.epochs} epochs)...")
    generate_images(
        traj,
        args.resolution_nm,
        args.width,
        args.height,
        args.epochs,
        args.dataset_size,
        distance=cv_array,
        batch_size=args.batch_size,
        min_z=0.0,
        noise_nm=args.noise_nm,
        max_tip_radius=args.max_tip_radius,
        min_tip_radius=args.min_tip_radius,
        max_tip_angle=args.max_tip_angle,
        min_tip_angle=args.min_tip_angle,
        ref_images=None,
        is_tqdm=True,
        match_histgram=False,
        save_dir=str(args.output_dir),
        device=args.device,
    )

    # Summary
    img_files = sorted(args.output_dir.glob("image_*.npy"))
    label_files = sorted(args.output_dir.glob("label_*.npy"))
    print(f"\nDone. Generated {len(img_files)} image batches, "
          f"{len(label_files)} label batches in {args.output_dir}/")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
