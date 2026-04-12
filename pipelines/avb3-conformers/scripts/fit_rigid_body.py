#!/usr/bin/env python3
"""Rigid-body fitting: find optimal SO(3) orientation per HS-AFM frame.

Uses AFMFold's pseudo-AFM generation + correlation maximization to find the
3D rotation that best matches each conformer to its real AFM frame.

This replaces the sparse 50-rotation approach with exhaustive SO(3) search
(2048+ orientations per frame), eliminating head-tail orientation ambiguity.

Outputs:
  - fitted_rotations.npy: (N_frames, 3, 3) best rotation per frame
  - fitted_coords.npy: (N_frames, N_atoms, 3) fitted coordinates
  - fitted_correlations.npy: (N_frames,) correlation at best orientation
  - fitted_pseudo_afm.npy: (N_frames, H, W) pseudo-AFM at best orientation
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import mdtraj as md
from PIL import Image


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gif", type=Path, required=True,
                   help="Input HS-AFM GIF")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--frames-dir", type=Path, required=True,
                   help="Directory of protein-only PDB conformer frames")
    p.add_argument("--matched-indices", type=Path, required=True,
                   help="matched_indices.npy from correlation matching")
    p.add_argument("--afmfold-root", type=Path, required=True,
                   help="Root of the afmfold repository")
    p.add_argument("--n-conformers", type=int, default=264)
    p.add_argument("--width", type=int, default=35)
    p.add_argument("--height", type=int, default=35)
    p.add_argument("--resolution-nm", type=float, default=0.98)
    p.add_argument("--tip-radius", type=float, default=3.57,
                   help="Tip radius in pixels (3.5nm at 0.98nm/px = 3.57)")
    p.add_argument("--tip-angle", type=float, default=20.0)
    p.add_argument("--n-rotations", type=int, default=2048,
                   help="Total SO(3) rotation samples per frame")
    p.add_argument("--rot-batch", type=int, default=32,
                   help="Rotations per GPU step")
    p.add_argument("--frame-batch", type=int, default=20,
                   help="Frames to fit simultaneously")
    p.add_argument("--device", default="cuda")
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


def resize_to_grid(frames: list[np.ndarray], h: int, w: int) -> np.ndarray:
    """Crop to square, resize to (h, w)."""
    from skimage.transform import resize as sk_resize
    out = np.zeros((len(frames), h, w), dtype=np.float32)
    for i, f in enumerate(frames):
        fh, fw = f.shape
        if fh != fw:
            s = min(fh, fw)
            y0, x0 = (fh - s) // 2, (fw - s) // 2
            f = f[y0:y0+s, x0:x0+s]
        out[i] = sk_resize(f, (h, w), anti_aliasing=True, preserve_range=True)
    return out


def fit_orientations(
    target_images: np.ndarray,
    xyz: "torch.Tensor",
    center: "torch.Tensor",
    xedges: "torch.Tensor",
    yedges: "torch.Tensor",
    tip: "torch.Tensor",
    n_rotations: int,
    rot_batch: int,
    height: int,
    width: int,
    device: "torch.device",
):
    """Find optimal SO(3) rotation for each frame via exhaustive search.

    Args:
        target_images: (N, H, W) real AFM frames (resized to grid)
        xyz: (N, natoms, 3) tensor — conformer coordinates per frame
        center: (N, 3) tensor — target image COM in nm
        xedges, yedges: grid edges for pseudo-AFM generation
        tip: tip shape tensor
        n_rotations: total rotation samples
        rot_batch: rotations per GPU step
        height, width: grid dimensions
        device: torch device

    Returns:
        best_rot, best_cc, best_coords, best_pseudo — numpy arrays
    """
    import torch
    from afmfold.images import sample_uniform_so3, apply_rotations, generate_landscape, idilation
    from afmfold.rigid_body_fitting import compute_correlation_coefficient

    N = xyz.shape[0]
    target = torch.from_numpy(target_images).to(device)

    # Target image statistics for normalization
    tmed = torch.median(target)
    tmax = torch.max(target)

    # Track best per frame
    best_cc = torch.full((N,), -1.0, device=device)
    best_rot = torch.zeros((N, 3, 3), device=device)
    best_coords = torch.zeros_like(xyz)
    best_pseudo = torch.zeros((N, height, width), device=device)

    n_steps = (n_rotations + rot_batch - 1) // rot_batch

    for step in range(n_steps):
        # Sample random SO(3) rotations
        rots = sample_uniform_so3(rot_batch, device=device)  # (rot_batch, 3, 3)

        # Apply to all frames: (N, rot_batch, natoms, 3)
        rotated = apply_rotations(xyz, rots)

        # Center on target image COM
        com = rotated.mean(dim=-2, keepdim=True)  # (N, rot_batch, 1, 3)
        centered = rotated - com + center.unsqueeze(1).unsqueeze(2)

        # Z-shift: minimum z to 0
        z_min = centered[..., 2].min(dim=-1, keepdim=True).values
        centered = centered.clone()
        centered[..., 2] = centered[..., 2] - z_min

        # Flatten for landscape generation
        flat = centered.reshape(-1, centered.shape[-2], 3)

        # Generate pseudo-AFM images
        landscape, _, _ = generate_landscape(flat, xedges, yedges)
        landscape = landscape.reshape(-1, height, width)
        pseudo = idilation(landscape, tip)

        # Normalize pseudo images to match target statistics
        pmin = pseudo.reshape(len(pseudo), -1).min(dim=1).values.view(-1, 1, 1)
        pmax = pseudo.reshape(len(pseudo), -1).max(dim=1).values.view(-1, 1, 1)
        pseudo_norm = (pseudo - pmin) / (pmax - pmin + 1e-6) * (tmax - tmed) + tmed

        # Reshape: (N, rot_batch, H, W)
        pseudo_norm = pseudo_norm.reshape(N, rot_batch, height, width)

        # Correlation: target (N, 1, H, W) vs pseudo (N, rot_batch, H, W)
        # → (N, 1, rot_batch) → squeeze → (N, rot_batch)
        cc = compute_correlation_coefficient(
            target.unsqueeze(1), pseudo_norm
        ).squeeze(1)

        # Find best rotation per frame in this batch
        batch_best_cc, batch_best_idx = cc.max(dim=-1)

        # Update global best where improved
        improved = batch_best_cc > best_cc
        if improved.any():
            imp_idx = torch.where(improved)[0]
            best_cc[imp_idx] = batch_best_cc[imp_idx]
            for f in imp_idx:
                bi = batch_best_idx[f]
                best_rot[f] = rots[bi]
                best_coords[f] = centered[f, bi]
                best_pseudo[f] = pseudo_norm[f, bi]

        if (step + 1) % 16 == 0 or step == n_steps - 1:
            print(f"    step {step+1}/{n_steps}: "
                  f"mean_cc={best_cc.mean():.4f} "
                  f"min={best_cc.min():.4f} max={best_cc.max():.4f}")

    return (best_rot.cpu().numpy(), best_cc.cpu().numpy(),
            best_coords.cpu().numpy(), best_pseudo.cpu().numpy())


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Add afmfold to path
    afmfold_src = str(args.afmfold_root / "src")
    if afmfold_src not in sys.path:
        sys.path.insert(0, afmfold_src)

    import torch
    from afmfold.images import generate_tip_shape
    from afmfold.rigid_body_fitting import threshold_and_mask
    from scipy.ndimage import center_of_mass

    device = torch.device(args.device)

    # Load matched indices
    matched_indices = np.load(str(args.matched_indices))
    n_frames = len(matched_indices)
    conformer_indices = (matched_indices % args.n_conformers).astype(int)
    print(f"Loaded {n_frames} matched frame indices "
          f"({len(np.unique(conformer_indices))} unique conformers)")

    # Extract and resize GIF frames to pseudo-AFM grid
    print("Extracting GIF frames...")
    real_frames = extract_gif_frames(args.gif)
    afm_grid = resize_to_grid(real_frames, args.height, args.width)
    print(f"  Resized {len(afm_grid)} frames to {args.height}x{args.width}")

    # Load conformer PDB coordinates
    print("Loading conformer PDBs...")
    pdb_files = sorted(args.frames_dir.glob("*.pdb"))
    ref = md.load(str(pdb_files[0]))
    keep = ref.topology.select("protein")
    ref = ref.atom_slice(keep)
    n_atoms = ref.n_atoms
    topology = ref.topology

    conformer_xyz = {}
    for pf in pdb_files:
        idx = int(pf.stem.split("_")[-1])
        t = md.load(str(pf))
        t = t.atom_slice(t.topology.select("protein"))
        if t.n_atoms == n_atoms:
            conformer_xyz[idx] = t.xyz[0]
    print(f"  Loaded {len(conformer_xyz)} conformers ({n_atoms} atoms each)")

    # Build per-frame xyz array
    all_xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    for i in range(n_frames):
        cidx = conformer_indices[i]
        if cidx in conformer_xyz:
            all_xyz[i] = conformer_xyz[cidx]
        else:
            nearest = min(conformer_xyz.keys(), key=lambda k: abs(k - cidx))
            all_xyz[i] = conformer_xyz[nearest]

    # Setup pseudo-AFM grid and tip
    xedges = args.resolution_nm * torch.arange(
        -int(args.width / 2) - 1, args.width - int(args.width / 2),
        device=device)
    yedges = args.resolution_nm * torch.arange(
        -int(args.height / 2) - 1, args.height - int(args.height / 2),
        device=device)
    tip = generate_tip_shape(args.tip_radius, args.tip_angle, device=device)

    # Compute target image COMs for centering
    cleaned = threshold_and_mask(afm_grid)
    coms = np.zeros((n_frames, 3), dtype=np.float32)
    for i in range(n_frames):
        cy, cx = center_of_mass(cleaned[i])
        coms[i] = [cx * args.resolution_nm, cy * args.resolution_nm, 0.0]
    center_all = torch.from_numpy(coms).to(device)

    # Allocate output arrays
    fitted_rotations = np.zeros((n_frames, 3, 3), dtype=np.float32)
    fitted_coords = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    fitted_correlations = np.zeros(n_frames, dtype=np.float32)
    fitted_pseudo = np.zeros((n_frames, args.height, args.width), dtype=np.float32)

    # Run fitting in batches
    print(f"\nRigid-body fitting: {args.n_rotations} rotations/frame, "
          f"rot_batch={args.rot_batch}, frame_batch={args.frame_batch}")

    for batch_start in range(0, n_frames, args.frame_batch):
        batch_end = min(batch_start + args.frame_batch, n_frames)
        n_batch = batch_end - batch_start

        print(f"\nBatch [{batch_start}:{batch_end}] "
              f"({batch_start // args.frame_batch + 1}/"
              f"{(n_frames + args.frame_batch - 1) // args.frame_batch})")

        batch_images = afm_grid[batch_start:batch_end]
        batch_xyz = torch.from_numpy(all_xyz[batch_start:batch_end]).to(device)
        batch_center = center_all[batch_start:batch_end]

        rot, cc, coords, pseudo = fit_orientations(
            batch_images, batch_xyz, batch_center,
            xedges, yedges, tip,
            args.n_rotations, args.rot_batch,
            args.height, args.width, device,
        )

        fitted_rotations[batch_start:batch_end] = rot
        fitted_coords[batch_start:batch_end] = coords
        fitted_correlations[batch_start:batch_end] = cc
        fitted_pseudo[batch_start:batch_end] = pseudo

        print(f"  → corr: mean={cc.mean():.4f} min={cc.min():.4f} max={cc.max():.4f}")

    # Save results
    np.save(str(args.output_dir / "fitted_rotations.npy"), fitted_rotations)
    np.save(str(args.output_dir / "fitted_coords.npy"), fitted_coords)
    np.save(str(args.output_dir / "fitted_correlations.npy"), fitted_correlations)
    np.save(str(args.output_dir / "fitted_pseudo_afm.npy"), fitted_pseudo)

    # Save topology PDB for rendering
    ref.save_pdb(str(args.output_dir / "topology.pdb"))

    print(f"\n{'='*60}")
    print(f"RIGID-BODY FITTING COMPLETE")
    print(f"  Frames: {n_frames}")
    print(f"  Rotations sampled: {args.n_rotations}")
    print(f"  Correlation: mean={fitted_correlations.mean():.4f} "
          f"min={fitted_correlations.min():.4f} max={fitted_correlations.max():.4f}")
    print(f"  Output: {args.output_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
