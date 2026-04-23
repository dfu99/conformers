#!/usr/bin/env python3
"""Sliding-window temporal smoothing for per-frame PDB fits.

The per-frame SO(3) fit is noisy — adjacent frames with almost-identical
AFM images can get assigned different conformers with different rotations,
producing jittery overlays. This script applies temporal regularization
by taking the mode (or median) assignment in a sliding window, weighted
by AFM image similarity.

Strategy:
  1. Load per-frame fitted coords + matched_indices + head positions.
  2. Compute AFM-frame similarity matrix (rolling).
  3. For each frame, find frames in window ±W with AFM similarity >= T.
  4. Use the median-CV conformer from those frames — i.e., if the AFM
     hasn't changed, don't change the conformer.
  5. For each frame, apply the median conformer's CA backbone with the
     per-frame rotation + head-tracked translation.

Output: fitted_coords_temporal.npy that smooths out jitter in stable regions.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import mdtraj as md
from PIL import Image
from scipy.ndimage import gaussian_filter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fitted-dir', type=Path, required=True)
    p.add_argument('--gif', type=Path, required=True)
    p.add_argument('--frames-dir', type=Path, required=True,
                   help='Library of conformer PDBs.')
    p.add_argument('--output', type=Path, default=None,
                   help='Output npy path (default: fitted-dir/fitted_coords_temporal.npy)')
    p.add_argument('--window', type=int, default=15,
                   help='Rolling window size (frames).')
    p.add_argument('--similarity-threshold', type=float, default=0.95,
                   help='AFM frame-to-frame correlation threshold. '
                        'Below this, frames are considered "changing" and '
                        'per-frame fit is kept. Above this, median in window used.')
    p.add_argument('--image-size', type=int, default=35)
    return p.parse_args()


def extract_gif_frames(gif_path, target_size):
    gif = Image.open(gif_path)
    frames = []
    try:
        while True:
            arr = np.array(gif.convert('L'), dtype=np.float32)
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

    # Load existing fit
    fitted_coords = np.load(str(args.fitted_dir / 'fitted_coords_smooth.npy'))
    matched_indices = np.load(str(args.fitted_dir / 'matched_indices.npy'))
    matched_cvs = np.load(str(args.fitted_dir / 'matched_cvs.npy'))
    fitted_rotations = np.load(str(args.fitted_dir / 'fitted_rotations_smooth.npy'))
    head_positions_px = np.load(str(args.fitted_dir / 'head_positions_px.npy'))
    import json
    with open(args.fitted_dir / 'fitting_metadata.json') as f:
        meta = json.load(f)
    skip_frames = meta.get('skip_frames', 0)

    n_fitted = fitted_coords.shape[0]
    print(f'Loaded {n_fitted} fitted frames (skip_frames={skip_frames})')

    # Load AFM frames and align to fitted frame count
    afm = extract_gif_frames(args.gif, args.image_size)
    print(f'Loaded {len(afm)} AFM frames')
    afm_fitted = afm[skip_frames:skip_frames + n_fitted]

    # Compute AFM self-similarity (frame i vs frame i+1)
    afm_flat = afm_fitted.reshape(n_fitted, -1)
    afm_flat = afm_flat - afm_flat.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(afm_flat, axis=1, keepdims=True)
    norms[norms < 1e-9] = 1.0
    afm_norm = afm_flat / norms

    # Frame-to-frame cosine similarity (diagonal + nearby bands)
    # For efficient windowed median we only need within-window similarities.
    W = args.window
    half = W // 2
    print(f'Running sliding-window smoothing (W={W}, thresh={args.similarity_threshold})')

    # For each frame, find neighbors with high AFM similarity and take their CV median
    new_matched_cvs = matched_cvs[:n_fitted].copy()
    new_matched_idx = matched_indices[:n_fitted].copy()
    stability = np.zeros(n_fitted)
    changed = np.zeros(n_fitted, dtype=bool)
    # AFM frames have first 30 skipped
    afm_fitted_matched = afm_norm[:n_fitted]
    assert len(afm_fitted_matched) == n_fitted, f"{len(afm_fitted_matched)} vs {n_fitted}"

    for i in range(n_fitted):
        lo = max(0, i - half)
        hi = min(n_fitted, i + half + 1)
        # cosine sim between frame i and window
        sims = afm_fitted_matched[lo:hi] @ afm_fitted_matched[i]
        mask = sims >= args.similarity_threshold
        stability[i] = mask.sum() / len(sims)
        if mask.sum() >= 3:
            # Stable region: take window median of the CVs
            win_idx = np.arange(lo, hi)[mask]
            new_matched_cvs[i] = np.median(matched_cvs[win_idx], axis=0)
            # Find library conformer with CV closest to the smoothed median
            library_cvs = matched_cvs[:n_fitted]  # proxy — real library CVs
            # Actually use the library CVs (from process_frames_to_afm)
            # For simplicity here: pick the window-mode conformer
            unique, counts = np.unique(matched_indices[win_idx], return_counts=True)
            mode_idx = unique[counts.argmax()]
            if mode_idx != matched_indices[i]:
                new_matched_idx[i] = mode_idx
                changed[i] = True

    n_changed = int(changed.sum())
    print(f'Changed assignments: {n_changed}/{n_fitted} ({100*n_changed/n_fitted:.1f}%)')
    print(f'Mean stability: {stability.mean():.3f} (frames where >= 3 nbrs are '
          f'cosine-sim >= {args.similarity_threshold})')

    # Now construct new fitted_coords using the temporally-smoothed indices
    # Load library conformer CAs once
    pdb_files = sorted(args.frames_dir.glob('*.pdb'))
    print(f'Loading {len(pdb_files)} library conformers for re-rendering...')
    library_xyz = []
    for pf in pdb_files:
        t = md.load(str(pf))
        library_xyz.append(t.xyz[0].astype(np.float32))  # nm, [N, 3]
    library_xyz = np.stack(library_xyz)  # [L, N, 3]
    print(f'Library shape: {library_xyz.shape}')

    # Compute library CVs to do CV-based lookup (matched_indices is into
    # the 25k-image pseudo-AFM library, not the conformer library).
    import sys as _sys
    _sys.path.insert(0, '/home2/Documents/code/afmfold/src')
    from afmfold.domain import compute_domain_distance
    from afmfold.domain import AVB3_DOMAINS
    # Compute library conformer CVs
    lib_cvs = np.zeros((len(pdb_files), 3), dtype=np.float32)
    for i, pf in enumerate(pdb_files):
        t = md.load(str(pf))
        d0 = compute_domain_distance(t, AVB3_DOMAINS['alpha_head_thigh'], AVB3_DOMAINS['alpha_calf'])
        d1 = compute_domain_distance(t, AVB3_DOMAINS['beta_head_hybrid_egf1'], AVB3_DOMAINS['beta_tail_egf2_3_4_btail'])
        d2 = compute_domain_distance(t, AVB3_DOMAINS['alpha_head_thigh'], AVB3_DOMAINS['beta_head_hybrid_egf1'])
        lib_cvs[i] = [d0[0,0]*10, d1[0,0]*10, d2[0,0]*10]  # nm→Å
        if (i+1) % 100 == 0:
            print(f'  computed CVs for {i+1}/{len(pdb_files)}')
    print(f'Library CV ranges: CV0=[{lib_cvs[:,0].min():.1f},{lib_cvs[:,0].max():.1f}], '
          f'CV1=[{lib_cvs[:,1].min():.1f},{lib_cvs[:,1].max():.1f}]')

    # Topology for head atom indices
    topo = md.load(str(args.fitted_dir / 'topology.pdb')).topology
    head_a = topo.select(f"chainid 0 and resSeq 1 to 440")
    head_b = topo.select(f"chainid 1 and resSeq 1 to 350")
    head_idx = np.concatenate([head_a, head_b])

    # Build new coords: take the mode conformer's atoms, apply the per-frame
    # rotation, translate head centroid to the tracked AFM head position.
    # Use the same approach as fit_with_head_tracking.resolve_flips_head_anchored.
    new_fitted = np.zeros_like(fitted_coords)
    resolution_nm = 0.98  # from fitting_metadata
    # head_positions_px is in pixels; convert to nm using resolution_nm
    # Same coordinate system as fitted_coords (which is in nm centered around head)
    for i in range(n_fitted):
        # Use CV-based lookup: find library conformer with CV closest to smoothed CV
        cv_target = new_matched_cvs[i]
        cv_dists = np.linalg.norm(lib_cvs - cv_target, axis=1)
        idx = int(cv_dists.argmin())
        conf = library_xyz[idx].copy()  # [N, 3] in nm
        # Apply per-frame rotation
        R = fitted_rotations[i]  # [3, 3]
        # Center before rotation
        conf = conf - conf.mean(axis=0)
        conf = conf @ R.T
        # Place head centroid at tracked AFM head position
        head_com = conf[head_idx].mean(axis=0)
        target_xy = head_positions_px[i, :2] * resolution_nm  # convert px → nm
        conf[:, 0] += target_xy[0] - head_com[0]
        conf[:, 1] += target_xy[1] - head_com[1]
        conf[:, 2] -= conf[:, 2].min()
        new_fitted[i] = conf

    out = args.output or (args.fitted_dir / 'fitted_coords_temporal.npy')
    np.save(str(out), new_fitted)
    print(f'Saved: {out}')

    # Also save the diagnostic info
    np.save(str(args.fitted_dir / 'temporal_stability.npy'), stability)
    np.save(str(args.fitted_dir / 'temporal_matched_idx.npy'), new_matched_idx)


if __name__ == '__main__':
    main()
