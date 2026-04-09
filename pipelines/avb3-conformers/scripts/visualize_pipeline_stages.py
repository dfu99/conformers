#!/usr/bin/env python3
"""Visualize each stage of the AFMFold pipeline with artifacts.

Mirrors Section 5 of Protenix.ipynb:
  Stage 1: Load PDB frames (MD steering output)
  Stage 2: Compute inter-domain CVs
  Stage 3: Generate pseudo-AFM images
  Stage 4: Training summary

Generates publication-quality plots at each stage for PI review.
"""
import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ──
FRAMES_DIR = Path("data/runs/avb3/domain_steering/cv_distance_extend_frames/frames_protein")
TRAINING_DATA = Path("data/runs/avb3/afmfold_training/training_data")
MODEL_DIR = Path("data/runs/avb3/afmfold_training/model")
OUT_DIR = Path("data/runs/avb3/afmfold_pipeline")
AFMFOLD_ROOT = Path("/home2/Documents/code/afmfold")

sys.path.insert(0, str(AFMFOLD_ROOT / "src"))


def stage1_load_frames():
    """Stage 1: Load PDB frames and report statistics."""
    import mdtraj as md

    s1_dir = OUT_DIR / "stage1_frames"
    s1_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = sorted(FRAMES_DIR.glob("*.pdb"))
    print(f"[Stage 1] Found {len(pdb_files)} PDB frames in {FRAMES_DIR}")

    # Load all frames
    ref = md.load(str(pdb_files[0]))
    frames = [ref]
    skipped = 0
    for pf in pdb_files[1:]:
        try:
            t = md.load(str(pf))
            if t.n_atoms != ref.n_atoms:
                skipped += 1
                continue
            t.unitcell_lengths = None
            t.unitcell_angles = None
            frames.append(t)
        except Exception:
            skipped += 1

    frames[0].unitcell_lengths = None
    frames[0].unitcell_angles = None
    traj = md.join(frames)
    print(f"  Loaded: {traj.n_frames} frames, {traj.n_atoms} atoms ({skipped} skipped)")

    # Save combined trajectory
    traj[0].save_pdb(str(s1_dir / "reference_topology.pdb"))
    traj.save_dcd(str(s1_dir / "combined_trajectory.dcd"))
    print(f"  Saved: reference_topology.pdb, combined_trajectory.dcd")

    # Plot: RMSD over frames (shows conformational progression)
    rmsd = md.rmsd(traj, traj, frame=0) * 10.0  # nm -> Å
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(rmsd, linewidth=1.5, color="#2196F3")
    ax.set_xlabel("Frame Index")
    ax.set_ylabel("RMSD to Frame 0 (Å)")
    ax.set_title(f"Stage 1: Conformational Progression ({traj.n_frames} frames, {traj.n_atoms} atoms)")
    ax.grid(alpha=0.3)
    fig.savefig(s1_dir / "rmsd_progression.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: rmsd_progression.png")

    # Plot: Radius of gyration
    rg = md.compute_rg(traj) * 10.0  # nm -> Å
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(rg, linewidth=1.5, color="#FF9800")
    ax.set_xlabel("Frame Index")
    ax.set_ylabel("Radius of Gyration (Å)")
    ax.set_title(f"Stage 1: Rg Progression (bent → extended)")
    ax.grid(alpha=0.3)
    fig.savefig(s1_dir / "rg_progression.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: rg_progression.png")

    return traj


def stage2_compute_cvs(traj):
    """Stage 2: Compute inter-domain collective variables."""
    from afmfold.domain import get_domain_pairs, compute_domain_distance

    s2_dir = OUT_DIR / "stage2_cvs"
    s2_dir.mkdir(parents=True, exist_ok=True)

    domain_pairs = get_domain_pairs("avb3")
    print(f"\n[Stage 2] Computing {len(domain_pairs)} inter-domain distances...")

    cv_array = np.zeros((len(traj), len(domain_pairs)))
    for i, (d1, d2) in enumerate(domain_pairs):
        cv_array[:, i] = compute_domain_distance(traj, d1, d2).ravel()
    cv_array *= 10.0  # nm -> Å

    np.save(str(s2_dir / "domain_distances.npy"), cv_array)
    print(f"  CV shape: {cv_array.shape}, range: [{cv_array.min():.1f}, {cv_array.max():.1f}] Å")
    print(f"  Saved: domain_distances.npy")

    # Domain pair labels
    pair_labels = [
        "α-head/thigh ↔ β-tail",
        "β-head ↔ α-tail",
        "α-head ↔ β-head",
    ]

    # Plot 1: CVs over frames (matches notebook cell 44)
    fig, axes = plt.subplots(cv_array.shape[1], 1, figsize=(10, 2.5 * cv_array.shape[1]), sharex=True)
    colors = ["#2196F3", "#FF5722", "#4CAF50"]
    for i in range(cv_array.shape[1]):
        axes[i].plot(cv_array[:, i], alpha=0.8, linewidth=1.5, color=colors[i],
                     label=pair_labels[i])
        axes[i].set_ylabel("Distance (Å)")
        axes[i].legend(loc="upper right")
        axes[i].grid(alpha=0.3)
    axes[-1].set_xlabel("Frame Index")
    plt.suptitle("Stage 2: Inter-Domain CVs Over Frames", fontsize=13, y=0.98)
    fig.savefig(s2_dir / "cvs_over_frames.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: cvs_over_frames.png")

    # Plot 2: CV scatter (CV0 vs CV1, matches notebook)
    if cv_array.shape[1] >= 2:
        fig, ax = plt.subplots(figsize=(7, 6))
        sc = ax.scatter(cv_array[:, 0], cv_array[:, 1], c=np.arange(len(cv_array)),
                        cmap="viridis", s=15, alpha=0.7)
        ax.set_xlabel(f"CV 0: {pair_labels[0]} (Å)")
        ax.set_ylabel(f"CV 1: {pair_labels[1]} (Å)")
        ax.set_title("Stage 2: CV Space (color = frame index)")
        plt.colorbar(sc, label="Frame Index")
        ax.grid(alpha=0.3)
        fig.savefig(s2_dir / "cv_scatter.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: cv_scatter.png")

    # Plot 3: All CV pairs
    if cv_array.shape[1] == 3:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        pairs = [(0, 1), (0, 2), (1, 2)]
        for ax, (a, b) in zip(axes, pairs):
            sc = ax.scatter(cv_array[:, a], cv_array[:, b], c=np.arange(len(cv_array)),
                            cmap="viridis", s=10, alpha=0.6)
            ax.set_xlabel(f"{pair_labels[a]} (Å)")
            ax.set_ylabel(f"{pair_labels[b]} (Å)")
            ax.grid(alpha=0.3)
        plt.suptitle("Stage 2: All CV Pair Combinations", fontsize=13)
        fig.savefig(s2_dir / "cv_all_pairs.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: cv_all_pairs.png")

    # Summary stats
    with open(s2_dir / "cv_summary.txt", "w") as f:
        f.write(f"CV shape: {cv_array.shape}\n")
        for i, lbl in enumerate(pair_labels):
            f.write(f"CV{i} ({lbl}): min={cv_array[:, i].min():.1f}, "
                    f"max={cv_array[:, i].max():.1f}, "
                    f"mean={cv_array[:, i].mean():.1f}, "
                    f"std={cv_array[:, i].std():.1f} Å\n")
    print(f"  Saved: cv_summary.txt")

    return cv_array


def stage3_pseudo_afm(traj, cv_array):
    """Stage 3: Generate pseudo-AFM images and visualize grid."""
    from afmfold.images import generate_images

    s3_dir = OUT_DIR / "stage3_pseudo_afm"
    s3_dir.mkdir(parents=True, exist_ok=True)

    # Check if we already have training data
    existing_imgs = sorted(TRAINING_DATA.glob("image_*.npy"))
    existing_lbls = sorted(TRAINING_DATA.glob("label_*.npy"))

    if existing_imgs and existing_lbls:
        print(f"\n[Stage 3] Loading existing pseudo-AFM data ({len(existing_imgs)} batches)...")
        all_images = np.concatenate([np.load(str(f)) for f in existing_imgs], axis=0)
        all_labels = np.concatenate([np.load(str(f)) for f in existing_lbls], axis=0)
    else:
        print(f"\n[Stage 3] Generating pseudo-AFM images...")
        import torch
        img_save_dir = str(s3_dir / "generated")
        os.makedirs(img_save_dir, exist_ok=True)
        generate_images(
            traj, 0.98, 35, 35, 5, 66,
            distance=cv_array, batch_size=1, min_z=0.0, noise_nm=0.1,
            max_tip_radius=12.0, min_tip_radius=6.0,
            max_tip_angle=30.0, min_tip_angle=10.0,
            ref_images=None, is_tqdm=True, match_histgram=False,
            save_dir=img_save_dir, device="cpu",
        )
        img_files = sorted(glob.glob(os.path.join(img_save_dir, "image_*.npy")))
        lbl_files = sorted(glob.glob(os.path.join(img_save_dir, "label_*.npy")))
        all_images = np.concatenate([np.load(f) for f in img_files], axis=0)
        all_labels = np.concatenate([np.load(f) for f in lbl_files], axis=0)

    num_samples = len(all_images)
    print(f"  Images shape: {all_images.shape}")
    print(f"  Labels shape: {all_labels.shape}")
    print(f"  Total samples: {num_samples}")

    # Save consolidated arrays
    np.save(str(s3_dir / "all_images.npy"), all_images)
    np.save(str(s3_dir / "all_labels.npy"), all_labels)

    # Plot 1: Sample Grid of Pseudo-AFM Images (matches notebook cell 53)
    grid_size = 10
    n_grid = min(grid_size * grid_size, num_samples)
    indices = np.linspace(0, num_samples - 1, n_grid, dtype=int)
    cols = grid_size
    rows = (n_grid + cols - 1) // cols

    fig = plt.figure(figsize=(14, 14))
    for i, idx in enumerate(indices):
        plt.subplot(rows, cols, i + 1)
        plt.imshow(all_images[idx], cmap="copper")
        plt.axis("off")
    plt.suptitle(f"Stage 3: Sample Grid of Pseudo-AFM Images (Total: {num_samples})", y=0.92, fontsize=14)
    fig.savefig(s3_dir / "pseudo_afm_grid.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: pseudo_afm_grid.png")

    # Plot 2: CV distribution of pseudo-AFM labels (matches notebook)
    pair_labels = ["α-head/thigh ↔ β-tail", "β-head ↔ α-tail", "α-head ↔ β-head"]
    num_cvs = all_labels.shape[1]
    fig, axes = plt.subplots(num_cvs, 1, figsize=(10, 2.5 * num_cvs), sharex=True)
    if num_cvs == 1:
        axes = [axes]
    colors = ["#2196F3", "#FF5722", "#4CAF50"]
    for i in range(num_cvs):
        axes[i].plot(all_labels[:, i], alpha=0.5, linewidth=0.5, color=colors[i],
                     label=pair_labels[i] if i < len(pair_labels) else f"CV {i}")
        axes[i].set_ylabel("Distance (Å)")
        axes[i].legend(loc="upper right")
        axes[i].grid(alpha=0.3)
    axes[-1].set_xlabel("Sample Index")
    plt.suptitle(f"Stage 3: CV Labels Across {num_samples} Pseudo-AFM Samples", fontsize=13, y=0.98)
    fig.savefig(s3_dir / "pseudo_afm_cv_distribution.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: pseudo_afm_cv_distribution.png")

    # Plot 3: CV scatter for pseudo-AFM labels
    if num_cvs >= 2:
        fig, ax = plt.subplots(figsize=(7, 6))
        sc = ax.scatter(all_labels[:, 0], all_labels[:, 1], c=np.arange(num_samples),
                        cmap="viridis", s=5, alpha=0.5)
        ax.set_xlabel(f"{pair_labels[0]} (Å)")
        ax.set_ylabel(f"{pair_labels[1]} (Å)")
        ax.set_title("Stage 3: Pseudo-AFM Label Space")
        plt.colorbar(sc, label="Sample Index")
        ax.grid(alpha=0.3)
        fig.savefig(s3_dir / "pseudo_afm_label_scatter.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: pseudo_afm_label_scatter.png")

    # Plot 4: Histogram of pixel values
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(all_images.ravel(), bins=100, color="#9C27B0", alpha=0.7)
    ax.set_xlabel("Pixel Value (height, nm)")
    ax.set_ylabel("Count")
    ax.set_title("Stage 3: Pseudo-AFM Pixel Distribution")
    ax.grid(alpha=0.3)
    fig.savefig(s3_dir / "pixel_distribution.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: pixel_distribution.png")

    # Summary
    with open(s3_dir / "pseudo_afm_summary.txt", "w") as f:
        f.write(f"Total samples: {num_samples}\n")
        f.write(f"Image shape: {all_images.shape}\n")
        f.write(f"Label shape: {all_labels.shape}\n")
        f.write(f"Pixel range: [{all_images.min():.3f}, {all_images.max():.3f}]\n")
        for i in range(num_cvs):
            lbl = pair_labels[i] if i < len(pair_labels) else f"CV {i}"
            f.write(f"Label CV{i} ({lbl}): [{all_labels[:, i].min():.1f}, {all_labels[:, i].max():.1f}] Å\n")
    print(f"  Saved: pseudo_afm_summary.txt")

    return all_images, all_labels


def stage4_training():
    """Stage 4: Training summary and model artifacts."""
    s4_dir = OUT_DIR / "stage4_training"
    s4_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[Stage 4] Training artifacts...")

    # Copy training curves if exists
    src_curves = Path("data/runs/avb3/afmfold_training/training_curves.png")
    if src_curves.exists():
        import shutil
        shutil.copy2(src_curves, s4_dir / "training_curves.png")
        print(f"  Copied: training_curves.png")

    # Load and report loss data
    train_path = MODEL_DIR / "train_losses.npy"
    val_path = MODEL_DIR / "val_losses.npy"
    if train_path.exists() and val_path.exists():
        train_losses = np.load(str(train_path))
        val_losses = np.load(str(val_path))
        best_epoch = np.argmin(val_losses) + 1
        best_val = val_losses.min()

        with open(s4_dir / "training_summary.txt", "w") as f:
            f.write(f"Model: CNSteerableCNN (C8-equivariant)\n")
            f.write(f"Parameters: 2,076,524\n")
            f.write(f"Image size: 35x35\n")
            f.write(f"Output dim: 3 (inter-domain distances in Å)\n")
            f.write(f"Epochs: {len(train_losses)}\n")
            f.write(f"Best val loss: {best_val:.3f} Å (epoch {best_epoch})\n")
            f.write(f"Final train loss: {train_losses[-1]:.3f} Å\n")
            f.write(f"Final val loss: {val_losses[-1]:.3f} Å\n")
        print(f"  Best val loss: {best_val:.3f} Å at epoch {best_epoch}")
        print(f"  Saved: training_summary.txt")

    # Verify model exists
    model_path = MODEL_DIR / "best_model.pt"
    if model_path.exists():
        size_mb = model_path.stat().st_size / 1e6
        print(f"  Model: {model_path} ({size_mb:.1f} MB)")
    else:
        print(f"  WARNING: No model found at {model_path}")


def print_manifest():
    """Print consolidated asset manifest."""
    print("\n" + "=" * 60)
    print("CONSOLIDATED PIPELINE ASSETS")
    print("=" * 60)
    for stage_dir in sorted(OUT_DIR.iterdir()):
        if stage_dir.is_dir():
            print(f"\n{stage_dir.name}/")
            for f in sorted(stage_dir.iterdir()):
                if f.is_file():
                    size = f.stat().st_size
                    if size > 1e6:
                        print(f"  {f.name:40s} {size/1e6:.1f} MB")
                    else:
                        print(f"  {f.name:40s} {size/1e3:.1f} KB")


if __name__ == "__main__":
    traj = stage1_load_frames()
    cv_array = stage2_compute_cvs(traj)
    all_images, all_labels = stage3_pseudo_afm(traj, cv_array)
    stage4_training()
    print_manifest()
