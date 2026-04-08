#!/usr/bin/env python3
"""Train AFMFold CNN on AVB3 integrin steering trajectories.

End-to-end pipeline:
1. Load PDB frames from domain steering output
2. Compute inter-domain CVs (centroid distances)
3. Generate simulated AFM images via afmfold
4. Train C8-equivariant CNN to predict CVs from AFM images
5. Save trained model checkpoint

Usage:
    python train_afmfold_cnn.py \
        --frames-dir data/runs/avb3/domain_steering/cv_distance_extend_frames/frames/ \
        --output-dir data/runs/avb3/afmfold_training/ \
        --epochs 50
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--frames-dir", type=Path, required=True,
                   help="Directory of PDB frames from steering (protein-only or solvated).")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Output directory for training data and model.")
    p.add_argument("--afmfold-root", type=Path,
                   default=Path("/home2/Documents/code/afmfold"),
                   help="Root of afmfold repository.")
    p.add_argument("--protein-name", default="avb3")
    # Image generation params
    p.add_argument("--image-size", type=int, default=35)
    p.add_argument("--resolution-nm", type=float, default=0.98)
    p.add_argument("--noise-nm", type=float, default=0.1)
    p.add_argument("--gen-epochs", type=int, default=500,
                   help="Image generation epochs (augmentation rounds).")
    p.add_argument("--dataset-size", type=int, default=100,
                   help="Samples per generation epoch.")
    # Training params
    p.add_argument("--train-epochs", type=int, default=50)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--device", default="cpu")
    p.add_argument("--skip-gen", action="store_true",
                   help="Skip image generation (use existing data).")
    return p.parse_args()


def load_protein_frames(frames_dir: Path) -> "mdtraj.Trajectory":
    """Load PDB frames and strip water/ions."""
    import mdtraj as md

    pdb_files = sorted(frames_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files in {frames_dir}")

    print(f"Loading {len(pdb_files)} PDB frames from {frames_dir}...")
    ref = md.load(str(pdb_files[0]))
    # Strip to protein only
    protein_idx = ref.topology.select("protein")
    if len(protein_idx) < ref.n_atoms:
        print(f"  Stripping to protein: {len(protein_idx)}/{ref.n_atoms} atoms")
        ref = ref.atom_slice(protein_idx)

    frames = [ref]
    skipped = 0
    for pf in pdb_files[1:]:
        try:
            t = md.load(str(pf))
            prot = t.topology.select("protein")
            if len(prot) < t.n_atoms:
                t = t.atom_slice(prot)
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
    print(f"  Loaded trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms "
          f"({skipped} skipped)")
    return traj


def compute_domain_cvs(traj: "mdtraj.Trajectory", protein_name: str,
                       afmfold_root: Path) -> np.ndarray:
    """Compute inter-domain centroid distances."""
    src = str(afmfold_root / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    from afmfold.domain import get_domain_pairs, compute_domain_distance

    pairs = get_domain_pairs(protein_name)
    print(f"Computing {len(pairs)} inter-domain distances...")
    dists = []
    for d1, d2 in pairs:
        d = compute_domain_distance(traj, d1, d2)
        dists.append(d)
    cv = np.concatenate(dists, axis=1) * 10.0  # nm → Å
    print(f"  CV shape: {cv.shape}, range: [{cv.min():.1f}, {cv.max():.1f}] Å")
    return cv


def generate_afm_images(traj, cv_array, args):
    """Generate simulated AFM images with augmentation."""
    src = str(args.afmfold_root / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    from afmfold.images import generate_images

    img_dir = args.output_dir / "training_data"
    img_dir.mkdir(parents=True, exist_ok=True)

    print(f"Generating AFM images: {args.gen_epochs} epochs × "
          f"{args.dataset_size} samples = {args.gen_epochs * args.dataset_size} total")

    generate_images(
        traj,
        args.resolution_nm,
        args.image_size,
        args.image_size,
        args.gen_epochs,
        args.dataset_size,
        distance=cv_array,
        batch_size=args.batch_size,
        min_z=0.0,
        noise_nm=args.noise_nm,
        max_tip_radius=12.0,
        min_tip_radius=6.0,
        max_tip_angle=30.0,
        min_tip_angle=10.0,
        ref_images=None,
        is_tqdm=True,
        match_histgram=False,
        save_dir=str(img_dir),
        device=args.device,
    )

    n_img = len(list(img_dir.glob("image_*.npy")))
    n_lbl = len(list(img_dir.glob("label_*.npy")))
    print(f"  Generated {n_img} image batches, {n_lbl} label batches")
    return img_dir


def load_training_data(data_dir: Path):
    """Load generated image/label .npy files into tensors."""
    import torch

    img_files = sorted(data_dir.glob("image_*.npy"))
    lbl_files = sorted(data_dir.glob("label_*.npy"))

    images = np.concatenate([np.load(str(f)) for f in img_files], axis=0)
    labels = np.concatenate([np.load(str(f)) for f in lbl_files], axis=0)

    print(f"Training data: {images.shape[0]} samples, "
          f"images={images.shape}, labels={labels.shape}")

    # Convert to tensors
    images_t = torch.from_numpy(images).float().unsqueeze(1)  # (N,1,H,W)
    labels_t = torch.from_numpy(labels).float()  # (N,D)

    return images_t, labels_t


def train_cnn(images, labels, args):
    """Train the C8-equivariant CNN."""
    import torch
    from torch.utils.data import TensorDataset, DataLoader, random_split

    src = str(args.afmfold_root / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    from afmfold.cnn import CNSteerableCNN

    device = torch.device(args.device)
    output_dim = labels.shape[1]

    # Create model
    model = CNSteerableCNN(
        N=8,
        image_shape=(1, args.image_size, args.image_size),
        channels=[3, 6, 6, 12, 12, 8],
        kernel_size=[7, 5, 5, 5, 5, 5],
        num_hidden_layers=3,
        hidden_dim=64,
        output_dim=output_dim,
    ).to(device)

    print(f"\nCNN: {sum(p.numel() for p in model.parameters())} parameters, "
          f"output_dim={output_dim}")

    # Split data
    dataset = TensorDataset(images, labels)
    n_val = max(1, int(0.2 * len(dataset)))
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(dataset, [n_train, n_val])

    train_loader = DataLoader(train_set, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_set, batch_size=args.batch_size)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    # Training loop
    model_dir = args.output_dir / "model"
    model_dir.mkdir(parents=True, exist_ok=True)
    best_val_loss = float("inf")
    train_losses, val_losses = [], []

    print(f"Training: {args.train_epochs} epochs, {n_train} train / {n_val} val samples")

    for epoch in range(1, args.train_epochs + 1):
        # Train
        model.train()
        epoch_loss = 0.0
        for imgs, lbls in train_loader:
            imgs, lbls = imgs.to(device), lbls.to(device)
            pred = model(imgs)
            loss = torch.mean(torch.norm(pred - lbls, dim=1))
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * imgs.size(0)
        train_loss = epoch_loss / n_train
        train_losses.append(train_loss)

        # Validate
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for imgs, lbls in val_loader:
                imgs, lbls = imgs.to(device), lbls.to(device)
                pred = model(imgs)
                loss = torch.mean(torch.norm(pred - lbls, dim=1))
                val_loss += loss.item() * imgs.size(0)
        val_loss = val_loss / n_val
        val_losses.append(val_loss)

        if epoch % 5 == 0 or epoch == 1:
            print(f"  Epoch {epoch:3d}: train={train_loss:.3f} Å, val={val_loss:.3f} Å")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            ckpt = {
                "epoch": epoch,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "val_loss": val_loss,
                "config": {
                    "N": 8,
                    "image_shape": (1, args.image_size, args.image_size),
                    "channels": [3, 6, 6, 12, 12, 8],
                    "kernel_size": [7, 5, 5, 5, 5, 5],
                    "num_hidden_layers": 3,
                    "hidden_dim": 64,
                    "output_dim": output_dim,
                    "protein_name": args.protein_name,
                },
            }
            torch.save(ckpt, model_dir / "best_model.pt")

    # Save final model and loss curves
    torch.save(ckpt, model_dir / "final_model.pt")
    np.save(str(model_dir / "train_losses.npy"), np.array(train_losses))
    np.save(str(model_dir / "val_losses.npy"), np.array(val_losses))

    print(f"\nBest val loss: {best_val_loss:.3f} Å (saved to {model_dir / 'best_model.pt'})")
    return model, train_losses, val_losses


def plot_training(train_losses, val_losses, output_dir):
    """Plot training curves."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(train_losses, label="Train", linewidth=1.5)
    ax.plot(val_losses, label="Validation", linewidth=1.5)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Mean Pairwise Distance Error (Å)")
    ax.set_title("AVB3 Integrin AFMFold CNN Training")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(output_dir / "training_curves.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved training curves to {output_dir / 'training_curves.png'}")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_gen:
        # Step 1: Load frames
        traj = load_protein_frames(args.frames_dir)

        # Step 2: Compute CVs
        cv_array = compute_domain_cvs(traj, args.protein_name, args.afmfold_root)
        np.save(str(args.output_dir / "cv_distances.npy"), cv_array)

        # Step 3: Generate AFM images
        img_dir = generate_afm_images(traj, cv_array, args)
    else:
        img_dir = args.output_dir / "training_data"
        print(f"Skipping generation, using existing data in {img_dir}")

    # Step 4: Load and train
    images, labels = load_training_data(img_dir)
    model, train_losses, val_losses = train_cnn(images, labels, args)

    # Step 5: Plot
    plot_training(train_losses, val_losses, args.output_dir)

    print("\n=== Pipeline complete ===")
    print(f"Model: {args.output_dir / 'model' / 'best_model.pt'}")
    print(f"Ready for inference with predict_from_afm_gif.py")


if __name__ == "__main__":
    main()
