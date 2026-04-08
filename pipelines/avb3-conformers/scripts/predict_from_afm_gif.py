#!/usr/bin/env python3
"""Predict integrin conformational states from HS-AFM GIF.

Takes an experimental HS-AFM GIF, extracts frames, preprocesses them to
match the training image format, runs the trained CNN to predict inter-domain
CVs per frame, and matches each to the nearest steering PDB conformer.

Usage:
    python predict_from_afm_gif.py \
        --gif input.gif \
        --model data/runs/avb3/afmfold_training/model/best_model.pt \
        --frames-dir data/runs/avb3/domain_steering/cv_distance_extend_frames/frames/ \
        --output-dir results/afm_predictions/
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gif", type=Path, required=True,
                   help="Input HS-AFM GIF file.")
    p.add_argument("--model", type=Path, required=True,
                   help="Trained CNN checkpoint (.pt).")
    p.add_argument("--frames-dir", type=Path, default=None,
                   help="Steering PDB frames for nearest-neighbor matching.")
    p.add_argument("--cv-distances", type=Path, default=None,
                   help="Pre-computed CV distances (.npy) for the steering frames.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--image-size", type=int, default=35,
                   help="CNN input image size (must match training).")
    p.add_argument("--device", default="cpu")
    p.add_argument("--afmfold-root", type=Path,
                   default=Path("/home2/Documents/code/afmfold"))
    return p.parse_args()


def extract_gif_frames(gif_path: Path) -> np.ndarray:
    """Extract frames from GIF as grayscale float arrays."""
    from PIL import Image

    gif = Image.open(gif_path)
    frames = []
    try:
        while True:
            frame = gif.convert("L")  # grayscale
            frames.append(np.array(frame, dtype=np.float32))
            gif.seek(gif.tell() + 1)
    except EOFError:
        pass

    print(f"Extracted {len(frames)} frames from {gif_path.name}")
    print(f"  Frame size: {frames[0].shape}")
    return frames


def preprocess_frame(frame: np.ndarray, target_size: int) -> np.ndarray:
    """Preprocess a single AFM frame for CNN input.

    Steps:
    1. Crop to square (center crop)
    2. Resize to target_size × target_size
    3. Normalize to [0, 1] range (height map normalization)
    """
    from PIL import Image

    h, w = frame.shape
    # Center crop to square
    if h != w:
        s = min(h, w)
        y0 = (h - s) // 2
        x0 = (w - s) // 2
        frame = frame[y0:y0+s, x0:x0+s]

    # Resize to target
    img = Image.fromarray(frame)
    img = img.resize((target_size, target_size), Image.BILINEAR)
    arr = np.array(img, dtype=np.float32)

    # Normalize: shift min to 0, scale max to ~reasonable AFM height range
    arr = arr - arr.min()
    if arr.max() > 0:
        arr = arr / arr.max()

    return arr


def load_cnn(model_path: Path, afmfold_root: Path, device: str):
    """Load trained CNN from checkpoint."""
    import torch
    src = str(afmfold_root / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    from afmfold.cnn import CNSteerableCNN

    ckpt = torch.load(str(model_path), map_location=device, weights_only=False)
    config = ckpt["config"]

    model = CNSteerableCNN(
        N=config["N"],
        image_shape=tuple(config["image_shape"]),
        channels=config["channels"],
        kernel_size=config["kernel_size"],
        num_hidden_layers=config["num_hidden_layers"],
        hidden_dim=config["hidden_dim"],
        output_dim=config["output_dim"],
    )
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    model.to(device)

    print(f"Loaded CNN: {config['output_dim']} CV dims, "
          f"val_loss={ckpt.get('val_loss', '?'):.3f} Å")
    return model, config


def predict_cvs(model, frames: list[np.ndarray], image_size: int,
                device: str) -> np.ndarray:
    """Run CNN inference on preprocessed frames."""
    import torch

    processed = np.stack([preprocess_frame(f, image_size) for f in frames])
    # (N, H, W) → (N, 1, H, W)
    tensor = torch.from_numpy(processed).float().unsqueeze(1).to(device)

    with torch.no_grad():
        preds = model(tensor).cpu().numpy()

    print(f"Predicted CVs for {len(frames)} frames, shape={preds.shape}")
    return preds


def match_to_conformers(predicted_cvs: np.ndarray, ref_cvs: np.ndarray,
                        frames_dir: Path) -> list[dict]:
    """Match predicted CVs to nearest steering frame by L2 distance."""
    pdb_files = sorted(frames_dir.glob("*.pdb"))

    if len(pdb_files) != ref_cvs.shape[0]:
        print(f"  WARNING: {len(pdb_files)} PDBs vs {ref_cvs.shape[0]} CV entries")
        n = min(len(pdb_files), ref_cvs.shape[0])
        pdb_files = pdb_files[:n]
        ref_cvs = ref_cvs[:n]

    matches = []
    for i, pred in enumerate(predicted_cvs):
        dists = np.linalg.norm(ref_cvs - pred, axis=1)
        best_idx = int(np.argmin(dists))
        matches.append({
            "afm_frame": i,
            "predicted_cv": pred.tolist(),
            "matched_pdb": str(pdb_files[best_idx].name),
            "matched_cv": ref_cvs[best_idx].tolist(),
            "cv_error_A": float(dists[best_idx]),
            "conformer_index": best_idx,
        })

    return matches


def plot_predictions(predicted_cvs, ref_cvs, matches, output_dir):
    """Generate visualization of predictions."""
    n_frames = len(predicted_cvs)
    n_cvs = predicted_cvs.shape[1]
    cv_names = [
        "αHead↔αCalf",
        "βHead↔βTail",
        "αHead↔βHead",
    ]
    if n_cvs > len(cv_names):
        cv_names.extend([f"CV{i}" for i in range(len(cv_names), n_cvs)])

    fig, axes = plt.subplots(n_cvs + 1, 1, figsize=(12, 3 * (n_cvs + 1)),
                             sharex=True)

    frame_idx = np.arange(n_frames)

    # Plot each CV over time
    for j in range(n_cvs):
        ax = axes[j]
        ax.plot(frame_idx, predicted_cvs[:, j], 'b-', linewidth=1.5,
                label="Predicted", alpha=0.8)
        if ref_cvs is not None:
            # Show range of reference CVs
            ax.axhspan(ref_cvs[:, j].min(), ref_cvs[:, j].max(),
                       alpha=0.1, color="green", label="Training range")
        ax.set_ylabel(f"{cv_names[j]} (Å)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Plot matching error
    errors = [m["cv_error_A"] for m in matches]
    axes[-1].bar(frame_idx, errors, color="coral", alpha=0.7)
    axes[-1].set_ylabel("Match Error (Å)")
    axes[-1].set_xlabel("AFM Frame")
    axes[-1].grid(alpha=0.3)

    fig.suptitle("HS-AFM → Integrin Conformational State Prediction",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_dir / "afm_predictions.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved prediction plot: {output_dir / 'afm_predictions.png'}")

    # Conformer index over time
    fig2, ax2 = plt.subplots(figsize=(10, 4))
    conf_idx = [m["conformer_index"] for m in matches]
    ax2.plot(frame_idx, conf_idx, 'ko-', markersize=3, linewidth=1)
    ax2.set_xlabel("AFM Frame")
    ax2.set_ylabel("Matched Conformer Index\n(0=bent → N=extended)")
    ax2.set_title("Conformational Trajectory from HS-AFM")
    ax2.grid(alpha=0.3)
    fig2.savefig(output_dir / "conformer_trajectory.png", dpi=150,
                 bbox_inches="tight")
    plt.close()
    print(f"Saved conformer trajectory: {output_dir / 'conformer_trajectory.png'}")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Extract GIF frames
    raw_frames = extract_gif_frames(args.gif)

    # 2. Load CNN
    model, config = load_cnn(args.model, args.afmfold_root, args.device)

    # 3. Predict CVs
    predicted_cvs = predict_cvs(model, raw_frames, args.image_size, args.device)

    # 4. Match to steering conformers (if frames available)
    ref_cvs = None
    matches = None
    if args.cv_distances and args.cv_distances.exists():
        ref_cvs = np.load(str(args.cv_distances))
    elif args.frames_dir and args.frames_dir.exists():
        # Compute CVs from frames on the fly
        print("Computing reference CVs from steering frames...")
        import mdtraj as md
        pdb_files = sorted(args.frames_dir.glob("*.pdb"))
        if pdb_files:
            from train_afmfold_cnn import load_protein_frames, compute_domain_cvs
            traj = load_protein_frames(args.frames_dir)
            ref_cvs = compute_domain_cvs(traj, config.get("protein_name", "avb3"),
                                         args.afmfold_root)

    if ref_cvs is not None and args.frames_dir:
        matches = match_to_conformers(predicted_cvs, ref_cvs, args.frames_dir)
        # Save matches
        with open(args.output_dir / "matches.json", "w") as f:
            json.dump(matches, f, indent=2)
        print(f"Saved {len(matches)} frame-conformer matches")

    # 5. Save predictions
    np.save(str(args.output_dir / "predicted_cvs.npy"), predicted_cvs)

    # 6. Plot
    plot_predictions(predicted_cvs, ref_cvs, matches or [], args.output_dir)

    # 7. Summary
    print(f"\n=== Results ===")
    print(f"Frames processed: {len(raw_frames)}")
    print(f"Predicted CVs: {args.output_dir / 'predicted_cvs.npy'}")
    if matches:
        errors = [m["cv_error_A"] for m in matches]
        print(f"Mean match error: {np.mean(errors):.2f} Å")
        print(f"Matches: {args.output_dir / 'matches.json'}")
    print(f"Plots: {args.output_dir / 'afm_predictions.png'}")


if __name__ == "__main__":
    main()
