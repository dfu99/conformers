#!/usr/bin/env python3
"""Predict integrin conformational states from HS-AFM GIF.

AFMFold pipeline (corrected):
  1. Extract frames from HS-AFM GIF
  2. CNN predicts inter-domain CVs (distances) per frame
  3. Cluster similar CVs to find unique conformational states
  4. For each unique state: run MD steering with those CV targets → PDB
  5. Map each HS-AFM frame to its generated PDB

Step 4 uses OpenMM domain steering (cv_distance_extend-style) with the
CNN-predicted distances as targets, NOT nearest-neighbor lookup.

Usage:
    python predict_from_afm_gif.py \
        --gif input.gif \
        --model data/runs/avb3/afmfold_training/model/best_model.pt \
        --reference-pdb data/avb3/AVB3_clean.pdb \
        --output-dir results/afm_predictions/
"""
from __future__ import annotations

import argparse
import json
import subprocess
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
    p.add_argument("--reference-pdb", type=Path, required=True,
                   help="Reference bent PDB for steering (starting structure).")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--image-size", type=int, default=35,
                   help="CNN input image size (must match training).")
    p.add_argument("--device", default="cpu")
    p.add_argument("--afmfold-root", type=Path,
                   default=Path("/home2/Documents/code/afmfold"))
    # MD generation options
    p.add_argument("--production-time", type=float, default=200.0,
                   help="Production time per conformer (ps). Short is fine since "
                        "we only need the endpoint structure.")
    p.add_argument("--n-clusters", type=int, default=None,
                   help="Number of CV clusters. If None, auto-select.")
    p.add_argument("--max-conformers", type=int, default=20,
                   help="Maximum unique conformers to generate via MD.")
    p.add_argument("--skip-md", action="store_true",
                   help="Skip MD generation (only predict CVs, useful for testing).")
    return p.parse_args()


# ── Step 1: Extract GIF frames ──────────────────────────────────────────────

def extract_gif_frames(gif_path: Path) -> list[np.ndarray]:
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


# ── Step 2: Preprocess + CNN inference ───────────────────────────────────────

def preprocess_frame(frame: np.ndarray, target_size: int) -> np.ndarray:
    """Preprocess a single AFM frame for CNN input.

    Center-crop to square, resize, normalize to [0, 1].
    """
    from PIL import Image

    h, w = frame.shape
    if h != w:
        s = min(h, w)
        y0 = (h - s) // 2
        x0 = (w - s) // 2
        frame = frame[y0:y0+s, x0:x0+s]

    img = Image.fromarray(frame)
    img = img.resize((target_size, target_size), Image.BILINEAR)
    arr = np.array(img, dtype=np.float32)

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
    """Run CNN inference on preprocessed frames → predicted CVs (Angstroms)."""
    import torch

    processed = np.stack([preprocess_frame(f, image_size) for f in frames])
    tensor = torch.from_numpy(processed).float().unsqueeze(1).to(device)

    with torch.no_grad():
        preds = model(tensor).cpu().numpy()

    print(f"Predicted CVs for {len(frames)} frames, shape={preds.shape}")
    return preds


# ── Step 3: Cluster CVs ─────────────────────────────────────────────────────

def cluster_cvs(predicted_cvs: np.ndarray, n_clusters: int | None = None,
                max_clusters: int = 20) -> tuple[np.ndarray, np.ndarray]:
    """Cluster predicted CVs to find unique conformational states.

    Returns:
        (cluster_centers, labels): centers in Angstroms, per-frame labels
    """
    from sklearn.cluster import KMeans

    n_frames = len(predicted_cvs)

    if n_clusters is None:
        # Auto-select: use silhouette or simple heuristic
        n_clusters = min(max_clusters, max(2, n_frames // 5))

    n_clusters = min(n_clusters, n_frames)

    km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = km.fit_predict(predicted_cvs)
    centers = km.cluster_centers_

    # Report
    print(f"Clustered {n_frames} frames into {n_clusters} conformational states")
    for i in range(n_clusters):
        count = (labels == i).sum()
        print(f"  Cluster {i}: {count} frames, "
              f"CVs = [{', '.join(f'{v:.1f}' for v in centers[i])}] Å")

    return centers, labels


# ── Step 4: Generate PDBs via MD steering ────────────────────────────────────

def generate_conformer_pdb(
    cluster_id: int,
    target_distances_angstrom: np.ndarray,
    reference_pdb: Path,
    output_dir: Path,
    production_time: float = 200.0,
) -> Path | None:
    """Run OpenMM MD steering to generate a PDB for the given CV targets.

    Converts CNN-predicted distances (Angstroms) to nm and runs
    run_domain_steering.py with --target-distances.
    """
    conformer_dir = output_dir / f"conformer_{cluster_id:03d}"
    conformer_dir.mkdir(parents=True, exist_ok=True)

    # CNN outputs Angstroms, steering expects nm
    target_nm = target_distances_angstrom / 10.0
    target_str = ",".join(f"{d:.2f}" for d in target_nm)

    steering_script = Path(__file__).resolve().parent / "run_domain_steering.py"

    cmd = [
        sys.executable, str(steering_script),
        "--input-pdb", str(reference_pdb),
        "--output-dir", str(conformer_dir),
        "--target-distances", target_str,
        "--production-time", str(production_time),
    ]

    print(f"\n--- Generating conformer {cluster_id} ---")
    print(f"  Target CVs (Å): {target_distances_angstrom}")
    print(f"  Target CVs (nm): {target_nm}")

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=7200,  # 2h max
        )
        if result.returncode != 0:
            print(f"  ERROR: steering failed for cluster {cluster_id}")
            print(f"  stderr: {result.stderr[-500:]}")
            return None

        final_pdb = conformer_dir / "final.pdb"
        if final_pdb.exists():
            print(f"  Generated: {final_pdb}")
            return final_pdb
        else:
            print(f"  WARNING: no final.pdb produced")
            return None

    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT: conformer {cluster_id} exceeded 2h")
        return None


def generate_all_conformers(
    cluster_centers: np.ndarray,
    reference_pdb: Path,
    output_dir: Path,
    production_time: float = 200.0,
) -> dict[int, Path]:
    """Generate PDBs for all cluster centers via MD steering."""
    md_dir = output_dir / "md_conformers"
    md_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for i, center in enumerate(cluster_centers):
        pdb_path = generate_conformer_pdb(
            i, center, reference_pdb, md_dir, production_time,
        )
        if pdb_path is not None:
            results[i] = pdb_path

    print(f"\nGenerated {len(results)}/{len(cluster_centers)} conformers")
    return results


# ── Step 5: Map frames to conformers ────────────────────────────────────────

def map_frames_to_conformers(
    labels: np.ndarray,
    conformer_pdbs: dict[int, Path],
    predicted_cvs: np.ndarray,
    cluster_centers: np.ndarray,
) -> list[dict]:
    """Map each AFM frame to its generated conformer PDB."""
    mappings = []
    for i in range(len(labels)):
        cluster_id = int(labels[i])
        pdb_path = conformer_pdbs.get(cluster_id)
        mappings.append({
            "afm_frame": i,
            "predicted_cv_A": predicted_cvs[i].tolist(),
            "cluster_id": cluster_id,
            "cluster_center_A": cluster_centers[cluster_id].tolist(),
            "conformer_pdb": str(pdb_path) if pdb_path else None,
            "md_generated": pdb_path is not None,
        })
    return mappings


# ── Visualization ────────────────────────────────────────────────────────────

def plot_predictions(predicted_cvs, cluster_centers, labels, mappings,
                     output_dir):
    """Generate visualization of predictions and conformer assignments."""
    n_frames = len(predicted_cvs)
    n_cvs = predicted_cvs.shape[1]
    cv_names = ["αHead↔αTail", "βHead↔αTail", "αHead↔βTail"]
    if n_cvs > len(cv_names):
        cv_names.extend([f"CV{i}" for i in range(len(cv_names), n_cvs)])

    fig, axes = plt.subplots(n_cvs + 1, 1, figsize=(12, 3 * (n_cvs + 1)),
                             sharex=True)

    frame_idx = np.arange(n_frames)
    colors = plt.cm.viridis(np.linspace(0, 1, len(cluster_centers)))

    # Plot each CV over time, colored by cluster
    for j in range(n_cvs):
        ax = axes[j]
        for k in range(len(cluster_centers)):
            mask = labels == k
            ax.scatter(frame_idx[mask], predicted_cvs[mask, j],
                       c=[colors[k]], s=8, alpha=0.7, label=f"State {k}")
            ax.axhline(cluster_centers[k, j], color=colors[k],
                       linestyle="--", alpha=0.4, linewidth=1)
        ax.set_ylabel(f"{cv_names[j]} (Å)")
        if j == 0:
            ax.legend(fontsize=6, ncol=4, loc="upper right")
        ax.grid(alpha=0.3)

    # Cluster assignment over time
    ax = axes[-1]
    ax.scatter(frame_idx, labels, c=[colors[l] for l in labels], s=10)
    ax.set_ylabel("Conformer State")
    ax.set_xlabel("AFM Frame")
    ax.set_yticks(range(len(cluster_centers)))
    ax.grid(alpha=0.3)

    fig.suptitle("HS-AFM → Integrin Conformational State Prediction\n"
                 "(CVs predicted by CNN, PDBs generated by MD steering)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_dir / "afm_predictions.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_dir / 'afm_predictions.png'}")

    # CV trajectory in 2D projection (first two CVs)
    if n_cvs >= 2:
        fig2, ax2 = plt.subplots(figsize=(8, 6))
        sc = ax2.scatter(predicted_cvs[:, 0], predicted_cvs[:, 1],
                         c=frame_idx, cmap="coolwarm", s=15, alpha=0.7)
        ax2.scatter(cluster_centers[:, 0], cluster_centers[:, 1],
                    c="black", marker="X", s=200, edgecolors="white",
                    linewidths=2, zorder=10, label="Cluster centers")
        ax2.set_xlabel(f"{cv_names[0]} (Å)")
        ax2.set_ylabel(f"{cv_names[1]} (Å)")
        ax2.set_title("CV Space: Frame Progression")
        plt.colorbar(sc, label="Frame index")
        ax2.legend()
        ax2.grid(alpha=0.3)
        fig2.savefig(output_dir / "cv_space.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Saved: {output_dir / 'cv_space.png'}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Extract GIF frames
    print("=" * 60)
    print("STEP 1: Extract HS-AFM frames")
    print("=" * 60)
    raw_frames = extract_gif_frames(args.gif)

    # Step 2: CNN inference → predicted CVs
    print("\n" + "=" * 60)
    print("STEP 2: CNN inference → inter-domain CVs")
    print("=" * 60)
    model, config = load_cnn(args.model, args.afmfold_root, args.device)
    predicted_cvs = predict_cvs(model, raw_frames, args.image_size, args.device)
    np.save(str(args.output_dir / "predicted_cvs.npy"), predicted_cvs)

    # Step 3: Cluster CVs
    print("\n" + "=" * 60)
    print("STEP 3: Cluster CVs → unique conformational states")
    print("=" * 60)
    cluster_centers, labels = cluster_cvs(
        predicted_cvs, n_clusters=args.n_clusters,
        max_clusters=args.max_conformers,
    )
    np.save(str(args.output_dir / "cluster_centers.npy"), cluster_centers)
    np.save(str(args.output_dir / "cluster_labels.npy"), labels)

    # Step 4: Generate PDBs via MD steering
    conformer_pdbs = {}
    if not args.skip_md:
        print("\n" + "=" * 60)
        print("STEP 4: MD steering → PDB conformers")
        print("=" * 60)
        conformer_pdbs = generate_all_conformers(
            cluster_centers, args.reference_pdb, args.output_dir,
            production_time=args.production_time,
        )

    # Step 5: Map frames to conformers
    print("\n" + "=" * 60)
    print("STEP 5: Map frames → conformers")
    print("=" * 60)
    mappings = map_frames_to_conformers(
        labels, conformer_pdbs, predicted_cvs, cluster_centers,
    )
    with open(args.output_dir / "frame_conformer_map.json", "w") as f:
        json.dump(mappings, f, indent=2)

    # Step 6: Visualize
    print("\n" + "=" * 60)
    print("STEP 6: Visualization")
    print("=" * 60)
    plot_predictions(predicted_cvs, cluster_centers, labels, mappings,
                     args.output_dir)

    # Summary
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    n_generated = len(conformer_pdbs)
    n_total = len(cluster_centers)
    print(f"HS-AFM frames:      {len(raw_frames)}")
    print(f"Unique states:       {n_total}")
    print(f"PDBs generated:      {n_generated}/{n_total}")
    print(f"Predicted CVs:       {args.output_dir / 'predicted_cvs.npy'}")
    print(f"Cluster centers:     {args.output_dir / 'cluster_centers.npy'}")
    print(f"Frame→conformer map: {args.output_dir / 'frame_conformer_map.json'}")
    if conformer_pdbs:
        print(f"Conformer PDBs:      {args.output_dir / 'md_conformers/'}")
    print(f"Plots:               {args.output_dir / 'afm_predictions.png'}")


if __name__ == "__main__":
    main()
