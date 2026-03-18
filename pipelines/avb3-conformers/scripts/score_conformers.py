#!/usr/bin/env python3
"""Score pulled conformer frames against MSA-subsampled Protenix predictions.

For each pulled frame, computes structural similarity (CA-RMSD after optimal
superposition, and a TM-score proxy) against Protenix predictions at each
MSA depth. Outputs a CSV matrix and a heatmap plot.

The key hypothesis: valid conformers should have some MSA depth where the
predicted structure is similar (low RMSD, high TM-score). Invalid (over-
stretched) conformers should have poor scores at ALL MSA depths.

Usage:
    python score_conformers.py \
        --frames-dir ../RoyalMD/results-avb3_04/splits/ \
        --predictions-dir data/runs/avb3/msa_validation/predictions/ \
        --output-dir data/runs/avb3/msa_validation/scores/
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np


def parse_pdb_ca_coords(pdb_path: Path) -> np.ndarray:
    """Extract CA atom coordinates from a PDB file. Returns (N, 3) array."""
    coords = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
    return np.array(coords, dtype=np.float64)


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Compute RMSD after optimal superposition (Kabsch algorithm).

    Args:
        P, Q: (N, 3) coordinate arrays (must have same N).

    Returns:
        RMSD in same units as input coordinates.
    """
    assert P.shape == Q.shape, f"Shape mismatch: {P.shape} vs {Q.shape}"
    # Center
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    # Covariance matrix
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    # Correct for reflection
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    # Optimal rotation
    R = Vt.T @ sign_matrix @ U.T
    P_rot = P_c @ R.T
    rmsd = np.sqrt(np.mean(np.sum((P_rot - Q_c) ** 2, axis=1)))
    return float(rmsd)


def tm_score_approx(P: np.ndarray, Q: np.ndarray) -> float:
    """Approximate TM-score after Kabsch superposition.

    Uses the Zhang & Skolnick (2004) normalization:
        TM = (1/L) * sum_i 1/(1 + (d_i/d0)^2)
    where d0 = 1.24 * (L - 15)^(1/3) - 1.8
    """
    assert P.shape == Q.shape
    L = len(P)
    if L < 20:
        return 0.0
    d0 = 1.24 * (L - 15) ** (1.0 / 3) - 1.8
    d0 = max(d0, 0.5)

    # Superpose
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T
    P_rot = P_c @ R.T

    distances = np.sqrt(np.sum((P_rot - Q_c) ** 2, axis=1))
    scores = 1.0 / (1.0 + (distances / d0) ** 2)
    return float(np.sum(scores) / L)


def find_prediction_dirs(predictions_base: Path) -> dict[str, list[Path]]:
    """Find Protenix prediction CIF/PDB files organized by MSA depth."""
    results: dict[str, list[Path]] = {}

    for depth_dir in sorted(predictions_base.iterdir()):
        if not depth_dir.is_dir():
            continue
        label = depth_dir.name  # e.g., "depth_1.00", "depth_0.50"

        # Look for prediction CIF/PDB files (Protenix output structure)
        pred_files = []
        for pattern in ["**/*sample_*.cif", "**/*sample_*.pdb",
                        "**/predictions/*.cif", "**/predictions/*.pdb"]:
            pred_files.extend(depth_dir.glob(pattern))

        # Deduplicate and prefer CIF
        seen_stems = set()
        unique_preds = []
        for f in sorted(pred_files):
            stem = f.stem.replace("_summary_confidence", "")
            if stem not in seen_stems and not f.name.endswith("_summary_confidence_sample_0.json"):
                seen_stems.add(stem)
                unique_preds.append(f)

        if unique_preds:
            results[label] = unique_preds

    return results


def extract_frame_number(filename: str) -> int | None:
    """Extract frame number from filename like 'production_trajectory_frame_042.pdb'."""
    m = re.search(r'(\d{3,})', filename)
    return int(m.group(1)) if m else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--frames-dir", type=Path, required=True,
                        help="Directory of pulled PDB frames.")
    parser.add_argument("--predictions-dir", type=Path, required=True,
                        help="Directory of Protenix predictions per MSA depth.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for scores CSV and plots.")
    parser.add_argument("--frame-indices", default=None,
                        help="Comma-separated frame indices to score (default: all).")
    parser.add_argument("--max-ca-mismatch", type=int, default=50,
                        help="Max allowed CA count difference before skipping.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load pulled frames
    frame_files = sorted(args.frames_dir.glob("*.pdb"))
    if args.frame_indices:
        wanted = set(int(x) for x in args.frame_indices.split(","))
        frame_files = [f for f in frame_files
                       if extract_frame_number(f.name) in wanted]

    if not frame_files:
        print("ERROR: No frame PDB files found.")
        return 1
    print(f"Loaded {len(frame_files)} pulled frames")

    # Find predictions
    pred_dirs = find_prediction_dirs(args.predictions_dir)
    if not pred_dirs:
        print("ERROR: No prediction directories found.")
        print(f"Looked in: {args.predictions_dir}")
        print("Expected structure: predictions_dir/depth_X.XX/.../*.cif")
        return 1
    print(f"Found predictions for {len(pred_dirs)} MSA depths: "
          f"{', '.join(pred_dirs.keys())}")

    # Score matrix: frames × depths
    depth_labels = sorted(pred_dirs.keys())
    frame_indices = []
    rmsd_matrix = []
    tmscore_matrix = []

    for frame_pdb in frame_files:
        frame_idx = extract_frame_number(frame_pdb.name)
        if frame_idx is None:
            continue
        frame_indices.append(frame_idx)
        frame_ca = parse_pdb_ca_coords(frame_pdb)

        rmsd_row = []
        tm_row = []
        for depth_label in depth_labels:
            pred_files = pred_dirs[depth_label]
            best_rmsd = float("inf")
            best_tm = 0.0

            for pred_file in pred_files:
                try:
                    pred_ca = parse_pdb_ca_coords(pred_file)
                except Exception:
                    continue

                # Handle CA count mismatch (truncate to shorter)
                n_frame, n_pred = len(frame_ca), len(pred_ca)
                if abs(n_frame - n_pred) > args.max_ca_mismatch:
                    continue
                n_min = min(n_frame, n_pred)
                rmsd = kabsch_rmsd(frame_ca[:n_min], pred_ca[:n_min])
                tm = tm_score_approx(frame_ca[:n_min], pred_ca[:n_min])
                best_rmsd = min(best_rmsd, rmsd)
                best_tm = max(best_tm, tm)

            rmsd_row.append(best_rmsd if best_rmsd < float("inf") else np.nan)
            tm_row.append(best_tm)

        rmsd_matrix.append(rmsd_row)
        tmscore_matrix.append(tm_row)

    rmsd_matrix = np.array(rmsd_matrix)
    tmscore_matrix = np.array(tmscore_matrix)

    # Save CSV
    csv_path = args.output_dir / "conformer_scores.csv"
    with open(csv_path, "w") as f:
        header = "frame," + ",".join(f"rmsd_{d}" for d in depth_labels) + \
                 "," + ",".join(f"tm_{d}" for d in depth_labels)
        f.write(header + "\n")
        for i, idx in enumerate(frame_indices):
            rmsd_vals = ",".join(f"{v:.2f}" for v in rmsd_matrix[i])
            tm_vals = ",".join(f"{v:.4f}" for v in tmscore_matrix[i])
            f.write(f"{idx},{rmsd_vals},{tm_vals}\n")
    print(f"\nScores CSV: {csv_path}")

    # Summary: best TM-score across all depths for each frame
    summary_path = args.output_dir / "conformer_summary.json"
    summary = []
    for i, idx in enumerate(frame_indices):
        best_depth_idx = int(np.argmax(tmscore_matrix[i]))
        summary.append({
            "frame": idx,
            "best_tm_score": float(tmscore_matrix[i, best_depth_idx]),
            "best_tm_depth": depth_labels[best_depth_idx],
            "best_rmsd": float(rmsd_matrix[i, best_depth_idx]),
            "min_rmsd": float(np.nanmin(rmsd_matrix[i])),
            "min_rmsd_depth": depth_labels[int(np.nanargmin(rmsd_matrix[i]))],
        })
    Path(summary_path).write_text(json.dumps(summary, indent=2))
    print(f"Summary: {summary_path}")

    # Generate heatmap plot
    try:
        _plot_heatmaps(frame_indices, depth_labels, rmsd_matrix, tmscore_matrix,
                       args.output_dir)
    except ImportError as e:
        print(f"Skipping plots (missing dependency: {e})")

    return 0


def _plot_heatmaps(frame_indices, depth_labels, rmsd_matrix, tmscore_matrix,
                   output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Subsample frame labels for readability
    n_frames = len(frame_indices)
    tick_step = max(1, n_frames // 20)

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    # TM-score heatmap
    ax = axes[0]
    im = ax.imshow(tmscore_matrix.T, aspect="auto", cmap="RdYlGn",
                   vmin=0, vmax=1)
    ax.set_ylabel("MSA Depth")
    ax.set_yticks(range(len(depth_labels)))
    ax.set_yticklabels(depth_labels)
    ax.set_title("TM-score: Pulled Frame vs Protenix Prediction\n"
                 "(higher = more similar → conformer is MSA-supported)")
    plt.colorbar(im, ax=ax, label="TM-score")

    # RMSD heatmap
    ax = axes[1]
    rmsd_clipped = np.clip(rmsd_matrix, 0, 30)  # clip for visualization
    im = ax.imshow(rmsd_clipped.T, aspect="auto", cmap="RdYlGn_r",
                   vmin=0, vmax=30)
    ax.set_ylabel("MSA Depth")
    ax.set_yticks(range(len(depth_labels)))
    ax.set_yticklabels(depth_labels)
    ax.set_xlabel("Pulled Frame Index")
    ax.set_xticks(range(0, n_frames, tick_step))
    ax.set_xticklabels([frame_indices[i] for i in range(0, n_frames, tick_step)],
                       rotation=45)
    ax.set_title("CA-RMSD (Å): Pulled Frame vs Protenix Prediction\n"
                 "(lower = more similar)")
    plt.colorbar(im, ax=ax, label="RMSD (Å)")

    plt.tight_layout()
    plot_path = output_dir / "conformer_validation_heatmap.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Heatmap: {plot_path}")
    plt.close()

    # Line plot: best TM-score per frame
    fig, ax = plt.subplots(figsize=(12, 5))
    best_tm = np.max(tmscore_matrix, axis=1)
    ax.plot(frame_indices, best_tm, "o-", markersize=3, color="#2ecc71",
            label="Best TM-score (any MSA depth)")
    ax.axhline(y=0.5, color="#e74c3c", ls="--", label="TM=0.5 (same fold)")
    ax.axhline(y=0.3, color="#f39c12", ls="--", label="TM=0.3 (marginal)")
    ax.set_xlabel("Pulled Frame Index")
    ax.set_ylabel("Best TM-score")
    ax.set_title("Conformer Validity: Best TM-score Across All MSA Depths\n"
                 "Expect: decreasing as frames become more stretched")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    line_path = output_dir / "conformer_validity_curve.png"
    plt.savefig(line_path, dpi=150, bbox_inches="tight")
    print(f"Validity curve: {line_path}")
    plt.close()


if __name__ == "__main__":
    raise SystemExit(main())
