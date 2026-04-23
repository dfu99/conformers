#!/usr/bin/env python3
"""Pairwise CA-RMSD across the 615-frame conformer library.

For each pair (i, j) of conformers, optimally superpose on CA atoms and
record the residual RMSD. Gives a 615x615 symmetric matrix. Also derive
per-residue variance: how much does this residue's position change
across the library (after superposition)?

This is intrinsic mechanical flexibility — not confounded by the rigid
rotations in the HS-AFM overlay.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frames-dir", type=Path, required=True,
                   help="Library of conformer PDBs (615 frames)")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--stride", type=int, default=1,
                   help="Subsample library by this stride (speeds up O(N²))")
    return p.parse_args()


def kabsch_align(P, Q):
    """Optimally align P onto Q using Kabsch, return rotated P."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return Pc @ R.T + Q.mean(axis=0)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(args.frames_dir.glob("*.pdb"))[::args.stride]
    n = len(pdbs)
    print(f"Loading {n} conformers (stride={args.stride})...")
    ca_list = []
    for i, p in enumerate(pdbs):
        t = md.load(str(p))
        ca = t.topology.select("name CA")
        ca_list.append(t.xyz[0, ca] * 10)  # nm → Å
        if (i+1) % 100 == 0:
            print(f"  {i+1}/{n}")
    X = np.stack(ca_list)  # [n, N_ca, 3]
    N_ca = X.shape[1]
    print(f"Loaded {n} conformers × {N_ca} CAs")

    # Superpose all to the first conformer
    print("Aligning all conformers to conformer 0...")
    ref = X[0]
    aligned = np.zeros_like(X)
    aligned[0] = ref
    for i in range(1, n):
        aligned[i] = kabsch_align(X[i], ref)

    # Per-residue variance across the library (std of position across n frames)
    pos_mean = aligned.mean(axis=0)
    pos_var = ((aligned - pos_mean) ** 2).sum(axis=-1).mean(axis=0)
    pos_std = np.sqrt(pos_var)  # Å, per residue
    print(f"Per-residue std: mean={pos_std.mean():.2f}Å, "
          f"min={pos_std.min():.2f}, max={pos_std.max():.2f}")
    np.save(str(args.output_dir / "per_residue_std_library.npy"), pos_std)

    # Pairwise RMSD matrix (O(n^2) but cheap — just broadcast)
    # Use aligned coords
    print(f"Computing {n}x{n} pairwise RMSD...")
    rmsd_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            d = aligned[i] - aligned[j]
            rmsd_mat[i, j] = np.sqrt((d ** 2).sum(axis=-1).mean())
            rmsd_mat[j, i] = rmsd_mat[i, j]
        if (i+1) % 50 == 0:
            print(f"  row {i+1}/{n}")
    np.save(str(args.output_dir / "rmsd_matrix.npy"), rmsd_mat)

    # Figure 1: RMSD matrix heatmap
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(rmsd_mat, cmap="viridis", aspect="equal")
    ax.set_xlabel("Conformer j")
    ax.set_ylabel("Conformer i")
    ax.set_title(f"Pairwise CA-RMSD matrix ({n} conformers, max {rmsd_mat.max():.1f}Å)")
    plt.colorbar(im, label="RMSD (Å)")
    plt.tight_layout()
    plt.savefig(str(args.output_dir / "cross_conformer_rmsd.png"), dpi=150)
    plt.close()

    # Figure 2: Per-residue variance line plot
    fig, ax = plt.subplots(figsize=(14, 4.5))
    x = np.arange(N_ca)
    ax.plot(x, pos_std, color="#2b8cbe", linewidth=0.8)
    # Highlight top 20
    top20 = np.argsort(pos_std)[-20:]
    ax.scatter(top20, pos_std[top20], color="red", s=30, zorder=5,
               label=f"Top 20 most variable (mean {pos_std[top20].mean():.1f}Å)")
    ax.set_xlabel("CA index")
    ax.set_ylabel("Per-residue std across library (Å)")
    ax.set_title(f"Per-residue flexibility across {n}-conformer library")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(str(args.output_dir / "per_residue_variance.png"), dpi=150)
    plt.close()

    # Annotate top-20 residues with chain/resSeq for README
    ref_top = md.load(str(pdbs[0]))
    ca_idx = ref_top.topology.select("name CA")
    ca_chain = np.array([ref_top.topology.atom(i).residue.chain.index for i in ca_idx])
    ca_res = np.array([ref_top.topology.atom(i).residue.resSeq for i in ca_idx])

    with open(args.output_dir / "top20_flexible_residues.md", "w") as f:
        f.write(f"# Top 20 most variable residues (CA std across {n} library conformers)\n\n")
        f.write("| rank | CA idx | chain | resSeq | std (Å) |\n")
        f.write("|------|--------|-------|--------|----------|\n")
        for rank, i in enumerate(reversed(top20), 1):
            f.write(f"| {rank} | {i} | {chr(ord('A') + ca_chain[i])} | "
                    f"{ca_res[i]} | {pos_std[i]:.2f} |\n")
    print(f"Saved results to {args.output_dir}/")


if __name__ == "__main__":
    main()
