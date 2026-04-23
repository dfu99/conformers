#!/usr/bin/env python3
"""Identify hinge residues via CA-CA-CA bond-angle angular variance.

For each residue i, compute the angle ∠(CA_{i-1}, CA_i, CA_{i+1}) at every
frame in the trajectory. Residues with the highest angular variance across
the trajectory are candidate hinges — places where the backbone flexes.

Pre-registered prediction: α-genu ~residue 600 (thigh→calf junction) and
β-knee ~residue 470 (β-head→β-tail junction) should be top hinges.
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
                   help="Library of conformer PDBs")
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def ca_triple_angles(ca_xyz):
    """Compute ∠(CA_{i-1}, CA_i, CA_{i+1}) for i in [1, N-2].

    ca_xyz: [T, N, 3] — T frames, N CA atoms.
    returns: [T, N-2] angles in degrees.
    """
    T, N, _ = ca_xyz.shape
    v1 = ca_xyz[:, :-2] - ca_xyz[:, 1:-1]  # CA_{i-1} - CA_i
    v2 = ca_xyz[:, 2:] - ca_xyz[:, 1:-1]   # CA_{i+1} - CA_i
    n1 = np.linalg.norm(v1, axis=-1, keepdims=True)
    n2 = np.linalg.norm(v2, axis=-1, keepdims=True)
    n1[n1 < 1e-9] = 1.0
    n2[n2 < 1e-9] = 1.0
    cos = (v1 * v2).sum(axis=-1) / (n1[..., 0] * n2[..., 0])
    cos = np.clip(cos, -1.0, 1.0)
    return np.degrees(np.arccos(cos))


DOMAIN_BOUNDARIES = {
    "α-genu (thigh→calf)": (0, 435),   # between head-thigh and calf
    "α-midcalf": (0, 600),              # within calf
    "β-knee (head→tail)": (1, 352),    # between β-head and β-tail
    "β-midtail": (1, 500),              # within β-tail
}


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(args.frames_dir.glob("*.pdb"))
    n = len(pdbs)
    print(f"Loading {n} conformers...")
    ca_list = []
    ca_chain = None
    ca_res = None
    for i, p in enumerate(pdbs):
        t = md.load(str(p))
        ca_idx = t.topology.select("name CA")
        ca_list.append(t.xyz[0, ca_idx])
        if ca_chain is None:
            ca_chain = np.array([t.topology.atom(j).residue.chain.index for j in ca_idx])
            ca_res = np.array([t.topology.atom(j).residue.resSeq for j in ca_idx])
        if (i+1) % 150 == 0:
            print(f"  {i+1}/{n}")
    X = np.stack(ca_list)  # [n, N_ca, 3]
    print(f"Loaded {n} conformers, {X.shape[1]} CAs")

    # Split by chain to avoid chain-boundary artifacts
    angles_all = np.zeros((n, X.shape[1] - 2))
    # Compute angles within each chain separately
    for ch_i in np.unique(ca_chain):
        mask = ca_chain == ch_i
        idx = np.where(mask)[0]
        if len(idx) < 3:
            continue
        chain_xyz = X[:, idx]  # [n, N_chain, 3]
        chain_angles = ca_triple_angles(chain_xyz)  # [n, N_chain-2]
        # Place into angles_all at positions idx[1:-1] (center residue of triples)
        angles_all[:, idx[1:-1] - 1] = chain_angles  # -1 because full array has len N_ca-2

    # Angular variance per residue
    angle_std = angles_all.std(axis=0)  # [N_ca-2]
    print(f"Angle std: mean={angle_std.mean():.2f}°, "
          f"min={angle_std.min():.2f}, max={angle_std.max():.2f}")
    np.save(str(args.output_dir / "angle_std_per_residue.npy"), angle_std)

    # Top 20 hinge candidates
    center_idx = np.arange(1, X.shape[1] - 1)  # center residue indices of triples
    top20 = np.argsort(angle_std)[-20:][::-1]

    # Figure: angular variance over CA index
    fig, ax = plt.subplots(figsize=(14, 4.5))
    x = center_idx
    ax.plot(x, angle_std, color="#2b8cbe", linewidth=0.7)
    ax.scatter(top20 + 1, angle_std[top20], color="red", s=40, zorder=5,
               label=f"Top 20 hinges (mean σ={angle_std[top20].mean():.1f}°)")
    # Annotate known integrin hinges
    # α-genu at the α-head-thigh / α-calf boundary (chain 0, resSeq ~435)
    # β-knee at β-head / β-tail boundary (chain 1, resSeq ~352)
    for label, (ch, rs) in DOMAIN_BOUNDARIES.items():
        mask = (ca_chain == ch) & (ca_res == rs)
        if mask.any():
            i = np.where(mask)[0][0]
            ax.axvline(i, color="green", linestyle="--", alpha=0.4)
            ax.text(i, angle_std.max() * 0.95, label, rotation=90,
                    fontsize=8, ha="right", va="top")
    ax.set_xlabel("CA index (α chain then β)")
    ax.set_ylabel("CA-CA-CA bond angle σ across library (°)")
    ax.set_title(f"Hinge candidates — angular variance across {n}-conformer library")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(str(args.output_dir / "hinge_candidates.png"), dpi=150)
    plt.close()

    # Write top-20 hinge candidates as markdown
    with open(args.output_dir / "hinge_candidates.md", "w") as f:
        f.write(f"# Top 20 hinge-angle candidates\n\n")
        f.write(f"σ computed as std of ∠(CA_{{i-1}}, CA_i, CA_{{i+1}}) across "
                f"{n} library conformers.\n\n")
        f.write("| rank | center CA idx | chain | resSeq | mean angle (°) | σ (°) |\n")
        f.write("|------|---------------|-------|--------|----------------|-------|\n")
        for rank, i_triple in enumerate(top20, 1):
            ca_i = i_triple + 1  # center residue index into ca_chain/ca_res
            mean_angle = angles_all[:, i_triple].mean()
            f.write(f"| {rank} | {ca_i} | {chr(ord('A') + ca_chain[ca_i])} | "
                    f"{ca_res[ca_i]} | {mean_angle:.1f} | {angle_std[i_triple]:.1f} |\n")

    print(f"Saved results to {args.output_dir}/")


if __name__ == "__main__":
    main()
