#!/usr/bin/env python3
"""Per-residue root-mean-square fluctuation (RMSF) across fitted trajectory.

After head-anchored alignment, compute per-CA RMSF — standard MD analysis.
Map to integrin domains (α-head-thigh, α-calf, α-coil, β-head, β-tail,
β-coil) and plot:
  (a) per-residue RMSF line plot colored by chain
  (b) heatmap over time windows
  (c) domain-averaged bar chart
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
    p.add_argument("--fitted-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--n-windows", type=int, default=8,
                   help="Number of time windows for heatmap")
    return p.parse_args()


DOMAIN_RANGES = {
    "α-head-thigh": (0, 1, 435),      # chain 0 (A), residues 1-435
    "α-calf":       (0, 436, 741),
    "α-coil":       (0, 742, 956),
    "β-head":       (1, 1, 352),
    "β-tail":       (1, 353, 692),
    # "β-coil": (1, 693, 999),  # often not in our topology
}
DOMAIN_COLORS = {
    "α-head-thigh": "#2c7bb6",
    "α-calf":       "#abd9e9",
    "α-coil":       "#d7191c",
    "β-head":       "#fdae61",
    "β-tail":       "#ffffbf",
    "β-coil":       "#2b83ba",
}


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fitted = np.load(str(args.fitted_dir / "fitted_coords_smooth.npy"))
    topo = md.load(str(args.fitted_dir / "topology.pdb")).topology
    ca_idx = np.array([a.index for a in topo.atoms if a.name == "CA"])
    n_frames = fitted.shape[0]
    print(f"{n_frames} frames, {len(ca_idx)} CA atoms")

    ca_traj = fitted[:, ca_idx, :]  # [T, N_ca, 3]
    # RMSF = sqrt(mean over t of |r(t) - r_mean|^2)
    mean_pos = ca_traj.mean(axis=0)
    delta = ca_traj - mean_pos
    rmsf_per_res = np.sqrt((delta ** 2).sum(axis=-1).mean(axis=0)) * 10  # nm→Å
    print(f"RMSF: mean={rmsf_per_res.mean():.2f}Å, "
          f"min={rmsf_per_res.min():.2f}, max={rmsf_per_res.max():.2f}")
    np.save(str(args.output_dir / "rmsf_per_residue.npy"), rmsf_per_res)

    # Tag each CA with its chain and residue
    ca_chain = np.array([topo.atom(i).residue.chain.index for i in ca_idx])
    ca_resseq = np.array([topo.atom(i).residue.resSeq for i in ca_idx])

    # Per-domain mean RMSF
    domain_means = {}
    for name, (ch, lo, hi) in DOMAIN_RANGES.items():
        mask = (ca_chain == ch) & (ca_resseq >= lo) & (ca_resseq <= hi)
        if mask.sum() > 0:
            domain_means[name] = rmsf_per_res[mask].mean()
            print(f"  {name}: {rmsf_per_res[mask].mean():.2f}Å (n={mask.sum()})")

    # Figure 1: per-residue RMSF line plot colored by chain
    fig, ax = plt.subplots(figsize=(14, 4.5))
    x = np.arange(len(rmsf_per_res))
    for ch_i, chain_name in [(0, "α (chain A)"), (1, "β (chain B)")]:
        mask = ca_chain == ch_i
        color = "#2c7bb6" if ch_i == 0 else "#d7191c"
        ax.plot(x[mask], rmsf_per_res[mask], label=chain_name, color=color, linewidth=0.8)
    # Shade domain regions
    for name, (ch, lo, hi) in DOMAIN_RANGES.items():
        mask = (ca_chain == ch) & (ca_resseq >= lo) & (ca_resseq <= hi)
        if mask.any():
            idx = np.where(mask)[0]
            ax.axvspan(idx[0], idx[-1], alpha=0.15,
                       color=DOMAIN_COLORS.get(name, "gray"),
                       label=name)
    ax.set_xlabel("CA index (α chain then β chain)")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title(f"Per-residue RMSF — {args.fitted_dir.name}")
    ax.legend(fontsize=7, ncol=4, loc="upper right")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(str(args.output_dir / "rmsf_per_residue.png"), dpi=150)
    plt.close()

    # Figure 2: heatmap over time windows
    win_size = n_frames // args.n_windows
    win_rmsf = np.zeros((args.n_windows, len(ca_idx)))
    for w in range(args.n_windows):
        lo = w * win_size
        hi = (w + 1) * win_size if w < args.n_windows - 1 else n_frames
        win_mean = ca_traj[lo:hi].mean(axis=0)
        win_delta = ca_traj[lo:hi] - win_mean
        win_rmsf[w] = np.sqrt((win_delta ** 2).sum(axis=-1).mean(axis=0)) * 10

    fig, ax = plt.subplots(figsize=(14, 5))
    im = ax.imshow(win_rmsf, aspect="auto", cmap="viridis",
                   extent=[0, len(ca_idx), args.n_windows, 0])
    ax.set_xlabel("CA index")
    ax.set_ylabel("Time window")
    ax.set_title(f"RMSF heatmap over time windows — {args.fitted_dir.name}")
    plt.colorbar(im, ax=ax, label="RMSF (Å)")
    plt.tight_layout()
    plt.savefig(str(args.output_dir / "rmsf_heatmap_timewindows.png"), dpi=150)
    plt.close()

    # Figure 3: domain-averaged bar chart
    fig, ax = plt.subplots(figsize=(9, 4.5))
    names = list(domain_means.keys())
    values = [domain_means[n] for n in names]
    colors = [DOMAIN_COLORS.get(n, "gray") for n in names]
    bars = ax.bar(names, values, color=colors, edgecolor="black")
    ax.set_ylabel("Mean RMSF (Å)")
    ax.set_title(f"Domain-averaged RMSF — {args.fitted_dir.name}")
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width()/2, v + 0.1, f"{v:.2f}",
                ha="center", fontsize=9)
    plt.xticks(rotation=20)
    plt.tight_layout()
    plt.savefig(str(args.output_dir / "rmsf_domain_averaged.png"), dpi=150)
    plt.close()

    # CSV table
    import csv
    with open(args.output_dir / "rmsf_table.csv", "w") as f:
        w = csv.writer(f)
        w.writerow(["ca_idx", "chain", "resSeq", "rmsf_A"])
        for i in range(len(rmsf_per_res)):
            w.writerow([i, ca_chain[i], ca_resseq[i], f"{rmsf_per_res[i]:.3f}"])
    print(f"Saved RMSF results to {args.output_dir}/")


if __name__ == "__main__":
    main()
