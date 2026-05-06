#!/usr/bin/env python3
"""Per-state domain-pair centroid distances — audit §37 v28 / obj-068.

Decomposes the BC → EC transition into pairwise domain-centroid
distances. The five canonical αVβ3 domains
(`pipelines/conformer-library/scripts/build_library_metadata.py`):
  - alpha_head_thigh  (A 1-435)   — αV β-propeller + thigh
  - alpha_calf        (A 436-741) — αV calf-1 + calf-2
  - alpha_tail        (A 742-956) — αV genu + foot
  - beta_head         (B 1-352)   — β3 βA + hybrid
  - beta_tail         (B 353-692) — β3 PSI + I-EGF + β-tail

10 pairwise distances per frame. Averaged per HMM Viterbi state
(obj-048) → identifies which domain pairs separate or stay together
during extension.

Output:
  - figures/domain_pair_distances_per_state_v1.png
  - results/afm_pipeline/domain_pair_distances/{summary.json, per_frame.npz}
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
V7 = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final"
HMM_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "hmm_states.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "domain_pair_distances"
FIG_DIR = ROOT / "figures"

DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 956),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}
STATE_NAMES = ["BC", "Inter", "EC"]
STATE_COLORS = ["#4a90e2", "#f5a623", "#d0021b"]


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    return cid if cid else chr(ord("A") + atom.residue.chain.index)


def domain_ca_mask(top, chain_letter: str, lo: int, hi: int) -> np.ndarray:
    mask = []
    for atom in top.atoms:
        if atom.name != "CA":
            continue
        ch = chain_id(atom)
        rs = atom.residue.resSeq
        mask.append(ch == chain_letter and lo <= rs <= hi)
    return np.array(mask)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading topology + computing per-domain CA masks ...")
    top = md.load(str(V7 / "video1" / "topology.pdb")).topology
    ca_idx = np.array([a.index for a in top.atoms if a.name == "CA"])
    n_ca = ca_idx.size
    print(f"{n_ca} CAs total")

    masks = {}
    for name, (chain, lo, hi) in DOMAINS.items():
        masks[name] = domain_ca_mask(top, chain, lo, hi)
        print(f"  {name:>20s}: {int(masks[name].sum()):>4d} CAs "
              f"({chain}:{lo}-{hi})")

    print("\nLoading V1 + V2 fitted coords ...")
    v1 = np.load(V7 / "video1" / "fitted_coords_median_reanchored.npy")
    v2 = np.load(V7 / "video2" / "fitted_coords_median_reanchored.npy")
    v1_ca = v1[:, ca_idx, :]
    v2_ca = v2[:, ca_idx, :]
    combined = np.concatenate([v1_ca, v2_ca], axis=0)
    n_frames = combined.shape[0]
    print(f"Combined: {n_frames} frames × {n_ca} CAs")

    print("Loading HMM Viterbi states ...")
    hmm = np.load(HMM_NPZ)
    states = np.concatenate([hmm["viterbi_v1"], hmm["viterbi_v2"]]).astype(int)
    assert states.size == n_frames

    print("\nComputing per-frame domain centroids ...")
    domain_names = list(DOMAINS.keys())
    centroids = np.zeros((n_frames, len(domain_names), 3))
    for d, name in enumerate(domain_names):
        m = masks[name]
        centroids[:, d, :] = combined[:, m, :].mean(axis=1)

    pairs = list(combinations(range(len(domain_names)), 2))
    pair_labels = [f"{domain_names[i]} ↔ {domain_names[j]}" for i, j in pairs]
    print(f"\n{len(pairs)} pairwise distances")

    print("Computing per-frame pair distances ...")
    distances = np.zeros((n_frames, len(pairs)))
    for k, (i, j) in enumerate(pairs):
        diff = centroids[:, i, :] - centroids[:, j, :]
        # nm → Å conversion
        distances[:, k] = np.linalg.norm(diff, axis=1) * 10.0

    print("\nPer-state pair-distance summary (Å):")
    summary = {"per_state": {STATE_NAMES[s]: {} for s in range(3)},
               "per_pair": {}}
    delta_per_pair = {}
    print(f"{'pair':>40s}  {'BC':>10s}  {'Inter':>10s}  {'EC':>10s}  "
          f"{'Δ EC-BC':>10s}")
    for k, label in enumerate(pair_labels):
        means = []
        stds = []
        for s in range(3):
            mask = states == s
            arr = distances[mask, k]
            means.append(float(arr.mean()))
            stds.append(float(arr.std(ddof=1)))
        delta = means[2] - means[0]
        delta_per_pair[label] = delta
        summary["per_pair"][label] = {
            "BC_mean_A": means[0],
            "BC_std_A": stds[0],
            "Inter_mean_A": means[1],
            "Inter_std_A": stds[1],
            "EC_mean_A": means[2],
            "EC_std_A": stds[2],
            "delta_EC_BC_A": delta,
            "delta_Inter_BC_A": means[1] - means[0],
            "delta_EC_Inter_A": means[2] - means[1],
        }
        print(f"{label:>40s}  "
              f"{means[0]:>5.1f}±{stds[0]:>3.1f}  "
              f"{means[1]:>5.1f}±{stds[1]:>3.1f}  "
              f"{means[2]:>5.1f}±{stds[2]:>3.1f}  "
              f"{delta:>+8.2f}")

    # Rank pairs by Δ EC-BC
    ranked_pairs = sorted(delta_per_pair.items(), key=lambda kv: -kv[1])
    print(f"\nPairs ranked by Δ EC-BC (largest separation first):")
    for label, delta in ranked_pairs:
        print(f"  {label:>40s}: Δ = {delta:+.2f} Å")

    # Per-state global summary
    for s in range(3):
        mask = states == s
        n_s = int(mask.sum())
        sub = distances[mask]
        summary["per_state"][STATE_NAMES[s]] = {
            "n": n_s,
            "mean_distance_A": float(sub.mean()),
            "max_distance_A": float(sub.max()),
            "min_distance_A": float(sub.min()),
        }

    summary["max_separation_pair"] = ranked_pairs[0][0]
    summary["max_separation_delta_A"] = float(ranked_pairs[0][1])
    summary["least_separation_pair"] = ranked_pairs[-1][0]
    summary["least_separation_delta_A"] = float(ranked_pairs[-1][1])

    with (OUT_DIR / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
    np.savez(OUT_DIR / "per_frame.npz",
             distances=distances,
             pair_labels=np.array(pair_labels),
             states=states)
    print(f"\nsaved {OUT_DIR / 'summary.json'} + per_frame.npz")

    plot(distances, states, pair_labels, summary, ranked_pairs, n_frames)
    return 0


def plot(distances, states, pair_labels, summary, ranked_pairs, n_frames):
    n_pairs = len(pair_labels)
    fig = plt.figure(figsize=(15.5, 12))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.4, 1.0],
                          hspace=0.55, wspace=0.20)

    # Top-left: per-state mean distance bar chart per pair
    ax = fig.add_subplot(gs[0, 0])
    short_labels = [lbl.replace("alpha_head_thigh", "α-head")
                       .replace("alpha_calf", "α-calf")
                       .replace("alpha_tail", "α-tail")
                       .replace("beta_head", "β-head")
                       .replace("beta_tail", "β-tail")
                       .replace(" ↔ ", "↔") for lbl in pair_labels]
    width = 0.27
    x = np.arange(n_pairs)
    for s in range(3):
        means = [summary["per_pair"][lbl][f"{STATE_NAMES[s]}_mean_A"]
                 for lbl in pair_labels]
        stds = [summary["per_pair"][lbl][f"{STATE_NAMES[s]}_std_A"]
                for lbl in pair_labels]
        ax.bar(x + (s - 1) * width, means, width=width,
               yerr=stds, capsize=2,
               color=STATE_COLORS[s], alpha=0.85, label=STATE_NAMES[s])
    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel("centroid distance (Å)")
    ax.set_title("(a) Per-state pairwise domain-centroid distances",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3, axis="y")

    # Top-right: Δ EC-BC ranked
    ax = fig.add_subplot(gs[0, 1])
    deltas = np.array([rp[1] for rp in ranked_pairs])
    labels = [rp[0].replace("alpha_head_thigh", "α-head")
                  .replace("alpha_calf", "α-calf")
                  .replace("alpha_tail", "α-tail")
                  .replace("beta_head", "β-head")
                  .replace("beta_tail", "β-tail")
                  .replace(" ↔ ", "↔") for rp in ranked_pairs]
    color_arr = np.where(deltas >= 0, "#d0021b", "#4a90e2")
    ax.barh(np.arange(len(deltas))[::-1], deltas, color=color_arr,
            alpha=0.85)
    ax.set_yticks(np.arange(len(deltas))[::-1])
    ax.set_yticklabels(labels, fontsize=9)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel("Δ centroid distance EC − BC (Å)")
    ax.set_title("(b) Pairs ranked by Δ EC-BC — top: separating, "
                 "bottom: stable",
                 fontsize=10, fontweight="bold")
    ax.grid(alpha=0.3, axis="x")
    for i, (lbl, dl) in enumerate(ranked_pairs):
        y = len(deltas) - 1 - i
        ax.text(dl + 0.5 if dl >= 0 else dl - 0.5, y,
                f"{dl:+.1f}", fontsize=8.5, va="center",
                ha="left" if dl >= 0 else "right")

    # Mid: per-pair distance histograms (3 most-separating pairs + α-head↔α-calf)
    top_pairs_to_plot = ranked_pairs[:4]
    for j, (label, delta) in enumerate(top_pairs_to_plot):
        ax = fig.add_subplot(gs[1, 0] if j < 2 else gs[1, 1])
        if j == 0 or j == 2:
            # Setup
            pass
        idx = pair_labels.index(label)
        bins = np.linspace(distances[:, idx].min(),
                           distances[:, idx].max(), 30)
        for s in range(3):
            mask = states == s
            ax.hist(distances[mask, idx], bins=bins, density=True,
                    alpha=0.45, color=STATE_COLORS[s],
                    label=f"{STATE_NAMES[s]} "
                          f"({summary['per_pair'][label][f'{STATE_NAMES[s]}_mean_A']:.1f} ± "
                          f"{summary['per_pair'][label][f'{STATE_NAMES[s]}_std_A']:.1f} Å)")
        # avoid stacking 4 plots in 2 axes
        if j > 1:
            continue
        short = label.replace("alpha_head_thigh", "α-head") \
                     .replace("alpha_calf", "α-calf") \
                     .replace("alpha_tail", "α-tail") \
                     .replace("beta_head", "β-head") \
                     .replace("beta_tail", "β-tail")
        ax.set_xlabel("centroid distance (Å)")
        ax.set_ylabel("density")
        ax.set_title(f"({chr(ord('c') + j)}) {short}  Δ = {delta:+.1f} Å",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=8.5, loc="upper right")
        ax.grid(alpha=0.3)

    # Bottom: summary table
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    rows = [["pair (short)", "BC mean ± std (Å)", "Inter mean ± std",
             "EC mean ± std", "Δ EC-BC", "Δ Inter-BC", "Δ EC-Inter",
             "rank"]]
    for rank_i, (label, delta) in enumerate(ranked_pairs):
        d = summary["per_pair"][label]
        short = label.replace("alpha_head_thigh", "α-head") \
                     .replace("alpha_calf", "α-calf") \
                     .replace("alpha_tail", "α-tail") \
                     .replace("beta_head", "β-head") \
                     .replace("beta_tail", "β-tail") \
                     .replace(" ↔ ", "↔")
        rows.append([
            short,
            f"{d['BC_mean_A']:.2f} ± {d['BC_std_A']:.2f}",
            f"{d['Inter_mean_A']:.2f} ± {d['Inter_std_A']:.2f}",
            f"{d['EC_mean_A']:.2f} ± {d['EC_std_A']:.2f}",
            f"{d['delta_EC_BC_A']:+.2f}",
            f"{d['delta_Inter_BC_A']:+.2f}",
            f"{d['delta_EC_Inter_A']:+.2f}",
            f"{rank_i + 1}",
        ])
    table = ax.table(cellText=rows, loc="center", cellLoc="center",
                     colWidths=[0.18, 0.13, 0.13, 0.13, 0.08, 0.08,
                                 0.08, 0.05])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.30)
    for k in range(len(rows[0])):
        table[(0, k)].set_text_props(weight="bold", color="white")
        table[(0, k)].set_facecolor("#444444")
    for i, (label, delta) in enumerate(ranked_pairs):
        if delta > 30:
            color = "#d0021b30"
        elif delta > 10:
            color = "#fcdedb"
        elif delta > 0:
            color = "#fff5e0"
        else:
            color = "#dde9f4"
        for k in range(len(rows[0])):
            table[(i + 1, k)].set_facecolor(color)
    ax.set_title("Domain-pair distances per state, ranked by Δ EC-BC. "
                 "(red = high separation, blue = stable / closer)",
                 fontsize=10, fontweight="bold")

    fig.suptitle(
        f"Per-state domain-pair distances (obj-068) — "
        f"V1+V2 fitted, n = {n_frames}, 5 αVβ3 domains × 10 pairs\n"
        f"Maximum separation: {summary['max_separation_pair']} "
        f"Δ = {summary['max_separation_delta_A']:+.2f} Å; "
        f"least: {summary['least_separation_pair']} "
        f"Δ = {summary['least_separation_delta_A']:+.2f} Å",
        fontsize=10.5, fontweight="bold", y=0.995,
    )

    out = OUT_DIR / "domain_pair_distances.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    fig.savefig(FIG_DIR / "domain_pair_distances_per_state_v1.png",
                dpi=140, bbox_inches="tight")
    print(f"saved {out}")
    print(f"copied to {FIG_DIR / 'domain_pair_distances_per_state_v1.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
