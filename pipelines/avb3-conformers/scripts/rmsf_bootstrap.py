#!/usr/bin/env python3
"""Bootstrap-confidence per-residue RMSF — audit §10.9 follow-up to obj-027.

Resamples the V1+V2 fitted-trajectory frames with replacement, recomputes
the rotation-corrected per-residue RMSF over the αV β-propeller +
β3 βA-domain headpiece (Kabsch-aligned per frame), and reports 95 %
confidence bands. Annotates the top-20 most-flexible-residue ranking
with whether ranking survives 95 % CI overlap.

Output:
- figures/rmsf_bootstrap_v1.png
- results/afm_pipeline/rmsf_bootstrap/{rmsf_ci.npz, top20_stable.json}
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
V7 = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "rmsf_bootstrap"
FIG_DIR = ROOT / "figures"

HEAD_A = (1, 440)   # αV β-propeller + thigh
HEAD_B = (1, 350)   # β3 βA-domain


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    return cid if cid else chr(ord("A") + atom.residue.chain.index)


def collect_ca_meta(top) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (ca_idx, ca_chain, ca_resseq) arrays."""
    ca_idx, ca_chain, ca_resseq = [], [], []
    for atom in top.atoms:
        if atom.name != "CA":
            continue
        ca_idx.append(atom.index)
        ca_chain.append(0 if chain_id(atom) == "A" else 1)
        ca_resseq.append(atom.residue.resSeq)
    return np.array(ca_idx), np.array(ca_chain), np.array(ca_resseq)


def head_mask(ca_chain: np.ndarray, ca_resseq: np.ndarray) -> np.ndarray:
    return (
        ((ca_chain == 0) & (ca_resseq >= HEAD_A[0]) & (ca_resseq <= HEAD_A[1]))
        | ((ca_chain == 1) & (ca_resseq >= HEAD_B[0]) & (ca_resseq <= HEAD_B[1]))
    )


def kabsch_R(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Optimal rotation R such that P @ R best matches Q.
    P, Q already centered (zero centroid)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return U @ D @ Vt


def head_align_ca_only(coords_ca_nm: np.ndarray,
                       head_mask_ca: np.ndarray) -> np.ndarray:
    """Kabsch-align all frames to frame 0 using only headpiece CAs.

    coords_ca_nm: [T, n_ca, 3] in nm (mdtraj convention).
    head_mask_ca: bool over n_ca picking the headpiece CAs.

    Returns rotated/translated CA coords in nm, shape [T, n_ca, 3].
    """
    ref_head = coords_ca_nm[0, head_mask_ca]
    ref_centroid = ref_head.mean(axis=0)
    ref_centered = ref_head - ref_centroid
    aligned = np.empty_like(coords_ca_nm)
    aligned[0] = coords_ca_nm[0] - ref_centroid + ref_centroid  # identity
    for t in range(1, coords_ca_nm.shape[0]):
        head = coords_ca_nm[t, head_mask_ca]
        head_centroid = head.mean(axis=0)
        moving_centered = head - head_centroid
        R = kabsch_R(moving_centered, ref_centered)
        translated = coords_ca_nm[t] - head_centroid
        aligned[t] = translated @ R + ref_centroid
    return aligned


def rmsf_per_residue(coords_ca_aligned_nm: np.ndarray,
                     idx: np.ndarray | None = None) -> np.ndarray:
    """RMSF in Å per residue.

    coords_ca_aligned_nm: [T, n_ca, 3] in nm.
    idx: optional frame indices for bootstrap.
    """
    if idx is not None:
        coords = coords_ca_aligned_nm[idx]
    else:
        coords = coords_ca_aligned_nm
    mean = coords.mean(axis=0)
    delta = coords - mean
    return np.sqrt((delta**2).sum(axis=-1).mean(axis=0)) * 10.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--top-k", type=int, default=20)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Both videos share the same topology by construction — the v7
    # pipeline fits the same library against both AFM streams.
    top1 = md.load(str(V7 / "video1" / "topology.pdb")).topology
    ca_idx, ca_chain, ca_resseq = collect_ca_meta(top1)
    n_ca = len(ca_idx)
    print(f"{n_ca} CAs total")
    head_mask_ca = head_mask(ca_chain, ca_resseq)
    print(f"{int(head_mask_ca.sum())} headpiece CAs (α 1-440 + β 1-350)")

    print("Loading + aligning V1 ...")
    v1 = np.load(V7 / "video1" / "fitted_coords_median_reanchored.npy")
    v1_ca = v1[:, ca_idx, :]
    v1_aligned = head_align_ca_only(v1_ca, head_mask_ca)

    print("Loading + aligning V2 ...")
    v2 = np.load(V7 / "video2" / "fitted_coords_median_reanchored.npy")
    v2_ca = v2[:, ca_idx, :]
    v2_aligned = head_align_ca_only(v2_ca, head_mask_ca)

    combined = np.concatenate([v1_aligned, v2_aligned], axis=0)
    print(f"  combined {combined.shape[0]} frames")

    rmsf_v1 = rmsf_per_residue(v1_aligned)
    rmsf_v2 = rmsf_per_residue(v2_aligned)
    rmsf_combined = rmsf_per_residue(combined)
    print(f"V1 RMSF mean={rmsf_v1.mean():.2f} Å  max={rmsf_v1.max():.2f}")
    print(f"V2 RMSF mean={rmsf_v2.mean():.2f} Å  max={rmsf_v2.max():.2f}")
    print(f"combined RMSF mean={rmsf_combined.mean():.2f} Å  max={rmsf_combined.max():.2f}")

    print(f"\nBootstrapping {args.n_boot} replicates ...")
    rng = np.random.default_rng(args.seed)
    n_frames = combined.shape[0]
    boot = np.empty((args.n_boot, n_ca), dtype=np.float64)
    for k in range(args.n_boot):
        idx = rng.integers(0, n_frames, size=n_frames)
        boot[k] = rmsf_per_residue(combined, idx)
        if (k + 1) % 100 == 0:
            print(f"  {k+1}/{args.n_boot}")

    rmsf_lo = np.quantile(boot, 0.025, axis=0)
    rmsf_hi = np.quantile(boot, 0.975, axis=0)
    rmsf_med = np.quantile(boot, 0.50, axis=0)
    ci_width = rmsf_hi - rmsf_lo

    # Top-K ranking stability
    rank_order = np.argsort(-rmsf_combined)
    top_k_idx = rank_order[:args.top_k]
    top_k_plus_idx = rank_order[args.top_k:args.top_k + 1]
    cutoff_lo_kth = rmsf_lo[top_k_idx[-1]]
    cutoff_hi_kplus_first = rmsf_hi[top_k_plus_idx[0]]
    rank_stable = bool(cutoff_lo_kth > cutoff_hi_kplus_first)

    # Per-residue stability inside top-K: does residue i's RMSF 95 %
    # CI dominate residue (i+1)'s in the same ranking?
    boot_rank = np.argsort(-boot, axis=1)  # [n_boot, n_ca]
    top_k_membership = np.zeros((args.n_boot, args.top_k), dtype=int)
    for k in range(args.top_k):
        top_k_membership[:, k] = boot_rank[:, k]
    consistency = []
    for k in range(args.top_k):
        # Fraction of bootstrap iterations where the point-estimate
        # k-th residue is in the top-K of that resample
        target = top_k_idx[k]
        in_top_k = (boot_rank[:, :args.top_k] == target).any(axis=1).mean()
        consistency.append(float(in_top_k))

    np.savez(OUT_DIR / "rmsf_ci.npz",
             rmsf_v1=rmsf_v1, rmsf_v2=rmsf_v2, rmsf_combined=rmsf_combined,
             rmsf_lo=rmsf_lo, rmsf_hi=rmsf_hi, rmsf_med=rmsf_med,
             ci_width=ci_width,
             ca_chain=ca_chain, ca_resseq=ca_resseq,
             rank_order=rank_order,
             top_k_consistency=np.array(consistency),
             n_boot=args.n_boot)
    print(f"saved {OUT_DIR / 'rmsf_ci.npz'}")

    top20 = []
    for k, target in enumerate(top_k_idx):
        top20.append({
            "rank": k + 1,
            "ca_idx": int(target),
            "chain": "A" if ca_chain[target] == 0 else "B",
            "resSeq": int(ca_resseq[target]),
            "rmsf_A": float(rmsf_combined[target]),
            "ci_lo_A": float(rmsf_lo[target]),
            "ci_hi_A": float(rmsf_hi[target]),
            "ci_width_A": float(ci_width[target]),
            "top_20_consistency": float(consistency[k]),
        })
    with (OUT_DIR / "top20_stable.json").open("w") as f:
        json.dump({"top20": top20,
                   "rank_separation_stable": rank_stable,
                   "cutoff_lo_kth": float(cutoff_lo_kth),
                   "cutoff_hi_kplus_first": float(cutoff_hi_kplus_first),
                   "n_boot": args.n_boot,
                   "method": "Kabsch head-aligned, V1+V2 combined, "
                             "bootstrap with replacement"},
                  f, indent=2)

    print(f"\nTop {args.top_k} most-flexible residues (V1+V2 combined):")
    print(f"{'rank':>4s} {'chain':>5s} {'resSeq':>6s} {'RMSF (Å)':>10s} {'95% CI (Å)':>16s} {'top-20 hold rate':>16s}")
    for r in top20[:20]:
        print(f"  {r['rank']:>2d}  {r['chain']:>4s} {r['resSeq']:>5d}  "
              f"{r['rmsf_A']:>9.2f}  [{r['ci_lo_A']:.2f}, {r['ci_hi_A']:.2f}]  "
              f"{r['top_20_consistency']:>15.2f}")
    print(f"\ntop-{args.top_k} cutoff vs (k+1) — separation stable: {rank_stable}")
    print(f"  RMSF lower bound of {args.top_k}-th = {cutoff_lo_kth:.2f} Å")
    print(f"  RMSF upper bound of ({args.top_k}+1)-th = {cutoff_hi_kplus_first:.2f} Å")

    plot(rmsf_combined, rmsf_lo, rmsf_hi, ci_width,
         rmsf_v1, rmsf_v2,
         ca_chain, ca_resseq, top20, args.n_boot,
         rank_stable, cutoff_lo_kth, cutoff_hi_kplus_first)
    return 0


def plot(rmsf, rmsf_lo, rmsf_hi, ci_width,
         rmsf_v1, rmsf_v2,
         ca_chain, ca_resseq, top20, n_boot,
         rank_stable, cutoff_lo_kth, cutoff_hi_kplus_first):
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 2,
                          height_ratios=[1.7, 1.0, 1.2],
                          hspace=0.50, wspace=0.22)

    # Top: per-residue RMSF with 95% band, separated by chain
    ax = fig.add_subplot(gs[0, :])
    n_ca = rmsf.size
    x = np.arange(n_ca)
    chain_a = ca_chain == 0
    chain_b = ca_chain == 1

    for mask, color, label in [(chain_a, "#2166ac", "αV (chain A)"),
                                (chain_b, "#b2182b", "β3 (chain B)")]:
        ax.fill_between(x[mask], rmsf_lo[mask], rmsf_hi[mask],
                        color=color, alpha=0.22)
        ax.plot(x[mask], rmsf[mask], color=color, linewidth=0.8,
                label=label)
    # Mark top-20 residues
    for r in top20:
        ax.scatter([r["ca_idx"]], [r["rmsf_A"]],
                   s=30, color="red", zorder=10, alpha=0.9,
                   edgecolor="black", linewidth=0.4)
    ax.set_xlabel("CA index (chain A residues 1-940 → chain B residues 941+)")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title(
        f"Bootstrap-confidence per-residue RMSF (audit F3 RMSF analogue) — "
        f"{n_boot} resamples, V1+V2 combined (n=1645)\n"
        f"Red dots = top-20 most-flexible residues (point estimate). "
        f"Shaded = 95 % central bootstrap band.",
        fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(alpha=0.3)

    # Mid-left: CI width across residues
    ax = fig.add_subplot(gs[1, 0])
    ax.fill_between(x, 0, ci_width, color="#888888", alpha=0.7)
    ax.axhline(1.0, color="#cc7a00", linestyle="--",
               label="1 Å reference")
    ax.set_xlabel("CA index")
    ax.set_ylabel("RMSF 95 % CI width (Å)")
    ax.set_title(f"Per-residue CI width — median {np.median(ci_width):.2f} Å, "
                 f"max {ci_width.max():.2f} Å",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Mid-right: V1 vs V2 RMSF (sanity / consistency)
    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(rmsf_v1, rmsf_v2, s=4, alpha=0.45, color="#1b9e77")
    lim = max(rmsf_v1.max(), rmsf_v2.max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="#444444", linewidth=1,
            label="y=x")
    # Pearson correlation
    r = float(np.corrcoef(rmsf_v1, rmsf_v2)[0, 1])
    ax.text(0.05, 0.95,
            f"Pearson r = {r:.3f}\n(V1 vs V2 per-residue agreement)",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round", facecolor="white",
                      edgecolor="#cccccc"))
    ax.set_xlabel("V1 RMSF (Å)")
    ax.set_ylabel("V2 RMSF (Å)")
    ax.set_title("V1 vs V2 per-residue RMSF — consistency check",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)

    # Bottom: top-20 table
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    rows = [["rank", "chain:resSeq", "RMSF (Å)", "95 % CI (Å)",
             "top-20 hold rate"]]
    for r in top20:
        rows.append([
            f"{r['rank']}",
            f"{r['chain']}:{r['resSeq']}",
            f"{r['rmsf_A']:.2f}",
            f"[{r['ci_lo_A']:.2f}, {r['ci_hi_A']:.2f}]",
            f"{r['top_20_consistency']:.2f}",
        ])
    table = ax.table(cellText=rows, loc="center", cellLoc="center",
                     colWidths=[0.06, 0.16, 0.10, 0.20, 0.16])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.25)
    for k in range(len(rows[0])):
        table[(0, k)].set_text_props(weight="bold", color="white")
        table[(0, k)].set_facecolor("#444444")
    # Highlight low-stability rows
    for i, r in enumerate(top20):
        if r["top_20_consistency"] < 0.8:
            for k in range(len(rows[0])):
                table[(i + 1, k)].set_facecolor("#fff4e0")
    rank_msg = ("✓ stable" if rank_stable
                else "⚠ unstable")
    ax.set_title(f"Top-20 most-flexible residue ranking — "
                 f"separation cutoff: lower bound of #20 = "
                 f"{cutoff_lo_kth:.2f} Å, upper bound of #21 = "
                 f"{cutoff_hi_kplus_first:.2f} Å  → {rank_msg}\n"
                 f"(Highlighted rows = residues that fall out of top-20 "
                 f"in > 20 % of bootstrap resamples)",
                 fontsize=10, fontweight="bold")

    out_path = OUT_DIR / "rmsf_bootstrap.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig.savefig(FIG_DIR / "rmsf_bootstrap_v1.png",
                dpi=140, bbox_inches="tight")
    print(f"saved {out_path}")
    print(f"copied to {FIG_DIR / 'rmsf_bootstrap_v1.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
