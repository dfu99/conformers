#!/usr/bin/env python3
"""Per-state shape anisotropy + Rg + asphericity — audit §36 v27 / obj-067.

Decomposes per-frame integrin shape into inertia-tensor moments per
HMM state. Tests whether bent (BC) and extended (EC) states differ
in global shape anisotropy: bent should be more "blob-like" (low
asphericity) and extended more "rod-like" (high asphericity).

Method:
  1. Load V1 + V2 fitted_coords_median_reanchored.npy (1645 frames).
  2. Center each frame on its centroid (using all 1654 CAs).
  3. Compute gyration tensor S_{ab} = (1/N) Σ_i Δr_{ia} Δr_{ib}.
  4. Diagonalize → eigenvalues λ1 ≥ λ2 ≥ λ3.
  5. Derived shape parameters per frame:
       Rg² = λ1 + λ2 + λ3
       Asphericity b = λ1 - (λ2 + λ3) / 2          (∈ [0, 3λ1/2])
       Acylindricity c = λ2 - λ3                    (∈ [0, λ2])
       Relative anisotropy κ² = (b² + 0.75c²) / Rg⁴ (∈ [0, 1])
  6. Per-state mean / std / KS pairwise tests.

Output:
  - figures/shape_anisotropy_per_state_v1.png
  - results/afm_pipeline/shape_anisotropy_per_state/{summary.json, per_frame.npz}
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
V7 = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final"
HMM_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "hmm_states.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "shape_anisotropy_per_state"
FIG_DIR = ROOT / "figures"

STATE_NAMES = ["BC", "Inter", "EC"]
STATE_COLORS = ["#4a90e2", "#f5a623", "#d0021b"]


def gyration_tensor(coords_centered: np.ndarray) -> np.ndarray:
    """coords_centered: [N, 3] in nm; return 3x3 gyration tensor."""
    return (coords_centered.T @ coords_centered) / coords_centered.shape[0]


def shape_params(coords_nm: np.ndarray) -> tuple[float, float, float, float]:
    """Return (Rg [Å], asphericity [Å²], acylindricity [Å²], κ²)."""
    centroid = coords_nm.mean(axis=0)
    centered = coords_nm - centroid
    S = gyration_tensor(centered)
    eig = np.linalg.eigvalsh(S)  # sorted ascending
    lam = eig[::-1]  # descending λ1 ≥ λ2 ≥ λ3
    rg2 = float(lam.sum())
    rg = np.sqrt(rg2) * 10.0  # nm → Å
    asphericity = float(lam[0] - (lam[1] + lam[2]) / 2) * 100.0  # nm² → Å²
    acylindricity = float(lam[1] - lam[2]) * 100.0
    kappa_sq = float(((lam[0] - (lam[1] + lam[2]) / 2) ** 2
                      + 0.75 * (lam[1] - lam[2]) ** 2)
                     / (rg2 ** 2)) if rg2 > 0 else 0.0
    return rg, asphericity, acylindricity, kappa_sq


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading V1 + V2 ...")
    v1 = np.load(V7 / "video1" / "fitted_coords_median_reanchored.npy")
    v2 = np.load(V7 / "video2" / "fitted_coords_median_reanchored.npy")
    print(f"V1 {v1.shape} + V2 {v2.shape}")

    # Use all CAs (fitted_coords stores all atoms) — strip non-CA later
    # actually fitted_coords are all atoms, 25139 each frame. We need CA-only.
    # Use V7 topology to find CAs.
    import mdtraj as md
    top = md.load(str(V7 / "video1" / "topology.pdb")).topology
    ca_idx = np.array([a.index for a in top.atoms if a.name == "CA"])
    print(f"{ca_idx.size} CAs")

    v1_ca = v1[:, ca_idx, :]
    v2_ca = v2[:, ca_idx, :]
    combined = np.concatenate([v1_ca, v2_ca], axis=0)
    n_frames = combined.shape[0]
    print(f"Combined: {n_frames} frames × {ca_idx.size} CAs")

    print("Loading HMM Viterbi states ...")
    hmm = np.load(HMM_NPZ)
    states = np.concatenate([hmm["viterbi_v1"], hmm["viterbi_v2"]]).astype(int)
    assert states.size == n_frames

    print("Computing per-frame shape params ...")
    rg = np.zeros(n_frames)
    asph = np.zeros(n_frames)
    acyl = np.zeros(n_frames)
    kappa = np.zeros(n_frames)
    for t in range(n_frames):
        rg[t], asph[t], acyl[t], kappa[t] = shape_params(combined[t])
        if (t + 1) % 200 == 0:
            print(f"  {t+1}/{n_frames}")

    print("\nPer-state shape summary:")
    print(f"{'state':>5s} {'n':>5s} {'Rg (Å)':>14s} "
          f"{'asphericity (Å²)':>20s} "
          f"{'acylindricity (Å²)':>22s} "
          f"{'κ²':>10s}")

    summary = {"per_state": {}, "pairwise_ks": {}}
    for s in (0, 1, 2):
        mask = states == s
        n_s = int(mask.sum())
        d = {
            "n": n_s,
            "Rg_mean_A": float(rg[mask].mean()),
            "Rg_std_A": float(rg[mask].std(ddof=1)),
            "asphericity_mean_A2": float(asph[mask].mean()),
            "asphericity_std_A2": float(asph[mask].std(ddof=1)),
            "acylindricity_mean_A2": float(acyl[mask].mean()),
            "acylindricity_std_A2": float(acyl[mask].std(ddof=1)),
            "kappa_sq_mean": float(kappa[mask].mean()),
            "kappa_sq_std": float(kappa[mask].std(ddof=1)),
        }
        summary["per_state"][STATE_NAMES[s]] = d
        print(f"  {STATE_NAMES[s]:>4s}  {n_s:>4d}  "
              f"{d['Rg_mean_A']:>5.2f} ± {d['Rg_std_A']:>5.2f}  "
              f"{d['asphericity_mean_A2']:>10.1f} ± {d['asphericity_std_A2']:>5.1f}  "
              f"{d['acylindricity_mean_A2']:>10.1f} ± {d['acylindricity_std_A2']:>5.1f}  "
              f"{d['kappa_sq_mean']:>5.3f} ± {d['kappa_sq_std']:>5.3f}")

    # Pairwise KS tests for each shape param
    print("\nPairwise KS tests:")
    for label, arr in [("Rg", rg), ("asphericity", asph),
                        ("acylindricity", acyl), ("kappa_sq", kappa)]:
        summary["pairwise_ks"][label] = {}
        for i, j in [(0, 1), (1, 2), (0, 2)]:
            a = arr[states == i]
            b = arr[states == j]
            r = stats.ks_2samp(a, b)
            key = f"{STATE_NAMES[i]}_vs_{STATE_NAMES[j]}"
            summary["pairwise_ks"][label][key] = {
                "ks_D": float(r.statistic),
                "p": float(r.pvalue),
            }
            print(f"  {label:>14s} {key:>14s}: D = {r.statistic:.4f}, "
                  f"p = {r.pvalue:.2e}")

    summary["bonferroni_alpha_per_param"] = 0.05 / 12  # 4 params × 3 pairs
    bonf = summary["bonferroni_alpha_per_param"]
    summary["interpretation"] = {
        "rg_increases_BC_to_EC": (
            summary["per_state"]["EC"]["Rg_mean_A"] >
            summary["per_state"]["BC"]["Rg_mean_A"]
        ),
        "asphericity_increases_BC_to_EC": (
            summary["per_state"]["EC"]["asphericity_mean_A2"] >
            summary["per_state"]["BC"]["asphericity_mean_A2"]
        ),
        "kappa_increases_BC_to_EC": (
            summary["per_state"]["EC"]["kappa_sq_mean"] >
            summary["per_state"]["BC"]["kappa_sq_mean"]
        ),
        "delta_Rg_BC_to_EC_A": (
            summary["per_state"]["EC"]["Rg_mean_A"]
            - summary["per_state"]["BC"]["Rg_mean_A"]
        ),
        "delta_asphericity_BC_to_EC_A2": (
            summary["per_state"]["EC"]["asphericity_mean_A2"]
            - summary["per_state"]["BC"]["asphericity_mean_A2"]
        ),
    }

    with (OUT_DIR / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
    np.savez(OUT_DIR / "per_frame.npz",
             rg=rg, asphericity=asph, acylindricity=acyl, kappa_sq=kappa,
             states=states)
    print(f"\nsaved {OUT_DIR / 'summary.json'} + per_frame.npz")

    print(f"\nVerdict — extension makes integrin more rod-like:")
    print(f"  Rg BC → EC: {summary['interpretation']['delta_Rg_BC_to_EC_A']:+.2f} Å")
    print(f"  asphericity BC → EC: "
          f"{summary['interpretation']['delta_asphericity_BC_to_EC_A2']:+.1f} Å²")
    print(f"  bonferroni-sig pairs at α = {bonf:.4f}: ", end="")
    n_sig = sum(1 for label in summary["pairwise_ks"]
                for k in summary["pairwise_ks"][label]
                if summary["pairwise_ks"][label][k]["p"] < bonf)
    print(f"{n_sig} / 12")

    plot(rg, asph, acyl, kappa, states, summary, n_frames)
    return 0


def plot(rg, asph, acyl, kappa, states, summary, n_frames):
    fig = plt.figure(figsize=(15.5, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.0, 1.0],
                          hspace=0.55, wspace=0.20)

    # Top-left: Rg per state histogram
    ax = fig.add_subplot(gs[0, 0])
    bins = np.linspace(rg.min(), rg.max(), 40)
    for s in range(3):
        sub = rg[states == s]
        d = summary["per_state"][STATE_NAMES[s]]
        ax.hist(sub, bins=bins, density=True, alpha=0.5,
                color=STATE_COLORS[s],
                label=f"{STATE_NAMES[s]} (n={d['n']}, "
                      f"μ = {d['Rg_mean_A']:.1f} ± {d['Rg_std_A']:.1f})")
    ax.set_xlabel("Rg (Å)")
    ax.set_ylabel("density")
    ax.set_title("(a) Radius of gyration per state",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Top-right: asphericity per state
    ax = fig.add_subplot(gs[0, 1])
    bins = np.linspace(asph.min(), asph.max(), 40)
    for s in range(3):
        sub = asph[states == s]
        d = summary["per_state"][STATE_NAMES[s]]
        ax.hist(sub, bins=bins, density=True, alpha=0.5,
                color=STATE_COLORS[s],
                label=f"{STATE_NAMES[s]} (μ = "
                      f"{d['asphericity_mean_A2']:.0f} ± "
                      f"{d['asphericity_std_A2']:.0f})")
    ax.set_xlabel(r"asphericity $b$ (Å²)")
    ax.set_ylabel("density")
    ax.set_title(f"(b) Asphericity per state — extension "
                 f"{'increases' if summary['interpretation']['asphericity_increases_BC_to_EC'] else 'decreases'} "
                 f"asphericity",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Mid-left: relative anisotropy κ²
    ax = fig.add_subplot(gs[1, 0])
    bins = np.linspace(0, kappa.max() * 1.05, 40)
    for s in range(3):
        sub = kappa[states == s]
        d = summary["per_state"][STATE_NAMES[s]]
        ax.hist(sub, bins=bins, density=True, alpha=0.5,
                color=STATE_COLORS[s],
                label=f"{STATE_NAMES[s]} (κ² = {d['kappa_sq_mean']:.3f})")
    ax.set_xlabel(r"relative anisotropy $\kappa^2$")
    ax.set_ylabel("density")
    ax.axvline(0.0, color="black", linewidth=0.6, label="κ²=0 (sphere)")
    ax.axvline(1.0, color="black", linewidth=0.6, linestyle="--",
               label="κ²=1 (rod)")
    ax.set_title("(c) Relative anisotropy κ² per state",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Mid-right: 2-D scatter Rg vs asphericity colored by state
    ax = fig.add_subplot(gs[1, 1])
    for s in range(3):
        mask = states == s
        ax.scatter(rg[mask], asph[mask], s=12, alpha=0.55,
                   color=STATE_COLORS[s], edgecolor="none",
                   label=f"{STATE_NAMES[s]}")
    ax.set_xlabel("Rg (Å)")
    ax.set_ylabel(r"asphericity $b$ (Å²)")
    ax.set_title("(d) Rg vs asphericity scatter by state",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)

    # Bottom: summary table
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    bonf = summary["bonferroni_alpha_per_param"]
    rows = [["param", "BC mean ± std", "Inter mean ± std",
             "EC mean ± std", "Δ EC-BC", "BC↔Inter p", "Inter↔EC p",
             "BC↔EC p", "sig"]]
    for label, key in [
        ("Rg (Å)", "Rg"),
        ("asphericity (Å²)", "asphericity"),
        ("acylindricity (Å²)", "acylindricity"),
        ("κ²", "kappa_sq"),
    ]:
        sk = "kappa_sq" if key == "kappa_sq" else key
        per = summary["per_state"]
        if sk == "kappa_sq":
            rs = [per[s]["kappa_sq_mean"] for s in STATE_NAMES]
            stds = [per[s]["kappa_sq_std"] for s in STATE_NAMES]
        elif sk == "Rg":
            rs = [per[s]["Rg_mean_A"] for s in STATE_NAMES]
            stds = [per[s]["Rg_std_A"] for s in STATE_NAMES]
        elif sk == "asphericity":
            rs = [per[s]["asphericity_mean_A2"] for s in STATE_NAMES]
            stds = [per[s]["asphericity_std_A2"] for s in STATE_NAMES]
        else:
            rs = [per[s]["acylindricity_mean_A2"] for s in STATE_NAMES]
            stds = [per[s]["acylindricity_std_A2"] for s in STATE_NAMES]
        ks = summary["pairwise_ks"][sk]
        n_sig = sum(1 for k in ks if ks[k]["p"] < bonf)
        rows.append([
            label,
            f"{rs[0]:.2f} ± {stds[0]:.2f}",
            f"{rs[1]:.2f} ± {stds[1]:.2f}",
            f"{rs[2]:.2f} ± {stds[2]:.2f}",
            f"{rs[2] - rs[0]:+.2f}",
            f"{ks['BC_vs_Inter']['p']:.1e}",
            f"{ks['Inter_vs_EC']['p']:.1e}",
            f"{ks['BC_vs_EC']['p']:.1e}",
            f"{n_sig}/3",
        ])
    table = ax.table(cellText=rows, loc="center", cellLoc="center",
                     colWidths=[0.10, 0.14, 0.14, 0.14, 0.08,
                                 0.10, 0.10, 0.10, 0.05])
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 1.40)
    for k in range(len(rows[0])):
        table[(0, k)].set_text_props(weight="bold", color="white")
        table[(0, k)].set_facecolor("#444444")
    # Highlight rows where extension increased the param
    for i, label_key in enumerate([("Rg",), ("asphericity",),
                                    ("acylindricity",), ("kappa_sq",)]):
        param = label_key[0]
        per = summary["per_state"]
        if param == "kappa_sq":
            increased = per["EC"]["kappa_sq_mean"] > per["BC"]["kappa_sq_mean"]
        elif param == "Rg":
            increased = per["EC"]["Rg_mean_A"] > per["BC"]["Rg_mean_A"]
        elif param == "asphericity":
            increased = (per["EC"]["asphericity_mean_A2"]
                         > per["BC"]["asphericity_mean_A2"])
        else:
            increased = (per["EC"]["acylindricity_mean_A2"]
                         > per["BC"]["acylindricity_mean_A2"])
        color = "#e8f4ea" if increased else "#fcdedb"
        for k in range(len(rows[0])):
            table[(i + 1, k)].set_facecolor(color)
    ax.set_title(
        f"Per-state shape summary — green = increases on extension, "
        f"orange = decreases  (Bonferroni α = {bonf:.4f}, 4 params × 3 pairs = 12 tests)",
        fontsize=10, fontweight="bold")

    fig.suptitle(
        f"Per-state shape anisotropy decomposition (obj-067) — "
        f"V1+V2 fitted, n = {n_frames}\n"
        f"Tests whether extension makes the integrin more rod-like "
        f"(higher asphericity / κ²)",
        fontsize=10.5, fontweight="bold", y=0.995,
    )

    out = OUT_DIR / "shape_anisotropy_per_state.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    fig.savefig(FIG_DIR / "shape_anisotropy_per_state_v1.png",
                dpi=140, bbox_inches="tight")
    print(f"saved {out}")
    print(f"copied to {FIG_DIR / 'shape_anisotropy_per_state_v1.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
