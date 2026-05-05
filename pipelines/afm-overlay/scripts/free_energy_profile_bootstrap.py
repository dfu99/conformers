#!/usr/bin/env python3
"""Bootstrap-confidence ΔG(CV0) — audit follow-up F3.

Resample the 1645 fitted-trajectory CV0 frames with replacement
(`n_boot=1000`), recompute ΔG(CV0) via the same KDE pipeline as
`free_energy_profile.py`, and plot a 95 % central band over the v1
curve. Validates obj-038 against finite-sample noise.

Output: figures/free_energy_profile_v2.png and
results/afm_pipeline/free_energy_profile/bootstrap.npz.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

T_KELVIN = 300.0
KT_KCAL = 0.001987 * T_KELVIN
LIB_MIN, LIB_MAX = 47.3, 85.0
BC_CV0_MAX = 65.0
EC_CV0_MIN = 78.0


def delta_g(cv0: np.ndarray, grid: np.ndarray, bw: float = 0.18) -> np.ndarray:
    kde = gaussian_kde(cv0, bw_method=bw)
    p = np.clip(kde(grid), 1e-9, None)
    g = -KT_KCAL * np.log(p)
    return g - g.min()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    data = np.load(SRC)
    cv0_v1 = data["cv0_v1"]
    cv0_v2 = data["cv0_v2"]
    cv0_all = np.concatenate([cv0_v1, cv0_v2])
    grid = data["grid"]
    g_point = np.array(data["delta_g_combined_kcal"])
    print(f"n_total={cv0_all.size}  n_v1={cv0_v1.size}  n_v2={cv0_v2.size}")
    print(f"running {args.n_boot} bootstrap replicates ...")

    rng = np.random.default_rng(args.seed)
    g_boot = np.empty((args.n_boot, grid.size), dtype=np.float64)
    for k in range(args.n_boot):
        sample = rng.choice(cv0_all, size=cv0_all.size, replace=True)
        g_boot[k] = delta_g(sample, grid)
        if (k + 1) % 100 == 0:
            print(f"  {k+1}/{args.n_boot}")

    lo = np.quantile(g_boot, 0.025, axis=0)
    hi = np.quantile(g_boot, 0.975, axis=0)
    median = np.quantile(g_boot, 0.50, axis=0)
    band_width = hi - lo

    np.savez(OUT_DIR / "bootstrap.npz",
             grid=grid, g_point=g_point, g_lo=lo, g_hi=hi,
             g_median=median, n_boot=args.n_boot, seed=args.seed,
             cv0_all=cv0_all)
    print(f"saved {OUT_DIR / 'bootstrap.npz'}")

    plot(grid, g_point, median, lo, hi, band_width, cv0_v1, cv0_v2,
         FIG_DIR, OUT_DIR, args.n_boot)
    return 0


def plot(grid, g_point, g_median, g_lo, g_hi, band_width,
         cv0_v1, cv0_v2, fig_dir, out_dir, n_boot):
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 7.0),
                             gridspec_kw={"height_ratios": [2.2, 1.0]},
                             sharex=True)

    ax = axes[0]
    ax.axvspan(40, BC_CV0_MAX, color="#2166ac", alpha=0.10, zorder=0)
    ax.axvspan(BC_CV0_MAX, EC_CV0_MIN, color="#fdae61", alpha=0.12, zorder=0)
    ax.axvspan(EC_CV0_MIN, LIB_MAX, color="#b2182b", alpha=0.10, zorder=0)
    ax.axvspan(LIB_MAX, 100, color="#7f7f7f", alpha=0.18,
               hatch="///", zorder=0)
    ax.text((40 + BC_CV0_MAX) / 2, 0.18, "BC", ha="center",
            color="#2166ac", fontsize=11, fontweight="bold")
    ax.text((BC_CV0_MAX + EC_CV0_MIN) / 2, 0.18, "Intermediate",
            ha="center", color="#cc7a00", fontsize=10, fontweight="bold")
    ax.text((EC_CV0_MIN + LIB_MAX) / 2, 0.18, "EC", ha="center",
            color="#b2182b", fontsize=11, fontweight="bold")
    ax.text((LIB_MAX + 100) / 2, 0.18,
            "EO (no library coverage)\n— ΔG undefined here —",
            ha="center", color="#222222", fontsize=9, fontweight="bold")

    ax.fill_between(grid, g_lo, g_hi, color="black", alpha=0.18,
                    label=f"95% bootstrap band (n={n_boot})")
    ax.plot(grid, g_median, color="black", linewidth=1.0, alpha=0.55,
            linestyle="--", label="bootstrap median")
    ax.plot(grid, g_point, color="black", linewidth=2.2,
            label=f"point estimate (n=1645)")

    ax.set_ylim(-0.3, 5.0)
    ax.set_ylabel("ΔG (kcal/mol)", fontsize=11)
    ax.set_title("Bootstrap-confidence ΔG(CV0) — audit F3 follow-up\n"
                 "1000 resamples (with replacement), KDE bw=0.18, T=300 K",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
    ax.grid(alpha=0.3)

    # Panel 2: band width across CV0
    ax2 = axes[1]
    ax2.fill_between(grid, 0, band_width, color="#888888", alpha=0.65)
    ax2.axvspan(LIB_MAX, 100, color="#7f7f7f", alpha=0.18, hatch="///")
    ax2.axhline(0.5, color="#cc7a00", linestyle="--", linewidth=1,
                label="0.5 kcal/mol reference")
    ax2.set_xlabel("CV0 — αV head ↔ αV calf centroid distance (Å)", fontsize=11)
    ax2.set_ylabel("ΔG 95 % CI width (kcal/mol)")
    ax2.set_xlim(40, 100)
    ax2.set_ylim(0, max(2.5, band_width[(grid > 50) & (grid < 90)].max() * 1.05))
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(alpha=0.3)

    # Annotate sample-size mass
    n_v1 = cv0_v1.size
    n_v2 = cv0_v2.size
    fig.text(0.01, 0.01, f"V1 n={n_v1}  •  V2 n={n_v2}  •  combined 1645",
             fontsize=8, alpha=0.6)

    fig.tight_layout()
    out_path = out_dir / "free_energy_bootstrap.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig.savefig(fig_dir / "free_energy_profile_v2.png",
                dpi=140, bbox_inches="tight")
    print(f"saved {out_path}")
    print(f"copied to {fig_dir / 'free_energy_profile_v2.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
