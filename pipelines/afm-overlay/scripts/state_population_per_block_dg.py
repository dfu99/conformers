#!/usr/bin/env python3
"""State-population v3: per-block ΔG(CV0) drift across HS-AFM time (obj-045).

Audit-2026-05-05 §14.4 follow-up to obj-044. obj-044 rejected the
stationary-Bernoulli null on band fractions (max |z| = 6.7/7.5 in
V1/V2). This v3 visualizes *how* the underlying free-energy
landscape drifts: per 50-frame block, KDE-smooth P(CV0), compute
ΔG = −kT log P normalized to 0 at the block minimum, and overlay
all per-block curves on top of the pooled ΔG.

Outputs:
  figures/state_populations_v3_per_block_dg.png
  results/afm_pipeline/free_energy_profile/per_block_dg.npz
  results/afm_pipeline/free_energy_profile/per_block_dg_summary.json

Per-block summary stats:
  - location of the FES minimum (CV0_min) per block
  - block-min displacement vs pooled minimum
  - per-block ΔG at fixed CV0 = (50, 70, 80, 85) Å for landscape
    drift visualization.
"""
from __future__ import annotations

import json
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
WINDOW = 50

BC_MAX = 65.0
EC_MIN = 78.0
EO_MIN = 85.0
LIB_MIN, LIB_MAX = 47.3, 85.0


def block_dg(cv0_block: np.ndarray, grid: np.ndarray, bw: float = 0.30) -> np.ndarray:
    """Compute ΔG(grid) from a block of CV0 samples via KDE.

    Use a slightly broader bandwidth than the pooled (0.18) since
    each block has only ~50 samples, where Silverman's rule is
    finite-sample noisy.
    """
    if cv0_block.size < 5 or np.std(cv0_block) < 1e-3:
        return np.full_like(grid, np.nan, dtype=float)
    kde = gaussian_kde(cv0_block, bw_method=bw)
    p = np.clip(kde(grid), 1e-9, None)
    g = -KT_KCAL * np.log(p)
    return g - g.min()


def block_min_location(g: np.ndarray, grid: np.ndarray) -> float:
    if not np.all(np.isfinite(g)):
        return float("nan")
    return float(grid[int(np.argmin(g))])


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    data = np.load(SRC)
    cv0_v1 = np.asarray(data["cv0_v1"])
    cv0_v2 = np.asarray(data["cv0_v2"])
    grid = np.linspace(45, 95, 501)
    pooled_v1_dg = block_dg(cv0_v1, grid, bw=0.18)
    pooled_v2_dg = block_dg(cv0_v2, grid, bw=0.18)
    pooled_all_dg = block_dg(np.concatenate([cv0_v1, cv0_v2]), grid, bw=0.18)
    print(f"V1 frames: {cv0_v1.size}  V2 frames: {cv0_v2.size}")
    print(f"Pooled all FES min at CV0 = {block_min_location(pooled_all_dg, grid):.2f} Å")

    # --- Per-block ΔG ---
    n_v1_blocks = cv0_v1.size // WINDOW
    n_v2_blocks = cv0_v2.size // WINDOW
    print(f"V1 blocks: {n_v1_blocks}  V2 blocks: {n_v2_blocks}")

    g_v1_blocks = np.zeros((n_v1_blocks, grid.size))
    g_v2_blocks = np.zeros((n_v2_blocks, grid.size))
    centers_v1 = np.zeros(n_v1_blocks, dtype=int)
    centers_v2 = np.zeros(n_v2_blocks, dtype=int)
    min_loc_v1 = np.zeros(n_v1_blocks)
    min_loc_v2 = np.zeros(n_v2_blocks)

    for b in range(n_v1_blocks):
        sl = slice(b * WINDOW, (b + 1) * WINDOW)
        g_v1_blocks[b] = block_dg(cv0_v1[sl], grid, bw=0.30)
        centers_v1[b] = (sl.start + sl.stop) // 2
        min_loc_v1[b] = block_min_location(g_v1_blocks[b], grid)

    for b in range(n_v2_blocks):
        sl = slice(b * WINDOW, (b + 1) * WINDOW)
        g_v2_blocks[b] = block_dg(cv0_v2[sl], grid, bw=0.30)
        centers_v2[b] = (sl.start + sl.stop) // 2
        min_loc_v2[b] = block_min_location(g_v2_blocks[b], grid)

    pooled_v1_min_loc = block_min_location(pooled_v1_dg, grid)
    pooled_v2_min_loc = block_min_location(pooled_v2_dg, grid)
    print(f"Pooled V1 FES min at CV0 = {pooled_v1_min_loc:.2f} Å")
    print(f"Pooled V2 FES min at CV0 = {pooled_v2_min_loc:.2f} Å")

    print()
    print("Per-block FES minimum (CV0_min):")
    print(f"  V1: {min_loc_v1.tolist()}")
    print(f"  V2: {min_loc_v2.tolist()}")
    print(f"  V1 range: [{np.nanmin(min_loc_v1):.2f}, {np.nanmax(min_loc_v1):.2f}] Å")
    print(f"  V2 range: [{np.nanmin(min_loc_v2):.2f}, {np.nanmax(min_loc_v2):.2f}] Å")
    print(f"  V1 std:   {np.nanstd(min_loc_v1):.2f} Å")
    print(f"  V2 std:   {np.nanstd(min_loc_v2):.2f} Å")

    # ΔG at fixed CV0 across blocks
    cv0_probes = np.array([50.0, 70.0, 80.0, 85.0])
    probe_idx = np.array([np.argmin(np.abs(grid - c)) for c in cv0_probes])

    g_v1_probes = g_v1_blocks[:, probe_idx]
    g_v2_probes = g_v2_blocks[:, probe_idx]

    print()
    print("Per-block ΔG at fixed CV0 (kcal/mol):")
    for k, c in enumerate(cv0_probes):
        v1_vals = g_v1_probes[:, k]
        v2_vals = g_v2_probes[:, k]
        print(f"  CV0={c:5.1f}: V1 mean={np.nanmean(v1_vals):4.2f}±{np.nanstd(v1_vals):4.2f}  "
              f"V2 mean={np.nanmean(v2_vals):4.2f}±{np.nanstd(v2_vals):4.2f}")

    # Compute pooled-vs-per-block deviation:
    # at each CV0 bin, std of per-block ΔG. Large = drift in landscape.
    drift_v1 = np.nanstd(g_v1_blocks, axis=0)
    drift_v2 = np.nanstd(g_v2_blocks, axis=0)
    drift_v1_lib = drift_v1[(grid >= LIB_MIN) & (grid <= LIB_MAX)]
    drift_v2_lib = drift_v2[(grid >= LIB_MIN) & (grid <= LIB_MAX)]
    print()
    print(f"Per-block ΔG drift std over CV0 ∈ [{LIB_MIN}, {LIB_MAX}] Å:")
    print(f"  V1: median {np.nanmedian(drift_v1_lib):.2f}, max {np.nanmax(drift_v1_lib):.2f} kcal/mol")
    print(f"  V2: median {np.nanmedian(drift_v2_lib):.2f}, max {np.nanmax(drift_v2_lib):.2f} kcal/mol")

    # --- Save ---
    out_npz = OUT_DIR / "per_block_dg.npz"
    np.savez(out_npz,
             grid=grid,
             pooled_v1=pooled_v1_dg,
             pooled_v2=pooled_v2_dg,
             pooled_all=pooled_all_dg,
             g_v1_blocks=g_v1_blocks,
             g_v2_blocks=g_v2_blocks,
             centers_v1=centers_v1,
             centers_v2=centers_v2,
             min_loc_v1=min_loc_v1,
             min_loc_v2=min_loc_v2,
             drift_v1=drift_v1,
             drift_v2=drift_v2,
             window=WINDOW,
             T_K=T_KELVIN,
             kT_kcal=KT_KCAL)
    print(f"\nSaved {out_npz}")

    summary = {
        "T_K": T_KELVIN,
        "window_frames": WINDOW,
        "n_v1_blocks": int(n_v1_blocks),
        "n_v2_blocks": int(n_v2_blocks),
        "pooled_v1_fes_min_cv0_A": float(pooled_v1_min_loc),
        "pooled_v2_fes_min_cv0_A": float(pooled_v2_min_loc),
        "v1_block_fes_min_range_A": [float(np.nanmin(min_loc_v1)),
                                      float(np.nanmax(min_loc_v1))],
        "v2_block_fes_min_range_A": [float(np.nanmin(min_loc_v2)),
                                      float(np.nanmax(min_loc_v2))],
        "v1_block_fes_min_std_A": float(np.nanstd(min_loc_v1)),
        "v2_block_fes_min_std_A": float(np.nanstd(min_loc_v2)),
        "v1_drift_std_lib_kcal_median": float(np.nanmedian(drift_v1_lib)),
        "v1_drift_std_lib_kcal_max": float(np.nanmax(drift_v1_lib)),
        "v2_drift_std_lib_kcal_median": float(np.nanmedian(drift_v2_lib)),
        "v2_drift_std_lib_kcal_max": float(np.nanmax(drift_v2_lib)),
        "verdict": "FES_drifts_significantly" if (
            np.nanmax(drift_v1_lib) > 1.0 or np.nanmax(drift_v2_lib) > 1.0
        ) else "FES_stable_within_block_noise",
    }
    out_json = OUT_DIR / "per_block_dg_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"Saved {out_json}")

    # --- Plot ---
    fig = plt.figure(figsize=(13.5, 9.5))
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.20,
                           height_ratios=[1.1, 1.0, 0.9])

    cmap_v1 = plt.colormaps["viridis"]
    cmap_v2 = plt.colormaps["viridis"]

    # (a) V1 per-block ΔG overlay
    ax_v1 = fig.add_subplot(gs[0, 0])
    for b in range(n_v1_blocks):
        col = cmap_v1(b / max(n_v1_blocks - 1, 1))
        ax_v1.plot(grid, g_v1_blocks[b], color=col, lw=0.9, alpha=0.6,
                   label=f"block {b} (frames {b*WINDOW}-{(b+1)*WINDOW})")
    ax_v1.plot(grid, pooled_v1_dg, color="black", lw=2.4,
               label=f"pooled (n={cv0_v1.size})")
    ax_v1.axvspan(0, BC_MAX, color="#999999", alpha=0.08)
    ax_v1.axvspan(BC_MAX, EC_MIN, color="#fee08b", alpha=0.12)
    ax_v1.axvspan(EC_MIN, EO_MIN, color="#fdae61", alpha=0.18)
    ax_v1.axvspan(EO_MIN, 100, color="#d73027", alpha=0.13)
    ax_v1.set_xlim(45, 95)
    ax_v1.set_ylim(0, 6)
    ax_v1.set_title(f"V1 — per-block ΔG(CV0)  ({n_v1_blocks} blocks of {WINDOW})",
                    fontsize=10.5)
    ax_v1.set_ylabel("ΔG (kcal/mol)")
    ax_v1.set_xlabel("CV0 (Å)")
    ax_v1.grid(alpha=0.3)

    # (b) V2 per-block ΔG overlay
    ax_v2 = fig.add_subplot(gs[0, 1])
    for b in range(n_v2_blocks):
        col = cmap_v2(b / max(n_v2_blocks - 1, 1))
        ax_v2.plot(grid, g_v2_blocks[b], color=col, lw=0.7, alpha=0.45)
    ax_v2.plot(grid, pooled_v2_dg, color="black", lw=2.4,
               label=f"pooled (n={cv0_v2.size})")
    ax_v2.axvspan(0, BC_MAX, color="#999999", alpha=0.08)
    ax_v2.axvspan(BC_MAX, EC_MIN, color="#fee08b", alpha=0.12)
    ax_v2.axvspan(EC_MIN, EO_MIN, color="#fdae61", alpha=0.18)
    ax_v2.axvspan(EO_MIN, 100, color="#d73027", alpha=0.13)
    ax_v2.set_xlim(45, 95)
    ax_v2.set_ylim(0, 6)
    ax_v2.set_title(f"V2 — per-block ΔG(CV0)  ({n_v2_blocks} blocks of {WINDOW})",
                    fontsize=10.5)
    ax_v2.set_xlabel("CV0 (Å)")
    ax_v2.legend(loc="upper right", fontsize=8)
    ax_v2.grid(alpha=0.3)

    # (c) Drift std over blocks per CV0 bin
    ax_d = fig.add_subplot(gs[1, :])
    ax_d.fill_between(grid, 0, drift_v1, color="#2ca02c", alpha=0.4,
                      label="V1 std across blocks")
    ax_d.fill_between(grid, 0, drift_v2, color="#d62728", alpha=0.4,
                      label="V2 std across blocks")
    ax_d.axvspan(LIB_MIN, LIB_MAX, color="#aaaaaa", alpha=0.15,
                 label=f"library support [{LIB_MIN}, {LIB_MAX}] Å")
    ax_d.axvline(BC_MAX, color="#333333", ls="--", lw=0.8, alpha=0.5)
    ax_d.axvline(EC_MIN, color="#333333", ls="--", lw=0.8, alpha=0.5)
    ax_d.axvline(EO_MIN, color="#333333", ls="--", lw=0.8, alpha=0.5)
    ax_d.set_xlim(45, 95)
    ax_d.set_xlabel("CV0 (Å)")
    ax_d.set_ylabel("std(ΔG) across blocks (kcal/mol)")
    ax_d.set_title(f"Per-block ΔG drift — V1 max {np.nanmax(drift_v1_lib):.2f}, "
                   f"V2 max {np.nanmax(drift_v2_lib):.2f} kcal/mol over library support",
                   fontsize=10.5)
    ax_d.legend(loc="upper right", fontsize=8)
    ax_d.grid(alpha=0.3)

    # (d) FES min location per block over time
    ax_m = fig.add_subplot(gs[2, 0])
    ax_m.plot(centers_v1, min_loc_v1, "o-", color="#2ca02c", label="V1")
    ax_m.axhline(pooled_v1_min_loc, color="#2ca02c", ls="--", alpha=0.5,
                 label=f"V1 pooled = {pooled_v1_min_loc:.1f} Å")
    ax_m.set_xlabel("HS-AFM frame index (block center)")
    ax_m.set_ylabel("CV0 of FES minimum (Å)")
    ax_m.set_title(f"V1 FES-min drift  std={np.nanstd(min_loc_v1):.2f} Å",
                   fontsize=10)
    ax_m.legend(fontsize=8)
    ax_m.grid(alpha=0.3)

    ax_m2 = fig.add_subplot(gs[2, 1])
    ax_m2.plot(centers_v2, min_loc_v2, "o-", color="#d62728", label="V2")
    ax_m2.axhline(pooled_v2_min_loc, color="#d62728", ls="--", alpha=0.5,
                  label=f"V2 pooled = {pooled_v2_min_loc:.1f} Å")
    ax_m2.set_xlabel("HS-AFM frame index (block center)")
    ax_m2.set_ylabel("CV0 of FES minimum (Å)")
    ax_m2.set_title(f"V2 FES-min drift  std={np.nanstd(min_loc_v2):.2f} Å",
                    fontsize=10)
    ax_m2.legend(fontsize=8)
    ax_m2.grid(alpha=0.3)

    fig.suptitle("State-population v3 (obj-045): per-block ΔG(CV0) drift across "
                 "HS-AFM acquisition timeline",
                 fontsize=12, fontweight="bold", y=0.995)

    out_png = FIG_DIR / "state_populations_v3_per_block_dg.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"Saved {out_png}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
