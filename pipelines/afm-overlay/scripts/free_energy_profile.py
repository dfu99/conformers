#!/usr/bin/env python3
"""Free-energy profile ΔG(CV0) from V1+V2 fitted trajectories.

The 615-frame conformer library is a *biased* MD ensemble (CV-distance
steering). The fitted trajectories from the AFM-overlay pipeline are
*unbiased* — each frame is the conformer best-matching the experimental
HS-AFM image at that timepoint, so the population P(CV0) it samples is
the experimental population (modulo library completeness).

We therefore estimate ΔG(CV0) = −kT ln P(CV0) directly from the V1+V2
fitted-trajectory CV0 series. T = 300 K (HS-AFM stage temperature),
kT = 0.596 kcal/mol = 2.494 kJ/mol.

Limitations: (1) cannot probe CV0 > 85 Å because the library has no EO
templates — the EO region of ΔG is invisible to this method. (2) the
estimate is conditioned on the library spanning the relevant range,
which is BC + Intermediate + EC for αVβ3 v7. (3) does not account for
HS-AFM contact mechanics biasing the population (any binding to the
mica surface is implicit).

Outputs:
- results/afm_pipeline/free_energy_profile/free_energy.png
- results/afm_pipeline/free_energy_profile/delta_g.npy   (CV0 grid + ΔG)
- results/afm_pipeline/free_energy_profile/cv_series.npy (per-frame CV0)
- figures/free_energy_profile_v1.png  (top-level audit copy)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parents[2].parent
V7 = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final"

DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 956),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}

CV_PAIRS = [
    ("alpha_head_thigh", "alpha_calf"),     # CV0
    ("beta_head",        "beta_tail"),      # CV1
    ("alpha_head_thigh", "beta_head"),      # CV2
]

BC_CV0_MAX = 65.0
EC_CV0_MIN = 78.0

T_KELVIN = 300.0
KT_KCAL = 0.001987 * T_KELVIN  # ≈ 0.596 kcal/mol


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    if cid:
        return cid
    return chr(ord("A") + atom.residue.chain.index)


def domain_ca_indices(top, chain: str, lo: int, hi: int) -> np.ndarray:
    out = []
    for atom in top.atoms:
        if atom.name != "CA":
            continue
        if chain_id(atom) == chain and lo <= atom.residue.resSeq <= hi:
            out.append(atom.index)
    return np.array(out, dtype=np.int64)


def cv_series(coords: np.ndarray, idx_a: np.ndarray,
              idx_b: np.ndarray) -> np.ndarray:
    """coords: (T, N_atoms, 3) in nm (mdtraj convention).
    Returns (T,) distance in Å (×10 to match the rest of the pipeline)."""
    ca = coords[:, idx_a].mean(axis=1)
    cb = coords[:, idx_b].mean(axis=1)
    return np.linalg.norm(ca - cb, axis=-1) * 10.0


def compute_video_cvs(coords: np.ndarray, top) -> dict:
    domain_idx = {
        name: domain_ca_indices(top, *spec)
        for name, spec in DOMAINS.items()
    }
    out = {}
    for i, (a, b) in enumerate(CV_PAIRS):
        ia, ib = domain_idx[a], domain_idx[b]
        if ia.size == 0 or ib.size == 0:
            raise RuntimeError(f"No CAs for {a} or {b}")
        out[f"cv{i}"] = cv_series(coords, ia, ib)
    return out


def delta_g_from_population(cv0: np.ndarray,
                            grid: np.ndarray) -> tuple[np.ndarray, gaussian_kde]:
    kde = gaussian_kde(cv0, bw_method=0.18)
    p = kde(grid)
    p = np.clip(p, 1e-9, None)
    g = -KT_KCAL * np.log(p)
    g = g - g.min()
    return g, kde


def main() -> int:
    out_dir = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading topologies + fitted coords ...")
    top1 = md.load(str(V7 / "video1" / "topology.pdb")).topology
    top2 = md.load(str(V7 / "video2" / "topology.pdb")).topology

    v1 = np.load(V7 / "video1" / "fitted_coords_median_reanchored.npy")
    v2 = np.load(V7 / "video2" / "fitted_coords_median_reanchored.npy")
    print(f"  V1 frames: {v1.shape[0]}   V2 frames: {v2.shape[0]}")

    # fitted_coords_*.npy stores coordinates in Angstroms already, per
    # the rest of the overlay pipeline. Verify by typical CV0 scale.
    cvs1 = compute_video_cvs(v1, top1)
    cvs2 = compute_video_cvs(v2, top2)
    cv0_1 = cvs1["cv0"]
    cv0_2 = cvs2["cv0"]
    cv1_1, cv1_2 = cvs1["cv1"], cvs2["cv1"]
    cv2_1, cv2_2 = cvs1["cv2"], cvs2["cv2"]

    print(f"  V1 CV0  mean={cv0_1.mean():.2f}  range=[{cv0_1.min():.2f},{cv0_1.max():.2f}]")
    print(f"  V2 CV0  mean={cv0_2.mean():.2f}  range=[{cv0_2.min():.2f},{cv0_2.max():.2f}]")

    cv0_all = np.concatenate([cv0_1, cv0_2])

    grid = np.linspace(40, 100, 601)
    g_combined, _ = delta_g_from_population(cv0_all, grid)
    g_v1, _ = delta_g_from_population(cv0_1, grid)
    g_v2, _ = delta_g_from_population(cv0_2, grid)

    lib_min, lib_max = 47.3, 85.0

    np.savez(out_dir / "cv_series.npz",
             cv0_v1=cv0_1, cv0_v2=cv0_2,
             cv1_v1=cv1_1, cv1_v2=cv1_2,
             cv2_v1=cv2_1, cv2_v2=cv2_2,
             grid=grid,
             delta_g_combined_kcal=g_combined,
             delta_g_v1_kcal=g_v1,
             delta_g_v2_kcal=g_v2,
             T_kelvin=T_KELVIN,
             kT_kcal=KT_KCAL)
    np.save(out_dir / "delta_g.npy",
            np.column_stack([grid, g_combined, g_v1, g_v2]))
    print(f"  saved arrays to {out_dir}/")

    plot(grid, g_combined, g_v1, g_v2,
         cv0_all, cv0_1, cv0_2,
         lib_min, lib_max, out_dir)

    figures_dir = ROOT / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copy(out_dir / "free_energy.png",
                figures_dir / "free_energy_profile_v1.png")
    print(f"  copied to {figures_dir}/free_energy_profile_v1.png")
    return 0


def plot(grid, g_combined, g_v1, g_v2,
         cv0_all, cv0_1, cv0_2,
         lib_min, lib_max, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 7),
                             gridspec_kw={"height_ratios": [2.0, 1.0]},
                             sharex=True)

    ax = axes[0]
    # State bands (shaded by CV0 alone, since this is a CV0 plot).
    ax.axvspan(40, BC_CV0_MAX, color="#2166ac", alpha=0.10, zorder=0)
    ax.axvspan(BC_CV0_MAX, EC_CV0_MIN, color="#fdae61", alpha=0.12, zorder=0)
    ax.axvspan(EC_CV0_MIN, lib_max, color="#b2182b", alpha=0.10, zorder=0)
    ax.axvspan(lib_max, 100, color="#7f7f7f", alpha=0.18,
               hatch="///", zorder=0)
    ax.text((40 + BC_CV0_MAX) / 2, 0.18, "BC",
            ha="center", color="#2166ac", fontsize=11, fontweight="bold")
    ax.text((BC_CV0_MAX + EC_CV0_MIN) / 2, 0.18, "Intermediate",
            ha="center", color="#cc7a00", fontsize=10, fontweight="bold")
    ax.text((EC_CV0_MIN + lib_max) / 2, 0.18, "EC",
            ha="center", color="#b2182b", fontsize=11, fontweight="bold")
    ax.text((lib_max + 100) / 2, 0.18,
            "EO (no library coverage)\n— ΔG undefined here —",
            ha="center", color="#222222", fontsize=9, fontweight="bold")

    # Curves
    ax.plot(grid, g_combined, color="black", linewidth=2.2,
            label="V1+V2 combined (n=1645)")
    ax.plot(grid, g_v1, color="#1f77b4", linewidth=1.4, alpha=0.9,
            linestyle="--", label="V1 (n=379)")
    ax.plot(grid, g_v2, color="#2ca02c", linewidth=1.4, alpha=0.9,
            linestyle="--", label="V2 (n=1266)")

    # Library-coverage lower edge
    ax.axvline(lib_min, color="#444444", linestyle=":", linewidth=1)
    ax.text(lib_min, ax.get_ylim()[1] * 0.95 if False else 4.0,
            f" lib min {lib_min} Å",
            color="#444444", fontsize=8, ha="left")

    ax.set_ylim(-0.3, 5.0)
    ax.set_ylabel("ΔG (kcal/mol)", fontsize=11)
    ax.set_title("Experimental ΔG(CV0) from unbiased AFM-fitted trajectories\n"
                 "CV0 = ‖centroid(αV head+thigh) − centroid(αV calf)‖    "
                 "T=300 K, kT=0.596 kcal/mol",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
    ax.grid(alpha=0.3)

    # Bottom panel: histograms + KDE
    ax2 = axes[1]
    ax2.hist(cv0_all, bins=np.linspace(40, 100, 121), color="#888888",
             alpha=0.55, density=True, label="V1+V2 histogram")
    kde_v1 = gaussian_kde(cv0_1, bw_method=0.18)
    kde_v2 = gaussian_kde(cv0_2, bw_method=0.18)
    ax2.plot(grid, kde_v1(grid), color="#1f77b4", linewidth=1.4,
             linestyle="--", label="V1 KDE")
    ax2.plot(grid, kde_v2(grid), color="#2ca02c", linewidth=1.4,
             linestyle="--", label="V2 KDE")
    ax2.axvspan(lib_max, 100, color="#7f7f7f", alpha=0.18, hatch="///")
    ax2.set_xlim(40, 100)
    ax2.set_xlabel("CV0 — αV head ↔ αV calf centroid distance (Å)",
                   fontsize=11)
    ax2.set_ylabel("P(CV0)")
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.95)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out_path = out_dir / "free_energy.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    print(f"  saved {out_path}")


if __name__ == "__main__":
    raise SystemExit(main())
