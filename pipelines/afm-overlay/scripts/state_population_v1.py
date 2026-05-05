#!/usr/bin/env python3
"""Boltzmann state-population breakdown + V1/V2 ΔG reproducibility (obj-043).

Audit-2026-05-05 §13 deepening pass v5. Pure-derivative analysis on the
existing 1645-frame fitted-trajectory cache (cv_series.npz from obj-038):

  1. Boltzmann-population fraction in BC, Intermediate, EC, EO bands
     under the unbiased experimental ΔG(CV0).
  2. V1 vs V2 ΔG cross-validation:
        - per-bin ΔG agreement (Pearson r over the supported region)
        - Wasserstein-1 distance between the V1 and V2 CV0 empirical
          distributions
        - Kolmogorov-Smirnov test
  3. EO-coverage gap quantification:
        - fraction of fitted CV0 frames in [85, 100] Å (audit §10.2 says 40)
        - Bayesian upper bound on P(CV0 ≥ 85 Å) given the observed counts
          (Jeffreys prior Beta(0.5, 0.5))
        - Boltzmann ΔG-implied EO population given the un-extrapolatable
          tail (i.e., what assumption is needed to claim ΔG_EO < X kcal/mol)

State definitions follow build_library_metadata.py + intuition.md:
  BC          : CV0 ≤ 65 Å
  Intermediate: 65 < CV0 < 78 Å
  EC          : CV0 ≥ 78 Å (still bent CV2)
  EO          : CV0 ≥ 85 Å + CV2 ≥ 50 Å (we don't have CV2 ≥ 50 — gap!)

For state populations from the ΔG profile the boundary is on CV0 only;
EO is reported as the CV0 ≥ 85 Å population since CV2 ≥ 50 Å is empty
in the fitted set.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp, wasserstein_distance, pearsonr, beta

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

T_KELVIN = 300.0
KT_KCAL = 0.001987 * T_KELVIN

BC_MAX = 65.0
EC_MIN = 78.0
EO_MIN = 85.0
LIB_MIN, LIB_MAX = 47.3, 85.0


def boltzmann_population(grid: np.ndarray, g_kcal: np.ndarray, lo: float, hi: float) -> float:
    mask = (grid >= lo) & (grid < hi)
    if not mask.any():
        return 0.0
    p = np.exp(-g_kcal / KT_KCAL)
    return float(p[mask].sum() / p.sum())


def empirical_fraction(cv0: np.ndarray, lo: float, hi: float) -> float:
    return float(((cv0 >= lo) & (cv0 < hi)).sum() / cv0.size)


def jeffreys_upper(k: int, n: int, alpha: float = 0.05) -> float:
    """Beta(0.5+k, 0.5+n−k) two-sided 95 % upper bound on the Bernoulli rate."""
    a = 0.5 + k
    b = 0.5 + (n - k)
    return float(beta.ppf(1 - alpha / 2, a, b))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    data = np.load(SRC)
    cv0_v1 = np.asarray(data["cv0_v1"])
    cv0_v2 = np.asarray(data["cv0_v2"])
    cv0_all = np.concatenate([cv0_v1, cv0_v2])
    grid = np.asarray(data["grid"])
    g_combined = np.asarray(data["delta_g_combined_kcal"])
    g_v1 = np.asarray(data["delta_g_v1_kcal"])
    g_v2 = np.asarray(data["delta_g_v2_kcal"])

    print(f"n_total={cv0_all.size}  n_v1={cv0_v1.size}  n_v2={cv0_v2.size}")
    print(f"grid range: [{grid.min()}, {grid.max()}]  step={grid[1]-grid[0]:.3f}")

    # --- 1. Boltzmann populations from the ΔG profile ---
    bc = boltzmann_population(grid, g_combined, 0, BC_MAX)
    inter = boltzmann_population(grid, g_combined, BC_MAX, EC_MIN)
    ec = boltzmann_population(grid, g_combined, EC_MIN, EO_MIN)
    eo = boltzmann_population(grid, g_combined, EO_MIN, 200)

    # Empirical for sanity
    e_bc = empirical_fraction(cv0_all, 0, BC_MAX)
    e_inter = empirical_fraction(cv0_all, BC_MAX, EC_MIN)
    e_ec = empirical_fraction(cv0_all, EC_MIN, EO_MIN)
    e_eo = empirical_fraction(cv0_all, EO_MIN, 200)

    print()
    print(f"{'band':<12} {'CV0 (Å)':<14} {'Boltz P':<10} {'Empirical':<10}")
    for name, lo, p_b, p_e in [
        ("BC", "≤65", bc, e_bc),
        ("Inter.", "65-78", inter, e_inter),
        ("EC", "78-85", ec, e_ec),
        ("EO*", "≥85", eo, e_eo),
    ]:
        print(f"{name:<12} {f'{lo}':<14} {p_b:>6.3f}    {p_e:>6.3f}")
    print("(* EO band uses CV0≥85 Å only; CV2≥50 Å is empty in fit set.)")

    # --- 2. V1 vs V2 ΔG cross-validation ---
    supp = (grid >= LIB_MIN) & (grid <= LIB_MAX)
    r_full, _ = pearsonr(g_v1[supp], g_v2[supp])
    r_kt, _ = pearsonr(np.minimum(g_v1[supp], 5.0), np.minimum(g_v2[supp], 5.0))
    w1 = float(wasserstein_distance(cv0_v1, cv0_v2))
    ks_stat, ks_p = ks_2samp(cv0_v1, cv0_v2, alternative="two-sided")

    print()
    print("V1 vs V2 reproducibility:")
    print(f"  ΔG Pearson r over CV0 ∈ [{LIB_MIN}, {LIB_MAX}] Å:  {r_full:.3f}")
    print(f"  ΔG Pearson r (clipped to ≤5 kcal/mol):           {r_kt:.3f}")
    print(f"  CV0-distribution Wasserstein-1 distance:          {w1:.2f} Å")
    print(f"  Kolmogorov-Smirnov:  D={ks_stat:.4f}  p={ks_p:.2e}")

    # --- 3. EO-gap quantification ---
    n_eo = int(((cv0_all >= EO_MIN)).sum())
    n_total = int(cv0_all.size)
    p_eo_emp = n_eo / n_total
    p_eo_upper = jeffreys_upper(n_eo, n_total, 0.05)
    n_eo_v1 = int((cv0_v1 >= EO_MIN).sum())
    n_eo_v2 = int((cv0_v2 >= EO_MIN).sum())
    print()
    print(f"EO-coverage gap:")
    print(f"  fitted CV0 ≥ 85 Å: {n_eo}/{n_total} ({100*p_eo_emp:.2f} %)")
    print(f"  V1: {n_eo_v1}/{cv0_v1.size}  V2: {n_eo_v2}/{cv0_v2.size}")
    print(f"  Jeffreys 97.5 % upper bound on P(CV0≥85 Å): {100*p_eo_upper:.2f} %")
    print(f"  -kT log P upper bound: ΔG_EO ≥ {-KT_KCAL*np.log(p_eo_upper):.2f} kcal/mol")

    # --- Save JSON ---
    out = {
        "T_K": T_KELVIN,
        "kT_kcal": KT_KCAL,
        "n_v1": cv0_v1.size,
        "n_v2": cv0_v2.size,
        "boltzmann_populations": {
            "BC": bc, "Intermediate": inter, "EC": ec, "EO_CV0_only": eo,
        },
        "empirical_populations": {
            "BC": e_bc, "Intermediate": e_inter, "EC": e_ec, "EO_CV0_only": e_eo,
        },
        "v1_v2_reproducibility": {
            "delta_g_pearson_r_lib_support": r_full,
            "delta_g_pearson_r_kT_clipped": r_kt,
            "cv0_wasserstein_A": w1,
            "ks_D": ks_stat,
            "ks_p": ks_p,
        },
        "eo_gap": {
            "n_fitted_eo": n_eo,
            "n_total": n_total,
            "p_eo_empirical": p_eo_emp,
            "p_eo_jeffreys_upper95": p_eo_upper,
            "delta_g_eo_lower_bound_kcal": -KT_KCAL * float(np.log(p_eo_upper)),
        },
    }
    out_json = OUT_DIR / "state_populations.json"
    out_json.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {out_json}")

    # --- Plot ---
    fig = plt.figure(figsize=(11, 7.5))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.30,
                          height_ratios=[1.0, 1.0])

    # (a) ΔG combined with state-population bars overlaid
    ax_g = fig.add_subplot(gs[0, :])
    ax_g.plot(grid, g_combined, color="#1f77b4", lw=1.8,
              label="ΔG combined (1645 frames)")
    ax_g.plot(grid, g_v1, color="#2ca02c", lw=1.0, alpha=0.7,
              label=f"ΔG V1 (n={cv0_v1.size})")
    ax_g.plot(grid, g_v2, color="#d62728", lw=1.0, alpha=0.7,
              label=f"ΔG V2 (n={cv0_v2.size})")
    ax_g.axvspan(0, BC_MAX, color="#999999", alpha=0.10, label="BC")
    ax_g.axvspan(BC_MAX, EC_MIN, color="#fee08b", alpha=0.20, label="Inter.")
    ax_g.axvspan(EC_MIN, EO_MIN, color="#fdae61", alpha=0.25, label="EC")
    ax_g.axvspan(EO_MIN, 100, color="#d73027", alpha=0.20, label="EO* (gap)")
    ax_g.set_xlim(40, 100)
    ax_g.set_ylim(0, 8)
    ax_g.set_xlabel("CV0 head-tail distance (Å)")
    ax_g.set_ylabel("ΔG (kcal/mol)")
    ax_g.set_title(f"ΔG(CV0) — V1 vs V2 reproducibility: Pearson r = {r_full:.3f}, "
                   f"Wasserstein-1 = {w1:.2f} Å",
                   fontsize=10.5)
    ax_g.legend(loc="upper right", fontsize=8, ncol=2)
    ax_g.grid(alpha=0.3)

    # (b) State populations bar
    ax_s = fig.add_subplot(gs[1, 0])
    states = ["BC", "Inter.", "EC", "EO*"]
    boltz = [bc, inter, ec, eo]
    emp = [e_bc, e_inter, e_ec, e_eo]
    x = np.arange(len(states))
    w = 0.36
    ax_s.bar(x - w/2, boltz, w, color="#1f77b4", label="Boltzmann (from ΔG)")
    ax_s.bar(x + w/2, emp, w, color="#ff7f0e", label="Empirical (count)")
    for i, (b_v, e_v) in enumerate(zip(boltz, emp)):
        ax_s.text(i - w/2, b_v + 0.005, f"{100*b_v:.1f}%", ha="center",
                  va="bottom", fontsize=8)
        ax_s.text(i + w/2, e_v + 0.005, f"{100*e_v:.1f}%", ha="center",
                  va="bottom", fontsize=8)
    ax_s.set_xticks(x)
    ax_s.set_xticklabels(states)
    ax_s.set_ylabel("population fraction")
    ax_s.set_ylim(0, max(max(boltz), max(emp)) * 1.20)
    ax_s.set_title("State populations (Boltzmann vs empirical) — "
                   "EO* is CV0 ≥ 85 Å only",
                   fontsize=10.5)
    ax_s.legend(fontsize=8)
    ax_s.grid(axis="y", alpha=0.3)

    # (c) V1 vs V2 CV0 histogram + EO Jeffreys bound
    ax_h = fig.add_subplot(gs[1, 1])
    bins = np.linspace(40, 100, 31).tolist()
    ax_h.hist(cv0_v1, bins=bins, alpha=0.55, color="#2ca02c",
              density=True, label=f"V1 (n={cv0_v1.size})")
    ax_h.hist(cv0_v2, bins=bins, alpha=0.55, color="#d62728",
              density=True, label=f"V2 (n={cv0_v2.size})")
    ax_h.axvline(EO_MIN, color="black", ls="--", lw=1,
                 label=f"EO threshold (85 Å)")
    ax_h.set_xlabel("CV0 (Å)")
    ax_h.set_ylabel("density")
    ax_h.set_title(f"EO gap: {n_eo}/{n_total} fitted at ≥85 Å "
                   f"({100*p_eo_emp:.2f}%); ΔG ≥ {-KT_KCAL*np.log(p_eo_upper):.2f} kcal/mol",
                   fontsize=10.5)
    ax_h.legend(fontsize=8)
    ax_h.grid(alpha=0.3)

    fig.suptitle("State populations + V1/V2 ΔG reproducibility (obj-043)",
                 fontsize=12, fontweight="bold", y=0.995)

    out_png = FIG_DIR / "state_populations_v1.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"Saved {out_png}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
