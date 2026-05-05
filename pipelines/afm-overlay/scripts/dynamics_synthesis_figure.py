#!/usr/bin/env python3
"""Synthesis publication figure (obj-053).

Audit-2026-05-05 §20.5 v12 P=3 follow-up. Collapses obj-038 / 043 /
047 / 048 / 049 / 050 / 051 into a single 4-panel manuscript figure:

  (a) ΔG(CV0) curve with 95 % bootstrap confidence band (obj-038 + F3)
      and state-band shading (BC / Intermediate / EC / EO).
  (b) State populations: Boltzmann (from ΔG) + empirical (count) +
      HMM stationary π — three side-by-side bars per state — with
      the Bayesian EO floor ΔG ≥ 2.02 kcal/mol annotated (obj-043).
  (c) HMM Viterbi state path on the V2 per-frame CV0 trace (1266
      frames; the longer of the two videos), with obj-047 BinSeg+BIC
      change-points overlaid as red dashed lines and an inset
      transition-matrix heatmap (obj-048).
  (d) Per-state Kaplan-Meier survival curves with single-exponential
      overlays + KS p-value annotation (obj-049).

Inputs (all already on disk from previous objectives):
  - results/afm_pipeline/free_energy_profile/cv_series.npz
  - results/afm_pipeline/free_energy_profile/bootstrap.npz (F3)
  - results/afm_pipeline/free_energy_profile/state_populations.json
  - results/afm_pipeline/free_energy_profile/hmm_states.{json,npz}
  - results/afm_pipeline/free_energy_profile/change_points.json
  - results/afm_pipeline/free_energy_profile/survival_times.json

Outputs:
  figures/dynamics_synthesis_v1.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

ROOT = Path(__file__).resolve().parents[3]
FE = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

T_KELVIN = 300.0
KT_KCAL = 0.001987 * T_KELVIN

STATE_LABELS = ["BC", "Intermediate", "EC"]
STATE_COLORS = ["#d62728", "#ff7f0e", "#2ca02c"]
EO_COLOR = "#9467bd"

BC_MAX = 65.0
EC_MIN = 78.0
EO_MIN = 85.0


def panel_a_dg(ax) -> None:
    """ΔG(CV0) with bootstrap band + state-band shading."""
    boot = np.load(FE / "bootstrap.npz")
    grid = boot["grid"]
    g_point = boot["g_point"]
    g_lo = boot["g_lo"]
    g_hi = boot["g_hi"]

    # Library support range
    lib_min, lib_max = 47.3, 85.0

    # Shade the four state bands
    ax.axvspan(grid.min(), BC_MAX, alpha=0.18, color=STATE_COLORS[0],
               label=f"BC ≤ {BC_MAX:.0f} Å")
    ax.axvspan(BC_MAX, EC_MIN, alpha=0.18, color=STATE_COLORS[1],
               label=f"Intermediate {BC_MAX:.0f}–{EC_MIN:.0f} Å")
    ax.axvspan(EC_MIN, EO_MIN, alpha=0.18, color=STATE_COLORS[2],
               label=f"EC {EC_MIN:.0f}–{EO_MIN:.0f} Å")
    ax.axvspan(EO_MIN, grid.max(), alpha=0.18, color=EO_COLOR,
               label=f"EO ≥ {EO_MIN:.0f} Å (Bayesian floor)")

    # Curve + bootstrap band
    mask = (grid >= lib_min) & (grid <= lib_max)
    ax.fill_between(grid[mask], g_lo[mask], g_hi[mask], color="grey",
                    alpha=0.4, label="95 % bootstrap CI")
    ax.plot(grid[mask], g_point[mask], color="black", lw=2,
            label="ΔG(CV0) point estimate")

    # EO Bayesian floor annotation
    with (FE / "state_populations.json").open() as fh:
        pop = json.load(fh)
    eo_floor = pop["eo_gap"]["delta_g_eo_lower_bound_kcal"]
    ax.axhline(eo_floor, color=EO_COLOR, ls="--", lw=1.3,
               label=f"Bayesian EO floor ≥ {eo_floor:.2f} kcal/mol")

    ax.set_xlim(45, 95)
    ax.set_ylim(-0.2, 6.5)
    ax.set_xlabel("CV0 (head–tail distance, Å)")
    ax.set_ylabel("ΔG (kcal/mol)")
    ax.set_title("(a) ΔG(CV0) profile with 95 % bootstrap CI + state bands "
                 "(obj-038 + F3 + obj-043)", fontsize=10.5)
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.92)
    ax.grid(alpha=0.3)


def panel_b_populations(ax) -> None:
    """Boltzmann + empirical + HMM stationary populations side by side."""
    with (FE / "state_populations.json").open() as fh:
        pop = json.load(fh)
    with (FE / "hmm_states.json").open() as fh:
        hmm = json.load(fh)

    boltz = pop["boltzmann_populations"]
    emp = pop["empirical_populations"]
    pi = hmm["stationary_distribution"]

    states = ["BC", "Intermediate", "EC", "EO"]
    bvals = [boltz["BC"], boltz["Intermediate"], boltz["EC"], boltz["EO_CV0_only"]]
    evals = [emp["BC"], emp["Intermediate"], emp["EC"], emp["EO_CV0_only"]]
    pvals = [pi[0], pi[1], pi[2], 0.0]

    x = np.arange(len(states))
    w = 0.27
    ax.bar(x - w, bvals, w, color="tab:blue", alpha=0.9,
           label="Boltzmann (from ΔG)")
    ax.bar(x, evals, w, color="tab:gray", alpha=0.9,
           label="Empirical (count)")
    ax.bar(x + w, pvals, w, color="tab:red", alpha=0.9,
           label="HMM stationary π")

    for xi, b in zip(x - w, bvals):
        ax.text(xi, b + 0.005, f"{b:.0%}", ha="center", fontsize=7.5)
    for xi, e in zip(x, evals):
        ax.text(xi, e + 0.005, f"{e:.0%}", ha="center", fontsize=7.5)
    for xi, p in zip(x + w, pvals):
        if p > 0:
            ax.text(xi, p + 0.005, f"{p:.0%}", ha="center", fontsize=7.5)

    eo_floor = pop["eo_gap"]["delta_g_eo_lower_bound_kcal"]
    ax.text(0.98, 0.95,
            f"Bayesian EO floor ΔG ≥ {eo_floor:.2f} kcal/mol\n"
            f"P_EO ≤ {pop['eo_gap']['p_eo_jeffreys_upper95']*100:.2f} % "
            f"(Jeffreys 95 % upper)\n"
            f"42 / 1645 frames at CV0 ≥ 85 Å",
            ha="right", va="top", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="grey"))

    ax.set_xticks(x)
    ax.set_xticklabels(states, fontsize=9)
    ax.set_ylabel("population fraction")
    ax.set_ylim(0, 0.55)
    ax.set_title("(b) State populations: Boltzmann + empirical + HMM "
                 "(obj-043 + obj-048)", fontsize=10.5)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", alpha=0.3)


def panel_c_hmm(ax_path, ax_inset) -> None:
    """HMM Viterbi state path + transition-matrix inset."""
    cv = np.load(FE / "cv_series.npz")
    cv0_v2 = np.asarray(cv["cv0_v2"], dtype=float)
    hmm_npz = np.load(FE / "hmm_states.npz")
    states_v2 = hmm_npz["viterbi_v2"]
    A = hmm_npz["transition_matrix"]

    # change-points
    with (FE / "change_points.json").open() as fh:
        cp = json.load(fh)
    cps_v2 = cp["results"]["v2_frame_cv0_smoothed"]["binseg_change_points_idx"]

    T = cv0_v2.size
    for k in range(3):
        mask = states_v2 == k
        ax_path.scatter(np.arange(T)[mask], cv0_v2[mask],
                        s=6, c=STATE_COLORS[k], alpha=0.7,
                        label=STATE_LABELS[k])
    for cpx in cps_v2:
        ax_path.axvline(cpx, color="red", ls="--", lw=1.3, alpha=0.85)

    # State band horizontal lines
    ax_path.axhline(BC_MAX, color="grey", lw=0.5, ls=":")
    ax_path.axhline(EC_MIN, color="grey", lw=0.5, ls=":")
    ax_path.axhline(EO_MIN, color="grey", lw=0.5, ls=":")

    ax_path.set_xlim(0, T - 1)
    ax_path.set_ylim(45, 95)
    ax_path.set_xlabel("frame (1 fps; ~ seconds)")
    ax_path.set_ylabel("CV0 (Å)")
    ax_path.set_title(
        f"(c) V2 HMM Viterbi state path + obj-047 change-points "
        f"({len(cps_v2)} red dashed)  —  obj-047 + obj-048",
        fontsize=10.5)
    ax_path.legend(loc="upper right", fontsize=7.5, ncol=3, framealpha=0.92)
    ax_path.grid(alpha=0.3)

    # Transition matrix inset
    im = ax_inset.imshow(A, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax_inset.set_xticks([0, 1, 2]); ax_inset.set_yticks([0, 1, 2])
    ax_inset.set_xticklabels(STATE_LABELS, fontsize=7)
    ax_inset.set_yticklabels(STATE_LABELS, fontsize=7)
    ax_inset.set_title("transition matrix A[from, to]", fontsize=8.5)
    for i in range(3):
        for j in range(3):
            ax_inset.text(j, i, f"{A[i, j]:.3f}", ha="center", va="center",
                          fontsize=7,
                          color="black" if A[i, j] < 0.5 else "white")
    cbar = plt.colorbar(im, ax=ax_inset, fraction=0.05)
    cbar.ax.tick_params(labelsize=6)


def panel_d_survival(ax) -> None:
    """Per-state Kaplan-Meier + exponential overlays + KS p-values."""
    with (FE / "survival_times.json").open() as fh:
        surv = json.load(fh)

    annotations = []
    for k, lab in enumerate(STATE_LABELS):
        if lab not in surv["fits"]:
            continue
        fit = surv["fits"][lab]
        if fit.get("insufficient_data"):
            continue
        ts, ss = surv["kaplan_meier_curves"][lab]
        ts = np.asarray(ts); ss = np.asarray(ss)
        ax.step(ts, ss, where="post", color=STATE_COLORS[k], lw=2,
                label=f"{lab} (n={fit['n']})")
        scale = fit["exp_scale_frames"]
        t_grid = np.linspace(0, ts.max() + 5, 200)
        ax.plot(t_grid, np.exp(-t_grid / scale),
                color=STATE_COLORS[k], lw=1.0, ls="--", alpha=0.7)
        annotations.append(
            f"{lab}:  μ = {scale:.1f} s, KS p = {fit['ks_pvalue_vs_exp']:.2g}"
        )

    ax.text(0.98, 0.95,
            "\n".join(annotations) + "\n\nAll three states are\n"
            "consistent with memoryless\n(Markovian) kinetics  (KS p > 0.13).",
            ha="right", va="top", transform=ax.transAxes, fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="grey"))

    ax.set_xlabel("dwell time (s @ 1 fps)")
    ax.set_ylabel("S(t) = P(dwell > t)")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, 95)
    ax.set_title("(d) Per-state survival curves + exponential fits  —  "
                 "obj-049 (V1+V2)", fontsize=10.5)
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.3)


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(15, 11.5))
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.27)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a_dg(ax_a)
    panel_b_populations(ax_b)

    # Inset for the transition matrix on panel C
    ax_c_inset = ax_c.inset_axes((0.62, 0.62, 0.35, 0.35))
    panel_c_hmm(ax_c, ax_c_inset)

    panel_d_survival(ax_d)

    fig.suptitle(
        "Conformers αVβ3 dynamics synthesis — V1+V2 HS-AFM (1645 frames @ 1 fps)\n"
        "ΔG profile + state populations + HMM kinetics + Markovian dwell-times",
        fontsize=12.5, fontweight="bold", y=0.995,
    )

    fig_path = FIG_DIR / "dynamics_synthesis_v1.png"
    fig.savefig(fig_path, dpi=160, bbox_inches="tight")
    print(f"Figure saved: {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
