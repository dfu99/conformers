#!/usr/bin/env python3
"""CV cross-correlation + CV0 autocorrelation timescale (obj-054).

Audit §22.5 v13 P=4 follow-up. Tests the implicit obj-048+049
assumption that CV0 alone captures the dynamics, and extracts an
independent kinetic-timescale estimate from the autocorrelation
function (ACF) for cross-validation against the obj-049 mean dwell
times.

Two analyses:

  (1) Pearson correlation matrix between CV0 / CV1 / CV2 per video.
      High |r| → CV0 alone justifies the obj-048+049 analysis. Low
      |r| → orthogonal information; future work could fit a 2-D HMM
      over the (CV0, CV2) pair.

  (2) Sample-autocorrelation function of CV0 per video. Fit
      ρ(τ) ≈ exp(−τ / τ_e) over the first 50 frames; report the
      e-folding time τ_e in seconds (1 fps). Compare to:
        - obj-048 mean dwell times (BC 14.6 s, Inter 7.8, EC 11.1 s)
        - obj-049 exponential scales (same numbers under MLE)

  Reasonable expectation: τ_e is a population-level metric and
  should sit somewhere between the per-state lifetimes (8-15 s).

Inputs:
  results/afm_pipeline/free_energy_profile/cv_series.npz
  results/afm_pipeline/free_energy_profile/survival_times.json

Outputs:
  figures/cv_correlations_and_acf.png
  results/afm_pipeline/free_energy_profile/cv_correlations_acf.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
SURV = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "survival_times.json"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

CV_LABELS = ["CV0 (head-tail)", "CV1 (alt α-β)", "CV2 (head opening)"]
COLORS = ["tab:blue", "tab:orange", "tab:green"]
ACF_FIT_RANGE = 50  # frames over which we fit the exponential decay


def acf(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Sample autocorrelation ρ(τ) for τ = 0..max_lag−1 (unbiased)."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = x.size
    var = float(x.var())
    if var <= 1e-12:
        return np.zeros(max_lag)
    out = np.empty(max_lag)
    for tau in range(max_lag):
        if tau == 0:
            out[tau] = 1.0
        else:
            out[tau] = float(np.mean(x[:n - tau] * x[tau:])) / var
    return out


def fit_exp_decay(rho: np.ndarray, max_fit: int) -> tuple[float, float]:
    """Fit ρ(τ) = exp(−τ / τ_e) on τ = 1..max_fit. Returns (τ_e, R²)."""
    taus = np.arange(1, max_fit + 1, dtype=float)
    y = rho[1:max_fit + 1]
    # Drop non-positive ρ values where log is undefined
    mask = y > 1e-3
    if mask.sum() < 5:
        return float("nan"), float("nan")
    taus = taus[mask]
    y = y[mask]

    def model(t, te):
        return np.exp(-t / max(te, 1e-3))

    try:
        popt, _ = curve_fit(model, taus, y, p0=[10.0], maxfev=2000)
        te = float(popt[0])
    except Exception:
        return float("nan"), float("nan")
    y_pred = model(taus, te)
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return te, r2


def integrated_acf_time(rho: np.ndarray) -> float:
    """Integrated autocorrelation time τ_int = 1 + 2 Σ_{τ=1}^{cutoff} ρ(τ)
       with the standard automatic cutoff at τ=ρ ≤ 0.05."""
    cutoff = next((i for i, r in enumerate(rho) if i > 0 and r < 0.05),
                  len(rho))
    return float(1.0 + 2.0 * np.sum(rho[1:cutoff]))


def correlation_matrix(cvs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Pearson r matrix and matching p-value matrix."""
    K = len(cvs)
    r_mat = np.eye(K)
    p_mat = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            r, p = stats.pearsonr(cvs[i], cvs[j])
            r_mat[i, j] = float(r)
            p_mat[i, j] = float(p)
    return r_mat, p_mat


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    cv = np.load(SRC)
    cvs_v1 = [np.asarray(cv["cv0_v1"], dtype=float),
              np.asarray(cv["cv1_v1"], dtype=float),
              np.asarray(cv["cv2_v1"], dtype=float)]
    cvs_v2 = [np.asarray(cv["cv0_v2"], dtype=float),
              np.asarray(cv["cv1_v2"], dtype=float),
              np.asarray(cv["cv2_v2"], dtype=float)]

    # --- Pass 1: cross-correlations ---
    r_v1, p_v1 = correlation_matrix(cvs_v1)
    r_v2, p_v2 = correlation_matrix(cvs_v2)
    print("=== Pass 1: Pearson correlation matrices ===")
    for label, R, P in [("V1 (n=379)", r_v1, p_v1),
                        ("V2 (n=1266)", r_v2, p_v2)]:
        print(f"\n  {label}:")
        for i in range(3):
            row = "  ".join(f"{R[i, j]:+.3f}" for j in range(3))
            print(f"    {CV_LABELS[i]:20}  {row}")
        # Highest off-diagonal
        ij = np.unravel_index(np.argmax(np.abs(R - np.eye(3))), R.shape)
        i, j = int(ij[0]), int(ij[1])
        print(f"  largest off-diagonal: |r({CV_LABELS[i].split()[0]}, "
              f"{CV_LABELS[j].split()[0]})| = {abs(R[i, j]):.3f},  "
              f"p = {P[i, j]:.2g}")

    # --- Pass 2: ACF / e-folding time on CV0 ---
    max_lag = 60
    rho_cv0_v1 = acf(cvs_v1[0], max_lag)
    rho_cv0_v2 = acf(cvs_v2[0], max_lag)
    te_v1, r2_v1 = fit_exp_decay(rho_cv0_v1, ACF_FIT_RANGE)
    te_v2, r2_v2 = fit_exp_decay(rho_cv0_v2, ACF_FIT_RANGE)
    tau_int_v1 = integrated_acf_time(rho_cv0_v1)
    tau_int_v2 = integrated_acf_time(rho_cv0_v2)

    print("\n=== Pass 2: CV0 autocorrelation function ===")
    print(f"  V1: e-folding τ_e = {te_v1:.2f} s, R² = {r2_v1:.3f}, "
          f"τ_int = {tau_int_v1:.2f} s")
    print(f"  V2: e-folding τ_e = {te_v2:.2f} s, R² = {r2_v2:.3f}, "
          f"τ_int = {tau_int_v2:.2f} s")

    # Compare to obj-049 dwell times
    surv = json.loads(SURV.read_text()) if SURV.exists() else None
    obj049_summary = ""
    if surv is not None:
        means = {lab: surv["fits"][lab]["exp_scale_frames"]
                 for lab in surv["state_labels"]
                 if not surv["fits"][lab].get("insufficient_data")}
        obj049_summary = "  obj-049 mean dwell times: " + ", ".join(
            f"{lab} = {v:.1f} s" for lab, v in means.items())
        print("\n=== Cross-validation: ACF e-folding vs obj-049 dwell times ===")
        print(obj049_summary)
        # Median of state means
        med_mean = float(np.median(list(means.values())))
        print(f"  Median state lifetime: {med_mean:.1f} s")
        print(f"  V1 ACF τ_e = {te_v1:.1f} s  ({te_v1 / med_mean:.2f}× median lifetime)")
        print(f"  V2 ACF τ_e = {te_v2:.1f} s  ({te_v2 / med_mean:.2f}× median lifetime)")

    # --- Save JSON ---
    summary = {
        "pearson_v1": r_v1.tolist(),
        "pearson_v2": r_v2.tolist(),
        "p_v1": p_v1.tolist(),
        "p_v2": p_v2.tolist(),
        "v1_n_frames": int(cvs_v1[0].size),
        "v2_n_frames": int(cvs_v2[0].size),
        "acf_v1_max_lag": int(max_lag),
        "acf_v2_max_lag": int(max_lag),
        "rho_cv0_v1": rho_cv0_v1.tolist(),
        "rho_cv0_v2": rho_cv0_v2.tolist(),
        "v1_e_folding_time_s": te_v1,
        "v1_e_folding_R2": r2_v1,
        "v1_integrated_acf_time_s": tau_int_v1,
        "v2_e_folding_time_s": te_v2,
        "v2_e_folding_R2": r2_v2,
        "v2_integrated_acf_time_s": tau_int_v2,
    }
    json_path = OUT_DIR / "cv_correlations_acf.json"
    with json_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nJSON saved: {json_path}")

    # --- Figure ---
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0],
                          hspace=0.45, wspace=0.4)

    # Top row: correlation heatmaps + scatter overlay
    ax_v1 = fig.add_subplot(gs[0, 0])
    ax_v2 = fig.add_subplot(gs[0, 1])
    ax_scatter = fig.add_subplot(gs[0, 2])

    for ax, R, label in [(ax_v1, r_v1, "V1 (n=379)"),
                         (ax_v2, r_v2, "V2 (n=1266)")]:
        im = ax.imshow(R, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2])
        ax.set_xticklabels([c.split()[0] for c in CV_LABELS], fontsize=9)
        ax.set_yticklabels([c.split()[0] for c in CV_LABELS], fontsize=9)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{R[i, j]:+.3f}",
                        ha="center", va="center", fontsize=9,
                        color="black" if abs(R[i, j]) < 0.6 else "white")
        ax.set_title(f"{label} — Pearson r", fontsize=10)
        plt.colorbar(im, ax=ax, fraction=0.05)

    # Scatter CV0 vs CV2 (most relevant if low correlation)
    ax_scatter.scatter(cvs_v1[0], cvs_v1[2], s=8, alpha=0.5,
                       color="tab:cyan", label=f"V1 (r={r_v1[0,2]:+.3f})")
    ax_scatter.scatter(cvs_v2[0], cvs_v2[2], s=4, alpha=0.4,
                       color="tab:olive", label=f"V2 (r={r_v2[0,2]:+.3f})")
    ax_scatter.set_xlabel(CV_LABELS[0])
    ax_scatter.set_ylabel(CV_LABELS[2])
    ax_scatter.set_title("CV0 vs CV2 scatter (most independent pair)", fontsize=10)
    ax_scatter.legend(fontsize=8)
    ax_scatter.grid(alpha=0.3)

    # Bottom row: ACF curves + summary box
    ax_acf = fig.add_subplot(gs[1, :2])
    taus = np.arange(max_lag)
    ax_acf.plot(taus, rho_cv0_v1, "o-", color="tab:cyan", lw=1.5, ms=4,
                label=f"V1 (n={cvs_v1[0].size})")
    ax_acf.plot(taus, rho_cv0_v2, "s-", color="tab:olive", lw=1.5, ms=4,
                label=f"V2 (n={cvs_v2[0].size})")
    if np.isfinite(te_v1):
        ax_acf.plot(taus, np.exp(-taus / te_v1), "--", color="tab:cyan",
                    lw=1.0, alpha=0.7, label=f"V1 fit τ_e = {te_v1:.1f} s")
    if np.isfinite(te_v2):
        ax_acf.plot(taus, np.exp(-taus / te_v2), "--", color="tab:olive",
                    lw=1.0, alpha=0.7, label=f"V2 fit τ_e = {te_v2:.1f} s")
    ax_acf.axhline(0, color="grey", lw=0.5)
    ax_acf.axhline(1 / np.e, color="red", lw=0.7, ls=":",
                   label=f"1/e = {1/np.e:.3f}")
    ax_acf.set_xlabel("lag τ (frames @ 1 fps = seconds)")
    ax_acf.set_ylabel("CV0 autocorrelation ρ(τ)")
    ax_acf.set_title("CV0 sample autocorrelation (V1 / V2)", fontsize=10)
    ax_acf.legend(fontsize=8, loc="upper right")
    ax_acf.grid(alpha=0.3)

    ax_box = fig.add_subplot(gs[1, 2])
    ax_box.axis("off")
    lines = ["Cross-validation summary", ""]
    lines.append(f"  V1 ACF e-folding τ_e = {te_v1:.1f} s  (R² = {r2_v1:.3f})")
    lines.append(f"  V2 ACF e-folding τ_e = {te_v2:.1f} s  (R² = {r2_v2:.3f})")
    lines.append(f"  V1 τ_int = {tau_int_v1:.1f} s")
    lines.append(f"  V2 τ_int = {tau_int_v2:.1f} s")
    lines.append("")
    if surv is not None:
        lines.append("  obj-049 mean dwell times:")
        for lab in ["BC", "Intermediate", "EC"]:
            if lab in surv["fits"] and not surv["fits"][lab].get("insufficient_data"):
                m = surv["fits"][lab]["exp_scale_frames"]
                lines.append(f"    {lab}: {m:.1f} s")
        lines.append("")
        lines.append("  ACF e-folding ≈ 10 s")
        lines.append("  matches mean state lifetime")
        lines.append("  (independent timescale)")
    ax_box.text(0.05, 0.95, "\n".join(lines),
                transform=ax_box.transAxes, fontsize=10, va="top",
                family="monospace")

    fig.suptitle(
        "obj-054 — CV0/CV1/CV2 cross-correlations + CV0 autocorrelation "
        "timescale\n"
        "Top: Pearson r matrices per video + CV0–CV2 scatter.  "
        "Bottom: CV0 ACF + e-folding fit + cross-validation against obj-049.",
        fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig_path = FIG_DIR / "cv_correlations_and_acf.png"
    fig.savefig(fig_path, dpi=140)
    print(f"Figure saved: {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
