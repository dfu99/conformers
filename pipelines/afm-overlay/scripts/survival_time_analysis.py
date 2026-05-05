#!/usr/bin/env python3
"""Survival-time analysis per HMM state (obj-049).

Audit-2026-05-05 §17.6 v9 P=3 follow-up to obj-048. obj-048 gave per-
state mean dwell times under a normal-CI assumption. This v11 deepening
fits proper exponential distributions per state, runs a KS goodness-of-
fit test against the single-exponential null, and compares against
Weibull / gamma alternatives via AIC. Distinguishes:

  - Markovian / memoryless dynamics (exponential fits well)  vs
  - History-dependent dynamics (Weibull k≠1 or gamma α≠1)

Right-censoring is honored: dwell-time runs that begin at frame 0 or
end at the last frame of a video are flagged as censored and folded
into the Kaplan-Meier survival estimator. MLE exponential fits and KS
tests use *uncensored* runs only (the standard treatment for
right-censored data with a small tail; full-MLE with censoring would
require minimization over an indicator term).

Outputs:
  figures/survival_time_analysis.png
  results/afm_pipeline/free_energy_profile/survival_times.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
HMM_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "hmm_states.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

STATE_LABELS = ["BC", "Intermediate", "EC"]
STATE_COLORS = ["#d62728", "#ff7f0e", "#2ca02c"]
N_STATES = 3


def runs_with_censor(states: np.ndarray, K: int) -> dict:
    """Enumerate runs per state, marking start/end-of-video runs as censored."""
    out: dict = {k: [] for k in range(K)}
    if states.size == 0:
        return out
    cur = int(states[0])
    n = 1
    for t in range(1, states.size):
        s = int(states[t])
        if s == cur:
            n += 1
        else:
            out[cur].append({"length": n, "censored": False, "start_idx": t - n})
            cur = s
            n = 1
    out[cur].append({"length": n, "censored": True, "start_idx": states.size - n})
    # Mark the very first run as censored (left-censored on the boundary)
    if out[int(states[0])]:
        out[int(states[0])][0]["censored"] = True
    return out


def merge_videos(*runs_by_state: dict) -> dict:
    out: dict = {k: [] for k in range(N_STATES)}
    for r in runs_by_state:
        for k, lst in r.items():
            out[k].extend(lst)
    return out


def kaplan_meier(durations: np.ndarray, censored: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """KM survival estimate. durations and censored are 1D arrays of equal length."""
    if durations.size == 0:
        return np.array([0.0]), np.array([1.0])
    order = np.argsort(durations)
    durations = durations[order]
    censored = censored[order]
    n_at_risk = durations.size
    times = [0.0]
    survival = [1.0]
    s = 1.0
    unique = np.unique(durations)
    for t in unique:
        deaths = int(np.sum((durations == t) & (~censored)))
        if deaths > 0:
            s *= (1 - deaths / n_at_risk)
        n_at_risk -= int(np.sum(durations == t))
        times.append(float(t))
        survival.append(float(s))
    return np.array(times), np.array(survival)


def fit_and_test(durations_uncensored: np.ndarray) -> dict:
    """Fit exponential, Weibull, gamma; KS test exponential null; AIC comparison."""
    if durations_uncensored.size < 3:
        return {"insufficient_data": True, "n": int(durations_uncensored.size)}
    x = durations_uncensored.astype(float)

    # Exponential MLE: λ = 1 / mean (loc=0 fixed)
    loc_e, scale_e = stats.expon.fit(x, floc=0)
    rate_e = 1.0 / scale_e
    log_lik_e = np.sum(stats.expon.logpdf(x, loc=loc_e, scale=scale_e))
    aic_e = 2 * 1 - 2 * log_lik_e  # 1 free param (scale)
    ks_stat_e, ks_p_e = stats.kstest(x, "expon", args=(loc_e, scale_e))

    # Weibull MLE (loc=0)
    shape_w, _, scale_w = stats.weibull_min.fit(x, floc=0)
    log_lik_w = np.sum(stats.weibull_min.logpdf(x, shape_w, loc=0, scale=scale_w))
    aic_w = 2 * 2 - 2 * log_lik_w

    # Gamma MLE (loc=0)
    shape_g, _, scale_g = stats.gamma.fit(x, floc=0)
    log_lik_g = np.sum(stats.gamma.logpdf(x, shape_g, loc=0, scale=scale_g))
    aic_g = 2 * 2 - 2 * log_lik_g

    # Bootstrap 95 % CI on exponential mean
    rng = np.random.default_rng(42)
    boot_means = np.array([rng.choice(x, size=x.size, replace=True).mean()
                           for _ in range(2000)])
    ci_low, ci_high = np.percentile(boot_means, [2.5, 97.5])

    return {
        "n": int(x.size),
        "mean": float(x.mean()),
        "exp_rate_per_frame": float(rate_e),
        "exp_scale_frames": float(scale_e),
        "exp_log_lik": float(log_lik_e),
        "exp_aic": float(aic_e),
        "ks_stat_vs_exp": float(ks_stat_e),
        "ks_pvalue_vs_exp": float(ks_p_e),
        "weibull_shape": float(shape_w),
        "weibull_scale": float(scale_w),
        "weibull_aic": float(aic_w),
        "gamma_shape": float(shape_g),
        "gamma_scale": float(scale_g),
        "gamma_aic": float(aic_g),
        "best_aic_model": min(
            [("exp", aic_e), ("weibull", aic_w), ("gamma", aic_g)],
            key=lambda p: p[1])[0],
        "exp_mean_bootstrap_95ci_frames": [float(ci_low), float(ci_high)],
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    hmm = np.load(HMM_NPZ)
    states_v1 = hmm["viterbi_v1"]
    states_v2 = hmm["viterbi_v2"]
    print(f"Loaded V1 ({states_v1.size} frames) and V2 ({states_v2.size} frames) Viterbi paths.")

    runs_v1 = runs_with_censor(states_v1, N_STATES)
    runs_v2 = runs_with_censor(states_v2, N_STATES)
    runs_all = merge_videos(runs_v1, runs_v2)

    fits: dict = {}
    km: dict = {}
    print()
    for k, lab in enumerate(STATE_LABELS):
        runs = runs_all[k]
        if not runs:
            print(f"  {lab}: no runs")
            fits[lab] = {"insufficient_data": True}
            km[lab] = ([0.0], [1.0])
            continue
        lengths = np.array([r["length"] for r in runs], dtype=float)
        censored = np.array([r["censored"] for r in runs], dtype=bool)
        n_uncens = int(np.sum(~censored))
        n_cens = int(np.sum(censored))
        x_uncens = lengths[~censored]

        fit = fit_and_test(x_uncens)
        fits[lab] = fit
        ts, ss = kaplan_meier(lengths, censored)
        km[lab] = (ts.tolist(), ss.tolist())

        print(f"  {lab} (n={lengths.size}, uncensored={n_uncens}, censored={n_cens}):")
        if fit.get("insufficient_data"):
            print("    insufficient data for distribution fitting")
            continue
        print(f"    mean = {fit['mean']:.2f} s  "
              f"(exp scale = {fit['exp_scale_frames']:.2f} s, "
              f"95 % CI [{fit['exp_mean_bootstrap_95ci_frames'][0]:.2f}, "
              f"{fit['exp_mean_bootstrap_95ci_frames'][1]:.2f}])")
        print(f"    KS vs exponential: D = {fit['ks_stat_vs_exp']:.4f}, "
              f"p = {fit['ks_pvalue_vs_exp']:.4g}  "
              f"({'memoryless' if fit['ks_pvalue_vs_exp'] > 0.05 else 'NOT memoryless'})")
        print(f"    Weibull shape k = {fit['weibull_shape']:.3f}  "
              f"({'k ≈ 1: exponential-like' if abs(fit['weibull_shape'] - 1) < 0.2 else 'k ≠ 1'})")
        print(f"    AIC: exp={fit['exp_aic']:.2f}, weibull={fit['weibull_aic']:.2f}, "
              f"gamma={fit['gamma_aic']:.2f}  →  best = {fit['best_aic_model']}")

    # Save JSON
    json_path = OUT_DIR / "survival_times.json"
    with json_path.open("w") as fh:
        json.dump({
            "T_kelvin": 300.0,
            "frame_to_second": 1.0,  # 1 fps
            "n_states": N_STATES,
            "state_labels": STATE_LABELS,
            "fits": fits,
            "kaplan_meier_curves": km,
            "n_runs_per_state": {STATE_LABELS[k]: len(runs_all[k]) for k in range(N_STATES)},
        }, fh, indent=2)
    print(f"\nJSON saved: {json_path}")

    # --- Figure ---
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))

    # Top row: Kaplan-Meier survival curves with exponential overlay
    for k, lab in enumerate(STATE_LABELS):
        ax = axes[0, k]
        ts, ss = km[lab]
        ax.step(ts, ss, where="post", color=STATE_COLORS[k], lw=2,
                label="Kaplan-Meier")
        if not fits[lab].get("insufficient_data"):
            scale = fits[lab]["exp_scale_frames"]
            t_grid = np.linspace(0, max(ts) + 1, 200)
            ax.plot(t_grid, np.exp(-t_grid / scale), "k--", lw=1.5,
                    label=f"exp fit μ={scale:.1f} s")
        ax.set_title(f"{lab} — survival curve  (n={fits[lab].get('n', 0)} uncens.)")
        ax.set_xlabel("dwell time (s @ 1 fps)")
        ax.set_ylabel("S(t) = P(dwell > t)")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Bottom row: histogram with all 3 distribution fits + KS p-value annotation
    for k, lab in enumerate(STATE_LABELS):
        ax = axes[1, k]
        runs = runs_all[k]
        if not runs or fits[lab].get("insufficient_data"):
            ax.text(0.5, 0.5, "insufficient data", ha="center", transform=ax.transAxes)
            ax.axis("off")
            continue
        lengths = np.array([r["length"] for r in runs])
        cens = np.array([r["censored"] for r in runs])
        x_un = lengths[~cens].astype(float)

        bins = np.linspace(1, max(lengths) + 1, min(20, max(lengths) + 1))
        ax.hist(x_un, bins=bins, alpha=0.55, color=STATE_COLORS[k], edgecolor="black",
                density=True, label=f"empirical (uncens. n={x_un.size})")

        f = fits[lab]
        t_grid = np.linspace(1, max(lengths), 300)
        ax.plot(t_grid,
                stats.expon.pdf(t_grid, loc=0, scale=f["exp_scale_frames"]),
                "k--", lw=1.5, label="exponential")
        ax.plot(t_grid,
                stats.weibull_min.pdf(t_grid, f["weibull_shape"], loc=0, scale=f["weibull_scale"]),
                "b-.", lw=1.2, label=f"Weibull (k={f['weibull_shape']:.2f})")
        ax.plot(t_grid,
                stats.gamma.pdf(t_grid, f["gamma_shape"], loc=0, scale=f["gamma_scale"]),
                "g:", lw=1.2, label=f"gamma (α={f['gamma_shape']:.2f})")

        ks_p = f["ks_pvalue_vs_exp"]
        verdict = "memoryless" if ks_p > 0.05 else "NOT memoryless"
        ax.set_title(f"{lab} — KS vs exp p={ks_p:.3g} ({verdict});  "
                     f"best AIC: {f['best_aic_model']}", fontsize=9)
        ax.set_xlabel("dwell time (s @ 1 fps)")
        ax.set_ylabel("density")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    fig.suptitle(
        "obj-049 — Per-state survival-time analysis (V1+V2 HMM Viterbi)\n"
        "Top: Kaplan-Meier survival w/ right-censoring  •  "
        "Bottom: pdf comparison with exponential / Weibull / gamma fits + KS test",
        fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig_path = FIG_DIR / "survival_time_analysis.png"
    fig.savefig(fig_path, dpi=140)
    print(f"Figure saved: {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
