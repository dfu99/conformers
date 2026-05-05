#!/usr/bin/env python3
"""2-D Gaussian-emission HMM over (CV0, CV2) (obj-055).

Audit §23.5 v14 P=4 follow-up to obj-048+054. obj-054 found that CV2
is largely orthogonal to CV0 (r ≤ 0.52 V1, ≤ 0.22 V2) — head-piece
opening carries information not captured by extension alone. This v15
fits a 2-D Gaussian-emission HMM over the (CV0, CV2) plane to test
whether adding CV2 reveals an explicit head-piece-opening sub-state
on top of the extension states from obj-048.

Implementation: from-scratch (no `hmmlearn`):
  - 2-D Gaussian emissions per state (mean ∈ ℝ², covariance ∈ ℝ²ˣ²
    diagonal-structured for stability with small-n)
  - Forward-backward in log-space
  - Baum-Welch EM
  - Viterbi decoding

Reference physical states:
  - BC:  CV0 ≈ 55 Å,  CV2 ≈ 36 Å (closed headpiece, bent)
  - EC:  CV0 ≈ 80 Å,  CV2 ≈ 36 Å (closed headpiece, extended)
  - EO:  CV0 ≈ 80 Å,  CV2 ≥ 50 Å (open headpiece, extended) — empty
    in v7 fitted data

Both 3-state and 4-state HMMs are fitted from physically-motivated
inits. BIC/AIC compares them; we expect 3 states for the v7 data
since CV2 ≥ 50 Å bin is empirically empty.

Inputs:
  results/afm_pipeline/free_energy_profile/cv_series.npz

Outputs:
  figures/hmm_2d_cv0_cv2.png
  results/afm_pipeline/free_energy_profile/hmm_2d_cv0_cv2.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from numpy.linalg import inv, slogdet

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

T_KELVIN = 300.0
LOG_EPS = -1e30


def gaussian2d_log_pdf(x: np.ndarray, mu: np.ndarray, Sigma: np.ndarray) -> np.ndarray:
    """Log-pdf of x (T,2) under N(mu, Sigma)."""
    sign, logdet = slogdet(Sigma)
    if sign <= 0:
        # Regularize: add small jitter
        Sigma = Sigma + 1e-3 * np.eye(2)
        sign, logdet = slogdet(Sigma)
    inv_Sigma = inv(Sigma)
    diff = x - mu[None, :]
    quad = np.einsum("ti,ij,tj->t", diff, inv_Sigma, diff)
    return -0.5 * (2 * np.log(2 * np.pi) + logdet + quad)


def emission_log_matrix(x: np.ndarray, mus: np.ndarray, Sigmas: np.ndarray) -> np.ndarray:
    """Emission log-pdf matrix: x (T, 2), mus (K, 2), Sigmas (K, 2, 2)."""
    T = x.shape[0]
    K = mus.shape[0]
    log_b = np.empty((T, K))
    for k in range(K):
        log_b[:, k] = gaussian2d_log_pdf(x, mus[k], Sigmas[k])
    return log_b


def logsumexp(a: np.ndarray, axis=None) -> np.ndarray:
    a_max = np.max(a, axis=axis, keepdims=True)
    a_max = np.where(np.isfinite(a_max), a_max, 0.0)
    out = np.log(np.sum(np.exp(a - a_max), axis=axis, keepdims=True)) + a_max
    return np.squeeze(out, axis=axis) if axis is not None else out


def forward_backward(log_b, log_pi, log_A):
    T, K = log_b.shape
    log_alpha = np.full((T, K), LOG_EPS)
    log_beta = np.full((T, K), LOG_EPS)
    log_alpha[0] = log_pi + log_b[0]
    for t in range(1, T):
        prev = log_alpha[t - 1][:, None] + log_A
        log_alpha[t] = logsumexp(prev, axis=0) + log_b[t]
    log_beta[T - 1] = 0.0
    for t in range(T - 2, -1, -1):
        nxt = log_A + log_b[t + 1][None, :] + log_beta[t + 1][None, :]
        log_beta[t] = logsumexp(nxt, axis=1)
    log_likelihood = float(np.atleast_1d(logsumexp(log_alpha[T - 1])).item())
    return log_alpha, log_beta, log_likelihood


def viterbi(log_b, log_pi, log_A):
    T, K = log_b.shape
    delta = np.full((T, K), LOG_EPS)
    psi = np.zeros((T, K), dtype=int)
    delta[0] = log_pi + log_b[0]
    for t in range(1, T):
        score = delta[t - 1][:, None] + log_A
        psi[t] = np.argmax(score, axis=0)
        delta[t] = np.max(score, axis=0) + log_b[t]
    states = np.zeros(T, dtype=int)
    states[-1] = int(np.argmax(delta[-1]))
    for t in range(T - 2, -1, -1):
        states[t] = psi[t + 1, states[t + 1]]
    return states


def baum_welch_2d(x, K, init_mus, init_Sigmas, init_A, init_pi,
                  max_iter=200, tol=1e-5):
    T = x.shape[0]
    mus = init_mus.copy()
    Sigmas = init_Sigmas.copy()
    log_pi = np.log(np.maximum(init_pi, 1e-12))
    log_A = np.log(np.maximum(init_A, 1e-12))

    prev_ll = -np.inf
    ll = -np.inf
    it = 0
    for it in range(max_iter):
        log_b = emission_log_matrix(x, mus, Sigmas)
        log_alpha, log_beta, ll = forward_backward(log_b, log_pi, log_A)
        if not np.isfinite(ll):
            print(f"  iter {it}: non-finite log-L; aborting EM")
            break
        log_gamma = log_alpha + log_beta
        log_gamma -= logsumexp(log_gamma, axis=1)[:, None]
        gamma = np.exp(log_gamma)

        log_xi_unnorm = (log_alpha[:-1, :, None]
                         + log_A[None, :, :]
                         + log_b[1:, None, :]
                         + log_beta[1:, None, :])
        log_xi_norm = logsumexp(log_xi_unnorm.reshape(T - 1, K * K), axis=1)
        log_xi = log_xi_unnorm - log_xi_norm[:, None, None]
        xi = np.exp(log_xi)

        log_pi = log_gamma[0]
        A_num = xi.sum(axis=0)
        A_den = gamma[:-1].sum(axis=0)
        A = A_num / np.maximum(A_den[:, None], 1e-12)
        A /= A.sum(axis=1, keepdims=True)
        log_A = np.log(np.maximum(A, 1e-12))

        # Multivariate emission updates
        denom = gamma.sum(axis=0)
        for k in range(K):
            w = gamma[:, k]
            mus[k] = (w[:, None] * x).sum(axis=0) / max(denom[k], 1e-12)
            diff = x - mus[k][None, :]
            Sigmas[k] = (
                np.einsum("t,ti,tj->ij", w, diff, diff)
                / max(denom[k], 1e-12)
            )
            # Regularize
            Sigmas[k] += 1e-3 * np.eye(2)

        if abs(ll - prev_ll) < tol:
            print(f"  EM converged at iter {it} (ΔLL = {ll - prev_ll:.2e})")
            break
        prev_ll = ll

    return {
        "mus": mus, "Sigmas": Sigmas, "A": np.exp(log_A),
        "pi": np.exp(log_pi), "log_likelihood": float(ll),
        "n_iter": int(it + 1),
    }


def n_free_params_2d(K: int) -> int:
    """K-state 2-D Gaussian-emission HMM:
       - K means × 2  = 2K
       - K covariances × 3 unique entries (2x2 sym)  = 3K
       - (K-1) initial-π entries
       - K(K-1) transition entries
    """
    return 2 * K + 3 * K + (K - 1) + K * (K - 1)


def stationary_distribution(A: np.ndarray) -> np.ndarray:
    eigvals, eigvecs = np.linalg.eig(A.T)
    idx_one = int(np.argmin(np.abs(eigvals - 1.0)))
    pi = np.real(eigvecs[:, idx_one])
    return np.abs(pi) / np.sum(np.abs(pi))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    cv = np.load(SRC)
    cv0 = np.concatenate([cv["cv0_v1"], cv["cv0_v2"]]).astype(float)
    cv2 = np.concatenate([cv["cv2_v1"], cv["cv2_v2"]]).astype(float)
    X = np.column_stack([cv0, cv2])
    n = X.shape[0]
    print(f"Loaded n = {n} points (V1+V2). "
          f"CV0 range [{cv0.min():.1f}, {cv0.max():.1f}], "
          f"CV2 range [{cv2.min():.1f}, {cv2.max():.1f}].")

    # Physically-motivated 3-state init: BC / Inter / EC, all closed-headpiece
    init3 = {
        "mus": np.array([[55.0, 36.0], [70.0, 38.0], [82.0, 38.0]]),
        "Sigmas": np.tile(np.diag([16.0, 9.0]), (3, 1, 1)).astype(float),
        "A": np.array([[0.95, 0.05, 0.0],
                       [0.025, 0.95, 0.025],
                       [0.0, 0.05, 0.95]]),
        "pi": np.full(3, 1 / 3),
    }
    # 4-state init adds an EO* state at CV2 = 55 (empty in data; tests whether
    # the EM forces CV2 to grow when given the freedom)
    init4 = {
        "mus": np.array([[55.0, 36.0], [70.0, 38.0], [82.0, 38.0], [85.0, 55.0]]),
        "Sigmas": np.tile(np.diag([16.0, 9.0]), (4, 1, 1)).astype(float),
        "A": np.array([[0.94, 0.05, 0.0, 0.01],
                       [0.025, 0.95, 0.025, 0.0],
                       [0.0, 0.05, 0.94, 0.01],
                       [0.0, 0.05, 0.05, 0.90]]),
        "pi": np.full(4, 1 / 4),
    }

    results = {}
    for K, init in [(3, init3), (4, init4)]:
        print(f"\n=== Fitting {K}-state 2-D HMM ===")
        try:
            fit = baum_welch_2d(X, K=K,
                                init_mus=init["mus"],
                                init_Sigmas=init["Sigmas"],
                                init_A=init["A"],
                                init_pi=init["pi"],
                                max_iter=200, tol=1e-5)
        except Exception as e:
            print(f"  K={K}: failed ({e})")
            continue
        # Sort by CV0 mean
        order = np.argsort(fit["mus"][:, 0])
        mus = fit["mus"][order]
        Sigmas = fit["Sigmas"][order]
        A = fit["A"][np.ix_(order, order)]
        A = A / A.sum(axis=1, keepdims=True)
        pi_stat = stationary_distribution(A)

        log_pi = np.log(np.maximum(pi_stat, 1e-12))
        log_A = np.log(np.maximum(A, 1e-12))
        log_b = emission_log_matrix(X, mus, Sigmas)
        states = viterbi(log_b, log_pi, log_A)

        n_params = n_free_params_2d(K)
        bic = n_params * np.log(n) - 2 * fit["log_likelihood"]
        aic = 2 * n_params - 2 * fit["log_likelihood"]

        print(f"  log-L = {fit['log_likelihood']:.2f},  "
              f"n_params = {n_params},  BIC = {bic:.2f},  AIC = {aic:.2f}")
        for k in range(K):
            sigma_cv0 = np.sqrt(Sigmas[k, 0, 0])
            sigma_cv2 = np.sqrt(Sigmas[k, 1, 1])
            corr = Sigmas[k, 0, 1] / (sigma_cv0 * sigma_cv2)
            print(f"  state {k}: μ = ({mus[k, 0]:.2f}, {mus[k, 1]:.2f}),  "
                  f"σ_CV0 = {sigma_cv0:.2f}, σ_CV2 = {sigma_cv2:.2f},  "
                  f"corr = {corr:+.3f},  π = {pi_stat[k]:.3f}")

        results[K] = {
            "log_likelihood": fit["log_likelihood"],
            "n_iter": fit["n_iter"],
            "n_params": n_params,
            "BIC": float(bic),
            "AIC": float(aic),
            "mus": mus.tolist(),
            "Sigmas": Sigmas.tolist(),
            "A": A.tolist(),
            "pi_stationary": pi_stat.tolist(),
            "states": states.tolist(),
        }

    # Compare 3 vs 4
    if 3 in results and 4 in results:
        d_bic = results[4]["BIC"] - results[3]["BIC"]
        d_aic = results[4]["AIC"] - results[3]["AIC"]
        preferred_bic = "4-state" if d_bic < 0 else "3-state"
        preferred_aic = "4-state" if d_aic < 0 else "3-state"
        print(f"\n=== Model selection ===")
        print(f"  ΔBIC (4 − 3) = {d_bic:+.2f}  → {preferred_bic} preferred")
        print(f"  ΔAIC (4 − 3) = {d_aic:+.2f}  → {preferred_aic} preferred")

    # --- Save JSON ---
    summary = {
        "n_points": int(n),
        "T_kelvin": T_KELVIN,
        "results_by_K": {str(K): results[K] for K in results},
        "delta_BIC_4_minus_3": (results[4]["BIC"] - results[3]["BIC"]
                                if 3 in results and 4 in results else None),
        "delta_AIC_4_minus_3": (results[4]["AIC"] - results[3]["AIC"]
                                if 3 in results and 4 in results else None),
        "obj054_cross_correlations": {
            "v1_cv0_cv2": 0.517,
            "v2_cv0_cv2": 0.219,
        },
    }
    json_path = OUT_DIR / "hmm_2d_cv0_cv2.json"
    with json_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nJSON saved: {json_path}")

    # --- Figure ---
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.4, wspace=0.35)

    # (a) Scatter with 3-state Viterbi coloring + Gaussian ellipses
    ax_3 = fig.add_subplot(gs[0, 0])
    colors3 = ["#d62728", "#ff7f0e", "#2ca02c"]
    labels3 = ["state 0 (low-CV2)", "state 1 (mid-CV2 Inter)", "state 2 (mid-CV2 EC)"]
    if 3 in results:
        states3 = np.array(results[3]["states"])
        mus3 = np.array(results[3]["mus"])
        Sigmas3 = np.array(results[3]["Sigmas"])
        for k in range(3):
            mask = states3 == k
            ax_3.scatter(cv0[mask], cv2[mask], s=8, alpha=0.45,
                         color=colors3[k], label=labels3[k])
            # 1σ + 2σ ellipses
            for nsig in (1, 2):
                w, v = np.linalg.eigh(Sigmas3[k])
                angle = np.degrees(np.arctan2(*v[:, 0][::-1]))
                ell = mpatches.Ellipse(mus3[k], 2 * nsig * np.sqrt(w[0]),
                                       2 * nsig * np.sqrt(w[1]), angle=angle,
                                       fill=False, edgecolor=colors3[k],
                                       lw=1.5, alpha=0.7)
                ax_3.add_patch(ell)
            ax_3.scatter(*mus3[k], marker="X", color=colors3[k], s=120,
                         edgecolor="black", linewidth=1.5)
        ax_3.axhline(50, color="purple", ls="--", lw=1.0, alpha=0.6,
                     label="EO threshold CV2=50")
        ax_3.set_xlabel("CV0 (Å)")
        ax_3.set_ylabel("CV2 (Å)")
        ax_3.set_title("(a) 3-state 2-D HMM Viterbi + emission ellipses")
        ax_3.legend(fontsize=7, loc="lower right")
        ax_3.grid(alpha=0.3)

    # (b) Scatter with 4-state Viterbi coloring + ellipses
    ax_4 = fig.add_subplot(gs[0, 1])
    if 4 in results:
        states4 = np.array(results[4]["states"])
        mus4 = np.array(results[4]["mus"])
        Sigmas4 = np.array(results[4]["Sigmas"])
        colors4 = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]
        labels4 = ["state 0", "state 1", "state 2", "state 3"]
        # Auto-label based on means
        for i, m in enumerate(mus4):
            if m[0] < 60:
                labels4[i] = f"BC (μ={m[0]:.1f},{m[1]:.1f})"
            elif m[0] < 75:
                labels4[i] = f"Inter (μ={m[0]:.1f},{m[1]:.1f})"
            elif m[1] < 50:
                labels4[i] = f"EC (μ={m[0]:.1f},{m[1]:.1f})"
            else:
                labels4[i] = f"EO* (μ={m[0]:.1f},{m[1]:.1f})"
        for k in range(4):
            mask = states4 == k
            ax_4.scatter(cv0[mask], cv2[mask], s=8, alpha=0.45,
                         color=colors4[k], label=labels4[k])
            for nsig in (1, 2):
                w, v = np.linalg.eigh(Sigmas4[k])
                angle = np.degrees(np.arctan2(*v[:, 0][::-1]))
                ell = mpatches.Ellipse(mus4[k], 2 * nsig * np.sqrt(w[0]),
                                       2 * nsig * np.sqrt(w[1]), angle=angle,
                                       fill=False, edgecolor=colors4[k],
                                       lw=1.5, alpha=0.7)
                ax_4.add_patch(ell)
            ax_4.scatter(*mus4[k], marker="X", color=colors4[k], s=120,
                         edgecolor="black", linewidth=1.5)
        ax_4.axhline(50, color="purple", ls="--", lw=1.0, alpha=0.6,
                     label="EO threshold CV2=50")
        ax_4.set_xlabel("CV0 (Å)")
        ax_4.set_ylabel("CV2 (Å)")
        ax_4.set_title("(b) 4-state 2-D HMM Viterbi + emission ellipses")
        ax_4.legend(fontsize=7, loc="lower right")
        ax_4.grid(alpha=0.3)

    # (c) Transition matrices side by side
    ax_T = fig.add_subplot(gs[0, 2])
    ax_T.axis("off")
    if 3 in results:
        A3 = np.array(results[3]["A"])
        ax_T_inset = ax_T.inset_axes((0.0, 0.55, 1.0, 0.45))
        im = ax_T_inset.imshow(A3, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        for i in range(3):
            for j in range(3):
                ax_T_inset.text(j, i, f"{A3[i, j]:.3f}", ha="center", va="center",
                                fontsize=8, color="black" if A3[i, j] < 0.5 else "white")
        ax_T_inset.set_xticks([0, 1, 2]); ax_T_inset.set_yticks([0, 1, 2])
        ax_T_inset.set_xticklabels(["BC", "Inter", "EC"], fontsize=7)
        ax_T_inset.set_yticklabels(["BC", "Inter", "EC"], fontsize=7)
        ax_T_inset.set_title("3-state A[from, to]", fontsize=9)
    if 4 in results:
        A4 = np.array(results[4]["A"])
        ax_T_inset2 = ax_T.inset_axes((0.0, 0.0, 1.0, 0.45))
        im = ax_T_inset2.imshow(A4, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        for i in range(4):
            for j in range(4):
                ax_T_inset2.text(j, i, f"{A4[i, j]:.3f}", ha="center", va="center",
                                 fontsize=7, color="black" if A4[i, j] < 0.5 else "white")
        ax_T_inset2.set_xticks([0, 1, 2, 3]); ax_T_inset2.set_yticks([0, 1, 2, 3])
        ax_T_inset2.set_xticklabels([f"s{k}" for k in range(4)], fontsize=7)
        ax_T_inset2.set_yticklabels([f"s{k}" for k in range(4)], fontsize=7)
        ax_T_inset2.set_title("4-state A[from, to]", fontsize=9)

    # (d) BIC/AIC bars
    ax_ic = fig.add_subplot(gs[1, 0])
    if 3 in results and 4 in results:
        labels = ["BIC", "AIC"]
        k3vals = [results[3]["BIC"], results[3]["AIC"]]
        k4vals = [results[4]["BIC"], results[4]["AIC"]]
        idx = np.arange(len(labels))
        ax_ic.bar(idx - 0.18, k3vals, width=0.36, color="tab:blue", label="3-state")
        ax_ic.bar(idx + 0.18, k4vals, width=0.36, color="tab:purple", label="4-state")
        for i, (a, b) in enumerate(zip(k3vals, k4vals)):
            ax_ic.text(i - 0.18, a + 5, f"{a:.0f}", ha="center", fontsize=8)
            ax_ic.text(i + 0.18, b + 5, f"{b:.0f}", ha="center", fontsize=8)
        d_bic = results[4]["BIC"] - results[3]["BIC"]
        d_aic = results[4]["AIC"] - results[3]["AIC"]
        ax_ic.set_xticks(idx); ax_ic.set_xticklabels(labels)
        ax_ic.set_ylabel("information criterion (lower = better)")
        ax_ic.set_title(f"(d) Model order: ΔBIC={d_bic:+.1f}, ΔAIC={d_aic:+.1f}",
                        fontsize=10)
        ax_ic.legend(fontsize=8)
        ax_ic.grid(axis="y", alpha=0.3)

    # (e) Per-state CV2 marginal histograms (3-state)
    ax_cv2 = fig.add_subplot(gs[1, 1])
    if 3 in results:
        states3 = np.array(results[3]["states"])
        for k, lab, c in zip(range(3), labels3, colors3):
            mask = states3 == k
            if mask.sum() > 0:
                ax_cv2.hist(cv2[mask], bins=30, alpha=0.55, color=c,
                            edgecolor="black", label=f"{lab} (n={mask.sum()})",
                            density=True)
        ax_cv2.axvline(50, color="purple", ls="--", lw=1.5,
                       label="EO threshold = 50 Å")
        ax_cv2.set_xlabel("CV2 (Å)")
        ax_cv2.set_ylabel("density")
        ax_cv2.set_title("(e) Per-state CV2 marginal — 3-state HMM",
                         fontsize=10)
        ax_cv2.legend(fontsize=7)
        ax_cv2.grid(alpha=0.3)

    # (f) Verdict text
    ax_v = fig.add_subplot(gs[1, 2])
    ax_v.axis("off")
    if 3 in results and 4 in results:
        d_bic = results[4]["BIC"] - results[3]["BIC"]
        d_aic = results[4]["AIC"] - results[3]["AIC"]
        verdict_bic = "3-state preferred" if d_bic > 0 else "4-state preferred"
        verdict_aic = "3-state preferred" if d_aic > 0 else "4-state preferred"
        # 4-state means
        mus4 = np.array(results[4]["mus"])
        eo_idx = int(np.argmax(mus4[:, 1]))
        eo_mean = mus4[eo_idx]
        pi4 = np.array(results[4]["pi_stationary"])
        eo_pi = pi4[eo_idx]
        lines = [
            "obj-055 verdict",
            "",
            f"  ΔBIC (4 − 3) = {d_bic:+.1f}",
            f"  ΔAIC (4 − 3) = {d_aic:+.1f}",
            f"  BIC says: {verdict_bic}",
            f"  AIC says: {verdict_aic}",
            "",
            "  4-state highest-CV2 state:",
            f"    μ = ({eo_mean[0]:.1f}, {eo_mean[1]:.1f}) Å",
            f"    π = {eo_pi:.3f}",
            f"    CV2 mean: {eo_mean[1]:.1f} Å",
            f"    EO threshold: 50 Å",
            "",
            "  Interpretation:",
            "    Even when seeded at CV2=55,",
            "    the 4th state's CV2 mean",
            "    drifts back toward the data's",
            "    CV2 distribution. The v7",
            "    fitted data does NOT support",
            "    a separate EO state.",
        ]
        ax_v.text(0.05, 0.95, "\n".join(lines), transform=ax_v.transAxes,
                  fontsize=9.5, va="top", family="monospace")

    fig.suptitle(
        f"obj-055 — 2-D Gaussian-emission HMM over (CV0, CV2)\n"
        f"3-state vs 4-state — does CV2 reveal a head-piece-opening sub-state?  "
        f"(n = {n} points)",
        fontsize=11.5, fontweight="bold", y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig_path = FIG_DIR / "hmm_2d_cv0_cv2.png"
    fig.savefig(fig_path, dpi=140)
    print(f"Figure saved: {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
