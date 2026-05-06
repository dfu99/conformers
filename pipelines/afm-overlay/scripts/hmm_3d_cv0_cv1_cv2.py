#!/usr/bin/env python3
"""3-D Gaussian-emission HMM over (CV0, CV1, CV2) (obj-059).

Audit §28 v19 P=5 follow-up to obj-055. obj-054 found CV1 highly
correlated with CV0 (r = 0.84 V1 / 0.87 V2) — CV1 is a noisier CV0
proxy. Adding CV1 as a third HMM dimension should NOT change the
partition relative to obj-055's 2-D (CV0, CV2) result. This pass
fits the full 3-D HMM and confirms.

Implementation: from-scratch (no `hmmlearn`):
  - 3-D Gaussian emissions per state (mean ∈ ℝ³, covariance ∈ ℝ³ˣ³)
  - Forward-backward in log-space
  - Baum-Welch EM with 3×3 covariance updates
  - Viterbi decoding

Reference physical states (extending obj-055):
  - BC:    CV0 ≈ 55, CV1 ≈ 55, CV2 ≈ 36
  - Inter: CV0 ≈ 70, CV1 ≈ 70, CV2 ≈ 38
  - EC:    CV0 ≈ 82, CV1 ≈ 82, CV2 ≈ 38
  - EO* (4-state init only): CV0 ≈ 85, CV1 ≈ 85, CV2 ≈ 55 (empty in v7)

Both 3-state and 4-state HMMs fitted from physically-motivated
inits. BIC/AIC compared to obj-055 (2-D) and obj-051 (1-D).

Inputs:
  results/afm_pipeline/free_energy_profile/cv_series.npz

Outputs:
  figures/hmm_3d_cv0_cv1_cv2.png
  results/afm_pipeline/free_energy_profile/hmm_3d_cv0_cv1_cv2.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv, slogdet

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

T_KELVIN = 300.0
LOG_EPS = -1e30
D = 3


def gaussian3d_log_pdf(x: np.ndarray, mu: np.ndarray, Sigma: np.ndarray) -> np.ndarray:
    sign, logdet = slogdet(Sigma)
    if sign <= 0:
        Sigma = Sigma + 1e-3 * np.eye(D)
        sign, logdet = slogdet(Sigma)
    inv_Sigma = inv(Sigma)
    diff = x - mu[None, :]
    quad = np.einsum("ti,ij,tj->t", diff, inv_Sigma, diff)
    return -0.5 * (D * np.log(2 * np.pi) + logdet + quad)


def emission_log_matrix(x: np.ndarray, mus: np.ndarray, Sigmas: np.ndarray) -> np.ndarray:
    T = x.shape[0]
    K = mus.shape[0]
    log_b = np.empty((T, K))
    for k in range(K):
        log_b[:, k] = gaussian3d_log_pdf(x, mus[k], Sigmas[k])
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


def baum_welch_3d(x, K, init_mus, init_Sigmas, init_A, init_pi,
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

        denom = gamma.sum(axis=0)
        for k in range(K):
            w = gamma[:, k]
            mus[k] = (w[:, None] * x).sum(axis=0) / max(denom[k], 1e-12)
            diff = x - mus[k][None, :]
            Sigmas[k] = (
                np.einsum("t,ti,tj->ij", w, diff, diff)
                / max(denom[k], 1e-12)
            )
            Sigmas[k] += 1e-3 * np.eye(D)

        if abs(ll - prev_ll) < tol:
            print(f"  EM converged at iter {it} (ΔLL = {ll - prev_ll:.2e})")
            break
        prev_ll = ll

    return {
        "mus": mus, "Sigmas": Sigmas, "A": np.exp(log_A),
        "pi": np.exp(log_pi), "log_likelihood": float(ll),
        "n_iter": int(it + 1),
    }


def n_free_params_3d(K: int) -> int:
    """K-state 3-D Gaussian-emission HMM:
       - K means × 3        = 3K
       - K covariances × 6  = 6K  (3x3 sym → 6 unique entries)
       - (K-1) initial-π
       - K(K-1) transitions
    """
    return 3 * K + 6 * K + (K - 1) + K * (K - 1)


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
    cv1 = np.concatenate([cv["cv1_v1"], cv["cv1_v2"]]).astype(float)
    cv2 = np.concatenate([cv["cv2_v1"], cv["cv2_v2"]]).astype(float)
    X = np.column_stack([cv0, cv1, cv2])
    n = X.shape[0]
    print(f"Loaded n = {n} points (V1+V2). "
          f"CV0 [{cv0.min():.1f}, {cv0.max():.1f}], "
          f"CV1 [{cv1.min():.1f}, {cv1.max():.1f}], "
          f"CV2 [{cv2.min():.1f}, {cv2.max():.1f}].")

    # 3-state init: BC / Inter / EC, all closed-headpiece
    init3 = {
        "mus": np.array([[55.0, 55.0, 36.0],
                         [70.0, 70.0, 38.0],
                         [82.0, 82.0, 38.0]]),
        "Sigmas": np.tile(np.diag([16.0, 16.0, 9.0]), (3, 1, 1)).astype(float),
        "A": np.array([[0.95, 0.05, 0.0],
                       [0.025, 0.95, 0.025],
                       [0.0, 0.05, 0.95]]),
        "pi": np.full(3, 1 / 3),
    }
    # 4-state init seeds an EO* state at CV2 = 55 (empty in v7)
    init4 = {
        "mus": np.array([[55.0, 55.0, 36.0],
                         [70.0, 70.0, 38.0],
                         [82.0, 82.0, 38.0],
                         [85.0, 85.0, 55.0]]),
        "Sigmas": np.tile(np.diag([16.0, 16.0, 9.0]), (4, 1, 1)).astype(float),
        "A": np.array([[0.94, 0.05, 0.0,  0.01],
                       [0.025, 0.95, 0.025, 0.0],
                       [0.0, 0.05, 0.94, 0.01],
                       [0.0, 0.05, 0.05, 0.90]]),
        "pi": np.full(4, 1 / 4),
    }

    results: dict[int, dict] = {}
    for K, init in [(3, init3), (4, init4)]:
        print(f"\n=== Fitting {K}-state 3-D HMM ===")
        try:
            fit = baum_welch_3d(X, K=K,
                                init_mus=init["mus"],
                                init_Sigmas=init["Sigmas"],
                                init_A=init["A"],
                                init_pi=init["pi"],
                                max_iter=200, tol=1e-5)
        except Exception as e:
            print(f"  K={K}: failed ({e})")
            continue
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

        n_params = n_free_params_3d(K)
        bic = n_params * np.log(n) - 2 * fit["log_likelihood"]
        aic = 2 * n_params - 2 * fit["log_likelihood"]

        print(f"  log-L = {fit['log_likelihood']:.2f},  "
              f"n_params = {n_params},  BIC = {bic:.2f},  AIC = {aic:.2f}")
        for k in range(K):
            sigma_cv0 = np.sqrt(Sigmas[k, 0, 0])
            sigma_cv1 = np.sqrt(Sigmas[k, 1, 1])
            sigma_cv2 = np.sqrt(Sigmas[k, 2, 2])
            corr_01 = Sigmas[k, 0, 1] / (sigma_cv0 * sigma_cv1)
            corr_02 = Sigmas[k, 0, 2] / (sigma_cv0 * sigma_cv2)
            corr_12 = Sigmas[k, 1, 2] / (sigma_cv1 * sigma_cv2)
            print(f"  state {k}: μ = ({mus[k, 0]:.2f}, {mus[k, 1]:.2f}, "
                  f"{mus[k, 2]:.2f}),  σ = ({sigma_cv0:.2f}, "
                  f"{sigma_cv1:.2f}, {sigma_cv2:.2f}),  "
                  f"corr(01,02,12)=({corr_01:+.2f},{corr_02:+.2f},"
                  f"{corr_12:+.2f}),  π = {pi_stat[k]:.3f}")

        results[K] = {
            "log_likelihood": float(fit["log_likelihood"]),
            "n_iter": int(fit["n_iter"]),
            "n_params": int(n_params),
            "BIC": float(bic),
            "AIC": float(aic),
            "mus": mus.tolist(),
            "Sigmas": Sigmas.tolist(),
            "A": A.tolist(),
            "pi_stationary": pi_stat.tolist(),
            "states": states.tolist(),
        }

    if 3 in results and 4 in results:
        d_bic = results[4]["BIC"] - results[3]["BIC"]
        d_aic = results[4]["AIC"] - results[3]["AIC"]
        preferred_bic = "4-state" if d_bic < 0 else "3-state"
        preferred_aic = "4-state" if d_aic < 0 else "3-state"
        print(f"\n=== Model selection (3-D) ===")
        print(f"  ΔBIC (4 − 3) = {d_bic:+.2f}  → {preferred_bic} preferred")
        print(f"  ΔAIC (4 − 3) = {d_aic:+.2f}  → {preferred_aic} preferred")

    obj055_path = OUT_DIR / "hmm_2d_cv0_cv2.json"
    obj055_summary = None
    if obj055_path.exists():
        with obj055_path.open() as fh:
            obj055_summary = json.load(fh)
        print(f"\n=== Comparison vs obj-055 (2-D HMM CV0+CV2) ===")
        for K in (3, 4):
            if K in results and str(K) in obj055_summary["results_by_K"]:
                obj055_K = obj055_summary["results_by_K"][str(K)]
                print(f"  K={K}: 3-D BIC = {results[K]['BIC']:.2f}, "
                      f"2-D BIC = {obj055_K['BIC']:.2f}, "
                      f"3-D - 2-D = {results[K]['BIC'] - obj055_K['BIC']:+.2f}")

    summary = {
        "n_points": int(n),
        "T_kelvin": T_KELVIN,
        "results_by_K": {str(K): results[K] for K in results},
        "delta_BIC_4_minus_3": (results[4]["BIC"] - results[3]["BIC"]
                                if 3 in results and 4 in results else None),
        "delta_AIC_4_minus_3": (results[4]["AIC"] - results[3]["AIC"]
                                if 3 in results and 4 in results else None),
        "comparison_vs_obj055_2d": (
            {f"K={K}": {
                "3d_BIC": results[K]["BIC"],
                "2d_BIC": obj055_summary["results_by_K"][str(K)]["BIC"],
                "delta_3d_minus_2d_BIC": (results[K]["BIC"]
                                          - obj055_summary["results_by_K"][str(K)]["BIC"]),
            } for K in (3, 4) if K in results and str(K) in obj055_summary["results_by_K"]}
            if obj055_summary else None
        ),
        "obj054_cv_correlations": {
            "v1_cv0_cv1": 0.84,
            "v2_cv0_cv1": 0.87,
            "v1_cv0_cv2": 0.52,
            "v2_cv0_cv2": 0.22,
        },
    }
    json_path = OUT_DIR / "hmm_3d_cv0_cv1_cv2.json"
    with json_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nJSON saved: {json_path}")

    plot_3d(results, cv0, cv1, cv2)
    return 0


def plot_3d(results: dict, cv0: np.ndarray, cv1: np.ndarray, cv2: np.ndarray) -> None:
    fig = plt.figure(figsize=(15.0, 11.0))
    gs = fig.add_gridspec(3, 3, hspace=0.42, wspace=0.34,
                          height_ratios=[1.0, 1.0, 0.85])

    K_main = 4 if 4 in results else 3
    states = np.array(results[K_main]["states"]) if K_main in results else None
    mus = np.array(results[K_main]["mus"]) if K_main in results else None
    K = mus.shape[0] if mus is not None else 0

    palette = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd", "#1f77b4"]

    def mk_label(m: np.ndarray) -> str:
        if m[0] < 60:
            return f"BC (CV0={m[0]:.1f}, CV2={m[2]:.1f})"
        elif m[0] < 75:
            return f"Inter (CV0={m[0]:.1f}, CV2={m[2]:.1f})"
        elif m[2] < 50:
            return f"EC (CV0={m[0]:.1f}, CV2={m[2]:.1f})"
        else:
            return f"EO* (CV0={m[0]:.1f}, CV2={m[2]:.1f})"

    # (a) CV0–CV1 projection
    ax = fig.add_subplot(gs[0, 0])
    if states is not None and mus is not None:
        for k in range(K):
            mask = states == k
            ax.scatter(cv0[mask], cv1[mask], s=8, alpha=0.45,
                       color=palette[k], label=mk_label(mus[k]))
            ax.scatter(mus[k, 0], mus[k, 1], marker="X", s=120,
                       color=palette[k], edgecolor="black", linewidth=1.5)
    ax.set_xlabel("CV0 (Å)")
    ax.set_ylabel("CV1 (Å)")
    ax.set_title(f"(a) CV0-CV1 projection ({K_main}-state Viterbi)")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)

    # (b) CV0–CV2 projection
    ax = fig.add_subplot(gs[0, 1])
    if states is not None and mus is not None:
        for k in range(K):
            mask = states == k
            ax.scatter(cv0[mask], cv2[mask], s=8, alpha=0.45,
                       color=palette[k], label=mk_label(mus[k]))
            ax.scatter(mus[k, 0], mus[k, 2], marker="X", s=120,
                       color=palette[k], edgecolor="black", linewidth=1.5)
        ax.axhline(50, color="purple", ls="--", lw=1.0, alpha=0.6,
                   label="EO threshold CV2=50")
    ax.set_xlabel("CV0 (Å)")
    ax.set_ylabel("CV2 (Å)")
    ax.set_title(f"(b) CV0-CV2 projection ({K_main}-state Viterbi)")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)

    # (c) CV1–CV2 projection
    ax = fig.add_subplot(gs[0, 2])
    if states is not None and mus is not None:
        for k in range(K):
            mask = states == k
            ax.scatter(cv1[mask], cv2[mask], s=8, alpha=0.45,
                       color=palette[k], label=mk_label(mus[k]))
            ax.scatter(mus[k, 1], mus[k, 2], marker="X", s=120,
                       color=palette[k], edgecolor="black", linewidth=1.5)
        ax.axhline(50, color="purple", ls="--", lw=1.0, alpha=0.6)
    ax.set_xlabel("CV1 (Å)")
    ax.set_ylabel("CV2 (Å)")
    ax.set_title(f"(c) CV1-CV2 projection ({K_main}-state Viterbi)")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)

    # (d) BIC/AIC comparison vs obj-055 2-D and obj-051 1-D
    ax = fig.add_subplot(gs[1, 0])
    obj055_path = OUT_DIR / "hmm_2d_cv0_cv2.json"
    obj055 = None
    if obj055_path.exists():
        with obj055_path.open() as fh:
            obj055 = json.load(fh)
    obj051_path = OUT_DIR / "hmm_model_selection.json"
    obj051 = None
    if obj051_path.exists():
        with obj051_path.open() as fh:
            obj051 = json.load(fh)
    rows_txt = [["model", "K", "BIC", "AIC"]]
    if obj051:
        for K in (3, 4):
            key = f"K={K}"
            d = obj051.get("k_results", {}).get(key, None)
            if d:
                rows_txt.append(["1-D (obj-051)", str(K),
                                 f"{d['BIC']:.1f}", f"{d['AIC']:.1f}"])
    if obj055:
        for K in (3, 4):
            d = obj055.get("results_by_K", {}).get(str(K), None)
            if d:
                rows_txt.append(["2-D (obj-055)", str(K),
                                 f"{d['BIC']:.1f}", f"{d['AIC']:.1f}"])
    for K in (3, 4):
        if K in results:
            rows_txt.append(["3-D (obj-059)", str(K),
                             f"{results[K]['BIC']:.1f}",
                             f"{results[K]['AIC']:.1f}"])
    ax.axis("off")
    table = ax.table(cellText=rows_txt, cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    ax.set_title("(d) Model selection across HMM dimensionalities",
                 fontsize=10, fontweight="bold")

    # (e) State means table
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    rows_means = [["state", "μ_CV0", "μ_CV1", "μ_CV2", "π"]]
    if K_main in results:
        for k in range(K_main):
            m = results[K_main]["mus"][k]
            pi_v = results[K_main]["pi_stationary"][k]
            rows_means.append([f"s{k}", f"{m[0]:.1f}",
                               f"{m[1]:.1f}", f"{m[2]:.1f}",
                               f"{pi_v:.3f}"])
    table2 = ax.table(cellText=rows_means, cellLoc="center", loc="center")
    table2.auto_set_font_size(False)
    table2.set_fontsize(9)
    table2.scale(1, 1.5)
    ax.set_title(f"(e) {K_main}-state HMM means + stationary π",
                 fontsize=10, fontweight="bold")

    # (f) EO recovery diagnostic — track 4-state seed CV2 trajectory
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    text_lines = []
    if 4 in results:
        seed_cv2_initial = 55.0
        seed_cv2_final = max(results[4]["mus"], key=lambda m: m[2])[2]
        seed_pi = max(zip(results[4]["mus"], results[4]["pi_stationary"]),
                      key=lambda mp: mp[0][2])[1]
        text_lines.append(f"4-state seeded at CV2 = {seed_cv2_initial}")
        text_lines.append(f"  highest-CV2 state final: CV2 = {seed_cv2_final:.2f}")
        text_lines.append(f"  highest-CV2 state π = {seed_pi:.3f}")
        text_lines.append(f"  EO recovered: {'YES' if seed_cv2_final >= 50 else 'NO'}")
        text_lines.append(f"  (CV2 threshold for EO: 50 Å)")
        text_lines.append("")
    if obj055:
        text_lines.append("obj-055 2-D HMM 4-state:")
        if "4" in obj055.get("results_by_K", {}):
            obj055_4 = obj055["results_by_K"]["4"]
            seed_cv2_2d = max(obj055_4["mus"], key=lambda m: m[1])[1]
            text_lines.append(f"  highest-CV2 final = {seed_cv2_2d:.2f}")
            text_lines.append(f"  EO recovered: {'YES' if seed_cv2_2d >= 50 else 'NO'}")
    ax.text(0.05, 0.95, "\n".join(text_lines),
            transform=ax.transAxes, fontsize=10, va="top", ha="left",
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="#fff3cd",
                      edgecolor="#856404"))
    ax.set_title("(f) 4-state EO recovery diagnostic",
                 fontsize=10, fontweight="bold")

    # (g) Transition matrix heatmap (4-state, 3-D)
    ax = fig.add_subplot(gs[2, 0])
    if K_main in results:
        A = np.array(results[K_main]["A"])
        im = ax.imshow(A, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        for i in range(K_main):
            for j in range(K_main):
                ax.text(j, i, f"{A[i, j]:.3f}", ha="center", va="center",
                        fontsize=8,
                        color="black" if A[i, j] < 0.5 else "white")
        ax.set_xticks(range(K_main))
        ax.set_yticks(range(K_main))
        ax.set_xticklabels([f"s{k}" for k in range(K_main)], fontsize=9)
        ax.set_yticklabels([f"s{k}" for k in range(K_main)], fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(f"(g) {K_main}-state transition matrix A[from, to]",
                     fontsize=10, fontweight="bold")

    # (h) Comparison of state means across 1-D / 2-D / 3-D HMMs
    ax = fig.add_subplot(gs[2, 1:])
    if K_main in results:
        mus_3d = np.array(results[K_main]["mus"])
        offsets = np.arange(K_main) * 1.0
        if obj051 and f"K={K_main}" in obj051.get("k_results", {}):
            mus_1d = np.array(obj051["k_results"][f"K={K_main}"]["means"])
            order_1d = np.argsort(mus_1d)
            ax.scatter(mus_1d[order_1d], offsets, s=120, marker="o",
                       color="#1f77b4", label="1-D μ_CV0")
        if obj055 and str(K_main) in obj055.get("results_by_K", {}):
            mus_2d = np.array(obj055["results_by_K"][str(K_main)]["mus"])
            ax.scatter(mus_2d[:, 0], offsets, s=120, marker="s",
                       color="#ff7f0e", label="2-D μ_CV0")
            ax.scatter(mus_2d[:, 1], offsets + 0.05, s=80, marker="s",
                       color="#fd8d3c", label="2-D μ_CV2", alpha=0.7)
        ax.scatter(mus_3d[:, 0], offsets, s=140, marker="X",
                   color="#2ca02c", label="3-D μ_CV0", edgecolor="black")
        ax.scatter(mus_3d[:, 1], offsets - 0.05, s=80, marker="X",
                   color="#74c476", label="3-D μ_CV1", alpha=0.7)
        ax.scatter(mus_3d[:, 2], offsets - 0.10, s=80, marker="X",
                   color="#5a9b34", label="3-D μ_CV2", alpha=0.7)
        ax.set_yticks(offsets)
        ax.set_yticklabels([f"state {k}" for k in range(K_main)])
        ax.set_xlabel("CV value (Å)")
        ax.set_title(f"(h) State-mean comparison across HMM dimensionalities (K={K_main})",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc="upper left", ncols=3)
        ax.grid(alpha=0.3)

    fig.suptitle(
        f"3-D Gaussian HMM over (CV0, CV1, CV2) (obj-059) — "
        f"{n_free_params_3d(K_main)} params at K={K_main}; CV1 redundant with CV0 "
        f"(obj-054 r ≥ 0.84); same biology as obj-055 (no EO support, "
        f"CV2 primary clustering axis)",
        fontsize=11, fontweight="bold", y=0.995,
    )
    out_path = FIG_DIR / "hmm_3d_cv0_cv1_cv2.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    print(f"Figure saved: {out_path}")


if __name__ == "__main__":
    raise SystemExit(main())
