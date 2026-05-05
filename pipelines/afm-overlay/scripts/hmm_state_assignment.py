#!/usr/bin/env python3
"""3-state Gaussian-emission HMM over per-frame CV0 (obj-048).

Audit-2026-05-05 §17.6 v9 follow-up. obj-047 BinSeg+BIC detected
3-4 change-points per video in the per-frame CV0 trace; this v10
formalizes the same picture as a 3-state Gaussian-emission HMM:

  state 0 = BC           (mean ≈ 55 Å)
  state 1 = Intermediate (mean ≈ 70 Å)
  state 2 = EC           (mean ≈ 82 Å)

Implementation is from-scratch (no `hmmlearn` dependency):
  - forward / backward algorithms in log-space for stability
  - Baum-Welch EM for parameter learning
  - Viterbi for state decoding

Train on the concatenated V1+V2 per-frame CV0 trace; decode V1 and
V2 separately. Output:
  - transition matrix and per-state Gaussian parameters
  - Viterbi-decoded state sequence per video
  - per-state dwell-time distribution + mean lifetime in seconds
    (HS-AFM at 1 fps)
  - co-localization of HMM state-transitions vs obj-047 BinSeg+BIC
    change-points (sanity check)

Outputs:
  figures/hmm_state_assignment.png
  results/afm_pipeline/free_energy_profile/hmm_states.json
  results/afm_pipeline/free_energy_profile/hmm_states.npz
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
CP_SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "change_points.json"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

T_KELVIN = 300.0
N_STATES = 3
LOG_EPS = -1e30


def gaussian_log_pdf(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    sigma = max(sigma, 1e-3)
    return -0.5 * np.log(2 * np.pi * sigma * sigma) - 0.5 * ((x - mu) / sigma) ** 2


def emission_log_matrix(x: np.ndarray, mus: np.ndarray, sigmas: np.ndarray) -> np.ndarray:
    """Return (T, K) matrix of log-emissions."""
    T = x.size
    K = mus.size
    log_b = np.empty((T, K))
    for k in range(K):
        log_b[:, k] = gaussian_log_pdf(x, mus[k], sigmas[k])
    return log_b


def logsumexp(a: np.ndarray, axis: int | None = None) -> np.ndarray:
    a_max = np.max(a, axis=axis, keepdims=True)
    a_max = np.where(np.isfinite(a_max), a_max, 0.0)
    out = np.log(np.sum(np.exp(a - a_max), axis=axis, keepdims=True)) + a_max
    return np.squeeze(out, axis=axis) if axis is not None else out


def forward_backward(log_b: np.ndarray, log_pi: np.ndarray, log_A: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    T, K = log_b.shape
    log_alpha = np.full((T, K), LOG_EPS)
    log_beta = np.full((T, K), LOG_EPS)

    log_alpha[0] = log_pi + log_b[0]
    for t in range(1, T):
        # log_alpha[t, j] = log sum_i alpha[t-1, i] * A[i, j] * b[t, j]
        prev = log_alpha[t - 1][:, None] + log_A           # (K, K)
        log_alpha[t] = logsumexp(prev, axis=0) + log_b[t]  # (K,)

    log_beta[T - 1] = 0.0
    for t in range(T - 2, -1, -1):
        nxt = log_A + log_b[t + 1][None, :] + log_beta[t + 1][None, :]
        log_beta[t] = logsumexp(nxt, axis=1)

    log_likelihood = float(logsumexp(log_alpha[T - 1]))
    return log_alpha, log_beta, log_likelihood


def viterbi(log_b: np.ndarray, log_pi: np.ndarray, log_A: np.ndarray) -> np.ndarray:
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


def baum_welch(x: np.ndarray, K: int, init_mus: np.ndarray,
               init_sigmas: np.ndarray, init_A: np.ndarray, init_pi: np.ndarray,
               max_iter: int = 100, tol: float = 1e-5,
               freeze_means: bool = False) -> dict:
    """EM training for a Gaussian-emission HMM.

    If freeze_means=True, only sigmas + A + pi are updated (means stay at init);
    useful for ensuring the BC/Intermediate/EC state interpretation is preserved.
    """
    T = x.size
    mus = init_mus.copy()
    sigmas = init_sigmas.copy()
    log_pi = np.log(np.maximum(init_pi, 1e-12))
    log_A = np.log(np.maximum(init_A, 1e-12))

    prev_ll = -np.inf
    ll = -np.inf
    it = 0
    for it in range(max_iter):
        log_b = emission_log_matrix(x, mus, sigmas)
        log_alpha, log_beta, ll = forward_backward(log_b, log_pi, log_A)
        if not np.isfinite(ll):
            print(f"  iter {it}: non-finite log-likelihood; aborting EM")
            break
        log_gamma = log_alpha + log_beta
        log_gamma -= logsumexp(log_gamma, axis=1)[:, None]
        gamma = np.exp(log_gamma)

        # ξ[t,i,j] = α[t,i] A[i,j] b[t+1,j] β[t+1,j] / Z
        log_xi_unnorm = (log_alpha[:-1, :, None]
                         + log_A[None, :, :]
                         + log_b[1:, None, :]
                         + log_beta[1:, None, :])
        log_xi_norm = logsumexp(
            log_xi_unnorm.reshape(T - 1, K * K), axis=1)
        log_xi = log_xi_unnorm - log_xi_norm[:, None, None]
        xi = np.exp(log_xi)

        # Updates
        log_pi = log_gamma[0]
        A_num = xi.sum(axis=0)            # (K, K)
        A_den = gamma[:-1].sum(axis=0)    # (K,)
        A = A_num / np.maximum(A_den[:, None], 1e-12)
        A /= A.sum(axis=1, keepdims=True)
        log_A = np.log(np.maximum(A, 1e-12))

        if not freeze_means:
            denom = gamma.sum(axis=0)
            mus = (gamma * x[:, None]).sum(axis=0) / np.maximum(denom, 1e-12)
        sq = (x[:, None] - mus[None, :]) ** 2
        sigmas = np.sqrt(np.maximum((gamma * sq).sum(axis=0)
                                    / np.maximum(gamma.sum(axis=0), 1e-12),
                                    1e-2))

        if abs(ll - prev_ll) < tol:
            print(f"  EM converged at iter {it} (ΔLL = {ll - prev_ll:.2e})")
            break
        prev_ll = ll

    return {
        "mus": mus,
        "sigmas": sigmas,
        "A": np.exp(log_A),
        "pi": np.exp(log_pi),
        "log_likelihood": float(ll),
        "n_iter": int(it + 1),
    }


def dwell_times(states: np.ndarray, K: int) -> dict:
    """Return list of dwell-time runs per state."""
    out = {k: [] for k in range(K)}
    if states.size == 0:
        return out
    cur = states[0]
    n = 1
    for s in states[1:]:
        if s == cur:
            n += 1
        else:
            out[int(cur)].append(int(n))
            cur = int(s)
            n = 1
    out[int(cur)].append(int(n))
    return out


def state_transitions(states: np.ndarray) -> list[int]:
    """Return frame indices where the state changes (the new-state index)."""
    return [int(i) for i in range(1, states.size) if states[i] != states[i - 1]]


def co_locate(hmm_transitions: list[int], cps: list[int], tol: int = 5) -> dict:
    matched = []
    for cp in cps:
        nearest = min(hmm_transitions, key=lambda t: abs(t - cp), default=None)
        if nearest is None:
            matched.append(None)
            continue
        matched.append(nearest if abs(nearest - cp) <= tol else None)
    n_matched = sum(1 for m in matched if m is not None)
    return {
        "binseg_change_points": list(cps),
        "matched_hmm_transition_frames": matched,
        "n_binseg_change_points": len(cps),
        "n_co_located": int(n_matched),
        "tolerance_frames": int(tol),
    }


def fit_exponential_lifetime(durations: list[int]) -> tuple[float, float]:
    """MLE estimator for an exponential dwell-time, plus its 95 % CI half-width."""
    arr = np.array(durations, dtype=float)
    if arr.size < 2:
        return float("nan"), float("nan")
    mean = float(arr.mean())
    se = mean / np.sqrt(arr.size)
    return mean, 1.96 * se


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    cv = np.load(SRC)
    cv0_v1 = np.asarray(cv["cv0_v1"], dtype=float)
    cv0_v2 = np.asarray(cv["cv0_v2"], dtype=float)
    cv0_all = np.concatenate([cv0_v1, cv0_v2])

    # State priors based on obj-043 / obj-038 boundaries
    init_mus = np.array([55.0, 70.0, 82.0])
    init_sigmas = np.array([4.0, 4.0, 4.0])
    init_A = np.array([
        [0.95, 0.05, 0.0],
        [0.025, 0.95, 0.025],
        [0.0, 0.05, 0.95],
    ])
    init_pi = np.array([1 / 3, 1 / 3, 1 / 3])

    print("Training Gaussian-emission HMM on V1+V2 (Baum-Welch)...")
    print(f"  data: {cv0_all.size} frames "
          f"(V1={cv0_v1.size}, V2={cv0_v2.size})")
    fit = baum_welch(cv0_all, K=N_STATES,
                     init_mus=init_mus, init_sigmas=init_sigmas,
                     init_A=init_A, init_pi=init_pi,
                     max_iter=200, tol=1e-6, freeze_means=False)

    # Sort states by mean so labels are stable: 0 = BC, 1 = Intermediate, 2 = EC
    order = np.argsort(fit["mus"])
    mus = fit["mus"][order]
    sigmas = fit["sigmas"][order]
    A = fit["A"][np.ix_(order, order)]
    A = A / A.sum(axis=1, keepdims=True)
    state_labels = ["BC", "Intermediate", "EC"]

    # Long-run stationary distribution of A (left eigenvector at λ=1)
    eigvals, eigvecs = np.linalg.eig(A.T)
    idx_one = int(np.argmin(np.abs(eigvals - 1.0)))
    pi = np.real(eigvecs[:, idx_one])
    pi = np.abs(pi) / np.sum(np.abs(pi))

    print("\n=== Learned HMM parameters ===")
    print(f"  log-likelihood = {fit['log_likelihood']:.2f},  n_iter = {fit['n_iter']}")
    for k, lab in enumerate(state_labels):
        print(f"  state {k} ({lab}):  μ = {mus[k]:.2f} Å,  σ = {sigmas[k]:.2f} Å")
    print("  Transition matrix:")
    print("              ", "    ".join(state_labels))
    for k, lab in enumerate(state_labels):
        row = "  ".join(f"{A[k, j]:.4f}" for j in range(N_STATES))
        print(f"   {lab:<13}{row}")
    print(f"  Stationary π:  BC={pi[0]:.3f}  Inter={pi[1]:.3f}  EC={pi[2]:.3f}")

    # --- Decode V1 and V2 separately with the trained HMM ---
    log_b_v1 = emission_log_matrix(cv0_v1, mus, sigmas)
    log_b_v2 = emission_log_matrix(cv0_v2, mus, sigmas)
    log_pi = np.log(np.maximum(pi, 1e-12))
    log_A = np.log(np.maximum(A, 1e-12))
    states_v1 = viterbi(log_b_v1, log_pi, log_A)
    states_v2 = viterbi(log_b_v2, log_pi, log_A)

    print("\n=== Viterbi state assignments ===")
    for label, s in [("V1", states_v1), ("V2", states_v2)]:
        bins = np.bincount(s, minlength=N_STATES)
        frac = bins / s.size
        bc, im, ec = bins
        print(f"  {label}: BC={bc} ({frac[0]:.1%})  "
              f"Intermediate={im} ({frac[1]:.1%})  "
              f"EC={ec} ({frac[2]:.1%})")

    # --- Dwell-time analysis ---
    dwell_v1 = dwell_times(states_v1, N_STATES)
    dwell_v2 = dwell_times(states_v2, N_STATES)
    dwell_all = {k: dwell_v1[k] + dwell_v2[k] for k in range(N_STATES)}
    print("\n=== Mean dwell time per state (frames @ 1 fps = seconds) ===")
    lifetime_summary: dict = {}
    for k, lab in enumerate(state_labels):
        mean, ci = fit_exponential_lifetime(dwell_all[k])
        n = len(dwell_all[k])
        max_run = max(dwell_all[k]) if dwell_all[k] else 0
        print(f"  {lab}:  n_visits={n}  mean lifetime = {mean:.1f} ± {ci:.1f} s  "
              f"(max run {max_run} frames)")
        lifetime_summary[lab] = {
            "n_visits": int(n),
            "mean_lifetime_frames": float(mean) if np.isfinite(mean) else None,
            "ci95_half_width_frames": float(ci) if np.isfinite(ci) else None,
            "max_run_frames": int(max_run),
            "lifetime_distribution_frames": list(map(int, dwell_all[k])),
        }

    # --- Co-localize HMM state-transitions with obj-047 BinSeg+BIC ---
    co_summary: dict = {"v1": None, "v2": None}
    if CP_SRC.exists():
        cp_data = json.loads(CP_SRC.read_text())
        cps_v1 = cp_data["results"]["v1_frame_cv0_smoothed"]["binseg_change_points_idx"]
        cps_v2 = cp_data["results"]["v2_frame_cv0_smoothed"]["binseg_change_points_idx"]
        trans_v1 = state_transitions(states_v1)
        trans_v2 = state_transitions(states_v2)
        co_v1 = co_locate(trans_v1, cps_v1, tol=10)
        co_v2 = co_locate(trans_v2, cps_v2, tol=10)
        co_summary["v1"] = co_v1
        co_summary["v2"] = co_v2
        print("\n=== Co-localization with obj-047 BinSeg+BIC change-points ===")
        for label, info in [("V1", co_v1), ("V2", co_v2)]:
            print(f"  {label}: {info['n_co_located']}/{info['n_binseg_change_points']} BinSeg "
                  f"change-points have an HMM transition within ±10 frames.")

    # --- Write JSON + NPZ ---
    json_path = OUT_DIR / "hmm_states.json"
    with json_path.open("w") as fh:
        json.dump({
            "n_states": N_STATES,
            "state_labels": state_labels,
            "T_kelvin": T_KELVIN,
            "log_likelihood": fit["log_likelihood"],
            "n_iter": fit["n_iter"],
            "means_A": mus.tolist(),
            "stds_A": sigmas.tolist(),
            "transition_matrix": A.tolist(),
            "stationary_distribution": pi.tolist(),
            "viterbi_state_fractions": {
                "v1": (np.bincount(states_v1, minlength=N_STATES) / states_v1.size).tolist(),
                "v2": (np.bincount(states_v2, minlength=N_STATES) / states_v2.size).tolist(),
            },
            "dwell_time_summary": lifetime_summary,
            "co_localization_with_obj047": co_summary,
        }, fh, indent=2)
    print(f"\nJSON saved:    {json_path}")

    npz_path = OUT_DIR / "hmm_states.npz"
    np.savez(npz_path,
             means=mus, stds=sigmas,
             transition_matrix=A, stationary=pi,
             viterbi_v1=states_v1, viterbi_v2=states_v2)
    print(f"NPZ saved:     {npz_path}")

    # --- Figure ---
    fig = plt.figure(figsize=(15, 12))
    gs = fig.add_gridspec(4, 3, height_ratios=[1.4, 1.4, 1.0, 1.0],
                          hspace=0.6, wspace=0.35)
    ax_v1 = fig.add_subplot(gs[0, :])
    ax_v2 = fig.add_subplot(gs[1, :])
    ax_A = fig.add_subplot(gs[2, 0])
    ax_pi = fig.add_subplot(gs[2, 1])
    ax_emiss = fig.add_subplot(gs[2, 2])
    ax_dt_bc = fig.add_subplot(gs[3, 0])
    ax_dt_im = fig.add_subplot(gs[3, 1])
    ax_dt_ec = fig.add_subplot(gs[3, 2])

    state_colors = ["#d62728", "#ff7f0e", "#2ca02c"]

    def plot_video(ax, x, states, label, cps_idx):
        T = x.size
        for k in range(N_STATES):
            mask = states == k
            ax.scatter(np.arange(T)[mask], x[mask],
                       s=8, c=state_colors[k], alpha=0.6,
                       label=state_labels[k])
        # HMM transitions
        for tr in state_transitions(states):
            ax.axvline(tr, color="black", lw=0.5, alpha=0.4)
        # obj-047 BinSeg+BIC change-points (vertical red dashed)
        for cp in cps_idx:
            ax.axvline(cp, color="red", ls="--", lw=1.4, alpha=0.85)
        ax.set_ylabel("CV0 (Å)")
        ax.set_xlim(0, T - 1)
        ax.set_title(f"{label} ({T} frames) — Viterbi state path "
                     "(red = obj-047 BinSeg+BIC change-points)")
        ax.legend(loc="upper right", fontsize=8, ncol=3, framealpha=0.9)
        ax.grid(alpha=0.3)
        ax.set_xlabel("frame")

    cps_v1_idx = co_summary["v1"]["binseg_change_points"] if co_summary["v1"] else []
    cps_v2_idx = co_summary["v2"]["binseg_change_points"] if co_summary["v2"] else []
    plot_video(ax_v1, cv0_v1, states_v1, "V1", cps_v1_idx)
    plot_video(ax_v2, cv0_v2, states_v2, "V2", cps_v2_idx)

    # Transition matrix heatmap
    im = ax_A.imshow(A, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax_A.set_xticks(range(N_STATES)); ax_A.set_yticks(range(N_STATES))
    ax_A.set_xticklabels(state_labels, fontsize=8)
    ax_A.set_yticklabels(state_labels, fontsize=8)
    for i in range(N_STATES):
        for j in range(N_STATES):
            ax_A.text(j, i, f"{A[i, j]:.3f}", ha="center", va="center",
                      fontsize=8, color="black" if A[i, j] < 0.5 else "white")
    ax_A.set_title("Transition matrix A[from, to]", fontsize=10)
    fig.colorbar(im, ax=ax_A, fraction=0.05)

    # Stationary distribution
    ax_pi.bar(range(N_STATES), pi, color=state_colors)
    ax_pi.set_xticks(range(N_STATES))
    ax_pi.set_xticklabels(state_labels, fontsize=8)
    ax_pi.set_ylim(0, 1)
    for i, p in enumerate(pi):
        ax_pi.text(i, p + 0.02, f"{p:.3f}", ha="center", fontsize=8)
    ax_pi.set_title("Stationary distribution π")
    ax_pi.set_ylabel("probability")
    ax_pi.grid(axis="y", alpha=0.3)

    # Emission Gaussians overlaid on CV0 histogram
    ax_emiss.hist(cv0_all, bins=80, density=True, alpha=0.4, color="gray",
                  label="empirical (V1+V2)")
    grid = np.linspace(cv0_all.min() - 5, cv0_all.max() + 5, 600)
    for k in range(N_STATES):
        pdf = pi[k] * np.exp(gaussian_log_pdf(grid, mus[k], sigmas[k]))
        ax_emiss.plot(grid, pdf, color=state_colors[k], lw=2,
                      label=f"{state_labels[k]} N({mus[k]:.1f},{sigmas[k]:.2f})")
    ax_emiss.set_xlabel("CV0 (Å)")
    ax_emiss.set_ylabel("density")
    ax_emiss.set_title("Emission distribution + π-weighted Gaussians")
    ax_emiss.legend(fontsize=7)
    ax_emiss.grid(alpha=0.3)

    # Dwell-time histograms
    for ax, k, lab in [(ax_dt_bc, 0, "BC"), (ax_dt_im, 1, "Intermediate"), (ax_dt_ec, 2, "EC")]:
        d = dwell_all[k]
        if not d:
            ax.text(0.5, 0.5, "no visits", ha="center", transform=ax.transAxes)
            ax.axis("off")
            continue
        bins = np.linspace(1, max(d) + 1, min(20, max(d) + 1))
        ax.hist(d, bins=bins, alpha=0.55, color=state_colors[k], edgecolor="black",
                density=True, label="empirical")
        mean = float(np.mean(d))
        ts = np.linspace(1, max(d), 200)
        ax.plot(ts, np.exp(-ts / mean) / mean, "k--", lw=1.5,
                label=f"exp fit, μ={mean:.1f} s")
        ax.set_xlabel("dwell time (s @ 1 fps)")
        ax.set_ylabel("density")
        ax.set_title(f"{lab} dwell-time (n={len(d)} visits)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    fig.suptitle(
        "obj-048 — 3-state Gaussian-emission HMM over per-frame CV0 (V1+V2)\n"
        f"BC μ={mus[0]:.1f} Å, Intermediate μ={mus[1]:.1f} Å, EC μ={mus[2]:.1f} Å  •  "
        f"log-L={fit['log_likelihood']:.0f}  •  "
        f"π = ({pi[0]:.2f}, {pi[1]:.2f}, {pi[2]:.2f})",
        fontsize=11.5, fontweight="bold", y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig_path = FIG_DIR / "hmm_state_assignment.png"
    fig.savefig(fig_path, dpi=140)
    print(f"Figure saved: {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
