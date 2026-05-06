#!/usr/bin/env python3
"""HMM rate-matrix consolidation across 1-D / 2-D / 3-D fits (obj-061).

Audit §30 v21 P=4 follow-up. obj-048 (1-D), obj-055 (2-D K=3 + K=4),
and obj-059 (3-D K=3) all share the same V1+V2 1645-frame trace at
1 fps. Each fit produced its own transition matrix A. This pass
compares rate matrices Q = (1/Δt)·ln(A) and slow-mode timescales
across all four HMMs.

For a Markov chain with transition matrix A at time step Δt:
  Q = ln(A) / Δt
  eigenvalues λ_i of A → timescales τ_i = -Δt / ln|Re(λ_i)|
  Largest eigenvalue is 1 (equilibrium).
  Slowest non-trivial mode corresponds to longest τ.

Cross-check: obj-049 mean dwell times (BC 14.6 s, Inter 6.9 s, EC
11.4 s) and obj-054 CV0 ACF e-folding (V1 5.1 s, V2 7.8 s).

Output: figures/hmm_rate_matrix_consolidation.png + JSON.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import logm

ROOT = Path(__file__).resolve().parents[3]
SRC_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"
OUT_JSON = SRC_DIR / "hmm_rate_matrix_consolidation.json"
FIG = FIG_DIR / "hmm_rate_matrix_consolidation.png"

DT_S = 1.0  # 1 fps HS-AFM


def load_a_matrices() -> dict:
    matrices: dict[str, dict] = {}
    with (SRC_DIR / "hmm_states.json").open() as f:
        d = json.load(f)
    matrices["1D K=3 (obj-048)"] = {
        "A": np.array(d["transition_matrix"]),
        "labels": d.get("state_labels",
                        [f"BC", "Inter", "EC"]),
        "means": d["means_A"],
        "stds": d["stds_A"],
        "dim": 1,
        "K": 3,
    }
    with (SRC_DIR / "hmm_2d_cv0_cv2.json").open() as f:
        d = json.load(f)
    for K in (3, 4):
        if str(K) in d.get("results_by_K", {}):
            r = d["results_by_K"][str(K)]
            matrices[f"2D K={K} (obj-055)"] = {
                "A": np.array(r["A"]),
                "labels": [f"s{k}" for k in range(K)],
                "means": r["mus"],
                "dim": 2,
                "K": K,
            }
    with (SRC_DIR / "hmm_3d_cv0_cv1_cv2.json").open() as f:
        d = json.load(f)
    for K in (3,):
        if str(K) in d.get("results_by_K", {}):
            r = d["results_by_K"][str(K)]
            matrices[f"3D K={K} (obj-059)"] = {
                "A": np.array(r["A"]),
                "labels": [f"s{k}" for k in range(K)],
                "means": r["mus"],
                "dim": 3,
                "K": K,
            }
    return matrices


def rate_matrix(A: np.ndarray, dt: float = DT_S) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    log_A = np.asarray(logm(A))
    Q_complex = log_A / dt
    Q = np.real(Q_complex)
    eigvals = np.linalg.eigvals(A)
    timescales = []
    for lam in eigvals:
        rl = float(np.real(lam))
        if rl >= 1.0 - 1e-9:
            timescales.append(np.inf)
        elif rl <= 0.0:
            timescales.append(np.nan)
        else:
            timescales.append(-dt / np.log(rl))
    return Q, np.array(eigvals), np.array(timescales)


def main() -> int:
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    FIG.parent.mkdir(parents=True, exist_ok=True)

    matrices = load_a_matrices()
    summary: dict[str, dict] = {}

    print("=== HMM rate-matrix consolidation ===\n")
    for name, m in matrices.items():
        A = m["A"]
        Q, eigvals, timescales = rate_matrix(A, DT_S)
        order = np.argsort(-np.real(eigvals))
        eigvals_sorted = eigvals[order]
        timescales_sorted = timescales[order]
        print(f"--- {name} ---")
        print("  A:\n", A)
        print(f"  eigenvalues: {[f'{ev:.4f}' for ev in eigvals_sorted]}")
        print(f"  timescales (s): "
              f"{[(f'{t:.2f}' if np.isfinite(t) else 'inf') for t in timescales_sorted]}")
        slow_finite = [t for t in timescales_sorted if np.isfinite(t)]
        slow_mode = max(slow_finite) if slow_finite else float("nan")
        print(f"  slowest non-trivial mode: τ = {slow_mode:.2f} s\n")

        summary[name] = {
            "dim": m["dim"],
            "K": m["K"],
            "A": A.tolist(),
            "Q_rate_matrix": Q.tolist(),
            "eigenvalues_complex": [(float(np.real(ev)), float(np.imag(ev)))
                                    for ev in eigvals_sorted],
            "timescales_s": [float(t) if np.isfinite(t) else None
                             for t in timescales_sorted],
            "slowest_nontrivial_timescale_s": float(slow_mode)
                if np.isfinite(slow_mode) else None,
        }

    benchmarks = {
        "obj-049 BC mean dwell s": 14.6,
        "obj-049 Inter mean dwell s": 6.9,
        "obj-049 EC mean dwell s": 11.4,
        "obj-054 V1 CV0 ACF τ_e s": 5.1,
        "obj-054 V2 CV0 ACF τ_e s": 7.8,
    }
    summary["__benchmarks_obj049_obj054"] = benchmarks
    summary["__params"] = {"dt_s": DT_S, "method": "Q = logm(A)/dt; τ_i = -dt/log(Re(λ_i))"}

    with OUT_JSON.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"JSON saved: {OUT_JSON}")

    plot(matrices, summary, benchmarks)
    print(f"Figure saved: {FIG}")
    return 0


def plot(matrices: dict, summary: dict, benchmarks: dict) -> None:
    n = len(matrices)
    fig = plt.figure(figsize=(15.5, 4 + 3.0 * n / 2))
    gs = fig.add_gridspec(n, 4, width_ratios=[1.0, 1.0, 1.2, 1.5],
                          hspace=0.55, wspace=0.45)

    cmap_blues = "Blues"
    cmap_diverging = "RdBu_r"

    for i, (name, m) in enumerate(matrices.items()):
        A = m["A"]
        K = A.shape[0]
        s = summary[name]
        Q = np.array(s["Q_rate_matrix"])

        ax = fig.add_subplot(gs[i, 0])
        im = ax.imshow(A, cmap=cmap_blues, vmin=0, vmax=1, aspect="auto")
        for r in range(K):
            for c in range(K):
                ax.text(c, r, f"{A[r, c]:.3f}", ha="center", va="center",
                        fontsize=8, color="black" if A[r, c] < 0.5 else "white")
        ax.set_xticks(range(K)); ax.set_yticks(range(K))
        ax.set_xticklabels([f"s{k}" for k in range(K)], fontsize=8)
        ax.set_yticklabels([f"s{k}" for k in range(K)], fontsize=8)
        ax.set_title(f"{name}\nA[from, to]", fontsize=9, fontweight="bold")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        ax = fig.add_subplot(gs[i, 1])
        Q_off = Q - np.diag(np.diag(Q))
        vmax = max(abs(Q_off.max()), abs(Q_off.min()), 1e-3)
        im = ax.imshow(Q_off, cmap=cmap_diverging, vmin=-vmax, vmax=vmax,
                       aspect="auto")
        for r in range(K):
            for c in range(K):
                if r != c:
                    ax.text(c, r, f"{Q[r, c]:+.3f}", ha="center", va="center",
                            fontsize=8, color="black")
        ax.set_xticks(range(K)); ax.set_yticks(range(K))
        ax.set_xticklabels([f"s{k}" for k in range(K)], fontsize=8)
        ax.set_yticklabels([f"s{k}" for k in range(K)], fontsize=8)
        ax.set_title(f"Q rate matrix\n(off-diag, s⁻¹)", fontsize=9,
                     fontweight="bold")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        ax = fig.add_subplot(gs[i, 2])
        eigvals = [complex(*ev) for ev in s["eigenvalues_complex"]]
        ax.scatter([np.real(e) for e in eigvals],
                   [np.imag(e) for e in eigvals],
                   s=140, color="#d62728", edgecolor="black", linewidth=0.8,
                   zorder=10)
        for j, e in enumerate(eigvals):
            ax.annotate(f"  λ{j} = {np.real(e):.3f}{'+' + str(np.imag(e)) + 'i' if abs(np.imag(e)) > 1e-6 else ''}",
                        (np.real(e), np.imag(e)),
                        fontsize=7, ha="left", va="center")
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(theta), np.sin(theta), color="#888",
                linewidth=0.8, alpha=0.5, label="unit circle")
        ax.axhline(0, color="black", linewidth=0.6, alpha=0.4)
        ax.axvline(0, color="black", linewidth=0.6, alpha=0.4)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.5, 0.5)
        ax.set_aspect("equal")
        ax.set_xlabel("Re(λ)")
        ax.set_ylabel("Im(λ)")
        ax.set_title("Eigenvalues of A", fontsize=9, fontweight="bold")
        ax.grid(alpha=0.3)

        ax = fig.add_subplot(gs[i, 3])
        timescales = [t for t in s["timescales_s"] if t is not None]
        if timescales:
            order = np.argsort(timescales)[::-1]
            sorted_t = [timescales[j] for j in order]
            ax.barh(np.arange(len(sorted_t)), sorted_t,
                    color=["#1f77b4" if t == max(sorted_t) else "#7570b3"
                           for t in sorted_t], edgecolor="black")
            for j, t in enumerate(sorted_t):
                ax.text(t + 0.5, j, f"{t:.2f} s", va="center", fontsize=8)
            ax.set_yticks(range(len(sorted_t)))
            ax.set_yticklabels([f"mode {j}" for j in range(len(sorted_t))],
                               fontsize=8)
        for label, val in benchmarks.items():
            ax.axvline(val, color="#fd8d3c", linestyle=":", linewidth=1.0,
                       alpha=0.6)
        ax.set_xlim(left=0, right=max([t for t in timescales] + list(benchmarks.values())) * 1.4)
        ax.set_xlabel("τ_i (s)")
        ax.set_title("Implied timescales", fontsize=9, fontweight="bold")
        ax.grid(alpha=0.3, axis="x")
        if i == 0:
            for label, val in benchmarks.items():
                ax.text(val, len(timescales) - 0.3, f" {label.split()[0]}",
                        rotation=90, fontsize=6.5, ha="left", va="top",
                        color="#fd8d3c")

    fig.suptitle(
        "HMM rate-matrix consolidation (1-D + 2-D + 3-D fits, obj-061) — "
        "Q = ln(A)/Δt at Δt = 1 s; benchmarks: obj-049 dwells (5-15 s) + obj-054 ACF τ (5-8 s)",
        fontsize=11, fontweight="bold", y=0.995,
    )
    fig.savefig(FIG, dpi=140, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
