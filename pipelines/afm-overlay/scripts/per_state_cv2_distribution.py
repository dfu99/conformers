#!/usr/bin/env python3
"""Per-state CV2 distribution × HMM state — audit §33 v24 / obj-064.

obj-063 stratified per-residue RMSF by HMM state. This pass asks the
complementary question: under the obj-048 1-D-CV0-only HMM, does the
state label carry information about CV2 (the orthogonal headpiece-
opening axis)?

Method:
  1. Load CV0/CV1/CV2 traces from cv_series.npz (V1+V2, n=1645).
  2. Load HMM Viterbi labels from hmm_states.npz (states 0/1/2 = BC/
     Inter/EC under obj-048's CV0-only 3-state HMM).
  3. Per state s: collect CV2 frame distribution.
  4. Compute mean / std / median / Q1 / Q3 per state.
  5. Pairwise KS test on CV2 distributions across states.
  6. R² (between-state-variance / total CV2 variance) — fraction of
     CV2 variance explained by 1-D HMM state assignment.
  7. Repeat for CV1 as a second sanity check (CV1 ≈ CV0, so per-state
     CV1 distributions should differ strongly).

Output:
  - figures/per_state_cv2_distribution.png
  - results/afm_pipeline/free_energy_profile/per_state_cv2_distribution.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
CV_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
HMM_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "hmm_states.npz"
OUT_JSON = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "per_state_cv2_distribution.json"
OUT_FIG = ROOT / "figures" / "per_state_cv2_distribution.png"

STATE_NAMES = ["BC", "Inter", "EC"]
STATE_COLORS = ["#4a90e2", "#f5a623", "#d0021b"]


def per_state_summary(cv: np.ndarray, states: np.ndarray) -> dict:
    out = {}
    overall_mean = float(cv.mean())
    overall_var = float(cv.var(ddof=1))
    n_total = cv.size
    between = 0.0
    for s in (0, 1, 2):
        sub = cv[states == s]
        n_s = sub.size
        between += n_s * (sub.mean() - overall_mean) ** 2
        out[f"{STATE_NAMES[s]}"] = {
            "n": int(n_s),
            "mean": float(sub.mean()),
            "std": float(sub.std(ddof=1)),
            "median": float(np.median(sub)),
            "q1": float(np.quantile(sub, 0.25)),
            "q3": float(np.quantile(sub, 0.75)),
            "min": float(sub.min()),
            "max": float(sub.max()),
        }
    between_var = between / (n_total - 1)
    out["overall_mean"] = overall_mean
    out["overall_var"] = overall_var
    out["between_state_var"] = float(between_var)
    out["R2_state_explains"] = float(between_var / overall_var)
    return out


def pairwise_ks(cv: np.ndarray, states: np.ndarray) -> dict:
    out = {}
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        a = cv[states == i]
        b = cv[states == j]
        result = stats.ks_2samp(a, b)
        out[f"{STATE_NAMES[i]}_vs_{STATE_NAMES[j]}"] = {
            "ks_D": float(result.statistic),
            "p": float(result.pvalue),
            "n_a": int(a.size),
            "n_b": int(b.size),
        }
    return out


def main() -> int:
    print("Loading CV series + HMM states ...")
    cvs = np.load(CV_NPZ)
    hmm = np.load(HMM_NPZ)

    cv0 = np.concatenate([cvs["cv0_v1"], cvs["cv0_v2"]])
    cv1 = np.concatenate([cvs["cv1_v1"], cvs["cv1_v2"]])
    cv2 = np.concatenate([cvs["cv2_v1"], cvs["cv2_v2"]])
    states = np.concatenate([hmm["viterbi_v1"], hmm["viterbi_v2"]]).astype(int)
    assert cv0.size == cv1.size == cv2.size == states.size

    print(f"n={cv0.size} frames; states {dict(zip(*np.unique(states, return_counts=True)))}")

    # CV0 sanity (the variable that defined the HMM)
    cv0_summary = per_state_summary(cv0, states)
    cv1_summary = per_state_summary(cv1, states)
    cv2_summary = per_state_summary(cv2, states)
    cv0_ks = pairwise_ks(cv0, states)
    cv1_ks = pairwise_ks(cv1, states)
    cv2_ks = pairwise_ks(cv2, states)

    print(f"CV0 R² state explains: {cv0_summary['R2_state_explains']:.4f}")
    print(f"CV1 R² state explains: {cv1_summary['R2_state_explains']:.4f}")
    print(f"CV2 R² state explains: {cv2_summary['R2_state_explains']:.4f}")

    print("\nCV2 per-state summary:")
    for s in (0, 1, 2):
        d = cv2_summary[STATE_NAMES[s]]
        print(f"  {STATE_NAMES[s]:>5s} (n={d['n']:>4d}): "
              f"mean {d['mean']:6.2f} ± {d['std']:5.2f} Å, "
              f"median {d['median']:.2f}, IQR [{d['q1']:.2f}, {d['q3']:.2f}]")

    print("\nCV2 pairwise KS:")
    for k, v in cv2_ks.items():
        print(f"  {k}: D = {v['ks_D']:.4f}, p = {v['p']:.2e}")

    # Bonferroni across 3 pairwise tests on each CV → 9 tests
    bonf_alpha = 0.05 / 9
    cv2_n_sig = sum(1 for v in cv2_ks.values() if v["p"] < bonf_alpha)
    print(f"\nCV2 Bonferroni-sig pairs: {cv2_n_sig} / 3 (α_corr = {bonf_alpha:.4f})")

    out = {
        "n_frames": int(cv0.size),
        "n_BC": int((states == 0).sum()),
        "n_Inter": int((states == 1).sum()),
        "n_EC": int((states == 2).sum()),
        "T_kelvin": float(cvs["T_kelvin"]),
        "cv0": {"per_state": cv0_summary, "pairwise_ks": cv0_ks},
        "cv1": {"per_state": cv1_summary, "pairwise_ks": cv1_ks},
        "cv2": {"per_state": cv2_summary, "pairwise_ks": cv2_ks},
        "bonferroni_alpha_per_cv": bonf_alpha,
        "summary": {
            "cv0_r2": cv0_summary["R2_state_explains"],
            "cv1_r2": cv1_summary["R2_state_explains"],
            "cv2_r2": cv2_summary["R2_state_explains"],
            "cv0_n_sig_pairs": sum(1 for v in cv0_ks.values() if v["p"] < bonf_alpha),
            "cv1_n_sig_pairs": sum(1 for v in cv1_ks.values() if v["p"] < bonf_alpha),
            "cv2_n_sig_pairs": cv2_n_sig,
        },
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSON.open("w") as f:
        json.dump(out, f, indent=2)
    print(f"\nsaved {OUT_JSON}")

    plot(cv0, cv1, cv2, states, cv0_summary, cv1_summary, cv2_summary,
         cv0_ks, cv1_ks, cv2_ks, bonf_alpha)
    return 0


def _hist_panel(ax, cv, states, summary, ks, bonf_alpha,
                cv_label, x_range, bins, title_suffix):
    for s in range(3):
        sub = cv[states == s]
        ax.hist(sub, bins=bins, range=x_range, density=True,
                color=STATE_COLORS[s], alpha=0.45,
                label=f"{STATE_NAMES[s]} (n={summary[STATE_NAMES[s]]['n']}, "
                      f"μ={summary[STATE_NAMES[s]]['mean']:.1f})")
    ax.set_xlabel(f"{cv_label} (Å)")
    ax.set_ylabel("density")
    n_sig = sum(1 for v in ks.values() if v["p"] < bonf_alpha)
    title = (f"{cv_label} per state — R² state explains "
             f"{summary['R2_state_explains']:.3f} "
             f"({n_sig}/3 pairs Bonferroni-sig){title_suffix}")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right", framealpha=0.95)
    ax.grid(alpha=0.3)


def plot(cv0, cv1, cv2, states, cv0_summary, cv1_summary, cv2_summary,
         cv0_ks, cv1_ks, cv2_ks, bonf_alpha):
    fig = plt.figure(figsize=(15.5, 10.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 1.0],
                          hspace=0.55, wspace=0.30)

    # Row 1: CV0/CV1/CV2 per-state histograms
    ax = fig.add_subplot(gs[0, 0])
    _hist_panel(ax, cv0, states, cv0_summary, cv0_ks, bonf_alpha,
                "CV0", (40, 100), bins=40,
                title_suffix="\n(HMM-defining axis — sanity)")
    ax = fig.add_subplot(gs[0, 1])
    _hist_panel(ax, cv1, states, cv1_summary, cv1_ks, bonf_alpha,
                "CV1", (10, 110), bins=40,
                title_suffix="")
    ax = fig.add_subplot(gs[0, 2])
    _hist_panel(ax, cv2, states, cv2_summary, cv2_ks, bonf_alpha,
                "CV2", (5, 55), bins=40,
                title_suffix="\n(headpiece-opening axis)")

    # Row 2: per-state CV0 vs CV2 scatter, colored by state — one
    # panel per video and combined
    cv0_v1 = cv0[: 379]
    cv0_v2 = cv0[379:]
    cv2_v1 = cv2[: 379]
    cv2_v2 = cv2[379:]
    states_v1 = states[: 379]
    states_v2 = states[379:]

    for j, (label, _cv0, _cv2, _states) in enumerate([
        ("V1 379 frames", cv0_v1, cv2_v1, states_v1),
        ("V2 1266 frames", cv0_v2, cv2_v2, states_v2),
        ("V1+V2 1645 frames", cv0, cv2, states),
    ]):
        ax = fig.add_subplot(gs[1, j])
        for s in range(3):
            mask = _states == s
            ax.scatter(_cv0[mask], _cv2[mask], s=10, alpha=0.55,
                       color=STATE_COLORS[s], edgecolor="none",
                       label=f"{STATE_NAMES[s]} (n={int(mask.sum())})")
        ax.axhline(50.0, color="#666", linestyle="--", linewidth=0.7,
                   label="EO threshold CV2 ≥ 50")
        ax.set_xlabel("CV0 (head-tail dist, Å)")
        ax.set_ylabel("CV2 (headpiece opening, Å)")
        ax.set_xlim(40, 100)
        ax.set_ylim(5, 55)
        ax.set_title(f"({chr(ord('a') + j)}) {label} — "
                     f"(CV0, CV2) by HMM state",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, loc="lower right", framealpha=0.92)
        ax.grid(alpha=0.3)

    # Row 3 (single span): summary table
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")

    rows = [["axis", "BC mean ± std", "Inter mean ± std", "EC mean ± std",
             "R²", "BC↔Inter p", "Inter↔EC p", "BC↔EC p", "sig pairs"]]
    for cv_name, summary, ks in [
        ("CV0", cv0_summary, cv0_ks),
        ("CV1", cv1_summary, cv1_ks),
        ("CV2", cv2_summary, cv2_ks),
    ]:
        rows.append([
            cv_name,
            f"{summary['BC']['mean']:.2f} ± {summary['BC']['std']:.2f}",
            f"{summary['Inter']['mean']:.2f} ± {summary['Inter']['std']:.2f}",
            f"{summary['EC']['mean']:.2f} ± {summary['EC']['std']:.2f}",
            f"{summary['R2_state_explains']:.3f}",
            f"{ks['BC_vs_Inter']['p']:.1e}",
            f"{ks['Inter_vs_EC']['p']:.1e}",
            f"{ks['BC_vs_EC']['p']:.1e}",
            f"{sum(1 for v in ks.values() if v['p'] < bonf_alpha)}/3",
        ])
    table = ax.table(cellText=rows, loc="center", cellLoc="center",
                     colWidths=[0.06, 0.13, 0.13, 0.13,
                                 0.06, 0.10, 0.10, 0.10, 0.07])
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 1.40)
    for k in range(len(rows[0])):
        table[(0, k)].set_text_props(weight="bold", color="white")
        table[(0, k)].set_facecolor("#444444")
    # Highlight CV2 row in green if its R² > 0.05 — i.e., 1-D HMM does
    # carry CV2 information
    cv2_r2 = cv2_summary["R2_state_explains"]
    if cv2_r2 > 0.05:
        for k in range(len(rows[0])):
            table[(3, k)].set_facecolor("#e8f4ea")

    ax.set_title(
        f"Per-state CV summary — does the obj-048 1-D-CV0-only HMM carry CV2 "
        f"(headpiece-opening) information?\n"
        f"CV0 R² = {cv0_summary['R2_state_explains']:.3f} (definitional); "
        f"CV1 R² = {cv1_summary['R2_state_explains']:.3f} "
        f"(propeller-A correlated; cf. obj-054 r ≥ 0.84); "
        f"CV2 R² = {cv2_summary['R2_state_explains']:.3f} "
        f"({'green' if cv2_r2 > 0.05 else 'no'}: state ≠ irrelevant for CV2)",
        fontsize=10, fontweight="bold")

    fig.suptitle(
        f"Per-state CV2 distribution × obj-048 HMM (obj-064) — does the 1-D "
        f"CV0-only state carry headpiece-opening info?\n"
        f"V1+V2 fitted (n=1645); {STATE_NAMES[0]} {(states == 0).sum()} / "
        f"{STATE_NAMES[1]} {(states == 1).sum()} / "
        f"{STATE_NAMES[2]} {(states == 2).sum()}",
        fontsize=11, fontweight="bold", y=0.995,
    )

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=140, bbox_inches="tight")
    print(f"saved {OUT_FIG}")


if __name__ == "__main__":
    raise SystemExit(main())
