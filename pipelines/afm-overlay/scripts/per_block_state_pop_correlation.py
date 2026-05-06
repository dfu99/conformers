#!/usr/bin/env python3
"""Per-block HMM state-population × FES-drift cross-tabulation (obj-062).

Audit §31 v22 P=4 follow-up. obj-046 already showed metadata signals
(fit correlation, head jitter, head drift, GIF intensity) do NOT
explain the per-block FES drift from obj-045 (0/10 hypothesis tests
reach p < 0.05). The natural follow-up: do per-block HMM Viterbi
state populations explain it instead?

Method:
  1. Load obj-048 HMM Viterbi state assignments (V1 379 frames, V2
     1266 frames).
  2. Load obj-045 per-block FES-min CV0 (V1 7 blocks of 50 frames,
     V2 25 blocks).
  3. Per block: empirical f_BC, f_Inter, f_EC from Viterbi.
  4. Per block: take obj-045 FES-min CV0.
  5. Pearson r between FES-min CV0 and each f_state, per video.
  6. Bonferroni correction across (3 states × 2 videos = 6 tests).
  7. Diagnostic: does state-population variation explain >50 % of
     FES-drift variance?

Output:
  figures/per_block_state_pop_correlation.png
  results/afm_pipeline/free_energy_profile/per_block_state_pop_correlation.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[3]
SRC_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG = ROOT / "figures" / "per_block_state_pop_correlation.png"
OUT_JSON = SRC_DIR / "per_block_state_pop_correlation.json"

WINDOW_FRAMES = 50
STATE_NAMES = ["BC", "Inter", "EC"]
STATE_COLORS = ["#d62728", "#ff7f0e", "#2ca02c"]


def per_block_populations(viterbi: np.ndarray, window: int) -> np.ndarray:
    n = len(viterbi)
    n_blocks = n // window
    pops = np.zeros((n_blocks, 3))
    for b in range(n_blocks):
        block = viterbi[b * window: (b + 1) * window]
        for k in range(3):
            pops[b, k] = (block == k).mean()
    return pops


def main() -> int:
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    FIG.parent.mkdir(parents=True, exist_ok=True)

    hmm = np.load(SRC_DIR / "hmm_states.npz")
    viterbi_v1 = hmm["viterbi_v1"]
    viterbi_v2 = hmm["viterbi_v2"]
    print(f"V1: {len(viterbi_v1)} frames; V2: {len(viterbi_v2)} frames")

    fes = np.load(SRC_DIR / "per_block_dg.npz")
    fes_min_v1 = fes["min_loc_v1"]
    fes_min_v2 = fes["min_loc_v2"]
    print(f"V1: {len(fes_min_v1)} FES blocks; V2: {len(fes_min_v2)} FES blocks")

    pops_v1 = per_block_populations(viterbi_v1, WINDOW_FRAMES)
    pops_v2 = per_block_populations(viterbi_v2, WINDOW_FRAMES)
    print(f"V1 populations: {pops_v1.shape}; V2 populations: {pops_v2.shape}")

    n_v1 = min(len(fes_min_v1), pops_v1.shape[0])
    n_v2 = min(len(fes_min_v2), pops_v2.shape[0])
    fes_min_v1 = fes_min_v1[:n_v1]
    pops_v1 = pops_v1[:n_v1]
    fes_min_v2 = fes_min_v2[:n_v2]
    pops_v2 = pops_v2[:n_v2]
    print(f"After alignment — V1: {n_v1} blocks, V2: {n_v2} blocks")

    results: dict[str, dict] = {}
    for video, pops, fes_min, n in [("V1", pops_v1, fes_min_v1, n_v1),
                                     ("V2", pops_v2, fes_min_v2, n_v2)]:
        for k, sn in enumerate(STATE_NAMES):
            r, p = pearsonr(pops[:, k], fes_min)
            r2 = r ** 2
            results[f"{video}_{sn}"] = {
                "video": video,
                "state": sn,
                "n_blocks": int(n),
                "pearson_r": float(r),
                "pearson_p": float(p),
                "r_squared": float(r2),
                "f_state_per_block": pops[:, k].tolist(),
                "fes_min_cv0_per_block": fes_min.tolist(),
            }

    n_tests = 6
    p_bonferroni_threshold = 0.05 / n_tests
    print(f"\nBonferroni-corrected α = {p_bonferroni_threshold:.4f} (6 tests)\n")
    print(f"{'video':<5} {'state':<7} {'n':<3} {'r':<8} {'r²':<7} {'p':<10} {'sig?'}")
    print("-" * 60)
    for key, d in results.items():
        sig_b = d["pearson_p"] < p_bonferroni_threshold
        sig_n = d["pearson_p"] < 0.05
        flag = "B" if sig_b else ("•" if sig_n else "")
        print(f"{d['video']:<5} {d['state']:<7} {d['n_blocks']:<3} "
              f"{d['pearson_r']:+.3f}   {d['r_squared']:.3f}   "
              f"{d['pearson_p']:.4f}    {flag}")

    n_b = sum(1 for d in results.values() if d["pearson_p"] < p_bonferroni_threshold)
    n_n = sum(1 for d in results.values() if d["pearson_p"] < 0.05)
    print(f"\n  n significant at α=0.05 (nominal): {n_n}/{n_tests}")
    print(f"  n significant after Bonferroni:    {n_b}/{n_tests}")

    max_r2 = max(d["r_squared"] for d in results.values())
    avg_r2 = np.mean([d["r_squared"] for d in results.values()])
    if max_r2 > 0.50 and n_b >= 1:
        verdict = ("STRONG — at least one HMM state population explains >50% of "
                   "FES-drift variance with Bonferroni-significant correlation. "
                   "FES drift IS reflected in state-level repopulation, not "
                   "intrinsic state stability.")
    elif max_r2 > 0.25 and n_n >= 2:
        verdict = ("MODERATE — HMM state populations explain 25-50% of FES-drift "
                   "variance with multiple nominally-significant correlations. "
                   "State repopulation contributes to FES drift but doesn't fully "
                   "explain it.")
    elif n_n >= 1:
        verdict = ("WEAK — only nominal significance, ≤25% variance explained. "
                   "FES drift is partially state-related.")
    else:
        verdict = ("NULL — no significant state-population × FES-drift correlation. "
                   "The FES drift may reflect intrinsic state-level free-energy "
                   "fluctuation rather than state repopulation. Combined with "
                   "obj-046, this means the FES drift has neither metadata nor "
                   "HMM-state cause — it is genuinely intrinsic biology.")

    summary = {
        "results": results,
        "verdict": verdict,
        "summary_stats": {
            "n_tests": n_tests,
            "n_sig_nominal_alpha_005": n_n,
            "n_sig_bonferroni_alpha_005": n_b,
            "max_r_squared": float(max_r2),
            "mean_r_squared": float(avg_r2),
            "bonferroni_threshold": float(p_bonferroni_threshold),
        },
        "interpretation": (
            f"Cross-tabulation of obj-048 HMM Viterbi state populations "
            f"with obj-045 FES-drift per block. Verdict: {verdict[:120]}"
        ),
        "context": {
            "obj_045_v1_drift_std_kcal": 4.08,
            "obj_045_v2_drift_std_kcal": 4.42,
            "obj_046_metadata_sig_count": 0,
            "obj_046_metadata_n_tests": 10,
        },
    }

    print(f"\n=== VERDICT: {verdict} ===\n")
    with OUT_JSON.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"JSON saved: {OUT_JSON}")

    plot(results, summary, pops_v1, pops_v2, fes_min_v1, fes_min_v2)
    print(f"Figure saved: {FIG}")
    return 0


def plot(results: dict, summary: dict, pops_v1: np.ndarray, pops_v2: np.ndarray,
         fes_min_v1: np.ndarray, fes_min_v2: np.ndarray) -> None:
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.85], hspace=0.42, wspace=0.36)

    for vi, (video, pops, fes_min) in enumerate([
        ("V1", pops_v1, fes_min_v1),
        ("V2", pops_v2, fes_min_v2),
    ]):
        for k in range(3):
            ax = fig.add_subplot(gs[0, vi * 2 + (k // 3)] if False else gs[vi, k])
            r = results[f"{video}_{STATE_NAMES[k]}"]
            ax.scatter(pops[:, k], fes_min, s=70, color=STATE_COLORS[k],
                       edgecolor="black", linewidth=0.6)
            if len(pops[:, k]) > 1:
                slope = np.polyfit(pops[:, k], fes_min, 1)
                xs = np.linspace(pops[:, k].min(), pops[:, k].max(), 50)
                ax.plot(xs, np.polyval(slope, xs), color=STATE_COLORS[k],
                        linewidth=1.5, alpha=0.6)
            ax.set_xlabel(f"f_{STATE_NAMES[k]} per block")
            ax.set_ylabel(f"FES-min CV0 (Å)")
            ax.set_title(f"{video} — {STATE_NAMES[k]} pop. vs FES-min CV0",
                         fontsize=10)
            ax.grid(alpha=0.3)
            sig_b = r["pearson_p"] < summary["summary_stats"]["bonferroni_threshold"]
            sig_n = r["pearson_p"] < 0.05
            flag = "B" if sig_b else ("•" if sig_n else "ns")
            color_box = ("#ccebc5" if sig_b
                         else "#fed9a6" if sig_n
                         else "#fbb4ae")
            ax.text(0.02, 0.97,
                    f"r = {r['pearson_r']:+.3f}\nr² = {r['r_squared']:.3f}\n"
                    f"p = {r['pearson_p']:.4f} ({flag})",
                    transform=ax.transAxes, va="top", ha="left", fontsize=8.5,
                    bbox=dict(boxstyle="round", facecolor=color_box,
                              edgecolor="#666"))

    ax = fig.add_subplot(gs[0, 3])
    ax.axis("off")
    rows = [["video", "state", "r", "r²", "p", "sig"]]
    for key, d in results.items():
        sig_b = d["pearson_p"] < summary["summary_stats"]["bonferroni_threshold"]
        sig_n = d["pearson_p"] < 0.05
        flag = "B" if sig_b else ("•" if sig_n else "ns")
        rows.append([d["video"], d["state"],
                     f"{d['pearson_r']:+.3f}", f"{d['r_squared']:.3f}",
                     f"{d['pearson_p']:.4f}", flag])
    table = ax.table(cellText=rows, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.4)
    ax.set_title("Per-test summary (B = Bonferroni-sig, • = nominal-sig, ns = NS)",
                 fontsize=9, fontweight="bold")

    ax = fig.add_subplot(gs[1, 3])
    ax.axis("off")
    summary_lines = [
        f"Tests: {summary['summary_stats']['n_tests']}",
        f"Bonferroni α: {summary['summary_stats']['bonferroni_threshold']:.4f}",
        "",
        f"Nominal sig (p<0.05): {summary['summary_stats']['n_sig_nominal_alpha_005']}",
        f"Bonferroni sig:       {summary['summary_stats']['n_sig_bonferroni_alpha_005']}",
        f"Max r²:  {summary['summary_stats']['max_r_squared']:.3f}",
        f"Mean r²: {summary['summary_stats']['mean_r_squared']:.3f}",
        "",
        "obj-045 FES drift std:",
        "  V1: 4.08 kcal/mol",
        "  V2: 4.42 kcal/mol",
        "obj-046 metadata sig:",
        "  0/10 reach p<0.05",
    ]
    ax.text(0.05, 0.95, "\n".join(summary_lines),
            transform=ax.transAxes, fontsize=10, va="top", ha="left",
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="#fff3cd",
                      edgecolor="#856404"))
    ax.set_title("Summary stats", fontsize=10, fontweight="bold")

    fig.suptitle(
        "Per-block HMM state-pop. × FES-drift correlation (obj-062) — "
        f"{summary['summary_stats']['n_sig_nominal_alpha_005']}/6 nominally sig, "
        f"{summary['summary_stats']['n_sig_bonferroni_alpha_005']}/6 Bonferroni-sig\n"
        f"Verdict: {summary['verdict'][:140]}",
        fontsize=10, fontweight="bold", y=0.995,
    )
    fig.savefig(FIG, dpi=140, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
