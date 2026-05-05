#!/usr/bin/env python3
"""V1 vs V2 breakpoint / dwell-time cross-validation (obj-050).

Audit-2026-05-05 §17.6 v9 P=4 follow-up to obj-047 + obj-048. The
question: are the V1 and V2 transition statistics drawn from the same
underlying stochastic process, or do the two HS-AFM acquisitions sample
distinguishably different kinetics?

Two passes:
  (1) Direct KS test on the obj-047 BinSeg+BIC inter-breakpoint
      intervals (V1: 4 intervals from 3 change-points + 2 boundaries;
      V2: 5 intervals from 4 change-points + 2 boundaries). Tiny
      samples but the simplest test.
  (2) Per-state KS test on the HMM Viterbi dwell-times (obj-048):
      richer (n ≈ 15-50 per state per video) and stratified by state.

Both passes report:
  - empirical KS statistic D and asymptotic p
  - bootstrap-permutation p-value (5000 perms): under H0 the V1/V2
    labels are exchangeable; resample labels and recompute D
  - bootstrap CI on the mean ratio V1/V2

Outputs:
  figures/breakpoint_cross_validation.png
  results/afm_pipeline/free_energy_profile/breakpoint_cross_validation.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
CP_SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "change_points.json"
HMM_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "hmm_states.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

STATE_LABELS = ["BC", "Intermediate", "EC"]
STATE_COLORS = ["#d62728", "#ff7f0e", "#2ca02c"]
N_STATES = 3
N_BOOTSTRAP = 5000


def inter_breakpoint_intervals(cps: list[int], n_frames: int) -> np.ndarray:
    """Return frame-counts between consecutive boundaries (start, cps..., end)."""
    boundaries = [0, *sorted(cps), n_frames]
    return np.array([b - a for a, b in zip(boundaries[:-1], boundaries[1:])], dtype=float)


def dwell_runs(states: np.ndarray) -> dict:
    """Per-state list of consecutive-frame run lengths (no censoring filter)."""
    out = {k: [] for k in range(N_STATES)}
    if states.size == 0:
        return out
    cur = int(states[0])
    n = 1
    for s in states[1:]:
        if int(s) == cur:
            n += 1
        else:
            out[cur].append(n)
            cur = int(s)
            n = 1
    out[cur].append(n)
    return out


def bootstrap_ks_p(a: np.ndarray, b: np.ndarray, n_bootstrap: int = N_BOOTSTRAP,
                   rng: np.random.Generator | None = None) -> tuple[float, float]:
    """Permutation null for the two-sample KS statistic.

    Pool a + b, randomly relabel into pieces of size |a| and |b|, recompute D.
    Returns (observed D, p-value)."""
    if rng is None:
        rng = np.random.default_rng(0)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan"), float("nan")
    observed_D, _ = stats.ks_2samp(a, b)
    pooled = np.concatenate([a, b])
    n_a = a.size
    null_D = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        rng.shuffle(pooled)
        x = pooled[:n_a]
        y = pooled[n_a:]
        d, _ = stats.ks_2samp(x, y)
        null_D[i] = d
    p = float(np.mean(null_D >= observed_D))
    return float(observed_D), p


def bootstrap_mean_ratio(a: np.ndarray, b: np.ndarray, n_bootstrap: int = N_BOOTSTRAP,
                         rng: np.random.Generator | None = None) -> dict:
    if rng is None:
        rng = np.random.default_rng(1)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return {"observed_ratio": float("nan"), "ci95": [float("nan"), float("nan")]}
    obs = float(a.mean() / b.mean())
    ratios = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        ai = rng.choice(a, size=a.size, replace=True)
        bi = rng.choice(b, size=b.size, replace=True)
        ratios[i] = ai.mean() / max(bi.mean(), 1e-12)
    lo, hi = np.percentile(ratios, [2.5, 97.5])
    return {
        "observed_ratio": obs,
        "ci95": [float(lo), float(hi)],
        "v1_mean": float(a.mean()),
        "v2_mean": float(b.mean()),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # --- Pass 1: obj-047 inter-breakpoint intervals ---
    cp_data = json.loads(CP_SRC.read_text())
    v1_info = cp_data["results"]["v1_frame_cv0_smoothed"]
    v2_info = cp_data["results"]["v2_frame_cv0_smoothed"]
    cps_v1 = list(v1_info["binseg_change_points_idx"])
    cps_v2 = list(v2_info["binseg_change_points_idx"])
    n_v1 = int(v1_info["n"])
    n_v2 = int(v2_info["n"])
    intervals_v1 = inter_breakpoint_intervals(cps_v1, n_v1)
    intervals_v2 = inter_breakpoint_intervals(cps_v2, n_v2)
    print(f"V1 inter-breakpoint intervals: {intervals_v1.tolist()}")
    print(f"V2 inter-breakpoint intervals: {intervals_v2.tolist()}")

    D_seg, p_seg = bootstrap_ks_p(intervals_v1, intervals_v2,
                                   rng=np.random.default_rng(42))
    ratio_seg = bootstrap_mean_ratio(intervals_v1, intervals_v2,
                                     rng=np.random.default_rng(43))
    print(f"\n=== Pass 1: obj-047 inter-breakpoint intervals ===")
    print(f"  V1 mean = {intervals_v1.mean():.1f} s,  V2 mean = {intervals_v2.mean():.1f} s")
    print(f"  KS D = {D_seg:.4f}, permutation p = {p_seg:.4f}  "
          f"({'reject H0' if p_seg < 0.05 else 'cannot reject H0 (same distribution)'})")
    print(f"  V1/V2 mean ratio = {ratio_seg['observed_ratio']:.3f}  "
          f"(95 % CI [{ratio_seg['ci95'][0]:.3f}, {ratio_seg['ci95'][1]:.3f}])")

    # --- Pass 2: per-state HMM dwell-times ---
    hmm = np.load(HMM_NPZ)
    states_v1 = hmm["viterbi_v1"]
    states_v2 = hmm["viterbi_v2"]
    dw_v1 = dwell_runs(states_v1)
    dw_v2 = dwell_runs(states_v2)

    per_state_results: dict = {}
    print("\n=== Pass 2: HMM Viterbi dwell-times per state ===")
    for k, lab in enumerate(STATE_LABELS):
        a = np.array(dw_v1[k], dtype=float)
        b = np.array(dw_v2[k], dtype=float)
        if a.size < 2 or b.size < 2:
            print(f"  {lab}: insufficient runs (V1 n={a.size}, V2 n={b.size})")
            per_state_results[lab] = {
                "v1_n": int(a.size), "v2_n": int(b.size),
                "insufficient_data": True,
            }
            continue
        D, p = bootstrap_ks_p(a, b, rng=np.random.default_rng(100 + k))
        rat = bootstrap_mean_ratio(a, b, rng=np.random.default_rng(200 + k))
        ks_async = stats.ks_2samp(a, b)
        per_state_results[lab] = {
            "v1_n": int(a.size),
            "v2_n": int(b.size),
            "v1_mean_s": float(a.mean()),
            "v2_mean_s": float(b.mean()),
            "ks_D": float(D),
            "ks_perm_pvalue": float(p),
            "ks_asymptotic_pvalue": float(ks_async.pvalue),
            "v1_v2_mean_ratio": rat["observed_ratio"],
            "v1_v2_mean_ratio_ci95": rat["ci95"],
        }
        verdict = ("reject H0" if p < 0.05
                   else "cannot reject H0 (same distribution)")
        print(f"  {lab} (V1 n={a.size}, V2 n={b.size}):")
        print(f"    means: V1={a.mean():.2f} s, V2={b.mean():.2f} s,"
              f"  ratio={rat['observed_ratio']:.3f}  CI [{rat['ci95'][0]:.3f}, {rat['ci95'][1]:.3f}]")
        print(f"    KS D={D:.4f}, perm p={p:.4f}  ({verdict})")

    # --- Save JSON ---
    summary = {
        "T_kelvin": 300.0,
        "frame_to_second": 1.0,
        "n_bootstrap": N_BOOTSTRAP,
        "obj047_intervals": {
            "v1_intervals_s": intervals_v1.tolist(),
            "v2_intervals_s": intervals_v2.tolist(),
            "v1_n": int(intervals_v1.size),
            "v2_n": int(intervals_v2.size),
            "v1_mean_s": float(intervals_v1.mean()),
            "v2_mean_s": float(intervals_v2.mean()),
            "ks_D": float(D_seg),
            "ks_perm_pvalue": float(p_seg),
            "v1_v2_mean_ratio": ratio_seg["observed_ratio"],
            "v1_v2_mean_ratio_ci95": ratio_seg["ci95"],
        },
        "obj048_dwell_per_state": per_state_results,
    }
    json_path = OUT_DIR / "breakpoint_cross_validation.json"
    with json_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nJSON saved: {json_path}")

    # --- Figure ---
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 4, hspace=0.4, wspace=0.35)

    # Top-left 2 panels: pass 1 (obj-047 intervals)
    ax_p1_hist = fig.add_subplot(gs[0, 0])
    ax_p1_ecdf = fig.add_subplot(gs[0, 1])

    bins = np.linspace(0, max(intervals_v1.max(), intervals_v2.max()) + 10, 12)
    ax_p1_hist.hist(intervals_v1, bins=bins, alpha=0.55, color="tab:cyan",
                    edgecolor="black", label=f"V1 (n={intervals_v1.size})")
    ax_p1_hist.hist(intervals_v2, bins=bins, alpha=0.55, color="tab:olive",
                    edgecolor="black", label=f"V2 (n={intervals_v2.size})")
    ax_p1_hist.set_xlabel("inter-breakpoint interval (s)")
    ax_p1_hist.set_ylabel("count")
    ax_p1_hist.set_title("Pass 1: obj-047 intervals")
    ax_p1_hist.legend(fontsize=8)
    ax_p1_hist.grid(alpha=0.3)

    for arr, color, lab in [(intervals_v1, "tab:cyan", "V1"),
                            (intervals_v2, "tab:olive", "V2")]:
        sx = np.sort(arr)
        sy = np.arange(1, sx.size + 1) / sx.size
        ax_p1_ecdf.step(sx, sy, where="post", color=color, lw=2, label=lab)
    ax_p1_ecdf.set_xlabel("inter-breakpoint interval (s)")
    ax_p1_ecdf.set_ylabel("ECDF")
    ax_p1_ecdf.set_title(f"Pass 1: KS D={D_seg:.3f}, p={p_seg:.3f}")
    ax_p1_ecdf.legend(fontsize=8)
    ax_p1_ecdf.grid(alpha=0.3)

    # Top-right 2 panels: ratio and summary
    ax_ratio = fig.add_subplot(gs[0, 2])
    ax_summary = fig.add_subplot(gs[0, 3])

    labels = ["obj-047\nintervals"] + STATE_LABELS
    ratios = [ratio_seg["observed_ratio"]]
    cis = [ratio_seg["ci95"]]
    for lab in STATE_LABELS:
        if per_state_results[lab].get("insufficient_data"):
            ratios.append(np.nan); cis.append([np.nan, np.nan])
        else:
            ratios.append(per_state_results[lab]["v1_v2_mean_ratio"])
            cis.append(per_state_results[lab]["v1_v2_mean_ratio_ci95"])
    cis = np.array(cis)
    xs = np.arange(len(labels))
    ax_ratio.errorbar(xs, ratios,
                      yerr=[np.array(ratios) - cis[:, 0], cis[:, 1] - np.array(ratios)],
                      fmt="o", color="black", capsize=4, lw=1.5)
    ax_ratio.axhline(1.0, color="red", ls="--", lw=1, alpha=0.7,
                     label="ratio = 1 (identical means)")
    ax_ratio.set_xticks(xs)
    ax_ratio.set_xticklabels(labels, fontsize=8)
    ax_ratio.set_ylabel("V1 / V2 mean ratio (95 % CI)")
    ax_ratio.set_title("V1/V2 mean ratio with bootstrap CI")
    ax_ratio.legend(fontsize=8)
    ax_ratio.grid(alpha=0.3)

    # Summary table
    ax_summary.axis("off")
    rows = [["sample", "V1 n", "V2 n", "KS D", "perm p", "ratio"]]
    rows.append(["obj-047 intervals",
                 str(intervals_v1.size), str(intervals_v2.size),
                 f"{D_seg:.3f}", f"{p_seg:.3f}",
                 f"{ratio_seg['observed_ratio']:.2f}"])
    for lab in STATE_LABELS:
        r = per_state_results[lab]
        if r.get("insufficient_data"):
            rows.append([lab, str(r["v1_n"]), str(r["v2_n"]), "—", "—", "—"])
        else:
            rows.append([lab, str(r["v1_n"]), str(r["v2_n"]),
                         f"{r['ks_D']:.3f}", f"{r['ks_perm_pvalue']:.3f}",
                         f"{r['v1_v2_mean_ratio']:.2f}"])
    table = ax_summary.table(cellText=rows[1:], colLabels=rows[0],
                             cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.05, 1.5)
    ax_summary.set_title("Summary table", fontsize=10)

    # Bottom row: per-state ECDF V1 vs V2
    for k, lab in enumerate(STATE_LABELS):
        ax = fig.add_subplot(gs[1, k])
        a = np.array(dw_v1[k], dtype=float)
        b = np.array(dw_v2[k], dtype=float)
        for arr, color, name in [(a, "tab:cyan", "V1"), (b, "tab:olive", "V2")]:
            if arr.size == 0:
                continue
            sx = np.sort(arr)
            sy = np.arange(1, sx.size + 1) / sx.size
            ax.step(sx, sy, where="post", color=color, lw=2, label=f"{name} (n={sx.size})")
        r = per_state_results[lab]
        if r.get("insufficient_data"):
            title = f"{lab} — insufficient data"
        else:
            title = f"{lab} — KS D={r['ks_D']:.3f}, p={r['ks_perm_pvalue']:.3f}"
        ax.set_title(title)
        ax.set_xlabel("dwell time (s)")
        ax.set_ylabel("ECDF")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Bottom-rightmost panel: forest plot of V1/V2 ratios for at-a-glance
    ax_forest = fig.add_subplot(gs[1, 3])
    ax_forest.axis("off")
    text_lines = ["V1 vs V2 mean-ratio summary",
                  "(95 % CI overlapping 1 ⇒ same)",
                  ""]
    for lab, r in [("obj-047 intervals", {"v1_v2_mean_ratio": ratio_seg["observed_ratio"],
                                          "v1_v2_mean_ratio_ci95": ratio_seg["ci95"]})] + \
                   [(lab, per_state_results[lab]) for lab in STATE_LABELS]:
        if r.get("insufficient_data"):
            text_lines.append(f"  {lab}:  —")
            continue
        ci = r["v1_v2_mean_ratio_ci95"]
        sym = "✓" if ci[0] <= 1 <= ci[1] else "✗"
        text_lines.append(f"  {lab}:  {r['v1_v2_mean_ratio']:.2f}  "
                          f"[{ci[0]:.2f}, {ci[1]:.2f}]  {sym}")
    ax_forest.text(0.05, 0.95, "\n".join(text_lines), transform=ax_forest.transAxes,
                   fontsize=10, va="top", family="monospace")

    fig.suptitle(
        "obj-050 — V1 vs V2 breakpoint / dwell-time cross-validation\n"
        "Top: obj-047 inter-breakpoint intervals.  "
        "Bottom: per-state HMM dwell-times.  "
        f"Bootstrap permutation null, n_perm = {N_BOOTSTRAP}",
        fontsize=11.5, fontweight="bold", y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig_path = FIG_DIR / "breakpoint_cross_validation.png"
    fig.savefig(fig_path, dpi=140)
    print(f"Figure saved: {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
