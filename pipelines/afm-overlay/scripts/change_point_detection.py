#!/usr/bin/env python3
"""Change-point detection on FES-min and CV0 trajectories (obj-047).

Audit-2026-05-05 §16.7 follow-up. obj-046 confirmed that the FES drift
identified in obj-045 is intrinsic biology, not metadata-driven. This
v9 deepening pass localizes the breakpoints in time using:

  (1) CUSUM: cumulative sum of standardized residuals against the
      sequence mean, max |CUSUM| pinpoints a single change.
  (2) Binary segmentation with BIC penalty: greedy recursive split
      using sum-of-squared-errors cost, BIC stops at the right
      number of change-points.

We run both detectors on:
  - per-block FES-min CV0 trajectories (V1: n=7, V2: n=25)
  - per-frame CV0 trajectories (V1: 379 frames, V2: 1266 frames)
    smoothed by a 25-frame moving average (~1/2 of obj-045 block
    size) so per-frame jitter does not dominate.

Outputs:
  figures/change_point_detection.png
  results/afm_pipeline/free_energy_profile/change_points.json
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SRC_BLOCK = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "per_block_dg.npz"
SRC_FRAME = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

WINDOW = 50  # obj-045 block size in frames; for axis labelling
SMOOTH = 25  # moving-average window for per-frame CV0


def cusum(x: np.ndarray) -> tuple[np.ndarray, int, float]:
    """Cumulative sum of (x - mean(x)). Returns (curve, argmax, max_abs)."""
    c = np.cumsum(x - np.mean(x))
    k = int(np.argmax(np.abs(c)))
    return c, k, float(np.abs(c[k]))


def cusum_pvalue(x: np.ndarray, n_bootstrap: int = 5000, rng: np.random.Generator | None = None) -> float:
    """Bootstrap p-value for CUSUM peak under permutation null."""
    if rng is None:
        rng = np.random.default_rng(0)
    _, _, observed = cusum(x)
    null = np.empty(n_bootstrap)
    x_local = np.asarray(x).copy()
    for i in range(n_bootstrap):
        rng.shuffle(x_local)
        _, _, peak = cusum(x_local)
        null[i] = peak
    p = float(np.mean(null >= observed))
    return p


def cost_sse(x: np.ndarray) -> float:
    """Sum of squared errors against segment mean."""
    if x.size == 0:
        return 0.0
    return float(np.sum((x - np.mean(x)) ** 2))


def best_split(x: np.ndarray, min_size: int = 2) -> tuple[int | None, float]:
    """Best single split point minimizing total within-segment SSE."""
    n = x.size
    if n < 2 * min_size:
        return None, float("inf")
    base = cost_sse(x)
    best_k, best_cost = None, base
    for k in range(min_size, n - min_size + 1):
        c = cost_sse(x[:k]) + cost_sse(x[k:])
        if c < best_cost:
            best_cost = c
            best_k = k
    return best_k, best_cost


def bic_score(x: np.ndarray, change_points: Iterable[int]) -> float:
    """BIC = n*log(SSE/n) + k*log(n) for k change-points (k+1 segments)."""
    cps = sorted(set(change_points))
    n = x.size
    boundaries = [0, *cps, n]
    sse_total = 0.0
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        sse_total += cost_sse(x[a:b])
    if sse_total <= 1e-12:
        return -np.inf
    k = len(cps)
    return n * np.log(sse_total / n) + k * np.log(n)


def binseg_bic(x: np.ndarray, max_cps: int = 4, min_size: int = 2) -> list[int]:
    """Binary segmentation with BIC stopping rule."""
    n = x.size
    cps: list[int] = []
    current_bic = bic_score(x, cps)
    for _ in range(max_cps):
        cps_sorted = sorted(cps)
        boundaries = [0, *cps_sorted, n]
        best_overall_k, best_overall_gain = None, 0.0
        for a, b in zip(boundaries[:-1], boundaries[1:]):
            seg = x[a:b]
            k_local, _ = best_split(seg, min_size=min_size)
            if k_local is None:
                continue
            k_global = a + k_local
            cps_trial = sorted([*cps, k_global])
            new_bic = bic_score(x, cps_trial)
            gain = current_bic - new_bic
            if gain > best_overall_gain:
                best_overall_gain = gain
                best_overall_k = k_global
        if best_overall_k is None:
            break
        cps.append(best_overall_k)
        current_bic = bic_score(x, sorted(cps))
    return sorted(cps)


def segment_means(x: np.ndarray, cps: list[int]) -> list[tuple[int, int, float]]:
    n = x.size
    boundaries = [0, *sorted(cps), n]
    out = []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        out.append((a, b, float(np.mean(x[a:b]))))
    return out


def moving_average(x: np.ndarray, w: int) -> np.ndarray:
    if w <= 1:
        return x.astype(float)
    pad = w // 2
    cum = np.cumsum(np.concatenate([[0.0], x.astype(float)]))
    out = (cum[w:] - cum[:-w]) / w
    return np.concatenate([np.full(pad, out[0]), out, np.full(x.size - out.size - pad, out[-1])])


def detect_all(label: str, x: np.ndarray, x_axis: np.ndarray, *,
               max_cps: int = 4, n_bootstrap: int = 5000) -> dict:
    """Run CUSUM + BinSeg+BIC on a sequence, return summary dict."""
    cusum_curve, cusum_k, cusum_peak = cusum(x)
    cusum_p = cusum_pvalue(x, n_bootstrap=n_bootstrap, rng=np.random.default_rng(42))
    cps = binseg_bic(x, max_cps=max_cps, min_size=max(2, x.size // 25))
    segs = segment_means(x, cps)
    return {
        "label": label,
        "n": int(x.size),
        "x": x.tolist(),
        "x_axis": x_axis.tolist(),
        "cusum_curve": cusum_curve.tolist(),
        "cusum_argmax": int(cusum_k),
        "cusum_argmax_x": float(x_axis[cusum_k]),
        "cusum_peak_abs": cusum_peak,
        "cusum_pvalue": cusum_p,
        "binseg_change_points_idx": list(cps),
        "binseg_change_points_x": [float(x_axis[k]) for k in cps],
        "segments": [
            {"start_idx": int(a), "end_idx": int(b), "mean": m,
             "start_x": float(x_axis[a]), "end_x": float(x_axis[min(b - 1, x.size - 1)])}
            for (a, b, m) in segs
        ],
        "max_jump_abs_A": (
            float(max(abs(segs[i + 1][2] - segs[i][2]) for i in range(len(segs) - 1)))
            if len(segs) > 1 else 0.0
        ),
    }


def plot_detection(ax_seq, ax_cusum, label: str, result: dict, *, ylab: str,
                   shade_color: str = "tab:blue") -> None:
    x = np.asarray(result["x"])
    x_axis = np.asarray(result["x_axis"])
    cusum_curve = np.asarray(result["cusum_curve"])
    cps_x = result["binseg_change_points_x"]

    # Top: sequence with segment means + change-point lines
    ax_seq.plot(x_axis, x, "o-", lw=1.5, ms=4, color=shade_color, label="data")
    for seg in result["segments"]:
        a, b = seg["start_idx"], seg["end_idx"]
        m = seg["mean"]
        if b > x.size:
            b = x.size
        if a >= b:
            continue
        ax_seq.hlines(m, x_axis[a], x_axis[min(b - 1, x.size - 1)],
                      color="black", lw=2.5, alpha=0.85)
    for cx in cps_x:
        ax_seq.axvline(cx, color="red", ls="--", lw=1.5, alpha=0.85)
    if cps_x:
        ax_seq.legend([f"change-points at x = {', '.join(f'{cx:.0f}' for cx in cps_x)}"],
                      loc="best", fontsize=8)
    ax_seq.set_ylabel(ylab)
    ax_seq.set_title(f"{label} — BinSeg+BIC: {len(cps_x)} change-point"
                     + ("s" if len(cps_x) != 1 else ""))
    ax_seq.grid(alpha=0.3)

    # Bottom: CUSUM curve
    ax_cusum.plot(x_axis, cusum_curve, color="tab:orange", lw=1.5)
    ax_cusum.axhline(0, color="grey", lw=0.5)
    peak_x = result["cusum_argmax_x"]
    ax_cusum.axvline(peak_x, color="red", ls="--", lw=1.2, alpha=0.7)
    ax_cusum.set_ylabel("CUSUM")
    ax_cusum.set_xlabel("frame")
    ax_cusum.set_title(f"CUSUM peak at x = {peak_x:.0f}, |peak| = {result['cusum_peak_abs']:.2f},"
                       f" perm. p = {result['cusum_pvalue']:.3g}")
    ax_cusum.grid(alpha=0.3)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    block = np.load(SRC_BLOCK)
    frame = np.load(SRC_FRAME)

    centers_v1 = np.asarray(block["centers_v1"], dtype=float)
    centers_v2 = np.asarray(block["centers_v2"], dtype=float)
    min_loc_v1 = np.asarray(block["min_loc_v1"], dtype=float)
    min_loc_v2 = np.asarray(block["min_loc_v2"], dtype=float)

    cv0_v1 = np.asarray(frame["cv0_v1"], dtype=float)
    cv0_v2 = np.asarray(frame["cv0_v2"], dtype=float)
    cv0_v1_smoothed = moving_average(cv0_v1, SMOOTH)
    cv0_v2_smoothed = moving_average(cv0_v2, SMOOTH)
    frames_v1 = np.arange(cv0_v1.size, dtype=float)
    frames_v2 = np.arange(cv0_v2.size, dtype=float)

    print("Running detectors...")
    res_block_v1 = detect_all("V1 per-block FES-min", min_loc_v1, centers_v1, max_cps=3)
    res_block_v2 = detect_all("V2 per-block FES-min", min_loc_v2, centers_v2, max_cps=4)
    res_frame_v1 = detect_all("V1 per-frame CV0 (25-MA)", cv0_v1_smoothed, frames_v1, max_cps=3,
                              n_bootstrap=2000)
    res_frame_v2 = detect_all("V2 per-frame CV0 (25-MA)", cv0_v2_smoothed, frames_v2, max_cps=4,
                              n_bootstrap=2000)

    summary = {
        "T_kelvin": 300.0,
        "block_window_frames": WINDOW,
        "smoothing_window_frames": SMOOTH,
        "results": {
            "v1_block_fes_min": res_block_v1,
            "v2_block_fes_min": res_block_v2,
            "v1_frame_cv0_smoothed": res_frame_v1,
            "v2_frame_cv0_smoothed": res_frame_v2,
        },
    }

    # Compact human-readable digest
    print("\n=== CUSUM peak locations ===")
    for r in [res_block_v1, res_block_v2, res_frame_v1, res_frame_v2]:
        print(f"  {r['label']}: peak at frame ≈ {r['cusum_argmax_x']:.1f},"
              f"  |peak| = {r['cusum_peak_abs']:.2f},  perm. p = {r['cusum_pvalue']:.4g}")

    print("\n=== BinSeg+BIC change-points ===")
    for r in [res_block_v1, res_block_v2, res_frame_v1, res_frame_v2]:
        cps_x = r['binseg_change_points_x']
        if not cps_x:
            print(f"  {r['label']}: 0 change-points")
        else:
            print(f"  {r['label']}: {len(cps_x)} change-point(s) at x ≈ "
                  + ", ".join(f"{cx:.1f}" for cx in cps_x)
                  + f"  (max segment-mean jump = {r['max_jump_abs_A']:.2f} Å)")

    # ---- Figure ----
    fig, axes = plt.subplots(4, 2, figsize=(14, 13.5),
                             gridspec_kw={"height_ratios": [1.6, 1.0, 1.6, 1.0]})
    plot_detection(axes[0, 0], axes[1, 0], "V1 per-block FES-min CV0",
                   res_block_v1, ylab="FES-min CV0 (Å)",
                   shade_color="tab:blue")
    plot_detection(axes[0, 1], axes[1, 1], "V2 per-block FES-min CV0",
                   res_block_v2, ylab="FES-min CV0 (Å)",
                   shade_color="tab:purple")
    plot_detection(axes[2, 0], axes[3, 0], "V1 per-frame CV0 (25-MA)",
                   res_frame_v1, ylab="CV0 (Å), MA-25",
                   shade_color="tab:cyan")
    plot_detection(axes[2, 1], axes[3, 1], "V2 per-frame CV0 (25-MA)",
                   res_frame_v2, ylab="CV0 (Å), MA-25",
                   shade_color="tab:olive")

    for ax in axes[1:2, :].flat:
        ax.set_xlabel("frame center")
    for ax in axes[3:4, :].flat:
        ax.set_xlabel("frame")

    fig.suptitle(
        "obj-047 — change-point detection on FES-min + CV0 trajectories\n"
        "BinSeg+BIC change-points (red dashed) + segment means (black) + CUSUM curves",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig_path = FIG_DIR / "change_point_detection.png"
    fig.savefig(fig_path, dpi=140)
    print(f"\nFigure saved: {fig_path}")

    json_path = OUT_DIR / "change_points.json"
    with json_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"JSON saved:    {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
