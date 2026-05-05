#!/usr/bin/env python3
"""State-population v2: windowed by HS-AFM time (obj-044).

Audit-2026-05-05 §13.5 P=3 follow-up to obj-043. The pooled v1 analysis
reports a single (BC, Inter, EC, EO*) breakdown, but it implicitly
assumes stationarity across the HS-AFM acquisition timeline. This v2
splits each video into 50-frame blocks (sliding non-overlapping
windows) and reports per-block fractions with Wilson 95 % CI.

Falsifying signal:
  * any block's BC/Inter/EC/EO fraction sits more than 2 standard
    errors outside the pooled-mean expectation under the
    stationary-Bernoulli null hypothesis. We aggregate by reporting
    the maximum z-score over all (block, state) pairs.

Inputs: results/afm_pipeline/free_energy_profile/cv_series.npz
Outputs:
  figures/state_populations_v2_windowed.png
  results/afm_pipeline/free_energy_profile/state_populations_windowed.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

WINDOW = 50  # frames per block

BC_MAX = 65.0
EC_MIN = 78.0
EO_MIN = 85.0

STATE_BANDS = [
    ("BC", 0.0, BC_MAX, "#7f7f7f"),
    ("Inter.", BC_MAX, EC_MIN, "#fee08b"),
    ("EC", EC_MIN, EO_MIN, "#fdae61"),
    ("EO*", EO_MIN, 1e3, "#d73027"),
]


def state_label(cv0: float) -> int:
    """Return integer index (0..3) into STATE_BANDS."""
    for i, (_, lo, hi, _) in enumerate(STATE_BANDS):
        if lo <= cv0 < hi:
            return i
    return len(STATE_BANDS) - 1


def windowed_fractions(cv0: np.ndarray, w: int = WINDOW) -> dict:
    """Return per-block state fractions + Wilson CI."""
    n = cv0.size
    n_blocks = max(1, n // w)
    frac = np.zeros((n_blocks, 4), dtype=float)
    ci_lo = np.zeros((n_blocks, 4), dtype=float)
    ci_hi = np.zeros((n_blocks, 4), dtype=float)
    centers = np.zeros(n_blocks, dtype=int)
    for b in range(n_blocks):
        sl = slice(b * w, min((b + 1) * w, n))
        block = cv0[sl]
        m = block.size
        centers[b] = (sl.start + sl.stop) // 2
        for s, (_, lo, hi, _) in enumerate(STATE_BANDS):
            k = int(((block >= lo) & (block < hi)).sum())
            p = k / m
            # Wilson 95 % CI (Beta(0.5+k, 0.5+m-k) Jeffreys form).
            ci_lo[b, s] = float(beta.ppf(0.025, 0.5 + k, 0.5 + m - k))
            ci_hi[b, s] = float(beta.ppf(0.975, 0.5 + k, 0.5 + m - k))
            frac[b, s] = p
    return dict(centers=centers, frac=frac, ci_lo=ci_lo, ci_hi=ci_hi,
                n_blocks=int(n_blocks), w=int(w))


def stationarity_z(frac: np.ndarray, pooled: np.ndarray, w: int) -> np.ndarray:
    """Per-(block, state) z-score under stationary-Bernoulli null."""
    n_blocks, n_state = frac.shape
    z = np.zeros_like(frac)
    for s in range(n_state):
        p = pooled[s]
        se = np.sqrt(p * (1 - p) / w)
        if se < 1e-12:
            z[:, s] = 0
        else:
            z[:, s] = (frac[:, s] - p) / se
    return z


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    data = np.load(SRC)
    cv0_v1 = np.asarray(data["cv0_v1"])
    cv0_v2 = np.asarray(data["cv0_v2"])
    print(f"V1 frames: {cv0_v1.size}  V2 frames: {cv0_v2.size}")

    pooled_v1 = np.array([
        ((cv0_v1 >= lo) & (cv0_v1 < hi)).mean()
        for _, lo, hi, _ in STATE_BANDS
    ])
    pooled_v2 = np.array([
        ((cv0_v2 >= lo) & (cv0_v2 < hi)).mean()
        for _, lo, hi, _ in STATE_BANDS
    ])

    print(f"Pooled V1: {dict(zip([s[0] for s in STATE_BANDS], np.round(pooled_v1, 3)))}")
    print(f"Pooled V2: {dict(zip([s[0] for s in STATE_BANDS], np.round(pooled_v2, 3)))}")

    out_v1 = windowed_fractions(cv0_v1, WINDOW)
    out_v2 = windowed_fractions(cv0_v2, WINDOW)
    print(f"V1 blocks: {out_v1['n_blocks']}  V2 blocks: {out_v2['n_blocks']}")

    z_v1 = stationarity_z(out_v1["frac"], pooled_v1, WINDOW)
    z_v2 = stationarity_z(out_v2["frac"], pooled_v2, WINDOW)
    z_max_v1 = float(np.abs(z_v1).max())
    z_max_v2 = float(np.abs(z_v2).max())
    n_outlier_v1 = int((np.abs(z_v1) > 2.0).sum())
    n_outlier_v2 = int((np.abs(z_v2) > 2.0).sum())
    n_total_v1 = int(z_v1.size)
    n_total_v2 = int(z_v2.size)
    print()
    print(f"V1 max |z| = {z_max_v1:.2f}  ({n_outlier_v1}/{n_total_v1} block-state pairs > 2 SE)")
    print(f"V2 max |z| = {z_max_v2:.2f}  ({n_outlier_v2}/{n_total_v2} block-state pairs > 2 SE)")
    # Bonferroni-adjusted z-threshold (normal-approx, two-sided)
    from scipy.stats import norm
    z_threshold = float(norm.ppf(1 - 0.05 / (2 * (z_v1.size + z_v2.size))))
    print(f"Bonferroni-corrected |z| threshold: {z_threshold:.2f}")
    n_significant_v1 = int((np.abs(z_v1) > z_threshold).sum())
    n_significant_v2 = int((np.abs(z_v2) > z_threshold).sum())
    print(f"V1 significant after Bonferroni: {n_significant_v1}/{n_total_v1}")
    print(f"V2 significant after Bonferroni: {n_significant_v2}/{n_total_v2}")

    # --- Save JSON ---
    out_json = {
        "window_frames": WINDOW,
        "state_bands": [
            {"name": s[0], "cv0_lo": s[1], "cv0_hi": s[2]}
            for s in STATE_BANDS
        ],
        "v1": {
            "n_frames": int(cv0_v1.size),
            "n_blocks": int(out_v1["n_blocks"]),
            "pooled_fractions": pooled_v1.tolist(),
            "max_abs_z": z_max_v1,
            "n_blocks_state_outside_2se": n_outlier_v1,
            "n_blocks_state_total": n_total_v1,
            "n_significant_bonferroni": n_significant_v1,
        },
        "v2": {
            "n_frames": int(cv0_v2.size),
            "n_blocks": int(out_v2["n_blocks"]),
            "pooled_fractions": pooled_v2.tolist(),
            "max_abs_z": z_max_v2,
            "n_blocks_state_outside_2se": n_outlier_v2,
            "n_blocks_state_total": n_total_v2,
            "n_significant_bonferroni": n_significant_v2,
        },
        "bonferroni_z_threshold": z_threshold,
        "stationarity_verdict": (
            "stationary"
            if (n_significant_v1 + n_significant_v2) == 0
            else "non_stationary_after_bonferroni"
        ),
    }
    out_path = OUT_DIR / "state_populations_windowed.json"
    out_path.write_text(json.dumps(out_json, indent=2))
    print(f"\nSaved {out_path}")

    # --- Plot ---
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5),
                              sharex="col", sharey="row",
                              gridspec_kw=dict(hspace=0.35, wspace=0.20,
                                               width_ratios=[1, 1]))

    for col, (out, z, label) in enumerate([
        (out_v1, z_v1, "V1 (379 frames)"),
        (out_v2, z_v2, "V2 (1266 frames)"),
    ]):
        ax_top = axes[0, col]
        ax_bot = axes[1, col]
        centers = out["centers"]
        frac = out["frac"]

        # Top: stacked area of state fractions per block
        ax_top.stackplot(
            centers,
            frac[:, 0], frac[:, 1], frac[:, 2], frac[:, 3],
            colors=[s[3] for s in STATE_BANDS],
            labels=[s[0] for s in STATE_BANDS],
            alpha=0.8,
        )
        ax_top.set_xlim(centers.min(), centers.max())
        ax_top.set_ylim(0, 1)
        ax_top.set_title(f"{label} — windowed state fractions ({WINDOW}-frame blocks)",
                         fontsize=10.5)
        ax_top.set_ylabel("cumulative fraction")
        if col == 1:
            ax_top.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8)
        ax_top.grid(alpha=0.3)

        # Bottom: per-block z-score per state
        for s, (name, _, _, color) in enumerate(STATE_BANDS):
            ax_bot.plot(centers, z[:, s], color=color, lw=1.2,
                        marker="o", ms=3, label=name)
        ax_bot.axhline(0, color="black", lw=0.6)
        ax_bot.axhline(2, color="red", ls="--", lw=0.7, alpha=0.6)
        ax_bot.axhline(-2, color="red", ls="--", lw=0.7, alpha=0.6)
        ax_bot.axhline(z_threshold, color="red", ls=":", lw=0.7, alpha=0.4)
        ax_bot.axhline(-z_threshold, color="red", ls=":", lw=0.7, alpha=0.4)
        verdict = (
            "stationary" if (np.abs(z) > z_threshold).sum() == 0
            else f"NON-STATIONARY (Bonf |z|>{z_threshold:.2f})"
        )
        ax_bot.set_title(f"per-block z-score ({verdict})  max|z|={np.abs(z).max():.2f}",
                         fontsize=10)
        ax_bot.set_xlabel("HS-AFM frame index (block center)")
        ax_bot.set_ylabel("z (block fraction − pooled)/SE")
        ax_bot.set_ylim(-8, 8)
        ax_bot.grid(alpha=0.3)
        if col == 1:
            ax_bot.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8)

    fig.suptitle("State-population v2 (obj-044): windowed BC/Inter/EC/EO* fractions "
                 "vs HS-AFM time — non-stationarity test",
                 fontsize=12, fontweight="bold", y=0.995)

    out_png = FIG_DIR / "state_populations_v2_windowed.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"Saved {out_png}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
