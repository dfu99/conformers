#!/usr/bin/env python3
"""Cross-correlate per-block FES drift with HS-AFM frame metadata (obj-046).

Audit §15.6 follow-up to obj-045. Five per-frame quality signals are
derived from existing artifacts (no GPU, no re-running the pipeline):

  1. fit correlation per frame — high = good fit (template registered
     well to the AFM image).
  2. head xy position (px) per frame — gross drift = molecule moving
     on the substrate (artifact) vs internal dynamics.
  3. inter-frame head jump (px) — frame-to-frame jitter; spikes
     indicate tracking failures or genuinely fast events.
  4. raw GIF per-frame mean intensity — proxy for tip-sample distance,
     imaging force, substrate brightness.
  5. raw GIF per-frame std — proxy for image contrast, tip condition.

For each video and each 50-frame block, summarize each signal and
compute the Pearson correlation across blocks between FES-min CV0
location (obj-045) and each metadata signal. Annotate the plot with
the strongest correlations and any aligned breakpoints.

Outputs:
  figures/fes_drift_metadata_correlation.png
  results/afm_pipeline/free_energy_profile/fes_drift_metadata.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[3]
V1_DIR = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final" / "video1"
V2_DIR = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final" / "video2"
GIF_V1 = ROOT / "inbox" / "linz_avb3_video1.gif"
GIF_V2 = ROOT / "inbox" / "linz_avb3_video2.gif"
FES_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "per_block_dg.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "free_energy_profile"
FIG_DIR = ROOT / "figures"

WINDOW = 50


def load_gif_intensities(gif_path: Path) -> tuple[np.ndarray, np.ndarray]:
    img = Image.open(gif_path)
    means = []
    stds = []
    try:
        while True:
            arr = np.asarray(img.convert("L"), dtype=np.float32)
            means.append(arr.mean())
            stds.append(arr.std())
            img.seek(img.tell() + 1)
    except EOFError:
        pass
    return np.array(means), np.array(stds)


def block_aggregate(arr: np.ndarray, w: int = WINDOW, fn=np.mean) -> np.ndarray:
    n_blocks = arr.size // w
    return np.array([fn(arr[b * w:(b + 1) * w]) for b in range(n_blocks)])


def head_drift_per_block(head_pos: np.ndarray, w: int = WINDOW) -> np.ndarray:
    n_blocks = head_pos.shape[0] // w
    return np.array([
        float(np.linalg.norm(head_pos[(b + 1) * w - 1] - head_pos[b * w]))
        for b in range(n_blocks)
    ])


def head_jitter_per_block(head_pos: np.ndarray, w: int = WINDOW) -> np.ndarray:
    n_blocks = head_pos.shape[0] // w
    return np.array([
        float(np.linalg.norm(np.diff(head_pos[b * w:(b + 1) * w], axis=0), axis=1).mean())
        for b in range(n_blocks)
    ])


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fes = np.load(FES_NPZ)
    min_loc_v1 = np.asarray(fes["min_loc_v1"])
    min_loc_v2 = np.asarray(fes["min_loc_v2"])
    centers_v1 = np.asarray(fes["centers_v1"])
    centers_v2 = np.asarray(fes["centers_v2"])

    print("Loading per-frame metadata...")
    v1_corr = np.load(V1_DIR / "fitted_correlations.npy")
    v2_corr = np.load(V2_DIR / "fitted_correlations.npy")
    v1_head = np.load(V1_DIR / "head_positions_px.npy")
    v2_head = np.load(V2_DIR / "head_positions_px.npy")
    print(f"  V1: corr {v1_corr.shape}, head {v1_head.shape}")
    print(f"  V2: corr {v2_corr.shape}, head {v2_head.shape}")

    print("Loading GIF intensities...")
    g1_mean, g1_std = load_gif_intensities(GIF_V1)
    g2_mean, g2_std = load_gif_intensities(GIF_V2)
    print(f"  V1 GIF: {g1_mean.size} frames, mean range [{g1_mean.min():.1f}, {g1_mean.max():.1f}]")
    print(f"  V2 GIF: {g2_mean.size} frames, mean range [{g2_mean.min():.1f}, {g2_mean.max():.1f}]")

    # The fitted set is GIF[skip_frames:], skip = 30
    SKIP = 30
    g1_mean_aligned = g1_mean[SKIP:]
    g1_std_aligned = g1_std[SKIP:]
    g2_mean_aligned = g2_mean[SKIP:]
    g2_std_aligned = g2_std[SKIP:]
    if g1_mean_aligned.size != v1_corr.size:
        # truncate to match
        n = min(g1_mean_aligned.size, v1_corr.size)
        g1_mean_aligned = g1_mean_aligned[:n]
        g1_std_aligned = g1_std_aligned[:n]
    if g2_mean_aligned.size != v2_corr.size:
        n = min(g2_mean_aligned.size, v2_corr.size)
        g2_mean_aligned = g2_mean_aligned[:n]
        g2_std_aligned = g2_std_aligned[:n]
    print(f"  Aligned V1: {g1_mean_aligned.size} frames; V2: {g2_mean_aligned.size}")

    # --- Per-block aggregates ---
    v1_block = {
        "fes_min": min_loc_v1,
        "fit_corr": block_aggregate(v1_corr, WINDOW),
        "head_drift": head_drift_per_block(v1_head, WINDOW),
        "head_jitter": head_jitter_per_block(v1_head, WINDOW),
        "intensity_mean": block_aggregate(g1_mean_aligned, WINDOW),
        "intensity_std": block_aggregate(g1_std_aligned, WINDOW),
        "centers": centers_v1,
    }
    v2_block = {
        "fes_min": min_loc_v2,
        "fit_corr": block_aggregate(v2_corr, WINDOW),
        "head_drift": head_drift_per_block(v2_head, WINDOW),
        "head_jitter": head_jitter_per_block(v2_head, WINDOW),
        "intensity_mean": block_aggregate(g2_mean_aligned, WINDOW),
        "intensity_std": block_aggregate(g2_std_aligned, WINDOW),
        "centers": centers_v2,
    }

    # --- Pearson correlations across blocks ---
    signals = ["fit_corr", "head_drift", "head_jitter", "intensity_mean", "intensity_std"]

    def safe_pearson(x, y):
        if np.std(x) < 1e-12 or np.std(y) < 1e-12 or x.size < 3:
            return float("nan"), float("nan")
        r, p = pearsonr(x, y)
        return float(r), float(p)

    correlations: dict = {"v1": {}, "v2": {}}
    print()
    print(f"{'signal':<18} {'V1 r':<10} {'V1 p':<10} {'V2 r':<10} {'V2 p':<10}")
    for sig in signals:
        r1, p1 = safe_pearson(v1_block["fes_min"], v1_block[sig])
        r2, p2 = safe_pearson(v2_block["fes_min"], v2_block[sig])
        correlations["v1"][sig] = {"r": r1, "p": p1}
        correlations["v2"][sig] = {"r": r2, "p": p2}
        print(f"{sig:<18} {r1:>+.3f}     {p1:.3f}     {r2:>+.3f}     {p2:.3f}")

    # Find which signals show |r| > 0.5 (notable) and p < 0.05 (significant)
    notable_v1 = [s for s in signals if abs(correlations["v1"][s]["r"]) > 0.5]
    notable_v2 = [s for s in signals if abs(correlations["v2"][s]["r"]) > 0.5]
    significant_v1 = [s for s in signals if correlations["v1"][s]["p"] < 0.05]
    significant_v2 = [s for s in signals if correlations["v2"][s]["p"] < 0.05]
    print()
    print(f"V1 notable (|r| > 0.5): {notable_v1}")
    print(f"V1 significant (p < 0.05): {significant_v1}")
    print(f"V2 notable (|r| > 0.5): {notable_v2}")
    print(f"V2 significant (p < 0.05): {significant_v2}")

    # --- V1 block 5 specifically: FES min jumps from 79.2 → 55.5 Å (block 4 → 5) ---
    v1_block5_change = {
        "fes_min_jump_A": float(v1_block["fes_min"][5] - v1_block["fes_min"][4]),
        "fit_corr_change": float(v1_block["fit_corr"][5] - v1_block["fit_corr"][4]),
        "head_drift_change": float(v1_block["head_drift"][5] - v1_block["head_drift"][4]),
        "intensity_mean_change": float(
            v1_block["intensity_mean"][5] - v1_block["intensity_mean"][4]
        ),
        "intensity_std_change": float(
            v1_block["intensity_std"][5] - v1_block["intensity_std"][4]
        ),
    }
    print()
    print("V1 block 4 → 5 transition (FES min jump 79.2 → 55.5 Å):")
    for k, v in v1_block5_change.items():
        print(f"  {k}: {v:+.4f}")

    # --- V2 block 17 specifically: lowest fit correlation (0.876) ---
    v2_block17 = {
        "fes_min_A": float(v2_block["fes_min"][17]),
        "fit_corr": float(v2_block["fit_corr"][17]),
        "head_drift": float(v2_block["head_drift"][17]),
        "intensity_mean": float(v2_block["intensity_mean"][17]),
    }
    print()
    print(f"V2 block 17 (lowest fit corr): {v2_block17}")

    # --- Save JSON summary ---
    summary = {
        "window_frames": WINDOW,
        "n_v1_blocks": int(min_loc_v1.size),
        "n_v2_blocks": int(min_loc_v2.size),
        "correlations": correlations,
        "v1_block5_transition": v1_block5_change,
        "v2_block17_lowest_fit": v2_block17,
        "v1_notable_signals": notable_v1,
        "v2_notable_signals": notable_v2,
        "v1_significant_signals": significant_v1,
        "v2_significant_signals": significant_v2,
        "verdict": (
            "FES_drift_correlates_with_metadata"
            if (significant_v1 or significant_v2)
            else "FES_drift_independent_of_metadata"
        ),
        "v1_block_data": {k: v.tolist() if hasattr(v, "tolist") else v
                          for k, v in v1_block.items()},
        "v2_block_data": {k: v.tolist() if hasattr(v, "tolist") else v
                          for k, v in v2_block.items()},
    }
    out_json = OUT_DIR / "fes_drift_metadata.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"\nSaved {out_json}")

    # --- Plot ---
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 2, hspace=0.40, wspace=0.20,
                           height_ratios=[1.0, 1.0, 1.1])

    pretty = {
        "fit_corr": "fit corr (mean)",
        "head_drift": "head xy drift (px)",
        "head_jitter": "head jitter (px/frame)",
        "intensity_mean": "GIF mean intensity",
        "intensity_std": "GIF intensity std",
    }
    sig_color = {
        "fit_corr": "#2ca02c",
        "head_drift": "#d62728",
        "head_jitter": "#9467bd",
        "intensity_mean": "#1f77b4",
        "intensity_std": "#ff7f0e",
    }

    # Top row: V1 timeline + signals
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(v1_block["centers"], v1_block["fes_min"], "ko-", lw=2,
             ms=8, label="FES min CV0 (Å)")
    ax1.set_ylim(50, 90)
    ax1.set_ylabel("FES-min CV0 (Å)")
    ax1.set_title("V1 — FES-min trajectory + metadata signals (z-scored)",
                  fontsize=10.5)
    ax1.grid(alpha=0.3)
    ax1b = ax1.twinx()
    for sig in signals:
        x = v1_block[sig]
        z = (x - x.mean()) / (x.std() + 1e-9)
        ax1b.plot(v1_block["centers"], z, "-", color=sig_color[sig],
                  lw=1.0, alpha=0.7,
                  label=f"{pretty[sig]} (r={correlations['v1'][sig]['r']:+.2f})")
    ax1b.axvline((4 + 5) * WINDOW / 2 + 25, color="red", ls=":", lw=1.0,
                 alpha=0.6, label="V1 block 4→5 break")
    ax1b.set_ylabel("metadata z-score")
    ax1b.legend(loc="upper right", fontsize=7, ncol=2)

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(v2_block["centers"], v2_block["fes_min"], "ko-", lw=2,
             ms=6, label="FES min CV0 (Å)")
    ax2.set_ylim(50, 90)
    ax2.set_ylabel("FES-min CV0 (Å)")
    ax2.set_title("V2 — FES-min trajectory + metadata signals (z-scored)",
                  fontsize=10.5)
    ax2.grid(alpha=0.3)
    ax2b = ax2.twinx()
    for sig in signals:
        x = v2_block[sig]
        z = (x - x.mean()) / (x.std() + 1e-9)
        ax2b.plot(v2_block["centers"], z, "-", color=sig_color[sig],
                  lw=1.0, alpha=0.6,
                  label=f"{pretty[sig]} (r={correlations['v2'][sig]['r']:+.2f})")
    ax2b.set_ylabel("metadata z-score")
    ax2b.legend(loc="upper right", fontsize=7, ncol=2)

    # Middle row: scatter of FES-min vs each signal
    ax_sv1 = fig.add_subplot(gs[1, 0])
    for sig in signals:
        x = v1_block[sig]
        x_z = (x - x.mean()) / (x.std() + 1e-9)
        y = v1_block["fes_min"]
        ax_sv1.scatter(x_z, y, color=sig_color[sig], s=40, alpha=0.7,
                       label=f"{pretty[sig]} r={correlations['v1'][sig]['r']:+.2f}")
    ax_sv1.set_xlabel("metadata signal (z-score)")
    ax_sv1.set_ylabel("FES-min CV0 (Å)")
    ax_sv1.set_title(
        f"V1 scatter — significant: {', '.join(significant_v1) if significant_v1 else 'none'}",
        fontsize=10,
    )
    ax_sv1.legend(loc="upper right", fontsize=7)
    ax_sv1.grid(alpha=0.3)

    ax_sv2 = fig.add_subplot(gs[1, 1])
    for sig in signals:
        x = v2_block[sig]
        x_z = (x - x.mean()) / (x.std() + 1e-9)
        y = v2_block["fes_min"]
        ax_sv2.scatter(x_z, y, color=sig_color[sig], s=30, alpha=0.6,
                       label=f"{pretty[sig]} r={correlations['v2'][sig]['r']:+.2f}")
    ax_sv2.set_xlabel("metadata signal (z-score)")
    ax_sv2.set_ylabel("FES-min CV0 (Å)")
    ax_sv2.set_title(
        f"V2 scatter — significant: {', '.join(significant_v2) if significant_v2 else 'none'}",
        fontsize=10,
    )
    ax_sv2.legend(loc="upper right", fontsize=7)
    ax_sv2.grid(alpha=0.3)

    # Bottom row: per-frame raw signals to look for non-block-aligned events
    ax_raw1 = fig.add_subplot(gs[2, 0])
    n_frames = v1_corr.size
    x = np.arange(n_frames)
    ax_raw1.plot(x, v1_corr, "-", color="#2ca02c", lw=0.6, alpha=0.7,
                 label="fit corr")
    ax_raw1.plot(x, g1_mean_aligned / 100, "-", color="#1f77b4", lw=0.6,
                 alpha=0.7, label="GIF mean / 100")
    ax_raw1.set_ylim(0.5, 1.5)
    ax_raw1.set_xlabel("HS-AFM frame index (after skip 30)")
    ax_raw1.set_ylabel("signal (normalized)")
    ax_raw1.set_title(
        f"V1 raw per-frame signals — V1 block 4→5 transition: "
        f"FES min {v1_block5_change['fes_min_jump_A']:+.2f} Å, "
        f"fit corr {v1_block5_change['fit_corr_change']:+.4f}",
        fontsize=10,
    )
    for b in range(min_loc_v1.size + 1):
        ax_raw1.axvline(b * WINDOW, color="black", ls="-", lw=0.3, alpha=0.3)
    ax_raw1.axvline(5 * WINDOW, color="red", ls=":", lw=1.5, alpha=0.7)
    ax_raw1.legend(loc="upper right", fontsize=8)
    ax_raw1.grid(alpha=0.3)

    ax_raw2 = fig.add_subplot(gs[2, 1])
    n_frames = v2_corr.size
    x = np.arange(n_frames)
    ax_raw2.plot(x, v2_corr, "-", color="#2ca02c", lw=0.4, alpha=0.7,
                 label="fit corr")
    ax_raw2.plot(x, g2_mean_aligned / 100, "-", color="#1f77b4", lw=0.4,
                 alpha=0.7, label="GIF mean / 100")
    ax_raw2.set_ylim(0.0, 1.5)
    ax_raw2.set_xlabel("HS-AFM frame index (after skip 30)")
    ax_raw2.set_ylabel("signal (normalized)")
    ax_raw2.set_title(
        f"V2 raw per-frame signals — V2 block 17 (lowest corr 0.876): "
        f"FES min = {v2_block17['fes_min_A']:.1f} Å",
        fontsize=10,
    )
    for b in range(min_loc_v2.size + 1):
        ax_raw2.axvline(b * WINDOW, color="black", ls="-", lw=0.3, alpha=0.2)
    ax_raw2.axvline(17 * WINDOW, color="red", ls=":", lw=1.5, alpha=0.7)
    ax_raw2.legend(loc="upper right", fontsize=8)
    ax_raw2.grid(alpha=0.3)

    fig.suptitle("FES-drift / HS-AFM-metadata cross-correlation (obj-046)",
                 fontsize=12, fontweight="bold", y=0.995)

    out_png = FIG_DIR / "fes_drift_metadata_correlation.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print(f"Saved {out_png}")
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
