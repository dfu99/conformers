#!/usr/bin/env python3
"""Side-by-side comparison of rotation-dominated vs rotation-corrected RMSF."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DOMAIN_RANGES = {
    "α-head-thigh": (0, 1, 435),
    "α-calf":       (0, 436, 741),
    "α-coil":       (0, 742, 956),
    "β-head":       (1, 1, 352),
    "β-tail":       (1, 353, 692),
}
DOMAIN_COLORS = {
    "α-head-thigh": "#2c7bb6",
    "α-calf":       "#abd9e9",
    "α-coil":       "#d7191c",
    "β-head":       "#fdae61",
    "β-tail":       "#fee090",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--orig-dir", type=Path, required=True)
    p.add_argument("--corrected-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--label", type=str, default="")
    return p.parse_args()


def load_rmsf(d: Path):
    arr = np.load(str(d / "rmsf_per_residue.npy"))
    chain = []
    resseq = []
    with (d / "rmsf_table.csv").open() as f:
        for row in csv.DictReader(f):
            chain.append(int(row["chain"]))
            resseq.append(int(row["resSeq"]))
    return arr, np.array(chain), np.array(resseq)


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    orig, ch_o, rs_o = load_rmsf(args.orig_dir)
    corr, ch_c, rs_c = load_rmsf(args.corrected_dir)
    assert len(orig) == len(corr), "Mismatched RMSF length"

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    x = np.arange(len(orig))
    axes[0].plot(x, orig, lw=0.7, color="tab:gray",
                 label=f"head-anchored only (mean={orig.mean():.1f} Å)")
    axes[0].set_title(f"RMSF — rotation-dominated (orig){' — ' + args.label if args.label else ''}")
    axes[0].set_ylabel("RMSF (Å)")
    axes[0].grid(alpha=0.3)
    axes[0].legend(loc="upper right", fontsize=9)

    axes[1].plot(x, corr, lw=0.7, color="tab:blue",
                 label=f"head-Kabsch aligned (mean={corr.mean():.1f} Å)")
    axes[1].set_title("RMSF — rotation-corrected (head-Kabsch alignment)")
    axes[1].set_ylabel("RMSF (Å)")
    axes[1].set_xlabel("CA index (α-chain then β-chain)")
    axes[1].grid(alpha=0.3)
    axes[1].legend(loc="upper right", fontsize=9)

    for name, (ch, lo, hi) in DOMAIN_RANGES.items():
        mask = (ch_o == ch) & (rs_o >= lo) & (rs_o <= hi)
        if mask.any():
            idx = np.where(mask)[0]
            for ax in axes:
                ax.axvspan(idx[0], idx[-1], alpha=0.15,
                           color=DOMAIN_COLORS.get(name, "gray"))

    fig.tight_layout()
    fig.savefig(str(args.output), dpi=140)
    plt.close(fig)
    print(f"Wrote {args.output}")

    # Domain-mean comparison table
    print("\nDomain-averaged comparison:")
    print(f"  {'domain':14s}  {'orig':>8s}  {'corrected':>10s}  {'ratio':>6s}")
    for name, (ch, lo, hi) in DOMAIN_RANGES.items():
        m = (ch_o == ch) & (rs_o >= lo) & (rs_o <= hi)
        if m.any():
            o = orig[m].mean()
            c = corr[m].mean()
            print(f"  {name:14s}  {o:8.2f}  {c:10.2f}  {c/o:6.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
