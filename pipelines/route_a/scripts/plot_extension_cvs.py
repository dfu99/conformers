#!/usr/bin/env python3
"""Visualize the real extension CVs vs the broken placeholder, plus the morph trajectory.

Left  : |% change| bent->extended per CV. The real CVs (genu angle, end-to-end, Rg,
        long-axis extent) move a lot; the placeholder cv0 (head-thigh) is flat -- the
        empirical reason it was replaced.
Right : the CVs traced along the incremental leg-swing morph (0 -> 144 deg), showing
        each is a smooth, monotonic reaction coordinate suitable for the string method.
"""
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", default="results/route_a/extension_cv_endpoints.json")
    ap.add_argument("--trajectory", default="results/route_a/stage2c_result.json")
    ap.add_argument("--out", default="figures/route_a_extension_cvs.png")
    args = ap.parse_args()

    ep = json.load(open(args.endpoints))["separation"]
    traj = json.load(open(args.trajectory))["trajectory"]

    order = ["genu_angle_deg", "end_to_end_A", "Rg_A", "long_extent_A", "cv0_placeholder_A"]
    labels = {
        "genu_angle_deg": "genu hinge angle\n(41 -> 174 deg)",
        "end_to_end_A": "head<->foot end-to-end\n(52 -> 155 A)",
        "Rg_A": "radius of gyration\n(39 -> 66 A)",
        "long_extent_A": "long-axis extent\n(129 -> 208 A)",
        "cv0_placeholder_A": "cv0 placeholder: head-thigh\n(54 -> 52 A)  [REPLACED]",
    }
    pct = [100.0 * ep[k]["delta"] / abs(ep[k]["bent"]) for k in order]
    colors = ["#1b7837", "#1b7837", "#1b7837", "#1b7837", "#b2182b"]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14.5, 5.6))

    # ---- Left: percent separation per CV ----
    y = np.arange(len(order))[::-1]
    bars = axL.barh(y, pct, color=colors, edgecolor="black", linewidth=0.6, height=0.62)
    axL.set_yticks(y)
    axL.set_yticklabels([labels[k] for k in order], fontsize=9)
    axL.axvline(0, color="black", lw=0.8)
    axL.set_xlabel("change bent (state A) -> extended (state B)   [% of bent value]", fontsize=10)
    axL.set_title("Real CVs separate the endpoints; the placeholder does not", fontsize=11, weight="bold")
    for bar, p in zip(bars, pct):
        x = bar.get_width()
        axL.text(x + (4 if x >= 0 else -4), bar.get_y() + bar.get_height() / 2,
                 f"{p:+.0f}%", va="center", ha="left" if x >= 0 else "right",
                 fontsize=9, weight="bold",
                 color="#1b7837" if x >= 0 else "#b2182b")
    axL.set_xlim(-30, max(pct) * 1.18)
    axL.grid(axis="x", alpha=0.25)

    # ---- Right: CVs along the morph (reaction coordinate = cumulative swing angle) ----
    deg = [t["cum_deg"] for t in traj]
    rg = [t["Rg_A"] for t in traj]
    ext = [t["long_axis_extent_A"] for t in traj]
    cv1 = [t["cv1_beta_head_tail_A"] for t in traj]

    axR.plot(deg, ext, "-o", ms=3.5, lw=1.8, color="#2166ac", label="long-axis extent (A)")
    axR.plot(deg, cv1, "-s", ms=3.5, lw=1.8, color="#762a83", label="beta head<->tail (A)")
    axR.plot(deg, rg, "-^", ms=3.5, lw=1.8, color="#1b7837", label="radius of gyration (A)")
    # genu angle endpoints (measured on the structures, deg on a twin axis)
    ax2 = axR.twinx()
    ax2.plot([0, 144], [40.98, 173.51], "--D", ms=6, lw=1.4, color="#d6604d",
             label="genu hinge angle (deg)")
    ax2.set_ylabel("genu hinge angle (deg)", color="#d6604d", fontsize=10)
    ax2.tick_params(axis="y", labelcolor="#d6604d")
    ax2.set_ylim(20, 195)

    axR.set_xlabel("cumulative leg swing about the genu hinge (deg)", fontsize=10)
    axR.set_ylabel("distance (A)", fontsize=10)
    axR.set_title("Each CV is smooth & monotonic along the bent->extended morph", fontsize=11, weight="bold")
    axR.grid(alpha=0.25)
    axR.annotate("bent\n(state A)", xy=(0, rg[0]), xytext=(8, 150),
                 fontsize=9, ha="left", color="gray")
    axR.annotate("extended\n(state B)", xy=(144, ext[-1]), xytext=(150, 150),
                 fontsize=9, ha="right", color="gray")
    l1, lab1 = axR.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    axR.legend(l1 + l2, lab1 + lab2, loc="center right", fontsize=8.5, framealpha=0.9)

    fig.suptitle("Route-A extension collective variables for the integrin alphaVbeta3 string method",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}  ({fig.get_size_inches()[0]*130:.0f}x{fig.get_size_inches()[1]*130:.0f})")


if __name__ == "__main__":
    main()
