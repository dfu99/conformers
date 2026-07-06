#!/usr/bin/env python3
"""Visualize the initial string (path of images) for the route-A string method.

Left  : the string in CV space -- genu hinge angle (primary reaction coordinate) and
        Rg vs normalized path position. Raw morph images (dots) are monotonic; the
        equal-arc reparametrized nodes (stars) are the evenly spaced starting images.
Right : per-segment Cα-RMSD (the string arc-length metric) -- the morph's equal
        angular steps already give near-uniform spacing (low coefficient of variation).
"""
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--string", default="results/route_a/initial_string.json")
    ap.add_argument("--out", default="figures/route_a_initial_string.png")
    args = ap.parse_args()

    d = json.load(open(args.string))
    img = d["per_image"]
    nodes = d["reparametrized_nodes"]

    frac = [p["arc_frac"] for p in img]
    genu = [p["genu_angle_deg"] for p in img]
    rg = [p["Rg_A"] for p in img]
    seg = [p["seg_rmsd_A"] for p in img][1:]  # segment i connects image i-1 -> i

    nfrac = [nd["arc_frac"] for nd in nodes]
    ngenu = [nd["genu_angle_deg"] for nd in nodes]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14.2, 5.5))

    # ---- Left: string in CV space ----
    axL.plot(frac, genu, "-o", color="#d6604d", ms=5, lw=1.6, label="genu angle, morph images")
    axL.plot(nfrac, ngenu, "*", color="#7f0000", ms=15, label=f"reparametrized nodes (n={len(nodes)})",
             zorder=5)
    axL.set_xlabel("normalized path position (Cα-RMSD arc length)", fontsize=10)
    axL.set_ylabel("genu hinge angle (deg)  — primary reaction coordinate",
                   color="#d6604d", fontsize=10)
    axL.tick_params(axis="y", labelcolor="#d6604d")
    axL.set_ylim(30, 190)
    ax2 = axL.twinx()
    ax2.plot(frac, rg, "-^", color="#1b7837", ms=4, lw=1.4, label="Rg")
    ax2.set_ylabel("radius of gyration (Å)", color="#1b7837", fontsize=10)
    ax2.tick_params(axis="y", labelcolor="#1b7837")
    axL.set_title("Initial string is monotonic in the genu-angle coordinate", fontsize=11, weight="bold")
    axL.grid(alpha=0.25)
    axL.annotate("bent (A)", xy=(0.0, genu[0]), xytext=(0.03, 55), fontsize=9, color="gray")
    axL.annotate("extended (B)", xy=(1.0, genu[-1]), xytext=(0.72, 165), fontsize=9, color="gray")
    h1, l1 = axL.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axL.legend(h1 + h2, l1 + l2, loc="center right", fontsize=8.5, framealpha=0.9)

    # ---- Right: per-segment arc length (spacing uniformity) ----
    x = np.arange(1, len(seg) + 1)
    mean = d["spacing_rmsd_A"]["mean"]
    cvsp = d["spacing_rmsd_A"]["cv_of_spacing"]
    axR.bar(x, seg, color="#4393c3", edgecolor="black", linewidth=0.5)
    axR.axhline(mean, color="#b2182b", lw=1.6, ls="--", label=f"mean {mean:.2f} Å")
    axR.set_xlabel("string segment (image i-1 → i)", fontsize=10)
    axR.set_ylabel("segment length  (Cα-RMSD, Å)", fontsize=10)
    axR.set_title(f"Even image spacing along the path  (spacing CV = {cvsp:.2f})",
                  fontsize=11, weight="bold")
    axR.set_xticks(x)
    axR.legend(fontsize=9)
    axR.grid(axis="y", alpha=0.25)

    fig.suptitle(f"Route-A initial string for the αVβ3 bent↔extended transition  "
                 f"({d['n_images']} images, path {d['total_path_rmsd_A']:.0f} Å)",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}  ({fig.get_size_inches()[0]*130:.0f}x{fig.get_size_inches()[1]*130:.0f})")


if __name__ == "__main__":
    main()
