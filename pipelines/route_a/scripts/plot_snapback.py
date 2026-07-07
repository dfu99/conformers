#!/usr/bin/env python3
"""Plot the snap-back test: knee angle (and Rg) vs time for WT vs linchpin mutants.

Reads results/route_a/snapback/<tag>/trajectory.csv for each system. The signal is the
WT-vs-mutant DIFFERENTIAL: WT should hold the extended knee (~175°) while a mutant that
released the genu salt-bridge lock should show the knee angle fall (re-bending).
"""
import argparse
import csv
import glob
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LABELS = {"wt": ("WT (control)", "#2166ac"), "k459a": ("K459A", "#d6604d"),
          "e598a": ("E598A", "#1b7837"), "double": ("K459A/E598A", "#762a83")}


def load(tag_dir):
    ps, knee, rg = [], [], []
    csvp = os.path.join(tag_dir, "trajectory.csv")
    if not os.path.exists(csvp):
        return None
    for row in csv.DictReader(open(csvp)):
        ps.append(float(row["ps"])); knee.append(float(row["knee_deg"])); rg.append(float(row["Rg_A"]))
    return ps, knee, rg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/route_a/snapback")
    ap.add_argument("--out", default="figures/route_a_snapback.png")
    args = ap.parse_args()

    fig, (axK, axR) = plt.subplots(1, 2, figsize=(13.5, 5.4))
    axK.axhline(175, ls=":", color="gray", lw=1); axK.text(0, 176, "extended (175°)", fontsize=8, color="gray")
    axK.axhline(41, ls=":", color="gray", lw=1); axK.text(0, 43, "bent (41°)", fontsize=8, color="gray")

    any_data = False
    for tag in ["wt", "k459a", "e598a", "double"]:
        d = load(os.path.join(args.dir, tag))
        if d is None:
            continue
        any_data = True
        ps, knee, rg = d
        lab, col = LABELS[tag]
        axK.plot(ps, knee, "-", color=col, lw=1.8, label=lab)
        axR.plot(ps, rg, "-", color=col, lw=1.8, label=lab)

    axK.set_xlabel("time (ps)"); axK.set_ylabel("genu knee angle (deg)")
    axK.set_title("Knee angle vs time — does the mutant snap back?", fontsize=11, weight="bold")
    axK.set_ylim(30, 185); axK.grid(alpha=0.25); axK.legend(fontsize=9)
    axR.set_xlabel("time (ps)"); axR.set_ylabel("radius of gyration (Å)")
    axR.set_title("Compaction (Rg) vs time", fontsize=11, weight="bold")
    axR.grid(alpha=0.25); axR.legend(fontsize=9)
    fig.suptitle("αVβ3 snap-back test: mutate the genu lock → MD → watch the knee",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    if not any_data:
        axK.text(0.5, 0.5, "no trajectory data yet", transform=axK.transAxes, ha="center")
    fig.savefig(args.out, dpi=125)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
