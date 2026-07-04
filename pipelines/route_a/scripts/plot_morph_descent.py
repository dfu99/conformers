#!/usr/bin/env python3
"""Plot the Stage-2c incremental leg-swing: extension + energy vs swing angle.

Reads stage2c_result.json (a trajectory of {cum_deg, E_kJ_mol, Rg_A,
long_axis_extent_A}) and shows that, as the legs swing open, the molecule
extends (Rg + long-axis grow bent->extended) while the energy stays physically
sane -- i.e. the incremental morph avoids the catastrophic clash that stalled
the one-shot rigid swing.
"""
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="results/route_a/stage2c_result.json")
    ap.add_argument("--out", default="figures/route_a_morph_descent.png")
    args = ap.parse_args()

    with open(args.json) as fh:
        data = json.load(fh)
    traj = data["trajectory"]
    deg = [t["cum_deg"] for t in traj]
    energy = [t["E_kJ_mol"] for t in traj]
    rg = [t["Rg_A"] for t in traj]
    extent = [t["long_axis_extent_A"] for t in traj]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    fig.suptitle("Route-A Stage 2c: bent → extended αVβ3 by incremental "
                 "leg-swing (vacuum, CPU)", fontsize=13, fontweight="bold")

    # extension metrics
    ax1.plot(deg, extent, "o-", color="#2b6cb0", lw=2, label="long-axis extent")
    ax1.plot(deg, rg, "s-", color="#dd6b20", lw=2, label="radius of gyration")
    ax1.set_xlabel("cumulative leg-swing (degrees)")
    ax1.set_ylabel("size (Å)")
    ax1.set_title("molecule extends as legs swing open")
    ax1.legend(frameon=False)
    ax1.grid(alpha=0.3)
    ax1.annotate("bent", (deg[0], extent[0]), textcoords="offset points",
                 xytext=(6, -14), fontsize=9, color="#2b6cb0")
    ax1.annotate("extended", (deg[-1], extent[-1]), textcoords="offset points",
                 xytext=(-30, 8), fontsize=9, color="#2b6cb0")

    # energy stays sane (no 1e13 blow-up)
    ax2.plot(deg, energy, "o-", color="#38a169", lw=2)
    ax2.set_xlabel("cumulative leg-swing (degrees)")
    ax2.set_ylabel("potential energy (kJ/mol)")
    ax2.set_title("energy stays physically sane every step")
    ax2.grid(alpha=0.3)

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    print(f"bent:     Rg={rg[0]} extent={extent[0]}  E={energy[0]:.3e}")
    print(f"extended: Rg={rg[-1]} extent={extent[-1]}  E={energy[-1]:.3e}")


if __name__ == "__main__":
    main()
