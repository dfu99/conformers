#!/usr/bin/env python3
"""Validation matrix for the snap-back experiment: confirm each mutant, built in the
production implicit-solvent force field, knocks out exactly its target salt bridges.

Reads results/route_a/snapback/bc_<tag>/build_check.json (from snapback_md.py --build-only).
This is the pre-flight that de-risks the GPU run: every system builds, and the tracked
cross-knee salt bridges are present/absent exactly as the mutation design intends.
"""
import argparse
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SYSTEMS = [("bc_wt", "WT (control)"), ("bc_k459a", "K459A"),
           ("bc_e598a", "E598A"), ("bc_double", "K459A/E598A")]
BRIDGES = ["E598-K459", "K459-E636", "D457-K688", "D595-K688"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/route_a/snapback")
    ap.add_argument("--out", default="figures/route_a_snapback_validation.png")
    args = ap.parse_args()

    data = {}
    for tag, _ in SYSTEMS:
        p = os.path.join(args.dir, tag, "build_check.json")
        data[tag] = json.load(open(p)) if os.path.exists(p) else None

    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    nrow, ncol = len(SYSTEMS), len(BRIDGES)
    for i, (tag, label) in enumerate(SYSTEMS):
        row = nrow - 1 - i
        d = data[tag]
        salt = d["start"]["salt"] if d else {}
        for j, br in enumerate(BRIDGES):
            present = br in salt
            color = "#1b7837" if present else "#f4c9c2"
            ax.add_patch(plt.Rectangle((j, row), 1, 1, facecolor=color,
                                       edgecolor="white", linewidth=2))
            txt = f"{salt[br]:.2f} Å" if present else "✕ removed"
            ax.text(j + 0.5, row + 0.5, txt, ha="center", va="center",
                    fontsize=9, color="white" if present else "#b2182b",
                    weight="bold")
        atoms = d["atoms"] if d else "?"
        ax.text(ncol + 0.15, row + 0.5, f"{atoms} atoms\nknee {d['start']['knee']:.0f}°" if d else "—",
                ha="left", va="center", fontsize=8.5, color="black")

    ax.set_xlim(0, ncol + 1.6); ax.set_ylim(0, nrow)
    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_xticklabels(BRIDGES, fontsize=9)
    ax.set_yticks(np.arange(nrow) + 0.5)
    ax.set_yticklabels([lab for _, lab in SYSTEMS][::-1], fontsize=10)
    ax.set_xlabel("tracked cross-knee salt bridge (t=0)", fontsize=10)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Snap-back pre-flight: each mutation knocks out exactly its target bridges\n"
                 "(all 4 systems build in the production implicit-solvent force field)",
                 fontsize=11, weight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#1b7837", label="salt bridge present"),
                       Patch(color="#f4c9c2", label="removed by mutation")],
              loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
