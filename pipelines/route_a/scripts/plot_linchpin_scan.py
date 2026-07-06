#!/usr/bin/env python3
"""Visualize the extended-state linchpin scan.

Left  : per-residue "extension-lock" score along the sequence. The residues that
        govern the state change cluster sharply at the genu (knee), in the αV
        thigh/calf-1 junction. Top candidates labeled.
Right : the key cross-knee salt-bridge locks -- charged-atom distance bent -> extended.
        The top three collapse from ~20-25 Å (apart, bent) to ~2.7 Å (locked, extended):
        electrostatic clasps that snap shut when the knee straightens. Mutating these is
        predicted to release the lock and let the extended state fall back to bent.
"""
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENU = {"A": 592, "B": 440}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", default="results/route_a/linchpin_scan.json")
    ap.add_argument("--out", default="figures/route_a_linchpin_scan.png")
    args = ap.parse_args()
    d = json.load(open(args.scan))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15.0, 6.0),
                                   gridspec_kw={"width_ratios": [1.35, 1.0]})

    # ---- Left: per-residue lock score along the sequence ----
    per = d["per_residue_score"]
    # x position: chain A at its resid; chain B appended after a gap
    gapA = 970
    xs, ys, cols = [], [], []
    for p in per:
        x = p["resid"] if p["chain"] == "A" else gapA + p["resid"]
        xs.append(x); ys.append(p["score"])
        cols.append("#2166ac" if p["chain"] == "A" else "#762a83")
    axL.vlines(xs, 0, ys, colors=cols, lw=0.8, alpha=0.85)
    axL.axvspan(GENU["A"] - 8, GENU["A"] + 8, color="orange", alpha=0.16, zorder=0)
    axL.axvspan(gapA + GENU["B"] - 8, gapA + GENU["B"] + 8, color="orange", alpha=0.16, zorder=0)
    axL.text(GENU["A"], axL.get_ylim()[1], " genu", color="darkorange", fontsize=8, va="top")

    # label top linchpins
    lab_pos = {(L["residue"]): (L["score"]) for L in d["extended_lock_linchpins"]}
    id2x = {}
    for p in per:
        rid = f'{p["chain"]}:{p["resname"]}{p["resid"]}'
        id2x[rid] = (p["resid"] if p["chain"] == "A" else gapA + p["resid"], p["score"])
    top = d["extended_lock_linchpins"][:9]
    for i, L in enumerate(top):
        if L["residue"] in id2x:
            x, y = id2x[L["residue"]]
            dy = 3 + (i % 3) * 3.2
            axL.annotate(L["residue"].split(":")[1], xy=(x, y), xytext=(x, y + dy),
                         fontsize=7.6, ha="center", color="black",
                         arrowprops=dict(arrowstyle="-", lw=0.5, color="gray"))
    axL.set_xlabel("residue  (αV chain A  |  β3 chain B)", fontsize=10)
    axL.set_ylabel("extension-lock score", fontsize=10)
    axL.set_title("Governing residues cluster at the genu (knee) hinge", fontsize=11, weight="bold")
    axL.set_xticks([1, 250, 500, 750, 956, gapA + 250, gapA + 550])
    axL.set_xticklabels(["A1", "250", "500", "750", "956", "B250", "B550"], fontsize=8)
    axL.set_ylim(0, max(ys) * 1.28)
    axL.grid(axis="y", alpha=0.2)
    from matplotlib.patches import Patch
    axL.legend(handles=[Patch(color="#2166ac", label="αV (chain A)"),
                        Patch(color="#762a83", label="β3 (chain B)"),
                        Patch(color="orange", alpha=0.3, label="genu ±8 res")],
               fontsize=8.5, loc="upper right")

    # ---- Right: cross-knee salt-bridge locks (dumbbell) ----
    locks = d["key_cross_knee_salt_bridges"][:8]
    locks = locks[::-1]  # largest change on top
    y = np.arange(len(locks))
    for i, k in enumerate(locks):
        axR.plot([k["bent_dist_A"], k["extended_dist_A"]], [i, i], "-",
                 color="gray", lw=1.4, zorder=1)
    axR.scatter([k["bent_dist_A"] for k in locks], y, s=70, color="#b2182b",
                zorder=3, label="bent (apart)", edgecolor="black", linewidth=0.5)
    axR.scatter([k["extended_dist_A"] for k in locks], y, s=70, color="#1b7837",
                zorder=3, label="extended (locked)", edgecolor="black", linewidth=0.5)
    axR.axvline(4.0, color="black", ls=":", lw=1, alpha=0.6)
    axR.text(4.1, len(locks) - 0.4, "4 Å salt-bridge\ncutoff", fontsize=7.5, color="black", va="top")
    axR.set_yticks(y)
    axR.set_yticklabels([k["pair"].replace(" — ", "\n↔ ") for k in locks], fontsize=8)
    axR.set_xlabel("charged-atom distance (Å)", fontsize=10)
    axR.set_title("Cross-knee salt bridges snap shut on extension", fontsize=11, weight="bold")
    axR.legend(fontsize=8.5, loc="lower right")
    axR.grid(axis="x", alpha=0.25)
    axR.set_xlim(0, max(k["bent_dist_A"] for k in locks) * 1.1)

    fig.suptitle("Which residues govern the αVβ3 bent↔extended switch — extended-state linchpins",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(args.out, dpi=125)
    print(f"wrote {args.out}  ({fig.get_size_inches()[0]*125:.0f}x{fig.get_size_inches()[1]*125:.0f})")


if __name__ == "__main__":
    main()
