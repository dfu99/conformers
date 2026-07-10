#!/usr/bin/env python3
"""In what ORDER do the genu salt-bridge locks engage as αVβ3 extends?

Complements the snap-back MD (pending A5000 slot) with a CPU-only mechanistic read from
the morph path: for every frame along bent→extended, measure the genu knee angle (the
extension reaction coordinate) and the four cross-knee salt-bridge distances. The knee
angle at which each bridge first closes (<4 Å) tells us the SEQUENCE of lock engagement
— i.e. which latch is the primary/earliest one to test first.

Residues E598/K459/D457/D595/E636/K688 are all αV chain A and identity-numbered across
the bent crystal, the morph frames, and the extended model, so they select directly.
"""
import argparse
import glob
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import define_extension_cvs as cv  # genu_angle + rank-matched loaders

# bridge: (name, (anion_resid, anion_atoms), (cation_resid, cation_atoms))  -- all chain A
BRIDGES = [
    ("E598–K459", (598, ["OE1", "OE2"]), (459, ["NZ"])),
    ("K459–E636", (636, ["OE1", "OE2"]), (459, ["NZ"])),
    ("D457–K688", (457, ["OD1", "OD2"]), (688, ["NZ"])),
    ("D595–K688", (595, ["OD1", "OD2"]), (688, ["NZ"])),
]
COLORS = ["#b2182b", "#d6604d", "#2166ac", "#1b7837"]


def charged_atoms(path, chain="A"):
    """Return {(resid, atomname): xyz} for the atoms we need (chain A)."""
    want = set()
    for _, (ar, aa), (cr, ca) in BRIDGES:
        for a in aa:
            want.add((ar, a))
        for a in ca:
            want.add((cr, a))
    out = {}
    for line in open(path):
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        resid = int(line[22:26]); nm = line[12:16].strip()
        if (resid, nm) in want:
            out[(resid, nm)] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return out


def bridge_dist(atoms, spec):
    _, (ar, aa), (cr, ca) = spec
    apts = [atoms[(ar, a)] for a in aa if (ar, a) in atoms]
    cpts = [atoms[(cr, a)] for a in ca if (cr, a) in atoms]
    if not apts or not cpts:
        return None
    return float(min(np.linalg.norm(x - y) for x in apts for y in cpts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bent", default="data/reference_pdbs/1jv2.pdb")
    ap.add_argument("--frames-dir", required=True)
    ap.add_argument("--out-json", default="results/route_a/linchpin_engagement.json")
    ap.add_argument("--out", default="figures/route_a_linchpin_engagement.png")
    args = ap.parse_args()

    bent_chains = cv.load_ca(args.bent)
    paths = [args.bent] + sorted(glob.glob(os.path.join(args.frames_dir, "frame_*.pdb")))

    rows = []
    for p in paths:
        # genu angle (rank-matched to original numbering)
        if p == args.bent:
            byorig = cv.orig_keyed_from_ref(bent_chains)
        else:
            byorig = cv.build_original_keyed(bent_chains, cv.load_ca(p))
        genu = cv.genu_angle(byorig)
        atoms = charged_atoms(p)
        dists = {name: bridge_dist(atoms, spec) for name, *_ in [(b[0], b) for b in BRIDGES]
                 for spec in [next(x for x in BRIDGES if x[0] == name)]}
        rows.append({"frame": os.path.basename(p), "genu_deg": round(genu, 2),
                     **{name: (round(dists[name], 2) if dists[name] else None) for name, *_ in BRIDGES}})
    rows.sort(key=lambda r: r["genu_deg"])

    # engagement knee angle: first genu at which the bridge closes (<4 Å)
    engage = {}
    for name, *_ in BRIDGES:
        got = None
        for r in rows:
            if r[name] is not None and r[name] < 4.0:
                got = r["genu_deg"]; break
        engage[name] = got
    order = sorted([b[0] for b in BRIDGES], key=lambda n: (engage[n] is None, engage[n] or 1e9))

    json.dump({"per_frame": rows, "engage_knee_deg": engage, "order": order},
              open(args.out_json, "w"), indent=2)

    # ---- figure: bridge distance vs knee angle ----
    fig, ax = plt.subplots(figsize=(9.6, 5.6))
    genu = [r["genu_deg"] for r in rows]
    for (name, *_), col in zip(BRIDGES, COLORS):
        y = [r[name] for r in rows]
        ax.plot(genu, y, "-o", ms=4, lw=1.8, color=col,
                label=f"{name}  (engages ~{engage[name]:.0f}°)" if engage[name] else name)
        if engage[name]:
            ax.plot(engage[name], 4.0, "*", ms=13, color=col, zorder=5)
    ax.axhline(4.0, ls=":", color="gray", lw=1)
    ax.text(43, 4.4, "salt-bridge closes (4 Å)", fontsize=8, color="gray")
    ax.set_xlabel("genu knee angle (deg)  —  bent (41°) → extended (175°)", fontsize=10)
    ax.set_ylabel("charged-atom distance (Å)", fontsize=10)
    ax.set_title("Order of genu-lock engagement as αVβ3 extends", fontsize=12, weight="bold")
    ax.legend(fontsize=9, title="cross-knee salt bridge", title_fontsize=9)
    ax.grid(alpha=0.25)
    ax.set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("engagement order (earliest first):")
    for n in order:
        print(f"  {n}: closes at knee ~{engage[n]}°" if engage[n] else f"  {n}: never < 4 Å along path")
    print(f"wrote {args.out} and {args.out_json}")


if __name__ == "__main__":
    main()
