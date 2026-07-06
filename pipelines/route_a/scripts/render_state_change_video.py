#!/usr/bin/env python3
"""Render the αVβ3 bent->extended leg-swing as an animated GIF (2D Cα projection).

All frames are superposed on the UPPER body (headpiece + thigh) of the bent state, so
the head stays fixed at the top and the lower legs visibly swing open from folded
(bent) to straight (extended) about the genu (knee). Colored upper vs lower so the
moving part is obvious. Pure numpy + matplotlib + imageio; no GPU.
"""
import argparse
import glob
import os
import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import imageio.v2 as imageio

import define_extension_cvs as cv

UPPER = {"A": (1, 592), "B": (1, 440)}
GENU = {"A": (588, 596), "B": (436, 444)}


def in_range(ch, resid, table):
    lo, hi = table[ch]
    return lo <= resid <= hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bent", default="data/reference_pdbs/1jv2.pdb")
    ap.add_argument("--frames-dir", required=True)
    ap.add_argument("--out", default="results/route_a/state_change_video.gif")
    args = ap.parse_args()

    bent_chains = cv.load_ca(args.bent)
    order = [(ch, r) for ch in sorted(bent_chains) for r, _ in bent_chains[ch]]
    upper_mask = np.array([in_range(ch, r, UPPER) for ch, r in order])
    genu_mask = np.array([in_range(ch, r, GENU) for ch, r in order])

    def mat(byorig):
        return np.array([byorig[ch][r] for ch, r in order])

    # image 0 = bent, images 1..16 = morph frames
    byorigs = [cv.orig_keyed_from_ref(bent_chains)]
    for f in sorted(glob.glob(os.path.join(args.frames_dir, "frame_*.pdb"))):
        byorigs.append(cv.build_original_keyed(bent_chains, cv.load_ca(f)))
    mats = [mat(b) for b in byorigs]
    genu_ang = [cv.genu_angle(b) for b in byorigs]

    # reference = bent; superpose every frame's UPPER onto bent UPPER (Kabsch)
    ref = mats[0]
    refU = ref[upper_mask]
    refUc = refU.mean(0)

    def superpose(P):
        Pu = P[upper_mask]
        Puc = Pu.mean(0)
        H = (Pu - Puc).T @ (refU - refUc)
        U, _, Vt = np.linalg.svd(H)
        dsign = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, dsign]) @ U.T
        return (P - Puc) @ R.T + refUc

    aligned = [superpose(m) for m in mats]

    # fixed 2D view: PCA plane over ALL aligned coords (long axis = PC1 vertical)
    allpts = np.vstack(aligned)
    c = allpts.mean(0)
    _, _, vt = np.linalg.svd(allpts - c, full_matrices=False)
    e1, e2 = vt[0], vt[1]

    def proj(P):
        q = P - c
        return q @ e1, q @ e2  # x=e1 (long axis), y=e2

    xs = np.concatenate([proj(a)[0] for a in aligned])
    ys = np.concatenate([proj(a)[1] for a in aligned])
    xlim = (xs.min() - 8, xs.max() + 8)
    ylim = (ys.min() - 8, ys.max() + 8)

    # boomerang sequence for a clean loop
    seq = list(range(len(aligned))) + list(range(len(aligned) - 2, 0, -1))
    imgs = []
    for k in seq:
        A = aligned[k]
        px, py = proj(A)
        fig, ax = plt.subplots(figsize=(4.0, 6.4))
        ax.scatter(py[upper_mask], px[upper_mask], s=7, c="#2166ac", label="headpiece + thigh (fixed)")
        ax.scatter(py[~upper_mask], px[~upper_mask], s=7, c="#c0392b", label="lower legs (swing)")
        ax.scatter(py[genu_mask], px[genu_mask], s=26, c="orange", edgecolor="black",
                   linewidth=0.4, zorder=5, label="genu (knee)")
        ax.set_xlim(ylim); ax.set_ylim(xlim)
        ax.set_aspect("equal"); ax.axis("off")
        state = "bent" if k == 0 else ("extended" if k == len(aligned) - 1 else "opening")
        ax.set_title(f"αVβ3 leg swing — {state}\nknee angle {genu_ang[k]:.0f}°",
                     fontsize=11, weight="bold")
        if k == 0:
            ax.legend(loc="lower center", fontsize=7.5, framealpha=0.9,
                      bbox_to_anchor=(0.5, -0.02))
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=95)
        plt.close(fig)
        buf.seek(0)
        imgs.append(imageio.imread(buf))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    # hold endpoints a bit longer
    durs = []
    for i, k in enumerate(seq):
        durs.append(0.5 if k in (0, len(aligned) - 1) else 0.14)
    imageio.mimsave(args.out, imgs, duration=durs, loop=0)
    print(f"wrote {args.out}  ({len(imgs)} frames, knee {genu_ang[0]:.0f}->{genu_ang[-1]:.0f} deg)")


if __name__ == "__main__":
    main()
