#!/usr/bin/env python3
"""Rebuild mechanical-sensitivity composite using rotation-corrected RMSF.

Combine three independent per-residue flexibility metrics into a single
triple-agreement score (rectified-z-score product) and overlay Matsumoto
2008 switch residues for cross-method validation.

Inputs:
- figures/rmsf_v7_corrected/{video1,video2}/rmsf_per_residue.npy
- figures/cross_conformer_v7/per_residue_std_library.npy
- figures/hinges_v7/angle_std_per_residue.npy
- figures/matsumoto_overlay/matsumoto_overlay.json
- One representative library PDB for chain/resSeq lookup

Output:
- figures/mechanical_sensitivity_composite_v2.png
- figures/mechanical_sensitivity_composite_v2.json (top-20 hotspots)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import mdtraj as md

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DOMAIN_RANGES = [
    ("α-head-thigh", 0, 1, 435, "#2c7bb6"),
    ("α-calf",       0, 436, 741, "#abd9e9"),
    ("α-coil",       0, 742, 956, "#d7191c"),
    ("β-head",       1, 1, 352, "#fdae61"),
    ("β-tail",       1, 353, 692, "#fee090"),
]

CLASSIFICATION_COLOR = {
    "direct_hit":       "tab:red",
    "near_hit":         "tab:orange",
    "neighborhood_hit": "tab:blue",
    "miss":             "tab:green",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rmsf-v1", type=Path,
                   default=Path("figures/rmsf_v7_corrected/video1/rmsf_per_residue.npy"))
    p.add_argument("--rmsf-v2", type=Path,
                   default=Path("figures/rmsf_v7_corrected/video2/rmsf_per_residue.npy"))
    p.add_argument("--cross-conformer", type=Path,
                   default=Path("figures/cross_conformer_v7/per_residue_std_library.npy"))
    p.add_argument("--angle-std", type=Path,
                   default=Path("figures/hinges_v7/angle_std_per_residue.npy"))
    p.add_argument("--matsumoto-json", type=Path,
                   default=Path("figures/matsumoto_overlay/matsumoto_overlay.json"))
    p.add_argument("--reference-pdb", type=Path,
                   default=Path("results/afm_pipeline/v7_smoothed_final/video1/topology.pdb"))
    p.add_argument("--output", type=Path,
                   default=Path("figures/mechanical_sensitivity_composite_v2.png"))
    return p.parse_args()


def rectified_z(arr):
    z = (arr - arr.mean()) / (arr.std() + 1e-9)
    return np.maximum(z, 0.0)


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rmsf_v1 = np.load(args.rmsf_v1)
    rmsf_v2 = np.load(args.rmsf_v2)
    rmsf = 0.5 * (rmsf_v1 + rmsf_v2)
    cross = np.load(args.cross_conformer)
    ang = np.load(args.angle_std)

    n = len(rmsf)
    # Pad ang to length n by 0 on the ends.
    if len(ang) == n - 2:
        ang_padded = np.zeros(n)
        ang_padded[1:-1] = ang
    else:
        ang_padded = ang
    assert len(cross) == n, "cross-conformer length mismatch"

    z_rmsf = rectified_z(rmsf)
    z_cross = rectified_z(cross)
    z_ang = rectified_z(ang_padded)
    composite = z_rmsf * z_cross * z_ang

    # Topology -> chain/resSeq.
    topo = md.load(str(args.reference_pdb)).topology
    ca_atoms = [a for a in topo.atoms if a.name == "CA"]
    chains = np.array([a.residue.chain.index for a in ca_atoms])
    resseqs = np.array([a.residue.resSeq for a in ca_atoms])
    assert len(ca_atoms) == n

    # Top hotspots.
    order = np.argsort(-composite)
    top_n = 20
    top_rows = []
    for i in order[:top_n]:
        ch_i = int(chains[i])
        chain_letter = "A" if ch_i == 0 else "B"
        top_rows.append(
            {
                "ca_idx": int(i),
                "chain": chain_letter,
                "resSeq": int(resseqs[i]),
                "composite": float(composite[i]),
                "z_rmsf": float(z_rmsf[i]),
                "z_cross": float(z_cross[i]),
                "z_angle": float(z_ang[i]),
            }
        )

    # Matsumoto residues with current classification.
    with args.matsumoto_json.open() as f:
        mat = json.load(f)
    matsumoto_rows = mat["rows"]

    # Plot composite + each metric panel.
    fig, axes = plt.subplots(4, 1, figsize=(15, 11), sharex=True)
    x = np.arange(n)

    panels = [
        (axes[0], rmsf, "RMSF (Å) — rotation-corrected, V1+V2 mean"),
        (axes[1], cross, "Cross-conformer CA std (Å) — 615-frame library"),
        (axes[2], ang_padded, "CA-CA-CA angular σ (°)"),
        (axes[3], composite, "Triple-agreement (rectified z-score product)"),
    ]
    for ax, arr, title in panels:
        ax.plot(x, arr, lw=0.7, color="tab:gray")
        # Domain bands.
        for name, ch, lo, hi, color in DOMAIN_RANGES:
            mask = (chains == ch) & (resseqs >= lo) & (resseqs <= hi)
            if mask.any():
                idx = np.where(mask)[0]
                ax.axvspan(idx[0], idx[-1], alpha=0.13, color=color)
        ax.grid(alpha=0.25)
        ax.set_title(title)
    axes[3].set_xlabel("CA index (chain A then chain B)")

    # Mark Matsumoto residues on composite panel.
    for r in matsumoto_rows:
        if r["ca_idx"] is None:
            continue
        ci = r["ca_idx"]
        color = CLASSIFICATION_COLOR.get(r["classification"], "k")
        axes[3].axvline(ci, color=color, alpha=0.45, lw=1.2)
        axes[3].annotate(
            r["name"],
            xy=(ci, composite[ci]),
            xytext=(ci, composite[ci] + max(composite) * 0.05),
            fontsize=7,
            ha="center",
            color=color,
            rotation=90,
        )
        axes[3].scatter([ci], [composite[ci]], color=color, s=24, zorder=5)

    # Mark top hotspots on the composite panel.
    for row in top_rows[:10]:
        axes[3].plot(row["ca_idx"], row["composite"], "o",
                     mfc="none", mec="black", ms=9, lw=1.2, zorder=6)
        axes[3].annotate(
            f"{row['chain']}:{row['resSeq']}",
            xy=(row["ca_idx"], row["composite"]),
            xytext=(row["ca_idx"], row["composite"] * 1.18 if row["composite"] > 0 else 1.0),
            fontsize=8, ha="center", color="black",
        )

    # Legend for Matsumoto classification.
    handles = []
    for cls, color in CLASSIFICATION_COLOR.items():
        handles.append(
            plt.Line2D([0], [0], color=color, lw=2, label=f"Matsumoto {cls.replace('_', ' ')}")
        )
    handles.append(
        plt.Line2D(
            [0], [0], marker="o", mfc="none", mec="black",
            ms=9, lw=0, label="Top-10 triple-agreement"
        )
    )
    axes[3].legend(handles=handles, fontsize=8, loc="upper right", ncol=2)

    fig.suptitle(
        "Mechanical sensitivity composite v2 — rotation-corrected RMSF + cross-conformer std + bond-angle σ\n"
        "Matsumoto 2008 ENM-NMA switch residues overlaid (color = classification)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(str(args.output), dpi=140)
    plt.close(fig)
    print(f"Wrote {args.output}")

    out_json = args.output.with_suffix(".json")
    with out_json.open("w") as f:
        json.dump(
            {"top_20_triple_agreement": top_rows, "matsumoto": matsumoto_rows},
            f, indent=2,
        )
    print(f"Wrote {out_json}")

    # Printed top-10 summary.
    print("\nTop 10 triple-agreement hotspots:")
    for row in top_rows[:10]:
        print(
            f"  {row['chain']}:{row['resSeq']:>4d}  composite={row['composite']:7.3f}  "
            f"(zRMSF={row['z_rmsf']:.2f}, zCross={row['z_cross']:.2f}, zAng={row['z_angle']:.2f})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
