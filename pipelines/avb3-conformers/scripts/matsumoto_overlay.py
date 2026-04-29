#!/usr/bin/env python3
"""Overlay Matsumoto 2008 ENM-NMA switch residues onto our v7 hinge-σ map.

Matsumoto et al., Biophys J 2008 (PMC2527288): elastic-network NMA on
bent αVβ3 ectodomain identified switch residues that stabilize the bent
state. Their snap residue β3:Arg633 plus surrounding interaction
networks are tabulated below.

Compare against `figures/hinges_v7/angle_std_per_residue.npy` and report:
1. σ at each Matsumoto residue (direct read).
2. Rank of each Matsumoto residue in the full σ distribution.
3. Nearest top-20 hinge candidate from our list (within ±10 resSeq).
4. Categorical hit/near-miss/miss classification.

Usage:
    python matsumoto_overlay.py \\
        --frames-dir data/runs/avb3/conformers/all_frames_bent_extended \\
        --angle-std figures/hinges_v7/angle_std_per_residue.npy \\
        --output-dir figures/matsumoto_overlay
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


MATSUMOTO_2008_RESIDUES = [
    {"chain": "B", "resSeq": 633, "name": "Arg633", "role": "primary_snap"},
    {"chain": "B", "resSeq": 375, "name": "Leu375", "role": "snap_sandwich"},
    {"chain": "B", "resSeq": 389, "name": "Leu389", "role": "snap_sandwich"},
    {"chain": "B", "resSeq": 374, "name": "Cys374", "role": "interaction_A"},
    {"chain": "B", "resSeq": 388, "name": "Gly388", "role": "interaction_A"},
    {"chain": "B", "resSeq": 404, "name": "Arg404", "role": "hybrid_egf"},
    {"chain": "B", "resSeq": 364, "name": "Glu364", "role": "hybrid_egf"},
    {"chain": "B", "resSeq": 543, "name": "Ser543", "role": "interaction_B"},
    {"chain": "B", "resSeq": 548, "name": "Leu548", "role": "interaction_B"},
    {"chain": "B", "resSeq": 549, "name": "Cys549", "role": "interaction_B"},
    {"chain": "B", "resSeq": 550, "name": "Asp550", "role": "interaction_B"},
    {"chain": "A", "resSeq": 305, "name": "Ser305", "role": "alpha_constraint"},
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frames-dir", type=Path, required=True)
    p.add_argument("--angle-std", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def load_residue_index(frames_dir: Path):
    """Load the first frame to map CA index -> (chain, resSeq)."""
    pdbs = sorted(frames_dir.glob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No PDBs in {frames_dir}")
    traj = md.load(str(pdbs[0]))
    top = traj.topology
    ca_atoms = [a for a in top.atoms if a.name == "CA"]
    ca_meta = []
    for a in ca_atoms:
        ca_meta.append(
            {
                "ca_idx": len(ca_meta),
                "chain": a.residue.chain.chain_id or chr(ord("A") + a.residue.chain.index),
                "resSeq": a.residue.resSeq,
                "resname": a.residue.name,
            }
        )
    return ca_meta


def find_ca_idx(ca_meta, chain: str, resSeq: int):
    for m in ca_meta:
        if m["chain"] == chain and m["resSeq"] == resSeq:
            return m["ca_idx"]
    return None


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sigma_per_ca = np.load(args.angle_std)
    print(f"Loaded σ array: {sigma_per_ca.shape}")
    ca_meta = load_residue_index(args.frames_dir)
    print(f"Loaded {len(ca_meta)} CA atoms")

    # σ array is for residues i with i in [1, N-2]; map to ca_idx => sigma_idx = ca_idx - 1.
    # Length of σ array should be N-2 if it followed `hinge_angles.py`.
    n_ca = len(ca_meta)
    if len(sigma_per_ca) == n_ca - 2:
        # σ[k] corresponds to ca_idx = k + 1
        def sigma_for_ca(idx: int) -> float | None:
            sig_idx = idx - 1
            if 0 <= sig_idx < len(sigma_per_ca):
                return float(sigma_per_ca[sig_idx])
            return None
    elif len(sigma_per_ca) == n_ca:
        def sigma_for_ca(idx: int) -> float | None:
            return float(sigma_per_ca[idx])
    else:
        raise ValueError(
            f"σ array length {len(sigma_per_ca)} doesn't match N_CA={n_ca}"
        )

    # Compute rank of each ca_idx by σ (descending).
    sigmas_indexed = []
    for m in ca_meta:
        s = sigma_for_ca(m["ca_idx"])
        if s is None:
            continue
        sigmas_indexed.append((m["ca_idx"], s))
    sigmas_indexed.sort(key=lambda t: -t[1])
    rank_lookup = {ca_idx: rank + 1 for rank, (ca_idx, _) in enumerate(sigmas_indexed)}
    total_with_sigma = len(sigmas_indexed)

    rows = []
    for matsumoto in MATSUMOTO_2008_RESIDUES:
        ca_idx = find_ca_idx(ca_meta, matsumoto["chain"], matsumoto["resSeq"])
        if ca_idx is None:
            rows.append(
                {
                    **matsumoto,
                    "ca_idx": None,
                    "sigma_deg": None,
                    "rank": None,
                    "percentile": None,
                    "classification": "absent",
                    "nearest_top20_in_window": None,
                }
            )
            continue
        s = sigma_for_ca(ca_idx)
        rank = rank_lookup.get(ca_idx)
        percentile = (
            100.0 * (1.0 - (rank - 1) / total_with_sigma) if rank else None
        )
        # Nearest top-20 candidate within ±10 resSeq.
        top20_meta = [ca_meta[ci] for (ci, _) in sigmas_indexed[:20]]
        nearest = None
        for cand in top20_meta:
            if cand["chain"] == matsumoto["chain"]:
                d = abs(cand["resSeq"] - matsumoto["resSeq"])
                if d <= 10 and (nearest is None or d < nearest["delta"]):
                    nearest = {
                        "chain": cand["chain"],
                        "resSeq": cand["resSeq"],
                        "delta": d,
                        "sigma_deg": sigma_for_ca(cand["ca_idx"]),
                    }
        if percentile is not None and percentile >= 95.0:
            classification = "direct_hit"
        elif percentile is not None and percentile >= 80.0:
            classification = "near_hit"
        elif nearest is not None:
            classification = "neighborhood_hit"
        else:
            classification = "miss"
        rows.append(
            {
                **matsumoto,
                "ca_idx": ca_idx,
                "sigma_deg": s,
                "rank": rank,
                "percentile": percentile,
                "classification": classification,
                "nearest_top20_in_window": nearest,
            }
        )

    overlay_json = args.output_dir / "matsumoto_overlay.json"
    with overlay_json.open("w") as f:
        json.dump(
            {"total_with_sigma": total_with_sigma, "rows": rows}, f, indent=2
        )
    print(f"Wrote {overlay_json}")

    # Markdown report.
    md_path = args.output_dir / "matsumoto_overlay.md"
    with md_path.open("w") as f:
        f.write("# Matsumoto 2008 switch residues vs our v7 hinge-σ\n\n")
        f.write(
            "σ = std of ∠(CA_{i-1}, CA_i, CA_{i+1}) across 615 v7 conformers.\n"
        )
        f.write(f"Total CA residues with σ measured: {total_with_sigma}.\n\n")
        f.write(
            "| Matsumoto residue | role | σ (°) | rank | percentile | "
            "nearest top-20 (±10) | classification |\n"
        )
        f.write("|---|---|---|---|---|---|---|\n")
        for r in rows:
            sig = f"{r['sigma_deg']:.1f}" if r["sigma_deg"] is not None else "—"
            rk = str(r["rank"]) if r["rank"] is not None else "—"
            pct = (
                f"{r['percentile']:.1f}%" if r["percentile"] is not None else "—"
            )
            near = r["nearest_top20_in_window"]
            near_str = (
                f"{near['chain']}:{near['resSeq']} (Δ={near['delta']}, σ={near['sigma_deg']:.1f}°)"
                if near
                else "—"
            )
            f.write(
                f"| {r['chain']}:{r['resSeq']} ({r['name']}) | {r['role']} | {sig} | {rk} | {pct} | {near_str} | {r['classification']} |\n"
            )
        # Summary counts.
        n_direct = sum(1 for r in rows if r["classification"] == "direct_hit")
        n_near = sum(1 for r in rows if r["classification"] == "near_hit")
        n_neigh = sum(1 for r in rows if r["classification"] == "neighborhood_hit")
        n_miss = sum(1 for r in rows if r["classification"] == "miss")
        n_absent = sum(1 for r in rows if r["classification"] == "absent")
        f.write("\n## Summary\n\n")
        f.write(
            f"- direct_hit (≥95th pct): {n_direct}\n"
            f"- near_hit (80–95th pct): {n_near}\n"
            f"- neighborhood_hit (top-20 within ±10): {n_neigh}\n"
            f"- miss: {n_miss}\n"
            f"- absent (residue not in our PDB): {n_absent}\n"
        )
    print(f"Wrote {md_path}")

    # Figure: σ trace over chain B with Matsumoto residues marked.
    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharey=False)
    for ax, chain_id in zip(axes, ("A", "B")):
        x_res, y_sig = [], []
        for m in ca_meta:
            if m["chain"] != chain_id:
                continue
            s = sigma_for_ca(m["ca_idx"])
            if s is not None:
                x_res.append(m["resSeq"])
                y_sig.append(s)
        ax.plot(x_res, y_sig, "-", color="tab:gray", alpha=0.6, lw=0.8)
        # Mark Matsumoto residues for this chain.
        for r in rows:
            if r["chain"] != chain_id or r["sigma_deg"] is None:
                continue
            color_map = {
                "direct_hit": "tab:red",
                "near_hit": "tab:orange",
                "neighborhood_hit": "tab:blue",
                "miss": "tab:green",
            }
            color = color_map.get(r["classification"], "k")
            ax.axvline(r["resSeq"], color=color, alpha=0.35, lw=1.2)
            ax.annotate(
                r["name"],
                xy=(r["resSeq"], r["sigma_deg"]),
                xytext=(r["resSeq"], r["sigma_deg"] + 2.0),
                fontsize=7,
                ha="center",
                color=color,
            )
            ax.scatter([r["resSeq"]], [r["sigma_deg"]], color=color, s=30, zorder=5)
        ax.set_title(f"Chain {chain_id} — angular σ (deg) with Matsumoto 2008 switch residues marked")
        ax.set_xlabel("resSeq")
        ax.set_ylabel("σ (deg)")
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig_path = args.output_dir / "matsumoto_overlay.png"
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)
    print(f"Wrote {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
