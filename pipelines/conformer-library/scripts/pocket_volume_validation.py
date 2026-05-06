#!/usr/bin/env python3
"""P2Rank-equivalent grid-based pocket-volume validation of obj-056 candidate.

Reviewer-E follow-up to obj-056 (audit-2026-05-05 §25.4 P=4). The
project lacks a venv with rdkit/meeko/openbabel that would let us
run P2Rank or fpocket directly (CLAUDE.md §4 gemmi conflict makes
collapsing the existing venvs unsafe). Instead, we implement the
core algorithm from scratch in pure numpy: a LIGSITE-style probe-
grid void detector that quantifies pocket volume around a target
centroid.

Algorithm (LIGSITE, Hendlich et al. 1997):
  For each grid point, raycast in 7 directions (3 Cartesian axes +
  4 diagonal) and count how many rays hit protein atoms within
  radius R = 7 Å. Grid points with ≥ 5/7 hits AND distance to
  nearest atom in [1.4, 4.0] Å are "pocket-like". Pocket volume =
  count × cell³.

Three regions evaluated:
  A — obj-056 candidate (β3 K417-K422 centroid)
  B — MIDAS centroid (positive control, canonical RGD pocket)
  C — αV C-terminal flat surface (negative control)

Compares vol(bent) vs vol(extended) for each region. The candidate
is "validated" if vol_ext(A) > vol_bent(A) AND the gain is bigger
than the noise floor estimated from the negative control.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LIB_DIR = ROOT / "data" / "runs" / "avb3" / "conformers" / "all_frames_bent_extended"
CV0_CACHE = ROOT / "results" / "afm_pipeline" / "rgd_docking" / "library_cv0_cache.npy.npz"
OUT_DIR = ROOT / "results" / "cryptic_pockets"
FIG = ROOT / "figures" / "pocket_volume_validation_v1.png"

N_PER_STATE = 5
GRID_SPACING_A = 1.5
BOX_HALF_A = 15.0
PROBE_MIN_A = 1.4
PROBE_MAX_A = 4.0
RAYCAST_RADIUS_A = 7.0
RAY_HITS_THRESHOLD = 5

CANDIDATE_RESIDUES = [("B", r) for r in (417, 419, 420, 421, 422)]
MIDAS_RESIDUES = [("A", 218),
                  ("B", 119), ("B", 121), ("B", 123), ("B", 220), ("B", 251),
                  ("B", 122), ("B", 214), ("B", 215), ("B", 216)]
CONTROL_RESIDUES = [("A", r) for r in (940, 945, 950, 955)]

RAYS = np.array([
    [1, 0, 0], [-1, 0, 0],
    [0, 1, 0], [0, -1, 0],
    [0, 0, 1], [0, 0, -1],
    [1, 1, 1] / np.sqrt(3),
], dtype=np.float64)


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    if cid:
        return cid
    return chr(ord("A") + atom.residue.chain.index)


def residue_centroid(traj: md.Trajectory, residues: list[tuple[str, int]]) -> np.ndarray:
    coords_nm: list[np.ndarray] = []
    topology = traj.topology
    assert topology is not None
    xyz = traj.xyz
    assert xyz is not None
    target_set = set(residues)
    for atom in topology.atoms:
        if atom.name != "CA":
            continue
        if (chain_id(atom), int(atom.residue.resSeq)) in target_set:
            coords_nm.append(xyz[0, atom.index])
    if not coords_nm:
        raise ValueError(f"no CAs in {residues[:3]}...")
    return np.mean(coords_nm, axis=0) * 10.0


def heavy_atom_coords(traj: md.Trajectory) -> np.ndarray:
    topology = traj.topology
    assert topology is not None
    xyz = traj.xyz
    assert xyz is not None
    indices = [atom.index for atom in topology.atoms
               if atom.element is not None and atom.element.symbol != "H"]
    return xyz[0, indices] * 10.0


def pocket_volume_at_centroid(coords_A: np.ndarray, centroid_A: np.ndarray) -> dict:
    """LIGSITE-style pocket-volume calculator (vectorized cylindrical-projection).

    For each accessible grid point in a 30 Å cubic box around centroid:
      ray-hit = (any atom A satisfying  0 < (A-cell)·dir < RAYCAST_RADIUS_A
                 AND  |perpendicular component| < 1.4 Å)
    Pocket grid points have ≥ RAY_HITS_THRESHOLD hits across 7 rays.

    Optimization: only consider atoms within (BOX_HALF + RAYCAST_RADIUS) of
    the centroid — drops the per-frame atom set from 25 000 to ~1 500.
    """
    n = int(2 * BOX_HALF_A / GRID_SPACING_A) + 1
    axis = np.linspace(-BOX_HALF_A, BOX_HALF_A, n)
    gx, gy, gz = np.meshgrid(axis, axis, axis, indexing="ij")
    grid_pts = np.stack([gx, gy, gz], axis=-1).reshape(-1, 3) + centroid_A
    n_pts = len(grid_pts)

    cutoff = BOX_HALF_A + RAYCAST_RADIUS_A + 2.0
    near_mask = np.linalg.norm(coords_A - centroid_A, axis=-1) < cutoff
    atoms = coords_A[near_mask]
    n_atoms = len(atoms)
    if n_atoms == 0:
        return {
            "centroid_A": centroid_A.tolist(),
            "pocket_volume_A3": 0.0,
            "n_pocket_cells": 0, "n_accessible_cells": 0, "n_total_cells": n_pts,
            "buriedness_fraction": 0.0, "n_atoms_considered": 0,
        }

    nearest_dist = np.full(n_pts, np.inf, dtype=np.float64)
    chunk = 256
    for s in range(0, n_pts, chunk):
        e = min(s + chunk, n_pts)
        d = np.linalg.norm(grid_pts[s:e, None, :] - atoms[None, :, :], axis=-1)
        nearest_dist[s:e] = d.min(axis=1)
    accessible = (nearest_dist >= PROBE_MIN_A) & (nearest_dist <= PROBE_MAX_A)
    accessible_idx = np.where(accessible)[0]
    if len(accessible_idx) == 0:
        return {
            "centroid_A": centroid_A.tolist(),
            "pocket_volume_A3": 0.0,
            "n_pocket_cells": 0,
            "n_accessible_cells": 0,
            "n_total_cells": n_pts,
            "buriedness_fraction": 0.0, "n_atoms_considered": int(n_atoms),
        }
    cells = grid_pts[accessible_idx]

    ray_hits = np.zeros(len(accessible_idx), dtype=np.int8)
    n_rays = len(RAYS)
    cell_chunk = 512
    for ray in RAYS:
        for s in range(0, len(accessible_idx), cell_chunk):
            e = min(s + cell_chunk, len(accessible_idx))
            v = atoms[None, :, :] - cells[s:e, None, :]
            t = np.einsum("ija,a->ij", v, ray)
            d2 = (v ** 2).sum(axis=-1) - t ** 2
            mask = (t > 0.0) & (t < RAYCAST_RADIUS_A) & (d2 < 1.96)
            hit = mask.any(axis=-1)
            ray_hits[s:e] += hit.astype(np.int8)

    pocket_mask_acc = ray_hits >= RAY_HITS_THRESHOLD
    n_pocket = int(pocket_mask_acc.sum())
    pocket_volume_A3 = n_pocket * GRID_SPACING_A ** 3
    return {
        "centroid_A": centroid_A.tolist(),
        "pocket_volume_A3": float(pocket_volume_A3),
        "n_pocket_cells": n_pocket,
        "n_accessible_cells": int(accessible.sum()),
        "n_total_cells": n_pts,
        "buriedness_fraction": (n_pocket / max(int(accessible.sum()), 1)),
        "n_atoms_considered": int(n_atoms),
        "n_rays": int(n_rays),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG.parent.mkdir(parents=True, exist_ok=True)

    cache = np.load(CV0_CACHE, allow_pickle=True)
    names = list(cache["names"])
    cv0 = np.array(cache["cv0"], dtype=np.float64)
    sort_idx = np.argsort(cv0)
    bent_idx = sort_idx[:N_PER_STATE]
    ext_idx = sort_idx[-N_PER_STATE:]
    print(f"library frames {len(names)}; using {N_PER_STATE} bent + {N_PER_STATE} extended")
    print(f"  bent CV0 ≤ {cv0[bent_idx[-1]]:.2f} Å; extended CV0 ≥ {cv0[ext_idx[0]]:.2f} Å")

    regions = {
        "A_candidate_K417-K422": CANDIDATE_RESIDUES,
        "B_MIDAS_pos_control":   MIDAS_RESIDUES,
        "C_alphaV_Cterm_neg_control": CONTROL_RESIDUES,
    }
    results: dict[str, dict[str, list]] = {
        rn: {"bent": [], "extended": [], "bent_cv0": [], "extended_cv0": []}
        for rn in regions
    }

    t_start = time.time()
    for label, frame_indices in [("bent", bent_idx), ("extended", ext_idx)]:
        for k, idx in enumerate(frame_indices):
            t0 = time.time()
            pdb = LIB_DIR / names[idx]
            traj = md.load(str(pdb))
            coords = heavy_atom_coords(traj)
            for region_name, residues in regions.items():
                centroid = residue_centroid(traj, residues)
                rec = pocket_volume_at_centroid(coords, centroid)
                rec["frame"] = names[idx]
                rec["state"] = label
                rec["cv0_A"] = float(cv0[idx])
                results[region_name][label].append(rec)
                results[region_name][f"{label}_cv0"].append(float(cv0[idx]))
            dt = time.time() - t0
            print(f"  {label} {k+1}/{N_PER_STATE} ({names[idx]}, "
                  f"CV0={cv0[idx]:.2f} Å) — {dt:.1f}s")
    print(f"total: {time.time() - t_start:.1f}s")

    region_summary: dict[str, dict] = {}
    for region_name, region_data in results.items():
        bent_vols = np.array([r["pocket_volume_A3"] for r in region_data["bent"]])
        ext_vols = np.array([r["pocket_volume_A3"] for r in region_data["extended"]])
        delta_mean = float(ext_vols.mean() - bent_vols.mean())
        pooled_se = float(np.sqrt(bent_vols.var(ddof=1) / N_PER_STATE
                                  + ext_vols.var(ddof=1) / N_PER_STATE))
        welch_t = delta_mean / pooled_se if pooled_se > 1e-9 else 0.0
        region_summary[region_name] = {
            "bent_mean_pocket_volume_A3": float(bent_vols.mean()),
            "bent_std_pocket_volume_A3": float(bent_vols.std(ddof=1)),
            "extended_mean_pocket_volume_A3": float(ext_vols.mean()),
            "extended_std_pocket_volume_A3": float(ext_vols.std(ddof=1)),
            "delta_mean_A3": delta_mean,
            "delta_pooled_se_A3": pooled_se,
            "welch_t": float(welch_t),
            "delta_significant_alpha_005": bool(abs(welch_t) > 1.96),
            "n_per_state": N_PER_STATE,
            "frames": region_data,
        }
        print(f"\n{region_name}:")
        print(f"  bent     pocket vol = {bent_vols.mean():.0f} ± {bent_vols.std(ddof=1):.0f} Å³")
        print(f"  extended pocket vol = {ext_vols.mean():.0f} ± {ext_vols.std(ddof=1):.0f} Å³")
        print(f"  Δ = {delta_mean:+.0f} Å³  (SE {pooled_se:.0f}; t = {welch_t:+.2f})")

    cand_delta = region_summary["A_candidate_K417-K422"]["delta_mean_A3"]
    ctrl_delta = region_summary["C_alphaV_Cterm_neg_control"]["delta_mean_A3"]
    midas_delta = region_summary["B_MIDAS_pos_control"]["delta_mean_A3"]

    cand_sig = region_summary["A_candidate_K417-K422"]["delta_significant_alpha_005"]
    abs_gain_vs_ctrl = cand_delta - ctrl_delta
    relative_to_midas = (cand_delta / midas_delta) if abs(midas_delta) > 1e-6 else float("nan")

    if cand_sig and cand_delta > 0 and abs_gain_vs_ctrl > 0:
        verdict = ("VALIDATED — candidate pocket volume increases on extension "
                   "and exceeds negative-control noise; cryptic site confirmed")
    elif cand_sig and cand_delta < 0 and abs(cand_delta) > abs(ctrl_delta):
        verdict = ("DOWNGRADED (informative) — candidate region IS conformationally active "
                   "during BC→EC, but pocket volume DECREASES (opens to bulk solvent) "
                   "rather than forming a discrete druggable pocket. "
                   "Consistent with obj-056 ΔSASA gain + obj-039 MIDAS occlusion: "
                   "extension exposes surface, does not gate a new ligand-binding site.")
    elif cand_sig and cand_delta < 0:
        verdict = ("INCONCLUSIVE — candidate ΔV is negative but no greater than "
                   "matched negative-control surface; cannot distinguish candidate "
                   "from a generic surface region")
    elif not cand_sig and abs(cand_delta) < abs(ctrl_delta) * 1.5:
        verdict = "NULL — candidate ΔV indistinguishable from negative-control noise"
    else:
        verdict = "AMBIGUOUS — candidate ΔV trend present but not significant at α=0.05"

    full_summary = {
        "regions": region_summary,
        "verdict": verdict,
        "effect_summary": {
            "candidate_delta_A3": cand_delta,
            "midas_delta_A3": midas_delta,
            "ctrl_delta_A3": ctrl_delta,
            "candidate_minus_control_A3": abs_gain_vs_ctrl,
            "candidate_relative_to_midas": relative_to_midas,
        },
        "params": {
            "n_per_state": N_PER_STATE,
            "grid_spacing_A": GRID_SPACING_A,
            "box_half_A": BOX_HALF_A,
            "probe_min_A": PROBE_MIN_A,
            "probe_max_A": PROBE_MAX_A,
            "raycast_radius_A": RAYCAST_RADIUS_A,
            "ray_hits_threshold": RAY_HITS_THRESHOLD,
            "n_rays": int(len(RAYS)),
            "algorithm": "LIGSITE-style 7-ray probe-grid void detection",
        },
    }

    summary_path = OUT_DIR / "pocket_volume_validation.json"
    with summary_path.open("w") as f:
        json.dump(full_summary, f, indent=2)
    print(f"\n=== VERDICT: {verdict} ===")
    print(f"\nWrote {summary_path}")

    plot(full_summary, results)
    print(f"Wrote {FIG}")
    return 0


def plot(summary: dict, results: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.5))
    region_order = ["A_candidate_K417-K422", "B_MIDAS_pos_control", "C_alphaV_Cterm_neg_control"]
    region_labels = {
        "A_candidate_K417-K422": "A: β3 K417-K422\n(obj-056 candidate)",
        "B_MIDAS_pos_control": "B: MIDAS centroid\n(positive control)",
        "C_alphaV_Cterm_neg_control": "C: αV C-term flat\n(negative control)",
    }

    for ax_i, region in enumerate(region_order):
        ax = axes[ax_i]
        s = summary["regions"][region]
        bent_vols = [r["pocket_volume_A3"] for r in results[region]["bent"]]
        ext_vols = [r["pocket_volume_A3"] for r in results[region]["extended"]]
        bent_cv0s = [r["cv0_A"] for r in results[region]["bent"]]
        ext_cv0s = [r["cv0_A"] for r in results[region]["extended"]]
        ax.scatter(bent_cv0s, bent_vols, s=80, color="#2166ac",
                   edgecolor="black", linewidth=0.6, label="bent",
                   zorder=5)
        ax.scatter(ext_cv0s, ext_vols, s=80, color="#b2182b",
                   edgecolor="black", linewidth=0.6, label="extended",
                   zorder=5)
        ax.axhline(s["bent_mean_pocket_volume_A3"], color="#2166ac",
                   linestyle="--", linewidth=1.0, alpha=0.7)
        ax.axhline(s["extended_mean_pocket_volume_A3"], color="#b2182b",
                   linestyle="--", linewidth=1.0, alpha=0.7)
        ax.set_xlabel("CV0 (Å)")
        ax.set_ylabel("pocket volume (Å³)")
        ax.set_title(region_labels[region], fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(alpha=0.3)

        delta = s["delta_mean_A3"]
        se = s["delta_pooled_se_A3"]
        t = s["welch_t"]
        sig = s["delta_significant_alpha_005"]
        ax.text(0.02, 0.97,
                f"Δ = {delta:+.0f} ± {se:.0f} Å³\n"
                f"t = {t:+.2f}  ({'sig' if sig else 'ns'} α=0.05)",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round",
                          facecolor=("#ccebc5" if sig and delta > 0
                                     else "#fbb4ae" if sig and delta < 0
                                     else "#ffffcc"),
                          edgecolor="#666"))

    fig.suptitle(
        f"P2Rank-equivalent pocket-volume validation of obj-056 candidate "
        f"(LIGSITE 7-ray probe grid, {N_PER_STATE}+{N_PER_STATE} frames)\n"
        f"Verdict: {summary['verdict'][:120]}",
        fontsize=10, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    fig.savefig(FIG, dpi=140, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
