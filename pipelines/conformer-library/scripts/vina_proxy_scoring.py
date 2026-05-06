#!/usr/bin/env python3
"""Vina-proxy ligand-fit scoring of obj-056 candidate (audit §27 v18).

Reviewer-E follow-up to obj-057 (audit-2026-05-05 §25.4 P=4). The
overseer queued blind Vina/smina docking on the obj-056 candidate
(β3 K417-K422) vs the MIDAS pocket. Real AutoDock-Vina docking is
gated on PACE A100 + a new docking venv with rdkit/meeko/openbabel.

Implements a simplified Vina-equivalent scoring function that mimics
the dominant terms of AutoDock-Vina's empirical scoring without
needing a real ligand library. For each candidate pocket centroid:
  • hydrophobic_score = count of (C/S) heavy atoms in 4-7 Å shell
                        × hydrophobic-fragment fit weight
  • hbond_score       = count of polar H-bond donors+acceptors in
                        3-5 Å shell × H-bond fit weight
  • clash_penalty     = count of heavy atoms in 0-3.5 Å shell
                        (would clash with a small-molecule fragment)
  • total_score       = hydrophobic + hbond - clash

Higher score → more "druggable" (the pocket can accommodate a
small-molecule fragment with favorable van-der-Waals + H-bond
interactions). Vina's actual scoring uses a torsional entropy term
+ Gaussians + electrostatics; this proxy captures the spatial
fragment-fit aspect that dominates buried-pocket scoring.

Three regions × 5+5 frames, same as obj-057:
  A — obj-056 candidate (β3 K417-K422 centroid)
  B — MIDAS centroid (positive control, canonical RGD pocket)
  C — αV C-terminal flat surface (negative control)

Validates whether obj-056 candidate has druggable atom-type
distribution AROUND its centroid (independent of pocket-volume
result from obj-057).
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
FIG = ROOT / "figures" / "vina_proxy_scoring_v1.png"

N_PER_STATE = 5

CANDIDATE_RESIDUES = [("B", r) for r in (417, 419, 420, 421, 422)]
MIDAS_RESIDUES = [("A", 218),
                  ("B", 119), ("B", 121), ("B", 123), ("B", 220), ("B", 251),
                  ("B", 122), ("B", 214), ("B", 215), ("B", 216)]
CONTROL_RESIDUES = [("A", r) for r in (940, 945, 950, 955)]

CLASH_R = 3.5
HYDROPH_INNER, HYDROPH_OUTER = 4.0, 7.0
HBOND_INNER, HBOND_OUTER = 3.0, 5.0

W_HYDROPH = 0.040
W_HBOND = 0.587
W_CLASH = 0.500


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    if cid:
        return cid
    return chr(ord("A") + atom.residue.chain.index)


def residue_centroid(traj: md.Trajectory, residues: list[tuple[str, int]]) -> np.ndarray:
    topology = traj.topology
    assert topology is not None
    xyz = traj.xyz
    assert xyz is not None
    coords_nm: list[np.ndarray] = []
    target_set = set(residues)
    for atom in topology.atoms:
        if atom.name != "CA":
            continue
        if (chain_id(atom), int(atom.residue.resSeq)) in target_set:
            coords_nm.append(xyz[0, atom.index])
    if not coords_nm:
        raise ValueError(f"no CAs in {residues[:3]}...")
    return np.mean(coords_nm, axis=0) * 10.0


def classify_atoms(traj: md.Trajectory) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (hydrophobic_coords, polar_donor_acceptor_coords, all_heavy_coords) in Å.

    Hydrophobic: C, S heavy atoms in V/I/L/M/F/Y/W/A residues.
    H-bond donor/acceptor: N, O, S in any residue (atom-level naming
    used by mdtraj is not Vina-style atomtype, so use element + residue
    context: side-chain N/O/S in non-Gly residues).
    All heavy: any non-H atom.
    """
    topology = traj.topology
    assert topology is not None
    xyz = traj.xyz
    assert xyz is not None
    coords = xyz[0] * 10.0
    hydroph_residues = {"VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "ALA"}
    hyd_idx, hbond_idx, heavy_idx = [], [], []
    for atom in topology.atoms:
        if atom.element is None or atom.element.symbol == "H":
            continue
        heavy_idx.append(atom.index)
        elem = atom.element.symbol
        rname = atom.residue.name
        is_sidechain = atom.name not in {"N", "C", "CA", "O", "OXT"}
        if elem == "C" and rname in hydroph_residues and is_sidechain:
            hyd_idx.append(atom.index)
        elif elem == "S":
            hyd_idx.append(atom.index)
        if elem in {"N", "O"}:
            hbond_idx.append(atom.index)
    return (coords[hyd_idx], coords[hbond_idx], coords[heavy_idx])


def vina_proxy_score(centroid: np.ndarray, hyd: np.ndarray,
                     hbond: np.ndarray, heavy: np.ndarray) -> dict:
    cutoff = max(HYDROPH_OUTER, HBOND_OUTER, CLASH_R) + 1.0
    near_heavy_mask = np.linalg.norm(heavy - centroid, axis=-1) < cutoff
    near_heavy = heavy[near_heavy_mask]
    near_hyd = hyd[np.linalg.norm(hyd - centroid, axis=-1) < cutoff]
    near_hbond = hbond[np.linalg.norm(hbond - centroid, axis=-1) < cutoff]

    d_heavy = np.linalg.norm(near_heavy - centroid, axis=-1)
    d_hyd = np.linalg.norm(near_hyd - centroid, axis=-1)
    d_hbond = np.linalg.norm(near_hbond - centroid, axis=-1)

    n_clash = int(((d_heavy > 0.5) & (d_heavy < CLASH_R)).sum())
    n_hyd = int(((d_hyd >= HYDROPH_INNER) & (d_hyd < HYDROPH_OUTER)).sum())
    n_hbond = int(((d_hbond >= HBOND_INNER) & (d_hbond < HBOND_OUTER)).sum())

    hydroph_score = n_hyd * W_HYDROPH
    hbond_score = n_hbond * W_HBOND
    clash_penalty = n_clash * W_CLASH
    total_score = hydroph_score + hbond_score - clash_penalty

    return {
        "centroid_A": centroid.tolist(),
        "n_clash": n_clash,
        "n_hydrophobic": n_hyd,
        "n_hbond": n_hbond,
        "hydroph_score": float(hydroph_score),
        "hbond_score": float(hbond_score),
        "clash_penalty": float(clash_penalty),
        "total_score": float(total_score),
        "n_heavy_in_shell": int(len(near_heavy)),
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

    regions = {
        "A_candidate_K417-K422": CANDIDATE_RESIDUES,
        "B_MIDAS_pos_control":   MIDAS_RESIDUES,
        "C_alphaV_Cterm_neg_control": CONTROL_RESIDUES,
    }
    results: dict[str, dict[str, list]] = {
        rn: {"bent": [], "extended": []} for rn in regions
    }

    t_start = time.time()
    for label, frame_indices in [("bent", bent_idx), ("extended", ext_idx)]:
        for k, idx in enumerate(frame_indices):
            traj = md.load(str(LIB_DIR / names[idx]))
            hyd, hbond, heavy = classify_atoms(traj)
            for region_name, residues in regions.items():
                centroid = residue_centroid(traj, residues)
                rec = vina_proxy_score(centroid, hyd, hbond, heavy)
                rec["frame"] = names[idx]
                rec["state"] = label
                rec["cv0_A"] = float(cv0[idx])
                results[region_name][label].append(rec)
            print(f"  {label} {k+1}/{N_PER_STATE} ({names[idx]}, "
                  f"CV0={cv0[idx]:.2f} Å)")
    print(f"total: {time.time() - t_start:.1f}s")

    region_summary: dict[str, dict] = {}
    for region_name, region_data in results.items():
        bent_scores = np.array([r["total_score"] for r in region_data["bent"]])
        ext_scores = np.array([r["total_score"] for r in region_data["extended"]])
        delta_mean = float(ext_scores.mean() - bent_scores.mean())
        pooled_se = float(np.sqrt(bent_scores.var(ddof=1) / N_PER_STATE
                                  + ext_scores.var(ddof=1) / N_PER_STATE))
        welch_t = delta_mean / pooled_se if pooled_se > 1e-9 else 0.0
        bent_hyd = np.mean([r["n_hydrophobic"] for r in region_data["bent"]])
        bent_hb = np.mean([r["n_hbond"] for r in region_data["bent"]])
        bent_clash = np.mean([r["n_clash"] for r in region_data["bent"]])
        ext_hyd = np.mean([r["n_hydrophobic"] for r in region_data["extended"]])
        ext_hb = np.mean([r["n_hbond"] for r in region_data["extended"]])
        ext_clash = np.mean([r["n_clash"] for r in region_data["extended"]])
        region_summary[region_name] = {
            "bent_mean_score": float(bent_scores.mean()),
            "bent_std_score": float(bent_scores.std(ddof=1)),
            "extended_mean_score": float(ext_scores.mean()),
            "extended_std_score": float(ext_scores.std(ddof=1)),
            "delta_mean": delta_mean,
            "delta_pooled_se": pooled_se,
            "welch_t": float(welch_t),
            "delta_significant_alpha_005": bool(abs(welch_t) > 1.96),
            "bent_components": {"n_hydroph": float(bent_hyd),
                                "n_hbond": float(bent_hb),
                                "n_clash": float(bent_clash)},
            "extended_components": {"n_hydroph": float(ext_hyd),
                                    "n_hbond": float(ext_hb),
                                    "n_clash": float(ext_clash)},
            "frames": region_data,
        }
        print(f"\n{region_name}:")
        print(f"  bent     score = {bent_scores.mean():+.2f} ± {bent_scores.std(ddof=1):.2f}  "
              f"(hyd={bent_hyd:.0f}, hb={bent_hb:.0f}, clash={bent_clash:.0f})")
        print(f"  extended score = {ext_scores.mean():+.2f} ± {ext_scores.std(ddof=1):.2f}  "
              f"(hyd={ext_hyd:.0f}, hb={ext_hb:.0f}, clash={ext_clash:.0f})")
        print(f"  Δ = {delta_mean:+.2f}  (SE {pooled_se:.2f}; t = {welch_t:+.2f})")

    cand = region_summary["A_candidate_K417-K422"]
    midas = region_summary["B_MIDAS_pos_control"]
    ctrl = region_summary["C_alphaV_Cterm_neg_control"]
    cand_score_ext = cand["extended_mean_score"]
    midas_score_ext = midas["extended_mean_score"]
    ctrl_score_ext = ctrl["extended_mean_score"]
    cand_vs_ctrl_ext = cand_score_ext - ctrl_score_ext
    cand_vs_midas_ext = cand_score_ext - midas_score_ext

    cand_delta_sig = cand["delta_significant_alpha_005"]
    cand_delta = cand["delta_mean"]
    cand_above_ctrl = cand_score_ext > ctrl_score_ext + 1.0
    cand_above_midas = cand_score_ext > midas_score_ext

    if cand_delta_sig and cand_delta > 0 and cand_above_ctrl:
        verdict = ("VALIDATED — candidate Vina-proxy score increases significantly "
                   "during BC→EC and exceeds matched negative-control surface; "
                   "consistent with a cryptic-pocket-opening event")
    elif cand_delta_sig and cand_delta < 0:
        verdict = ("DOWNGRADED — candidate Vina-proxy score DECREASES during BC→EC "
                   "(opposite of cryptic-opening expectation)")
    elif (not cand_delta_sig) and cand_above_ctrl and not cand_above_midas:
        verdict = ("NEGATIVE / INTERMEDIATE — candidate druggability is STABLE across "
                   "BC→EC (Δ ns) and absolute extended-state score sits between MIDAS "
                   "(positive control) and a flat-surface negative control. The region "
                   "has some druggable atom-type features but its druggability is not "
                   "gated by extension. Consistent with obj-057 LIGSITE: K417-K422 "
                   "is conformationally responsive (loses pocket volume) but its atom-"
                   "type composition stays moderate/similar — neither a cryptic-opening "
                   "site nor a discrete druggable pocket.")
    elif (not cand_delta_sig) and not cand_above_ctrl:
        verdict = ("NULL — candidate Vina-proxy score is indistinguishable from negative "
                   "control; not druggable in any state")
    else:
        verdict = "AMBIGUOUS — see per-component breakdown"

    full_summary = {
        "regions": region_summary,
        "verdict": verdict,
        "comparisons": {
            "candidate_ext_score": cand_score_ext,
            "midas_ext_score": midas_score_ext,
            "control_ext_score": ctrl_score_ext,
            "candidate_minus_control_ext": cand_vs_ctrl_ext,
            "candidate_minus_midas_ext": cand_vs_midas_ext,
        },
        "params": {
            "n_per_state": N_PER_STATE,
            "clash_radius_A": CLASH_R,
            "hydrophobic_shell_A": [HYDROPH_INNER, HYDROPH_OUTER],
            "hbond_shell_A": [HBOND_INNER, HBOND_OUTER],
            "weight_hydrophobic": W_HYDROPH,
            "weight_hbond": W_HBOND,
            "weight_clash": W_CLASH,
            "algorithm": "Vina-proxy: atom-type weighted shell counts (hydroph + hbond - clash)",
        },
    }

    summary_path = OUT_DIR / "vina_proxy_scoring.json"
    with summary_path.open("w") as f:
        json.dump(full_summary, f, indent=2)
    print(f"\n=== VERDICT: {verdict} ===")
    print(f"\nWrote {summary_path}")

    plot(full_summary, results)
    print(f"Wrote {FIG}")
    return 0


def plot(summary: dict, results: dict) -> None:
    fig = plt.figure(figsize=(16, 6.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.65], wspace=0.35)
    region_order = ["A_candidate_K417-K422", "B_MIDAS_pos_control",
                    "C_alphaV_Cterm_neg_control"]
    region_labels = {
        "A_candidate_K417-K422": "A: β3 K417-K422\n(obj-056 candidate)",
        "B_MIDAS_pos_control": "B: MIDAS centroid\n(positive control)",
        "C_alphaV_Cterm_neg_control": "C: αV C-term flat\n(negative control)",
    }

    for ax_i, region in enumerate(region_order):
        ax = fig.add_subplot(gs[0, ax_i])
        s = summary["regions"][region]
        bent_sc = [r["total_score"] for r in results[region]["bent"]]
        ext_sc = [r["total_score"] for r in results[region]["extended"]]
        bent_cv0s = [r["cv0_A"] for r in results[region]["bent"]]
        ext_cv0s = [r["cv0_A"] for r in results[region]["extended"]]
        ax.scatter(bent_cv0s, bent_sc, s=80, color="#2166ac",
                   edgecolor="black", linewidth=0.6, label="bent", zorder=5)
        ax.scatter(ext_cv0s, ext_sc, s=80, color="#b2182b",
                   edgecolor="black", linewidth=0.6, label="extended", zorder=5)
        ax.axhline(s["bent_mean_score"], color="#2166ac", linestyle="--",
                   linewidth=1.0, alpha=0.7)
        ax.axhline(s["extended_mean_score"], color="#b2182b", linestyle="--",
                   linewidth=1.0, alpha=0.7)
        ax.axhline(0, color="black", linewidth=0.6, alpha=0.4)
        ax.set_xlabel("CV0 (Å)")
        ax.set_ylabel("Vina-proxy score (a.u.)")
        ax.set_title(region_labels[region], fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.3)

        delta = s["delta_mean"]
        se = s["delta_pooled_se"]
        t = s["welch_t"]
        sig = s["delta_significant_alpha_005"]
        bent_c = s["bent_components"]
        ext_c = s["extended_components"]
        ax.text(0.02, 0.97,
                f"Δ = {delta:+.2f} ± {se:.2f}\n"
                f"t = {t:+.2f}  ({'sig' if sig else 'ns'})\n\n"
                f"bent: hyd={bent_c['n_hydroph']:.0f} "
                f"hb={bent_c['n_hbond']:.0f} cl={bent_c['n_clash']:.0f}\n"
                f"ext:  hyd={ext_c['n_hydroph']:.0f} "
                f"hb={ext_c['n_hbond']:.0f} cl={ext_c['n_clash']:.0f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=8.0,
                family="monospace",
                bbox=dict(boxstyle="round",
                          facecolor=("#ccebc5" if sig and delta > 0
                                     else "#fbb4ae" if sig and delta < 0
                                     else "#ffffcc"),
                          edgecolor="#666"))

    ax = fig.add_subplot(gs[0, 3])
    ax.axis("off")
    cmp = summary["comparisons"]
    rows = [
        ["", ""],
        ["Region (extended state)", "score"],
        ["A candidate K417-K422", f"{cmp['candidate_ext_score']:+.2f}"],
        ["B MIDAS pos control", f"{cmp['midas_ext_score']:+.2f}"],
        ["C αV-Cterm neg control", f"{cmp['control_ext_score']:+.2f}"],
        ["", ""],
        ["A − C", f"{cmp['candidate_minus_control_ext']:+.2f}"],
        ["A − B (MIDAS)", f"{cmp['candidate_minus_midas_ext']:+.2f}"],
    ]
    table = ax.table(cellText=rows, loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    ax.set_title("Score comparison", fontsize=10, fontweight="bold")

    fig.suptitle(
        f"Vina-proxy ligand-fit scoring of obj-056 candidate "
        f"({N_PER_STATE}+{N_PER_STATE} frames; atom-type-weighted shell counts)\n"
        f"Verdict: {summary['verdict'][:130]}",
        fontsize=10, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    fig.savefig(FIG, dpi=140, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
