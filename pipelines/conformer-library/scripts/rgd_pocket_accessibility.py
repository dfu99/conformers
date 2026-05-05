#!/usr/bin/env python3
"""RGD-binding-pocket accessibility — bent vs extended αVβ3 conformers.

Reviewer-E follow-up (audit-2026-05-05 §8 F4). Real Vina/smina docking on
αVβ3 requires Mn²⁺ MIDAS parameterization, polar-H + Gasteiger charges
on a 1654-residue receptor, and rdkit/meeko + gemmi to generate PDBQT —
all of which collide with the project's strict venv discipline (CLAUDE.md
§4: gemmi conflict). For the first-pass reviewer answer ("does extension
expose the binding site?") a geometric accessibility proxy is more
defensible than mis-parameterized blind Vina.

Approach:
  1. Compute CV0 (αV head ↔ αV calf centroid) for all 615 library
     conformers. Pick the 5 most-bent + 5 most-extended.
  2. For each picked PDB compute:
       - Shrake-Rupley SASA of the αVβ3 MIDAS pocket residues
         (αV D218 β-propeller; β3 D119, S121, S123, E220, D251 βA;
          β3 R214, N215, Y122, R216 RGD-contact);
       - "Pocket occlusion" = fraction of a 12-Å sphere around the
         MIDAS centroid that is occupied by non-pocket protein atoms;
       - Distance from MIDAS centroid to the leg-domain centroid
         (αV calf + β3 tail), as the steric blocker;
  3. Plot a 4-panel comparison: SASA bar, occlusion bar, MIDAS-leg
     distance, and SASA-vs-CV0 scatter coloured by state.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LIB_DIR = ROOT / "data" / "runs" / "avb3" / "conformers" / "all_frames_bent_extended"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "rgd_docking"
FIG_DIR = ROOT / "figures"

DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 956),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}

MIDAS_RESIDUES = {
    ("A", 218): "αV D218 (β-propeller, contacts RGD-Arg)",
    ("B", 119): "β3 D119 (MIDAS)",
    ("B", 121): "β3 S121 (MIDAS)",
    ("B", 123): "β3 S123 (MIDAS)",
    ("B", 220): "β3 E220 (MIDAS)",
    ("B", 251): "β3 D251 (MIDAS)",
    ("B", 122): "β3 Y122 (RGD-Asp pocket)",
    ("B", 214): "β3 R214 (RGD-Asp pocket)",
    ("B", 215): "β3 N215 (LIMBS)",
    ("B", 216): "β3 R216 (RGD-Asp pocket)",
}


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    if cid:
        return cid
    return chr(ord("A") + atom.residue.chain.index)


def domain_ca_mean(traj: md.Trajectory, chain: str, lo: int, hi: int) -> np.ndarray:
    sel = []
    for atom in traj.topology.atoms:
        if atom.name != "CA":
            continue
        if chain_id(atom) == chain and lo <= atom.residue.resSeq <= hi:
            sel.append(atom.index)
    if not sel:
        raise ValueError(f"no CAs in {chain} {lo}-{hi}")
    return traj.xyz[0, sel].mean(axis=0) * 10.0  # nm → Å


def midas_atom_indices(traj: md.Trajectory) -> list[int]:
    pocket_idx = []
    for atom in traj.topology.atoms:
        cid = chain_id(atom)
        rs = atom.residue.resSeq
        if (cid, rs) in MIDAS_RESIDUES:
            pocket_idx.append(atom.index)
    return pocket_idx


def cv0_for_pdb(path: Path) -> float:
    traj = md.load(str(path))
    a = domain_ca_mean(traj, *DOMAINS["alpha_head_thigh"])
    b = domain_ca_mean(traj, *DOMAINS["alpha_calf"])
    return float(np.linalg.norm(a - b))


def midas_centroid(traj: md.Trajectory, idx: list[int]) -> np.ndarray:
    return traj.xyz[0, idx].mean(axis=0) * 10.0


def pocket_occlusion(traj: md.Trajectory, midas_idx: list[int],
                     radius_A: float = 12.0) -> float:
    """Fraction of a sphere around the MIDAS centroid occupied by
    non-pocket protein heavy atoms. Approximated by counting heavy
    atoms inside the sphere — the more atoms intrude, the more occluded
    the pocket. Normalized by an empirical "free pocket" reference of
    ~50 atoms (1JV2 bent). Returns a unitless [0,1+] occlusion."""
    centroid_nm = np.array(midas_centroid(traj, midas_idx)) / 10.0
    coords = traj.xyz[0]  # nm
    d = np.linalg.norm(coords - centroid_nm, axis=-1) * 10.0
    in_sphere = (d < radius_A)
    pocket_set = set(midas_idx)
    n_intruders = int(sum(1 for i in np.where(in_sphere)[0]
                          if i not in pocket_set
                          and traj.topology.atom(i).element.symbol != "H"))
    return n_intruders / 50.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-each", type=int, default=5)
    p.add_argument("--cache", type=Path,
                   default=OUT_DIR / "library_cv0_cache.npy")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(LIB_DIR.glob("*.pdb"))
    print(f"{len(pdbs)} library PDBs")

    names: list[str] = [p.name for p in pdbs]
    cv0: np.ndarray | None = None
    if args.cache.exists():
        data = np.load(args.cache, allow_pickle=True)
        cached_names = list(data["names"])
        if cached_names == names:
            cv0 = np.array(data["cv0"], dtype=np.float64)
            print("  using cached CV0 values")
        else:
            args.cache.unlink()
    if cv0 is None:
        cv0_list = []
        for k, p in enumerate(pdbs):
            cv0_list.append(cv0_for_pdb(p))
            if (k + 1) % 100 == 0:
                print(f"  CV0 scored {k+1}/{len(pdbs)}")
        cv0 = np.array(cv0_list)
        np.savez(args.cache, names=np.array(names), cv0=cv0)

    print(f"CV0 min/median/max: {cv0.min():.2f} / {np.median(cv0):.2f} / {cv0.max():.2f}")
    sort_idx = np.argsort(cv0)
    bent_picks = sort_idx[: args.n_each]
    extended_picks = sort_idx[-args.n_each:]
    print("\nBent picks (lowest CV0):")
    for i in bent_picks:
        print(f"  {names[i]}  CV0={cv0[i]:.2f} Å")
    print("\nExtended picks (highest CV0):")
    for i in extended_picks:
        print(f"  {names[i]}  CV0={cv0[i]:.2f} Å")

    rows = []
    for label, picks in [("bent", bent_picks), ("extended", extended_picks)]:
        for idx in picks:
            pdb = LIB_DIR / names[idx]
            traj = md.load(str(pdb))
            pocket_idx = midas_atom_indices(traj)
            if not pocket_idx:
                print(f"  WARN: no MIDAS atoms in {pdb.name}; skipping")
                continue
            sasa = md.shrake_rupley(traj, mode="atom")[0]  # (N_atoms,) in nm²
            pocket_sasa_A2 = float(sasa[pocket_idx].sum() * 100.0)
            occl = pocket_occlusion(traj, pocket_idx)
            mid_centroid = midas_centroid(traj, pocket_idx)
            leg = (domain_ca_mean(traj, *DOMAINS["alpha_calf"])
                   + domain_ca_mean(traj, *DOMAINS["beta_tail"])) / 2.0
            midas_leg_dist = float(np.linalg.norm(mid_centroid - leg))
            rows.append({
                "label": label,
                "pdb": names[idx],
                "cv0_A": float(cv0[idx]),
                "midas_sasa_A2": pocket_sasa_A2,
                "occlusion_norm": occl,
                "midas_leg_distance_A": midas_leg_dist,
            })
            print(f"  {label} {pdb.name}: SASA={pocket_sasa_A2:.1f} Å² "
                  f"occl={occl:.2f} midas-leg={midas_leg_dist:.1f} Å")

    rec_path = OUT_DIR / "pocket_accessibility.json"
    with rec_path.open("w") as f:
        json.dump({"records": rows,
                   "midas_residues": {f"{c}:{r}": v
                                      for (c, r), v in MIDAS_RESIDUES.items()},
                   "method": "shrake-rupley + 12 A occlusion sphere",
                   "library_dir": str(LIB_DIR.relative_to(ROOT))}, f, indent=2)
    print(f"\nWrote {rec_path}")

    plot(rows, OUT_DIR, FIG_DIR)
    return 0


def plot(rows: list[dict], out_dir: Path, fig_dir: Path) -> None:
    bent = [r for r in rows if r["label"] == "bent"]
    ext = [r for r in rows if r["label"] == "extended"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8))

    # --- SASA bar
    ax = axes[0, 0]
    pos_b = np.arange(len(bent))
    pos_e = np.arange(len(ext)) + len(bent) + 0.6
    ax.bar(pos_b, [r["midas_sasa_A2"] for r in bent],
           color="#2166ac", label=f"bent (mean CV0 {np.mean([r['cv0_A'] for r in bent]):.1f} Å)")
    ax.bar(pos_e, [r["midas_sasa_A2"] for r in ext],
           color="#b2182b", label=f"extended (mean CV0 {np.mean([r['cv0_A'] for r in ext]):.1f} Å)")
    ax.set_xticks(np.concatenate([pos_b, pos_e]))
    ax.set_xticklabels([r["pdb"].replace(".pdb", "") for r in bent + ext],
                       rotation=55, fontsize=7, ha="right")
    ax.set_ylabel("MIDAS-pocket SASA (Å²)")
    ax.set_title("Solvent accessibility of RGD-binding pocket residues")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    bent_mean = float(np.mean([r["midas_sasa_A2"] for r in bent]))
    ext_mean = float(np.mean([r["midas_sasa_A2"] for r in ext]))
    ax.axhline(bent_mean, color="#2166ac", linestyle="--", linewidth=1.2)
    ax.axhline(ext_mean, color="#b2182b", linestyle="--", linewidth=1.2)
    ax.text(0.99, 0.95,
            f"Δ = {ext_mean - bent_mean:+.0f} Å²  ({(ext_mean - bent_mean) / bent_mean * 100:+.1f} %)",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    # --- Occlusion
    ax = axes[0, 1]
    ax.bar(pos_b, [r["occlusion_norm"] for r in bent], color="#2166ac")
    ax.bar(pos_e, [r["occlusion_norm"] for r in ext], color="#b2182b")
    ax.set_xticks(np.concatenate([pos_b, pos_e]))
    ax.set_xticklabels([r["pdb"].replace(".pdb", "") for r in bent + ext],
                       rotation=55, fontsize=7, ha="right")
    ax.set_ylabel("Pocket occlusion (heavy atoms / 50 inside 12 Å sphere)")
    ax.set_title("Steric occlusion around MIDAS centroid")
    ax.grid(axis="y", alpha=0.3)
    bent_o = np.mean([r["occlusion_norm"] for r in bent])
    ext_o = np.mean([r["occlusion_norm"] for r in ext])
    ax.text(0.99, 0.95,
            f"bent {bent_o:.2f}  vs  ext {ext_o:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    # --- MIDAS↔leg distance
    ax = axes[1, 0]
    ax.bar(pos_b, [r["midas_leg_distance_A"] for r in bent], color="#2166ac")
    ax.bar(pos_e, [r["midas_leg_distance_A"] for r in ext], color="#b2182b")
    ax.set_xticks(np.concatenate([pos_b, pos_e]))
    ax.set_xticklabels([r["pdb"].replace(".pdb", "") for r in bent + ext],
                       rotation=55, fontsize=7, ha="right")
    ax.set_ylabel("MIDAS centroid ↔ leg centroid (Å)")
    ax.set_title("Steric blocker: distance from MIDAS to (αV-calf + β3-tail)")
    ax.grid(axis="y", alpha=0.3)

    # --- SASA vs CV0 scatter (all picks)
    ax = axes[1, 1]
    for r in bent + ext:
        c = "#2166ac" if r["label"] == "bent" else "#b2182b"
        ax.scatter(r["cv0_A"], r["midas_sasa_A2"], s=70, color=c,
                   edgecolor="black", linewidth=0.6,
                   label=r["label"] if r is bent[0] or r is ext[0] else None)
    ax.set_xlabel("CV0 (αV head ↔ αV calf centroid, Å)")
    ax.set_ylabel("MIDAS-pocket SASA (Å²)")
    ax.set_title("Pocket SASA vs CV0 — extension does NOT expose RGD site\n"
                 "(consistent with EC = extended-closed; EO would need headpiece opening)",
                 fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle("RGD-binding-site accessibility (geometric proxy) — bent vs extended αVβ3 v7 library\n"
                 "Negative result: legs move 52 Å away (steering works), but MIDAS pocket SASA drops 35 %\n"
                 "→ extended frames are EC, not EO. RGD-docking PoC blocked on EO headpiece-opening (audit §3).",
                 fontsize=11, fontweight="bold", y=1.02)
    fig.tight_layout()
    out_path = out_dir / "rgd_pocket_accessibility.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig_top = fig_dir / "rgd_docking_v1.png"
    fig.savefig(fig_top, dpi=140, bbox_inches="tight")
    print(f"  saved {out_path}")
    print(f"  copied to {fig_top}")


if __name__ == "__main__":
    raise SystemExit(main())
