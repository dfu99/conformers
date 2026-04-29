#!/usr/bin/env python3
"""First-principles bent-fraction prediction across integrin heterodimers.

For each available bent ectodomain PDB, compute:
1. Genu contact area — buried-SASA between the headpiece and the
   leg domains. Larger contact ⇒ more stabilized bent ⇒ higher
   predicted-bent fraction.
2. Head-tail centroid distance (analog of αVβ3 CV0).
3. Maximum molecular dimension (long-axis radius of gyration).
4. Per-chain residue counts.

Predict: more bent contact + smaller head-tail distance ⇒ higher
bent-state stability ⇒ higher experimental bent population.

Headpiece-only structures (no legs) are still parsed but flagged as
"insufficient for bent-fraction prediction".

Usage:
    python multi_integrin_first_principles.py \\
        --pdb-dir data/multi_integrin/raw_pdbs \\
        --output-dir results/multi_integrin
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


HETERODIMER_REGISTRY = {
    "1JV2": {
        "name": "αVβ3 (bent crystal)",
        "alpha_chain": "A",
        "beta_chain": "B",
        "alpha_full_length": 957,
        "beta_full_length": 690,
        "expected_state": "bent",
        "is_full_ectodomain": True,
    },
    "3FCS": {
        "name": "αIIbβ3 (bent crystal)",
        "alpha_chain": "A",
        "beta_chain": "B",
        "alpha_full_length": 1009,
        "beta_full_length": 690,
        "expected_state": "bent",
        "is_full_ectodomain": True,
    },
    "3VI4": {
        "name": "α5β1 headpiece (bound peptide)",
        "alpha_chain": "A",
        "beta_chain": "B",
        "alpha_full_length": 991,
        "beta_full_length": 778,
        "expected_state": "headpiece-only",
        "is_full_ectodomain": False,
    },
    "4UM9": {
        "name": "αVβ6 head + peptide",
        "alpha_chain": "A",
        "beta_chain": "B",
        "alpha_full_length": 957,
        "beta_full_length": 707,
        "expected_state": "headpiece-only",
        "is_full_ectodomain": False,
    },
    "5FFG": {
        "name": "αVβ6 headpiece (apo)",
        "alpha_chain": "A",
        "beta_chain": "B",
        "alpha_full_length": 957,
        "beta_full_length": 707,
        "expected_state": "headpiece-only",
        "is_full_ectodomain": False,
    },
    "6OM2": {
        "name": "αVβ8 headpiece + LAP",
        "alpha_chain": "A",
        "beta_chain": "B",
        "alpha_full_length": 957,
        "beta_full_length": 690,
        "expected_state": "headpiece-only",
        "is_full_ectodomain": False,
    },
    "6UJC": {
        "name": "αVβ8 headpiece + Fab",
        "alpha_chain": "A",
        "beta_chain": "B",
        "alpha_full_length": 957,
        "beta_full_length": 690,
        "expected_state": "headpiece-only",
        "is_full_ectodomain": False,
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdb-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def get_chain_atoms(traj, chain_id: str):
    chain_indices = []
    for c in traj.topology.chains:
        cid = c.chain_id or chr(ord("A") + c.index)
        if cid == chain_id:
            for a in c.atoms:
                chain_indices.append(a.index)
    return np.array(chain_indices, dtype=int)


def get_first_protein_chain_with_aa_count(traj, expected_aa: int, tol: int = 50):
    """Find first chain whose AA count is within `tol` of expected_aa."""
    for c in traj.topology.chains:
        n_aa = sum(1 for r in c.residues if r.is_protein)
        if abs(n_aa - expected_aa) <= tol:
            return c
    return None


def compute_features(pdb_path: Path, meta: dict):
    print(f"\n=== {pdb_path.name} — {meta['name']} ===")
    traj = md.load(str(pdb_path))
    out = {"pdb_id": pdb_path.stem, **meta}

    # Identify the alpha & beta chains by position (first big chain is
    # usually α, second is β). Trust author chain IDs first.
    alpha_atoms = get_chain_atoms(traj, meta["alpha_chain"])
    beta_atoms = get_chain_atoms(traj, meta["beta_chain"])
    if len(alpha_atoms) == 0 or len(beta_atoms) == 0:
        # Fall back: find chains ranked by AA count.
        chain_aa = []
        for c in traj.topology.chains:
            n_aa = sum(1 for r in c.residues if r.is_protein)
            if n_aa > 100:
                chain_aa.append((c, n_aa))
        chain_aa.sort(key=lambda t: -t[1])
        if len(chain_aa) >= 2:
            alpha_chain_obj = chain_aa[0][0]
            beta_chain_obj = chain_aa[1][0]
            alpha_atoms = np.array([a.index for a in alpha_chain_obj.atoms], dtype=int)
            beta_atoms = np.array([a.index for a in beta_chain_obj.atoms], dtype=int)

    # Residue counts.
    alpha_resids = sorted(
        {traj.topology.atom(int(i)).residue.resSeq for i in alpha_atoms}
    )
    beta_resids = sorted(
        {traj.topology.atom(int(i)).residue.resSeq for i in beta_atoms}
    )
    alpha_n = len(alpha_resids)
    beta_n = len(beta_resids)
    alpha_coverage = alpha_n / meta["alpha_full_length"]
    beta_coverage = beta_n / meta["beta_full_length"]
    print(f"  α chain {meta['alpha_chain']}: {alpha_n} aa "
          f"(of {meta['alpha_full_length']}, coverage {alpha_coverage:.2f})")
    print(f"  β chain {meta['beta_chain']}: {beta_n} aa "
          f"(of {meta['beta_full_length']}, coverage {beta_coverage:.2f})")
    out.update(
        alpha_residues=alpha_n,
        beta_residues=beta_n,
        alpha_coverage=round(alpha_coverage, 3),
        beta_coverage=round(beta_coverage, 3),
    )

    # Centroid distances and Rg.
    ca = traj.topology.select("name CA")
    if len(ca) == 0:
        print("  No CA atoms; skipping geometry")
        return out
    ca_xyz = traj.xyz[0, ca] * 10.0  # nm -> Å
    ca_chain = np.array([traj.topology.atom(int(i)).residue.chain.index for i in ca])
    ca_alpha_mask = np.array(
        [traj.topology.atom(int(i)).residue.chain.chain_id ==
         meta["alpha_chain"] for i in ca]
    )
    ca_beta_mask = np.array(
        [traj.topology.atom(int(i)).residue.chain.chain_id ==
         meta["beta_chain"] for i in ca]
    )
    if not ca_alpha_mask.any():
        # fallback by chain-aa size
        idx_by_chain_count = {}
        for i_local, i in enumerate(ca):
            c = traj.topology.atom(int(i)).residue.chain
            idx_by_chain_count.setdefault(c.index, []).append(i_local)
        sizes = sorted(idx_by_chain_count.items(), key=lambda kv: -len(kv[1]))
        if len(sizes) >= 2:
            alpha_local = sizes[0][1]
            beta_local = sizes[1][1]
            ca_alpha_mask = np.zeros(len(ca), dtype=bool)
            ca_alpha_mask[alpha_local] = True
            ca_beta_mask = np.zeros(len(ca), dtype=bool)
            ca_beta_mask[beta_local] = True

    alpha_xyz = ca_xyz[ca_alpha_mask]
    beta_xyz = ca_xyz[ca_beta_mask]
    overall_xyz = np.concatenate([alpha_xyz, beta_xyz], axis=0)

    # Long-axis dimension via PCA.
    com = overall_xyz.mean(axis=0)
    rel = overall_xyz - com
    cov = rel.T @ rel / len(rel)
    eig_w, eig_v = np.linalg.eigh(cov)
    long_axis = eig_v[:, -1]  # largest eigenvalue
    proj = rel @ long_axis
    long_axis_extent = float(proj.max() - proj.min())
    rg = float(np.sqrt(((rel ** 2).sum(axis=1)).mean()))
    print(f"  Rg = {rg:.1f} Å, long-axis extent = {long_axis_extent:.1f} Å")
    out.update(rg_A=round(rg, 2), long_axis_extent_A=round(long_axis_extent, 2))

    # Head-tail centroid distance: head = first 1/3 of α + first 1/3 of β,
    # tail = last 1/3 of α + last 1/3 of β. This is a structure-agnostic
    # analog to our αVβ3 CV0.
    a_n = len(alpha_xyz)
    b_n = len(beta_xyz)
    head_alpha = alpha_xyz[: a_n // 3]
    head_beta = beta_xyz[: b_n // 3]
    tail_alpha = alpha_xyz[-(a_n // 3):]
    tail_beta = beta_xyz[-(b_n // 3):]
    if len(head_alpha) > 0 and len(tail_alpha) > 0:
        head_com_a = head_alpha.mean(axis=0)
        tail_com_a = tail_alpha.mean(axis=0)
        head_com_b = head_beta.mean(axis=0) if len(head_beta) > 0 else head_com_a
        tail_com_b = tail_beta.mean(axis=0) if len(tail_beta) > 0 else tail_com_a
        head_com = (head_com_a + head_com_b) / 2
        tail_com = (tail_com_a + tail_com_b) / 2
        head_tail_dist = float(np.linalg.norm(head_com - tail_com))
        alpha_head_to_alpha_tail = float(np.linalg.norm(head_com_a - tail_com_a))
        beta_head_to_beta_tail = float(np.linalg.norm(head_com_b - tail_com_b))
        head_alpha_to_head_beta = float(np.linalg.norm(head_com_a - head_com_b))
        print(f"  head-tail (combined) = {head_tail_dist:.1f} Å")
        print(f"  α head-tail = {alpha_head_to_alpha_tail:.1f} Å, "
              f"β head-tail = {beta_head_to_beta_tail:.1f} Å, "
              f"α-head ↔ β-head = {head_alpha_to_head_beta:.1f} Å")
        out.update(
            head_tail_distance_A=round(head_tail_dist, 2),
            alpha_head_to_alpha_tail_A=round(alpha_head_to_alpha_tail, 2),
            beta_head_to_beta_tail_A=round(beta_head_to_beta_tail, 2),
            head_alpha_to_head_beta_A=round(head_alpha_to_head_beta, 2),
        )

    # Genu contact area: only meaningful for full ectodomain. Compute as
    # buried-SASA between the headpiece (residue index < 1/3 of chain)
    # and the legs (residue index > 1/3). For headpiece-only structures
    # this is meaningless; flag and skip.
    if meta["is_full_ectodomain"]:
        try:
            # Build two sub-trajectories: head (residue idx < n/3 in each chain)
            # and the leg (residue idx >= 2n/3) — middle 1/3 is the genu
            # boundary and gets included in the head residue subset.
            head_atoms = []
            leg_atoms = []
            for c in traj.topology.chains:
                cid = c.chain_id or chr(ord("A") + c.index)
                if cid not in (meta["alpha_chain"], meta["beta_chain"]):
                    continue
                residues = [r for r in c.residues if r.is_protein]
                cut1 = len(residues) // 3
                cut2 = 2 * len(residues) // 3
                for r in residues[:cut1]:
                    for a in r.atoms:
                        head_atoms.append(a.index)
                for r in residues[cut2:]:
                    for a in r.atoms:
                        leg_atoms.append(a.index)
            if not head_atoms or not leg_atoms:
                raise RuntimeError("Empty head or leg atom set")
            # SASA of head alone (other atoms removed), leg alone, and combined.
            sasa_combined = md.shrake_rupley(traj, mode="atom")[0]
            head_atoms_a = np.array(head_atoms, dtype=int)
            leg_atoms_a = np.array(leg_atoms, dtype=int)
            sasa_head_part = float(sasa_combined[head_atoms_a].sum()) * 100.0  # nm² → Å²
            sasa_leg_part = float(sasa_combined[leg_atoms_a].sum()) * 100.0
            # SASA of head alone
            head_only = traj.atom_slice(head_atoms_a)
            leg_only = traj.atom_slice(leg_atoms_a)
            sasa_head_alone = float(md.shrake_rupley(head_only, mode="atom").sum()) * 100.0
            sasa_leg_alone = float(md.shrake_rupley(leg_only, mode="atom").sum()) * 100.0
            buried_total = max(0.0, sasa_head_alone + sasa_leg_alone - (sasa_head_part + sasa_leg_part))
            buried_head = max(0.0, sasa_head_alone - sasa_head_part)
            buried_leg = max(0.0, sasa_leg_alone - sasa_leg_part)
            print(f"  buried head SASA = {buried_head:.0f} Å², "
                  f"buried leg SASA = {buried_leg:.0f} Å², "
                  f"total head↔leg buried = {buried_total:.0f} Å²")
            out.update(
                head_only_sasa_A2=round(sasa_head_alone, 1),
                leg_only_sasa_A2=round(sasa_leg_alone, 1),
                head_in_complex_sasa_A2=round(sasa_head_part, 1),
                leg_in_complex_sasa_A2=round(sasa_leg_part, 1),
                head_leg_buried_sasa_A2=round(buried_total, 1),
            )
        except Exception as e:
            print(f"  SASA computation failed: {e}")

    return out


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for pdb_id, meta in HETERODIMER_REGISTRY.items():
        pdb_path = args.pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            print(f"Skipping {pdb_id}: not downloaded")
            continue
        try:
            row = compute_features(pdb_path, meta)
            rows.append(row)
        except Exception as e:
            print(f"  ERROR on {pdb_id}: {e}")

    out_json = args.output_dir / "first_principles.json"
    with out_json.open("w") as f:
        json.dump(rows, f, indent=2)
    print(f"\nWrote {out_json}")

    # Markdown summary
    md_path = args.output_dir / "first_principles.md"
    with md_path.open("w") as f:
        f.write("# Multi-integrin first-principles bent-state features\n\n")
        f.write(
            "Bent-state geometry from RCSB-deposited heterodimer ectodomain "
            "PDBs. Features are computed without reference to experimental "
            "activation kinetics or HS-AFM data.\n\n"
        )
        f.write(
            "| PDB | name | full ectodomain | α / β residues | Rg (Å) | "
            "long axis (Å) | head-tail (Å) | α head↔tail (Å) | β head↔tail (Å) | "
            "head↔head (Å) | head-leg buried SASA (Å²) |\n"
        )
        f.write("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for r in rows:
            full = "yes" if r.get("is_full_ectodomain") else "no"
            sasa = (
                f"{r['head_leg_buried_sasa_A2']:.0f}"
                if "head_leg_buried_sasa_A2" in r
                else "—"
            )
            f.write(
                f"| {r['pdb_id']} | {r['name']} | {full} | "
                f"{r.get('alpha_residues', '—')}/{r.get('beta_residues', '—')} | "
                f"{r.get('rg_A', '—')} | "
                f"{r.get('long_axis_extent_A', '—')} | "
                f"{r.get('head_tail_distance_A', '—')} | "
                f"{r.get('alpha_head_to_alpha_tail_A', '—')} | "
                f"{r.get('beta_head_to_beta_tail_A', '—')} | "
                f"{r.get('head_alpha_to_head_beta_A', '—')} | "
                f"{sasa} |\n"
            )
        f.write("\n## Predicted bent-fraction ranking\n\n")
        f.write(
            "Bent-fraction prediction uses head-leg buried SASA as the "
            "primary stabilizer (more contact ⇒ stable bent ⇒ higher "
            "experimental bent population), tied with smaller head-tail "
            "distance (more compact ⇒ bent). Headpiece-only structures are "
            "excluded since they cannot constrain bent-fraction.\n\n"
        )
        full_rows = [r for r in rows if r.get("is_full_ectodomain")
                     and "head_leg_buried_sasa_A2" in r]
        full_rows.sort(
            key=lambda r: (-r["head_leg_buried_sasa_A2"], r.get("head_tail_distance_A", 1e9))
        )
        f.write(
            "| rank | PDB | name | head-leg buried SASA (Å²) | head-tail (Å) | predicted bent fraction (relative) |\n"
        )
        f.write("|---:|---|---|---:|---:|---|\n")
        for i, r in enumerate(full_rows, 1):
            label = "highest" if i == 1 else f"#{i}"
            f.write(
                f"| {i} | {r['pdb_id']} | {r['name']} | "
                f"{r['head_leg_buried_sasa_A2']:.0f} | "
                f"{r.get('head_tail_distance_A', '—')} | {label} |\n"
            )
    print(f"Wrote {md_path}")

    # Bar chart of buried-SASA + head-tail distance side by side.
    full_rows = [r for r in rows if r.get("is_full_ectodomain")
                 and "head_leg_buried_sasa_A2" in r]
    if full_rows:
        labels = [f"{r['pdb_id']}\n{r['name'].split('(')[0].strip()}" for r in full_rows]
        sasa = [r["head_leg_buried_sasa_A2"] for r in full_rows]
        ht = [r.get("head_tail_distance_A", 0.0) for r in full_rows]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].bar(labels, sasa, color="tab:blue", edgecolor="black")
        axes[0].set_ylabel("Head-leg buried SASA (Å²)")
        axes[0].set_title("Bent-state head/leg interface contact")
        for x, v in zip(range(len(sasa)), sasa):
            axes[0].text(x, v + max(sasa) * 0.01, f"{v:.0f}", ha="center", fontsize=9)
        axes[0].grid(alpha=0.3, axis="y")
        axes[1].bar(labels, ht, color="tab:orange", edgecolor="black")
        axes[1].set_ylabel("Head-tail centroid distance (Å)")
        axes[1].set_title("Compactness")
        for x, v in zip(range(len(ht)), ht):
            axes[1].text(x, v + max(ht) * 0.01, f"{v:.1f}", ha="center", fontsize=9)
        axes[1].grid(alpha=0.3, axis="y")
        fig.suptitle("First-principles bent-state features — full ectodomain only")
        fig.tight_layout()
        fig_path = args.output_dir / "bent_state_features.png"
        fig.savefig(str(fig_path), dpi=140)
        plt.close(fig)
        print(f"Wrote {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
