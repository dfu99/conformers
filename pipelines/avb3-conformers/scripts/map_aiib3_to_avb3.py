#!/usr/bin/env python3
"""Map αIIbβ3 string method pathway structures onto αVβ3 via domain alignment.

Uses the 19 PDB structures from Dasetty et al. (bioRxiv 2025) that trace the
bent-closed → extended-open pathway of αIIbβ3 integrin. Maps each structure
onto αVβ3 by domain-wise CA superposition, producing αVβ3 conformers along
the same activation pathway.

Approach:
1. Identify matching domains between αIIbβ3 (chains A/B) and αVβ3 (chains A/B)
2. For each αIIbβ3 pathway image, align αVβ3 domains to the corresponding
   αIIbβ3 domains using Kabsch superposition
3. Output transformed αVβ3 structures

Usage:
    python map_aiib3_to_avb3.py \
        --aiib3-dir path/to/image_structures_relaxed/ \
        --avb3-pdb path/to/AVB3_clean.pdb \
        --output-dir path/to/avb3_pathway/
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np


def parse_pdb_atoms(path: Path) -> list[dict]:
    """Parse ATOM records from PDB file."""
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        atoms.append({
            "line": line,
            "name": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "chain": line[21],
            "resseq": int(line[22:26]),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
        })
    return atoms


def get_ca_coords(atoms: list[dict], chain: str) -> tuple[np.ndarray, list[int]]:
    """Get CA coordinates and residue numbers for a chain."""
    cas = [(a["resseq"], np.array([a["x"], a["y"], a["z"]]))
           for a in atoms if a["name"] == "CA" and a["chain"] == chain]
    cas.sort()
    resseqs = [r for r, _ in cas]
    coords = np.array([c for _, c in cas])
    return coords, resseqs


def kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute optimal rotation/translation to align P onto Q.

    Returns (R, t_P, t_Q) such that: P_aligned = (P - t_P) @ R.T + t_Q
    """
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Pc, Qc


def transform_pdb_atoms(atoms: list[dict], chain: str, R: np.ndarray,
                        t_src: np.ndarray, t_dst: np.ndarray) -> list[dict]:
    """Apply rotation/translation to all atoms of a chain."""
    result = []
    for a in atoms:
        if a["chain"] != chain:
            continue
        xyz = np.array([a["x"], a["y"], a["z"]])
        xyz_new = (xyz - t_src) @ R.T + t_dst
        new_line = (a["line"][:30] +
                    f"{xyz_new[0]:8.3f}{xyz_new[1]:8.3f}{xyz_new[2]:8.3f}" +
                    a["line"][54:])
        result.append({**a, "line": new_line,
                       "x": xyz_new[0], "y": xyz_new[1], "z": xyz_new[2]})
    return result


def find_common_residue_range(resseqs_a: list[int], resseqs_b: list[int],
                               min_overlap: int = 50) -> list[int]:
    """Find overlapping residue numbers between two chains."""
    common = sorted(set(resseqs_a) & set(resseqs_b))
    if len(common) < min_overlap:
        # Fall back to positional alignment (first N residues)
        n = min(len(resseqs_a), len(resseqs_b), 200)
        return list(range(n))
    return common


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--aiib3-dir", type=Path, required=True,
                        help="Directory with αIIbβ3 image PDBs (image0.pdb - image18.pdb)")
    parser.add_argument("--avb3-pdb", type=Path, required=True,
                        help="Reference αVβ3 PDB (bent state)")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for mapped αVβ3 pathway structures")
    parser.add_argument("--aiib3-alpha-chain", default="A",
                        help="αIIb chain ID in pathway PDBs (default: A)")
    parser.add_argument("--aiib3-beta-chain", default="B",
                        help="β3 chain ID in pathway PDBs (default: B)")
    parser.add_argument("--avb3-alpha-chain", default="A",
                        help="αV chain ID in reference PDB (default: A)")
    parser.add_argument("--avb3-beta-chain", default="B",
                        help="β3 chain ID in reference PDB (default: B)")
    # Domain ranges for alignment (approximate, based on integrin domain architecture)
    # These define which residues to use for superposition per domain
    parser.add_argument("--align-mode", default="per-chain",
                        choices=["per-chain", "head-only", "full"],
                        help="Alignment mode (default: per-chain)")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load αVβ3 reference
    avb3_atoms = parse_pdb_atoms(args.avb3_pdb)
    avb3_ca_a, avb3_res_a = get_ca_coords(avb3_atoms, args.avb3_alpha_chain)
    avb3_ca_b, avb3_res_b = get_ca_coords(avb3_atoms, args.avb3_beta_chain)
    print(f"αVβ3 reference: chain {args.avb3_alpha_chain}={len(avb3_ca_a)} CA, "
          f"chain {args.avb3_beta_chain}={len(avb3_ca_b)} CA")

    # Load αIIbβ3 bent state (image0) for initial alignment reference
    img0_path = args.aiib3_dir / "image0.pdb"
    if not img0_path.exists():
        print(f"ERROR: {img0_path} not found")
        return 1

    img0_atoms = parse_pdb_atoms(img0_path)
    img0_ca_a, img0_res_a = get_ca_coords(img0_atoms, args.aiib3_alpha_chain)
    img0_ca_b, img0_res_b = get_ca_coords(img0_atoms, args.aiib3_beta_chain)
    print(f"αIIbβ3 image0: chain {args.aiib3_alpha_chain}={len(img0_ca_a)} CA, "
          f"chain {args.aiib3_beta_chain}={len(img0_ca_b)} CA")

    # Use positional alignment (truncate to shorter chain)
    # αIIbβ3 and αVβ3 β3 subunits are identical, α subunits differ
    n_align_a = min(len(avb3_ca_a), len(img0_ca_a))
    n_align_b = min(len(avb3_ca_b), len(img0_ca_b))
    print(f"Alignment residues: α={n_align_a}, β={n_align_b}")

    # Process each pathway image
    summary = []
    image_files = sorted(args.aiib3_dir.glob("image*.pdb"),
                         key=lambda p: int(re.search(r'(\d+)', p.stem).group(1)))

    for img_path in image_files:
        img_idx = int(re.search(r'(\d+)', img_path.stem).group(1))
        img_atoms = parse_pdb_atoms(img_path)
        img_ca_a, _ = get_ca_coords(img_atoms, args.aiib3_alpha_chain)
        img_ca_b, _ = get_ca_coords(img_atoms, args.aiib3_beta_chain)

        if args.align_mode == "per-chain":
            # Align each chain independently
            R_a, t_a_src, t_a_dst = kabsch_rotation(
                avb3_ca_a[:n_align_a], img_ca_a[:n_align_a])
            R_b, t_b_src, t_b_dst = kabsch_rotation(
                avb3_ca_b[:n_align_b], img_ca_b[:n_align_b])

            # Transform αVβ3 atoms per chain
            mapped_a = transform_pdb_atoms(avb3_atoms, args.avb3_alpha_chain,
                                           R_a, t_a_src, t_a_dst)
            mapped_b = transform_pdb_atoms(avb3_atoms, args.avb3_beta_chain,
                                           R_b, t_b_src, t_b_dst)
        else:
            # Full alignment (concatenate both chains)
            src = np.vstack([avb3_ca_a[:n_align_a], avb3_ca_b[:n_align_b]])
            dst = np.vstack([img_ca_a[:n_align_a], img_ca_b[:n_align_b]])
            R, t_src, t_dst = kabsch_rotation(src, dst)
            mapped_a = transform_pdb_atoms(avb3_atoms, args.avb3_alpha_chain,
                                           R, t_src, t_dst)
            mapped_b = transform_pdb_atoms(avb3_atoms, args.avb3_beta_chain,
                                           R, t_src, t_dst)

        # Write output PDB
        out_path = args.output_dir / f"avb3_pathway_image{img_idx:02d}.pdb"
        with open(out_path, "w") as f:
            f.write(f"REMARK Mapped from αIIbβ3 image{img_idx} via {args.align_mode} alignment\n")
            f.write(f"REMARK Source: {img_path.name}\n")
            for a in mapped_a + mapped_b:
                f.write(a["line"] + "\n")
            f.write("END\n")

        # Compute alignment RMSD
        mapped_ca_a = np.array([[a["x"], a["y"], a["z"]] for a in mapped_a if a["name"] == "CA"])
        mapped_ca_b = np.array([[a["x"], a["y"], a["z"]] for a in mapped_b if a["name"] == "CA"])

        rmsd_a = np.sqrt(np.mean(np.sum(
            (mapped_ca_a[:n_align_a] - img_ca_a[:n_align_a])**2, axis=1)))
        rmsd_b = np.sqrt(np.mean(np.sum(
            (mapped_ca_b[:n_align_b] - img_ca_b[:n_align_b])**2, axis=1)))

        summary.append({
            "image": img_idx, "source": img_path.name,
            "output": out_path.name,
            "rmsd_alpha": float(rmsd_a), "rmsd_beta": float(rmsd_b),
        })
        print(f"  image{img_idx:2d} → {out_path.name}  "
              f"RMSD: α={rmsd_a:.1f}Å  β={rmsd_b:.1f}Å")

    # Write summary
    summary_path = args.output_dir / "mapping_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"\nMapped {len(summary)} structures to {args.output_dir}")
    print(f"Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
