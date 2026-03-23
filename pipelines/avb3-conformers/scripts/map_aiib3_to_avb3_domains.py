#!/usr/bin/env python3
"""Map αIIbβ3 pathway structures onto αVβ3 via domain-by-domain alignment.

Unlike the per-chain approach, this aligns each integrin domain independently
(head, calf, legs) so that internal domain structure is preserved while
inter-domain orientations follow the αIIbβ3 activation pathway.

Domains (αVβ3 residue ranges, 1-indexed PDB numbering):
  Chain A (αV):
    - head+thigh: 1-435
    - calf: 436-741
    - tail: 742-962
  Chain B (β3):
    - head+hybrid+EGF1: 1-352  (shared β3 subunit)
    - EGF2-4+tail: 353-692

For αIIbβ3, approximate domain boundaries (from structural homology):
  Chain A (αIIb):
    - head+thigh: ~1-450
    - calf: ~451-750
    - tail: ~751-999
  Chain B (β3):
    - head+hybrid+EGF1: 1-352  (identical subunit)
    - EGF2-4+tail: 353-762

Usage:
    python map_aiib3_to_avb3_domains.py \
        --aiib3-dir image_structures_relaxed/ \
        --avb3-pdb AVB3_clean.pdb \
        --output-dir avb3_domain_mapped/
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np


def parse_pdb_atoms(path: Path) -> list[dict]:
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


def get_ca_by_resrange(atoms: list[dict], chain: str, start: int, end: int) -> np.ndarray:
    """Get CA coords for a residue range (1-indexed, inclusive)."""
    cas = []
    for a in atoms:
        if a["name"] == "CA" and a["chain"] == chain and start <= a["resseq"] <= end:
            cas.append([a["x"], a["y"], a["z"]])
    return np.array(cas) if cas else np.empty((0, 3))


def kabsch_rotation(P: np.ndarray, Q: np.ndarray):
    """Optimal rotation to align P onto Q. Returns R, center_P, center_Q."""
    Pc, Qc = P.mean(0), Q.mean(0)
    P0, Q0 = P - Pc, Q - Qc
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    rmsd = float(np.sqrt(np.mean(np.sum((P0 @ R.T - Q0) ** 2, 1))))
    return R, Pc, Qc, rmsd


def transform_atoms(atoms: list[dict], chain: str, start: int, end: int,
                    R: np.ndarray, t_src: np.ndarray, t_dst: np.ndarray) -> list[dict]:
    """Apply rotation/translation to all atoms in a chain/residue range."""
    result = []
    for a in atoms:
        if a["chain"] != chain or not (start <= a["resseq"] <= end):
            continue
        xyz = np.array([a["x"], a["y"], a["z"]])
        xyz_new = (xyz - t_src) @ R.T + t_dst
        new_line = (a["line"][:30] +
                    f"{xyz_new[0]:8.3f}{xyz_new[1]:8.3f}{xyz_new[2]:8.3f}" +
                    a["line"][54:])
        result.append({**a, "line": new_line,
                       "x": xyz_new[0], "y": xyz_new[1], "z": xyz_new[2]})
    return result


# Domain definitions
AVB3_DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 962),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}

# αIIbβ3 approximate domain boundaries (from structural homology)
AIIB3_DOMAINS = {
    "alpha_head_thigh": ("A", 1, 450),
    "alpha_calf":       ("A", 451, 750),
    "alpha_tail":       ("A", 751, 999),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 762),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--aiib3-dir", type=Path, required=True)
    parser.add_argument("--avb3-pdb", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    avb3_atoms = parse_pdb_atoms(args.avb3_pdb)
    print(f"αVβ3 reference loaded: {sum(1 for a in avb3_atoms if a['name']=='CA')} CA atoms")

    image_files = sorted(args.aiib3_dir.glob("image*.pdb"),
                         key=lambda p: int(re.search(r'(\d+)', p.stem).group(1)))
    print(f"Found {len(image_files)} αIIbβ3 pathway images")

    summary = []
    for img_path in image_files:
        img_idx = int(re.search(r'(\d+)', img_path.stem).group(1))
        img_atoms = parse_pdb_atoms(img_path)

        all_transformed = []
        domain_rmsds = {}

        for domain_name in AVB3_DOMAINS:
            avb3_chain, avb3_start, avb3_end = AVB3_DOMAINS[domain_name]
            aiib3_chain, aiib3_start, aiib3_end = AIIB3_DOMAINS[domain_name]

            avb3_ca = get_ca_by_resrange(avb3_atoms, avb3_chain, avb3_start, avb3_end)
            aiib3_ca = get_ca_by_resrange(img_atoms, aiib3_chain, aiib3_start, aiib3_end)

            if len(avb3_ca) == 0 or len(aiib3_ca) == 0:
                print(f"  WARNING: empty domain {domain_name} for image{img_idx}")
                # Just keep original coords
                for a in avb3_atoms:
                    if (a["chain"] == avb3_chain and
                            avb3_start <= a["resseq"] <= avb3_end):
                        all_transformed.append(a)
                domain_rmsds[domain_name] = float("nan")
                continue

            # Truncate to shorter length for alignment
            n = min(len(avb3_ca), len(aiib3_ca))
            R, t_src, t_dst, rmsd = kabsch_rotation(avb3_ca[:n], aiib3_ca[:n])
            domain_rmsds[domain_name] = rmsd

            transformed = transform_atoms(avb3_atoms, avb3_chain,
                                          avb3_start, avb3_end,
                                          R, t_src, t_dst)
            all_transformed.extend(transformed)

        # Write output
        out_path = args.output_dir / f"avb3_domain_mapped_image{img_idx:02d}.pdb"
        with open(out_path, "w") as f:
            f.write(f"REMARK Domain-by-domain mapping from αIIbβ3 image{img_idx}\n")
            for dname, rmsd in domain_rmsds.items():
                f.write(f"REMARK   {dname}: alignment RMSD = {rmsd:.1f} Å\n")
            for a in sorted(all_transformed, key=lambda x: (x["chain"], x["resseq"])):
                f.write(a["line"] + "\n")
            f.write("END\n")

        rmsd_str = " ".join(f"{d}={domain_rmsds[d]:.1f}" for d in domain_rmsds)
        print(f"  image{img_idx:2d} → {out_path.name}  {rmsd_str}")

        summary.append({"image": img_idx, "output": out_path.name,
                        "domain_rmsds": {k: float(v) for k, v in domain_rmsds.items()}})

    summary_path = args.output_dir / "domain_mapping_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"\nMapped {len(summary)} structures (domain-by-domain) to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
