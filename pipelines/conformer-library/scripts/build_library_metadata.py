#!/usr/bin/env python3
"""Build library.json for a conformer library directory.

Walks one or more frame directories, computes per-frame CVs from the
heterodimer domain definitions, classifies the state, and writes a
unified library.json that pipelines/sim-afm-video consumes.

Usage:
    python build_library_metadata.py \\
        --frame-dirs run_bent/frames run_extend/frames \\
        --domains avb3 \\
        --output data/runs/avb3/conformers/all_frames_bent_extended
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import mdtraj as md


AVB3_DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 956),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}

DOMAIN_REGISTRY = {
    "avb3": AVB3_DOMAINS,
    # add aiib3, a5b1, etc. when those become available
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frame-dirs", nargs="+", type=Path, required=True,
                   help="One or more directories containing frame_*.pdb")
    p.add_argument("--domains", choices=list(DOMAIN_REGISTRY.keys()),
                   default="avb3")
    p.add_argument("--output", type=Path, required=True,
                   help="Library directory to write (PDBs copied + library.json)")
    p.add_argument("--cv-pairs", nargs="+", default=[
                   "alpha_head_thigh:alpha_calf",
                   "beta_head:beta_tail",
                   "alpha_head_thigh:beta_head"],
                   help="Domain pairs (colon-separated) for CV0, CV1, CV2, ...")
    return p.parse_args()


def domain_centroid(traj, chain_id: str, lo: int, hi: int) -> np.ndarray:
    sel = []
    for atom in traj.topology.atoms:
        if atom.name != "CA":
            continue
        cid = atom.residue.chain.chain_id or chr(ord("A") + atom.residue.chain.index)
        if cid == chain_id and lo <= atom.residue.resSeq <= hi:
            sel.append(atom.index)
    if not sel:
        raise ValueError(f"No CAs in chain {chain_id} {lo}-{hi}")
    xyz = traj.xyz[0, sel]
    return xyz.mean(axis=0) * 10.0  # nm → Å


def classify_state(cvs: dict) -> str:
    cv0 = cvs.get("cv0_A", -1)
    cv1 = cvs.get("cv1_A", -1)
    cv2 = cvs.get("cv2_A", 0)
    if cv0 < 60 and cv1 < 50:
        return "BC"
    if cv0 > 65 and cv1 > 55:
        if cv2 > 50:
            return "EO"
        return "EC"
    return "Intermediate"


def main() -> int:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    domains = DOMAIN_REGISTRY[args.domains]

    pairs = []
    for spec in args.cv_pairs:
        a, b = spec.split(":")
        if a not in domains or b not in domains:
            raise ValueError(f"Domain {a} or {b} not in registry {args.domains}")
        pairs.append((a, b))
    print(f"CV pairs: {pairs}")

    all_pdbs = []
    for d in args.frame_dirs:
        files = sorted(d.glob("*.pdb"))
        for f in files:
            all_pdbs.append((f, d.name))

    if not all_pdbs:
        raise FileNotFoundError("No PDBs found in --frame-dirs")
    print(f"Found {len(all_pdbs)} PDBs")

    conformers = []
    for k, (src, source_dir) in enumerate(all_pdbs):
        new_name = f"frame_{k:04d}.pdb"
        dst = args.output / new_name
        shutil.copy(str(src), str(dst))
        traj = md.load(str(dst))
        cv_record = {"file": new_name, "source_dir": source_dir}
        for ci, (a, b) in enumerate(pairs):
            ca = domain_centroid(traj, *domains[a])
            cb = domain_centroid(traj, *domains[b])
            cv = float(np.linalg.norm(ca - cb))
            cv_record[f"cv{ci}_A"] = round(cv, 2)
        cv_record["state"] = classify_state(cv_record)
        conformers.append(cv_record)
        if (k + 1) % 50 == 0:
            print(f"  scored {k+1}/{len(all_pdbs)}")

    # Trajectory order: just sequential within each source_dir, then
    # concatenate. Adjacency: pairs of consecutive entries.
    trajectory_order = list(range(len(conformers)))
    adjacency = [[i, i + 1] for i in range(len(conformers) - 1)]

    library = {
        "topology_pdb": conformers[0]["file"],
        "conformers": conformers,
        "trajectory_order": trajectory_order,
        "adjacency": adjacency,
        "domain_set": args.domains,
        "cv_pairs": [{"name": f"cv{i}_A", "domain_a": a, "domain_b": b}
                     for i, (a, b) in enumerate(pairs)],
        "states": {
            "BC": "bent closed (CV0 < 60 Å, CV1 < 50 Å)",
            "EC": "extended closed (CV0 > 65 Å, CV1 > 55 Å, CV2 < 40 Å)",
            "EO": "extended open (CV0 > 65 Å, CV1 > 55 Å, CV2 > 50 Å)",
            "Intermediate": "anything else",
        },
    }
    out_path = args.output / "library.json"
    with out_path.open("w") as f:
        json.dump(library, f, indent=2)
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
