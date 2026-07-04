#!/usr/bin/env python3
"""Route-A Stage 2 (seed): construct an EXTENDED αVβ3 pose (state B) from bent.

No extended αVβ3 crystal exists (obj-041), so the string-method endpoint B must
be built. Integrin extension is dominated by a rigid-body rotation of the lower
legs about the genu (knee) hinge: the headpiece (β-propeller + βI/hybrid) stays
together while the legs (αV thigh/calf-1/calf-2, β3 I-EGF/β-tail) swing from
folded-back (bent) to pointing-away (extended, ~straight molecule).

Construction (pure numpy, CPU — no GPU):
  upper body = headpiece + thigh   : αV 1-592,   β3 1-440
  lower body = legs                : αV 593-956, β3 441-690
  hinge = centroid of the genu Cα (αV 588-596 + β3 436-444)
  rotate the lower body about the hinge so the leg-centroid direction flips from
  ~parallel (bent) to ~antiparallel (extended) w.r.t. the head-centroid direction.

This is a SEED endpoint — a rigid rotation leaves a strained hinge; a short GPU
minimize+equilibration (next step) repairs local geometry, and the string method
relaxes the full path A→B. Output: extended-pose PDB + CV report.
"""
import argparse
import json
import numpy as np

# domain boundaries (αV chain A / β3 chain B), approximate genu split
UPPER = {"A": (1, 592), "B": (1, 440)}
LOWER = {"A": (593, 956), "B": (441, 690)}
GENU = {"A": (588, 596), "B": (436, 444)}

CV_DEFS = {
    "cv0_head_calf": (("A", 1, 435), ("A", 442, 605)),
    "cv1_beta_head_tail": (("B", 1, 352), ("B", 353, 690)),
    "cv2_head_open": (("A", 1, 435), ("B", 1, 352)),
}


def parse_atoms(pdb_path):
    """Return list of dicts for protein ATOM records (full line preserved)."""
    atoms = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atoms.append({
                "chain": line[21], "resseq": int(line[22:26]),
                "name": line[12:16].strip(),
                "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                "line": line.rstrip("\n"),
            })
    return atoms


def in_range(a, table):
    r = table.get(a["chain"])
    return r is not None and r[0] <= a["resseq"] <= r[1]


def ca_centroid(atoms, chain, lo, hi):
    pts = [a["xyz"] for a in atoms if a["chain"] == chain and a["name"] == "CA"
           and lo <= a["resseq"] <= hi]
    return np.mean(pts, axis=0) if pts else None


def cvs(atoms):
    out = {}
    for name, (a, b) in CV_DEFS.items():
        ca, cb = ca_centroid(atoms, *a), ca_centroid(atoms, *b)
        out[name] = round(float(np.linalg.norm(ca - cb)), 2) if ca is not None and cb is not None else None
    return out


def rot_matrix(axis, theta):
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = axis
    return np.array([
        [c + x*x*(1-c),   x*y*(1-c) - z*s, x*z*(1-c) + y*s],
        [y*x*(1-c) + z*s, c + y*y*(1-c),   y*z*(1-c) - x*s],
        [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)],
    ])


def write_pdb(atoms, path):
    with open(path, "w") as fh:
        for a in atoms:
            l = a["line"]
            x, y, z = a["xyz"]
            fh.write(f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}\n")
        fh.write("END\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out-pdb", default="results/route_a/extended_seed.pdb")
    ap.add_argument("--out-json", default="results/route_a/stage2_extended_seed.json")
    args = ap.parse_args()

    atoms = parse_atoms(args.pdb)
    cvs_bent = cvs(atoms)

    # hinge = genu Cα centroid across both chains
    genu_pts = [a["xyz"] for a in atoms if a["name"] == "CA"
                and GENU.get(a["chain"]) and GENU[a["chain"]][0] <= a["resseq"] <= GENU[a["chain"]][1]]
    hinge = np.mean(genu_pts, axis=0)

    upper_cen = np.mean([a["xyz"] for a in atoms if a["name"] == "CA" and in_range(a, UPPER)], axis=0)
    lower_cen = np.mean([a["xyz"] for a in atoms if a["name"] == "CA" and in_range(a, LOWER)], axis=0)

    v_up = upper_cen - hinge
    v_low = lower_cen - hinge
    u_low = v_low / np.linalg.norm(v_low)
    target = -v_up / np.linalg.norm(v_up)          # legs point opposite the head
    axis = np.cross(u_low, target)
    angle = np.arccos(np.clip(np.dot(u_low, target), -1, 1))
    R = rot_matrix(axis, angle)

    # rotate all lower-body atoms about the hinge
    n_moved = 0
    for a in atoms:
        if in_range(a, LOWER):
            a["xyz"] = hinge + R @ (a["xyz"] - hinge)
            n_moved += 1

    cvs_ext = cvs(atoms)
    write_pdb(atoms, args.out_pdb)

    report = {
        "source_pdb": args.pdb,
        "method": "rigid-body lower-leg rotation about genu hinge (seed for string method)",
        "hinge_xyz": [round(float(v), 2) for v in hinge],
        "rotation_angle_deg": round(float(np.degrees(angle)), 1),
        "n_atoms_moved": n_moved,
        "cvs_bent_A": cvs_bent,
        "cvs_extended_seed_A": cvs_ext,
        "cv0_change_A": round(cvs_ext["cv0_head_calf"] - cvs_bent["cv0_head_calf"], 2),
        "out_pdb": args.out_pdb,
        "note": "SEED: rigid rotation leaves a strained hinge; needs GPU minimize+equil before use as state B.",
    }
    print(json.dumps(report, indent=2))
    with open(args.out_json, "w") as fh:
        json.dump(report, fh, indent=2)


if __name__ == "__main__":
    main()
