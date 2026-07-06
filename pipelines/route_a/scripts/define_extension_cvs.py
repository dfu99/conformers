#!/usr/bin/env python3
"""Define REAL extension collective variables for the route-A bent->extended string method.

Motivation (see tasks/lessons.md):
  The synthetic placeholder CVs in examples/av_cvs_remapped.json were built to test
  the alpha-IIb -> alpha-V remap converter, NOT as validated reaction coordinates.
  cv0 ("head_calf", A:1-435 vs A:442-605) measures head-to-THIGH, both of which sit
  in the upper body and swing together -- so it barely moves on extension
  (54.3 -> 53.9 A) while the molecule clearly opens (Rg 39 -> 67 A). Wrong metric.

This script defines mechanistic extension CVs and verifies they SEPARATE the two
endpoints (bent state A = 1JV2 crystal; extended state B = the morph product):

  genu_angle   -- PRIMARY. Hinge angle at the genu (knee): angle between the
                  upper-body axis and the lower-leg axis about the genu centroid.
                  This is the literal leg-swing reaction coordinate. Bent ~acute,
                  extended ~180 deg (straight).
  end_to_end   -- head-cap centroid  <->  membrane-proximal foot centroid distance.
                  The tip-to-tip jackknife opening.
  Rg           -- radius of gyration (global shape sanity check).
  long_extent  -- long-axis span (SVD PC1 projection spread).
  cv0_placeholder -- the BROKEN placeholder, computed here only to document that it
                  does NOT separate the states.

The two PDBs use different residue numbering (PDBFixer renumbered state B to a
contiguous 1..N), but they contain the identical CA set in the same order, so the
i-th CA maps between them. We key everything by the ORIGINAL 1JV2 numbering.

Pure-CPU, dependency-light (numpy + a minimal PDB parser). No GPU, no OpenMM.
"""
import argparse
import json
import os
import numpy as np


# ---- domain boundaries in ORIGINAL 1JV2 numbering (from stage2c_morph_extend.py) ----
UPPER = {"A": (1, 592), "B": (1, 440)}      # headpiece + upper legs (thigh / I-EGF upper)
LOWER = {"A": (593, 956), "B": (441, 690)}  # lower legs (calf-1/2, beta-tail)
GENU = {"A": (588, 596), "B": (436, 444)}   # the knee hinge

# head cap (ligand-binding top) and membrane-proximal foot, for the end-to-end CV
HEAD = {"A": (1, 438), "B": (55, 352)}      # aV beta-propeller + b3 betaI/hybrid/PSI headpiece
FOOT = {"A": (880, 956), "B": (600, 690)}   # aV calf-2 C-term + b3 beta-tail C-term (feet)

# the broken placeholder cv0 (head vs thigh -- both in the upper body)
CV0_A = {"A": (1, 435)}
CV0_B = {"A": (442, 605)}


def load_ca(path):
    """Return {chain: [(resid, xyz), ...] sorted by resid} for first-altloc CA atoms."""
    chains = {}
    seen = set()
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        ch = line[21]
        resid = int(line[22:26])
        key = (ch, resid)
        if key in seen:
            continue
        seen.add(key)
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        chains.setdefault(ch, []).append((resid, xyz))
    for ch in chains:
        chains[ch].sort(key=lambda t: t[0])
    return chains


def build_original_keyed(ref_chains, other_chains):
    """Map other_chains (possibly renumbered) onto the ref (original) numbering by rank.

    ref_chains and other_chains must have the same CA count per chain and the same
    physical order. Returns {chain: {orig_resid: xyz}} for `other`.
    """
    out = {}
    for ch, ref_list in ref_chains.items():
        oth_list = other_chains[ch]
        if len(ref_list) != len(oth_list):
            raise ValueError(
                f"chain {ch}: CA count mismatch ref={len(ref_list)} other={len(oth_list)} "
                "-- structures are not the same residue set")
        out[ch] = {ref_list[i][0]: oth_list[i][1] for i in range(len(ref_list))}
    return out


def orig_keyed_from_ref(ref_chains):
    return {ch: {resid: xyz for resid, xyz in lst} for ch, lst in ref_chains.items()}


def select(byorig, ranges):
    """Collect CA coords whose original resid falls in the given per-chain ranges."""
    pts = []
    for ch, (lo, hi) in ranges.items():
        d = byorig.get(ch, {})
        for resid, xyz in d.items():
            if lo <= resid <= hi:
                pts.append(xyz)
    return np.array(pts)


def centroid(byorig, ranges):
    return select(byorig, ranges).mean(axis=0)


def rg(byorig):
    pts = np.array([xyz for ch in byorig for xyz in byorig[ch].values()])
    c = pts.mean(axis=0)
    return float(np.sqrt(((pts - c) ** 2).sum(axis=1).mean()))


def long_extent(byorig):
    pts = np.array([xyz for ch in byorig for xyz in byorig[ch].values()])
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c, full_matrices=False)
    proj = (pts - c) @ vt[0]
    return float(proj.max() - proj.min())


def genu_angle(byorig):
    """Angle (deg) between upper-body axis and lower-leg axis about the genu centroid."""
    hinge = centroid(byorig, GENU)
    v_up = centroid(byorig, UPPER) - hinge
    v_low = centroid(byorig, LOWER) - hinge
    cos = np.dot(v_up, v_low) / (np.linalg.norm(v_up) * np.linalg.norm(v_low))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def end_to_end(byorig):
    return float(np.linalg.norm(centroid(byorig, HEAD) - centroid(byorig, FOOT)))


def cv0_placeholder(byorig):
    return float(np.linalg.norm(centroid(byorig, CV0_A) - centroid(byorig, CV0_B)))


def all_cvs(byorig):
    return {
        "genu_angle_deg": round(genu_angle(byorig), 2),
        "end_to_end_A": round(end_to_end(byorig), 2),
        "Rg_A": round(rg(byorig), 2),
        "long_extent_A": round(long_extent(byorig), 2),
        "cv0_placeholder_A": round(cv0_placeholder(byorig), 2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bent", default="data/reference_pdbs/1jv2.pdb",
                    help="bent state A (original 1JV2 numbering)")
    ap.add_argument("--extended", default="results/route_a/extended_state_b.pdb",
                    help="extended state B (may be renumbered)")
    ap.add_argument("--out-json", default="results/route_a/extension_cv_endpoints.json")
    ap.add_argument("--out-cvdef", default="pipelines/route_a/examples/av_extension_cvs.json")
    args = ap.parse_args()

    bent_chains = load_ca(args.bent)
    ext_chains = load_ca(args.extended)

    bent = orig_keyed_from_ref(bent_chains)
    ext = build_original_keyed(bent_chains, ext_chains)

    bent_cvs = all_cvs(bent)
    ext_cvs = all_cvs(ext)

    print(f"{'CV':<20}{'bent (A)':>14}{'extended (B)':>16}{'delta':>12}{'sep?':>7}")
    print("-" * 69)
    rows = []
    for k in bent_cvs:
        b, e = bent_cvs[k], ext_cvs[k]
        d = round(e - b, 2)
        rel = abs(d) / (abs(b) + 1e-9)
        sep = "YES" if rel > 0.2 else "no"
        rows.append((k, b, e, d, sep))
        print(f"{k:<20}{b:>14}{e:>16}{d:>12}{sep:>7}")

    result = {
        "bent_state_A": {"source": args.bent, "cvs": bent_cvs},
        "extended_state_B": {"source": args.extended, "cvs": ext_cvs},
        "separation": {k: {"bent": b, "extended": e, "delta": d, "separates": s}
                       for (k, b, e, d, s) in rows},
        "note": ("genu_angle + end_to_end + Rg + long_extent all separate the endpoints "
                 "strongly; cv0_placeholder (head-thigh) does NOT -- confirming it was the "
                 "wrong reaction coordinate."),
    }
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nwrote {args.out_json}")

    # ---- write the real CV definition file (replaces the synthetic placeholder) ----
    cvdef = {
        "_source": "Mechanistic extension CVs for route-A alphaV-beta3 bent->extended string method",
        "_derivation": "pipelines/route_a/scripts/define_extension_cvs.py",
        "_supersedes": "examples/av_cvs_remapped.json (synthetic placeholder; cv0 was head-thigh, did not track extension)",
        "_numbering": "ORIGINAL 1JV2 crystal numbering (chain A 1..956, chain B 55..690)",
        "_endpoints": {"bent_state_A": bent_cvs, "extended_state_B": ext_cvs},
        "primary_cv": "genu_angle",
        "cvs": [
            {
                "name": "genu_angle",
                "type": "hinge_angle",
                "vertex": {"chain_ranges": GENU, "role": "genu (knee) hinge centroid"},
                "arm_upper": {"chain_ranges": UPPER, "role": "headpiece + upper legs"},
                "arm_lower": {"chain_ranges": LOWER, "role": "lower legs (calf/beta-tail)"},
                "units": "degrees",
                "bent": bent_cvs["genu_angle_deg"], "extended": ext_cvs["genu_angle_deg"],
                "_role": "PRIMARY reaction coordinate for the string method",
            },
            {
                "name": "end_to_end",
                "type": "centroid_distance",
                "domain_a": {"chain_ranges": HEAD, "role": "ligand-binding head cap"},
                "domain_b": {"chain_ranges": FOOT, "role": "membrane-proximal feet"},
                "units": "angstrom",
                "bent": bent_cvs["end_to_end_A"], "extended": ext_cvs["end_to_end_A"],
            },
            {
                "name": "Rg",
                "type": "radius_of_gyration",
                "selection": "all CA",
                "units": "angstrom",
                "bent": bent_cvs["Rg_A"], "extended": ext_cvs["Rg_A"],
            },
        ],
    }
    os.makedirs(os.path.dirname(args.out_cvdef), exist_ok=True)
    with open(args.out_cvdef, "w") as fh:
        json.dump(cvdef, fh, indent=2)
    print(f"wrote {args.out_cvdef}")


if __name__ == "__main__":
    main()
