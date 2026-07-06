#!/usr/bin/env python3
"""Stage 3 (initialization): build the initial STRING (path of images) for the
route-A bent<->extended transition from the incremental leg-swing morph frames.

The string method optimizes a discretized path (a "string" of images) connecting
two metastable states. It needs a good initial path; ours is the 16-frame morph
(obj-075) plus the bent crystal as image 0. This script:

  1. Loads image 0 = bent state A (1JV2 crystal) and images 1..16 = the morph frames
     (all share the identical Cα set; keyed to original 1JV2 numbering by rank).
  2. Computes the validated extension CVs per image (genu hinge angle = primary RC,
     plus head<->foot end-to-end, Rg, long-axis extent).  [see obj-076]
  3. Measures adjacent-image Cα-RMSD (Kabsch-aligned) -> cumulative arc length, the
     natural string metric. This reveals where images bunch up (uneven spacing).
  4. Reparametrizes: places N nodes at EQUAL arc length and reports the interpolated
     CV target at each -- the standard string reparametrization step, here applied
     once to the initial path so the string method starts from evenly spaced images.

Pure-CPU, numpy only. No GPU, no OpenMM, no pod contention. The per-image MD
evolution (the other half of a string iteration) is a separate GPU step, wired in
later; this establishes and validates the initial path in CV space.
"""
import argparse
import glob
import json
import os
import numpy as np

import define_extension_cvs as cv  # reuse the validated CV computation


def load_image_byorig(path, ref_chains):
    """Load a PDB's CA, rank-match onto the bent reference's original numbering."""
    chains = cv.load_ca(path)
    return cv.build_original_keyed(ref_chains, chains)


def ca_matrix(byorig, order):
    """Stack CA coords in a fixed (chain, resid) order -> (N,3) array."""
    return np.array([byorig[ch][r] for ch, r in order])


def kabsch_rmsd(P, Q):
    """RMSD between P and Q (N,3) after optimal rigid superposition."""
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    Pr = Pc @ R.T
    return float(np.sqrt(((Pr - Qc) ** 2).sum(1).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bent", default="data/reference_pdbs/1jv2.pdb")
    ap.add_argument("--frames-dir", required=True,
                    help="dir with frame_01.pdb .. frame_16.pdb (morph output)")
    ap.add_argument("--n-nodes", type=int, default=12,
                    help="number of equally-spaced images in the reparametrized string")
    ap.add_argument("--out-json", default="results/route_a/initial_string.json")
    args = ap.parse_args()

    bent_chains = cv.load_ca(args.bent)
    # canonical (chain, resid) order from the bent reference
    order = [(ch, r) for ch in sorted(bent_chains) for r, _ in bent_chains[ch]]

    images = [("bent_stateA", cv.orig_keyed_from_ref(bent_chains))]
    for f in sorted(glob.glob(os.path.join(args.frames_dir, "frame_*.pdb"))):
        images.append((os.path.basename(f).replace(".pdb", ""),
                       load_image_byorig(f, bent_chains)))
    n = len(images)
    print(f"loaded {n} images (bent + {n-1} morph frames)")

    # ---- per-image CVs ----
    per_image = []
    mats = []
    for name, byorig in images:
        c = cv.all_cvs(byorig)
        per_image.append({"name": name, **c})
        mats.append(ca_matrix(byorig, order))

    # ---- adjacent Cα-RMSD -> cumulative arc length ----
    seg = [0.0]
    for i in range(1, n):
        seg.append(kabsch_rmsd(mats[i - 1], mats[i]))
    arc = np.cumsum(seg)
    total = float(arc[-1])
    for i, p in enumerate(per_image):
        p["seg_rmsd_A"] = round(seg[i], 3)
        p["arc_A"] = round(float(arc[i]), 3)
        p["arc_frac"] = round(float(arc[i] / total), 4)

    # ---- reparametrize: N nodes at equal arc length, interpolate CVs ----
    genu = np.array([p["genu_angle_deg"] for p in per_image])
    e2e = np.array([p["end_to_end_A"] for p in per_image])
    rg = np.array([p["Rg_A"] for p in per_image])
    targets = np.linspace(0.0, total, args.n_nodes)
    repar = []
    for t in targets:
        repar.append({
            "arc_A": round(float(t), 3),
            "arc_frac": round(float(t / total), 4),
            "genu_angle_deg": round(float(np.interp(t, arc, genu)), 2),
            "end_to_end_A": round(float(np.interp(t, arc, e2e)), 2),
            "Rg_A": round(float(np.interp(t, arc, rg)), 2),
        })

    # ---- spacing diagnostics ----
    seg_arr = np.array(seg[1:])
    cv_gaps = np.abs(np.diff(genu))
    result = {
        "stage": "3_initial_string",
        "n_images": n,
        "primary_cv": "genu_angle",
        "total_path_rmsd_A": round(total, 2),
        "genu_span_deg": [round(float(genu.min()), 2), round(float(genu.max()), 2)],
        "monotonic_genu": bool(np.all(np.diff(genu) > 0)),
        "spacing_rmsd_A": {
            "mean": round(float(seg_arr.mean()), 3),
            "min": round(float(seg_arr.min()), 3),
            "max": round(float(seg_arr.max()), 3),
            "cv_of_spacing": round(float(seg_arr.std() / seg_arr.mean()), 3),
        },
        "genu_gap_deg": {
            "mean": round(float(cv_gaps.mean()), 2),
            "min": round(float(cv_gaps.min()), 2),
            "max": round(float(cv_gaps.max()), 2),
        },
        "per_image": per_image,
        "reparametrized_nodes": repar,
    }
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"genu angle span: {genu.min():.1f} -> {genu.max():.1f} deg  "
          f"(monotonic: {result['monotonic_genu']})")
    print(f"path length (Cα-RMSD): {total:.1f} Å over {n-1} segments")
    print(f"segment RMSD: mean {seg_arr.mean():.2f} Å, "
          f"min {seg_arr.min():.2f}, max {seg_arr.max():.2f} "
          f"(CV {seg_arr.std()/seg_arr.mean():.2f})")
    print(f"reparametrized to {args.n_nodes} equally-spaced nodes")
    print(f"wrote {args.out_json}")


if __name__ == "__main__":
    main()
