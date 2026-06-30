#!/usr/bin/env python3
"""Route-A Stage 1 (part a): bent-state CV baseline + obj-029 reproduction check.

Parses the 1JV2 αVβ3 bent-ectodomain crystal and computes:
  (1) obj-029's head-tail combined centroid distance (first/last third of CA
      atoms per chain) — should reproduce ~42.6 Å, validating the structure
      loaded correctly before any MD.
  (2) the three route-A steering CVs in αV numbering (per
      pipelines/route_a/examples/av_cvs_remapped.json): cv0_head_calf,
      cv1_beta_head_tail, cv2_head_open. These are the bent endpoint (state A)
      the string method steers FROM.

Pure stdlib (no numpy / OpenMM) so it runs before the conda env is built.
Protein CA only (ATOM records) — excludes calcium HETATM ions, whose atom
name is also "CA".
"""
import argparse
import json
import math
import os
from collections import OrderedDict


def parse_ca(pdb_path):
    """chain -> list of (resseq, (x,y,z)) for protein CA atoms, in file order."""
    chains = OrderedDict()
    seen = set()
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):  # ATOM only -> excludes Ca2+ HETATM
                continue
            if line[12:16].strip() != "CA":
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            chain = line[21]
            try:
                rs = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            key = (chain, rs, icode)
            if key in seen:
                continue
            seen.add(key)
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            chains.setdefault(chain, []).append((rs, xyz))
    return chains


def centroid(points):
    n = len(points)
    if n == 0:
        return None
    return tuple(sum(p[i] for p in points) / n for i in range(3))


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def sel_range(ca_list, lo, hi):
    return [xyz for (rs, xyz) in ca_list if lo <= rs <= hi]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", default="/workspace/route_a/stage1/cv_baseline.json")
    args = ap.parse_args()

    chains = parse_ca(args.pdb)
    A_ca = chains.get("A", [])
    B_ca = chains.get("B", [])
    A = [xyz for (_, xyz) in A_ca]
    B = [xyz for (_, xyz) in B_ca]

    res = OrderedDict()
    res["pdb"] = args.pdb
    res["n_ca"] = {"A": len(A), "B": len(B)}

    # (1) obj-029 head-tail combined (first/last third of CA per chain, file order)
    def third(lst):
        n = len(lst)
        return lst[: n // 3], lst[-(n // 3):]

    ha, ta = third(A)
    hb, tb = third(B)
    head_com = centroid(ha + hb)
    tail_com = centroid(ta + tb)
    res["obj029_head_tail_combined_A"] = round(dist(head_com, tail_com), 2)
    res["alpha_head_tail_A"] = round(dist(centroid(ha), centroid(ta)), 2)
    res["beta_head_tail_A"] = round(dist(centroid(hb), centroid(tb)), 2)
    res["obj029_reference_1jv2_head_tail_A"] = 42.6

    # (2) route-A steering CVs (αV numbering, by PDB resseq)
    cv0 = dist(centroid(sel_range(A_ca, 1, 435)), centroid(sel_range(A_ca, 442, 605)))
    cv1 = dist(centroid(sel_range(B_ca, 1, 352)), centroid(sel_range(B_ca, 353, 690)))
    cv2 = dist(centroid(sel_range(A_ca, 1, 435)), centroid(sel_range(B_ca, 1, 352)))
    res["cvs_bent_state_A"] = {
        "cv0_head_calf": round(cv0, 2),
        "cv1_beta_head_tail": round(cv1, 2),
        "cv2_head_open": round(cv2, 2),
    }
    res["resseq_span"] = {
        "A": [A_ca[0][0], A_ca[-1][0]] if A_ca else None,
        "B": [B_ca[0][0], B_ca[-1][0]] if B_ca else None,
    }

    print(json.dumps(res, indent=2))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=2)


if __name__ == "__main__":
    main()
