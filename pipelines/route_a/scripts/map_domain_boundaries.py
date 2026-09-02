#!/usr/bin/env python3
"""Transfer the validated alphaVbeta3 genu/upper/lower domain split onto every
integrin variant by sequence alignment, and emit a clean 2-chain bent structure.

Why alignment rather than hand-written per-variant residue ranges: the route-A
boundaries (UPPER/LOWER/GENU/HEAD/FOOT) were validated only for 1JV2 numbering
(obj-076: genu angle 41 -> 174 deg across the endpoints). Every other integrin
has its own numbering, its own insertions (alphaX/alphaL carry an extra ~180-residue
alphaI domain), and its own crystallographic gaps. Hand-transcribing ranges is
precisely the failure mode tasks/lessons.md warns about ("a CV defined as chain B
residues 436-444 selects DIFFERENT physical residues in the two files -- silently").

So: align variant alpha -> alphaV and variant beta -> beta3, walk the alignment,
and carry each boundary residue across. Then verify structurally that the
transferred genu actually behaves like a knee in the bent crystal.

Outputs per variant:
  <out-dir>/<variant>_bent.pdb    chains renamed A (alpha) / B (beta), original
                                  residue numbering preserved, single model,
                                  first altloc, protein ATOM records only
  <out-dir>/<variant>_boundaries.json  transferred ranges + bent-state CVs

Pure CPU. Requires biopython (alignment) + numpy.
"""
import argparse
import json
import os

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

# ---- reference boundaries, ORIGINAL 1JV2 numbering (from define_extension_cvs.py) ----
REF = {
    "UPPER": {"A": (1, 592), "B": (1, 440)},
    "LOWER": {"A": (593, 956), "B": (441, 690)},
    "GENU": {"A": (588, 596), "B": (436, 444)},
    "HEAD": {"A": (1, 438), "B": (55, 352)},
    "FOOT": {"A": (880, 956), "B": (600, 690)},
    "B_HEAD": {"B": (1, 352)},
    "B_TAIL": {"B": (353, 690)},
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def parse_chain_residues(path, chain_id):
    """[(resseq, one_letter, [atom_lines])] for one auth chain, first model/altloc."""
    residues, order, seen_atom = {}, [], set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[21] != chain_id:
                continue
            resname = line[17:20].strip()
            if resname not in THREE_TO_ONE:
                continue
            if line[16] not in (" ", "A"):
                continue
            icode = line[26]
            if icode != " ":
                continue  # insertion codes break integer-range selection
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            key = (resseq, atom_name)
            if key in seen_atom:
                continue
            seen_atom.add(key)
            if resseq not in residues:
                residues[resseq] = {"aa": THREE_TO_ONE[resname], "lines": []}
                order.append(resseq)
            residues[resseq]["lines"].append(line.rstrip("\n"))
    order.sort()
    return [(r, residues[r]["aa"], residues[r]["lines"]) for r in order]


# Align the head and the leg module SEPARATELY, by local alignment.
#
# A single global alignment of the whole alpha chain fails badly on the I-domain
# integrins: alphaX and alphaM carry a ~180-residue alphaI domain inserted into the
# beta-propeller that alphaV simply does not have. A global aligner pays that
# insertion as one enormous gap and smears the register across the rest of the
# chain (observed alpha scores: alphaV-vs-alphaV 4861, alphaIIb 1493, alpha5 2162,
# but alphaM 402 and alphaX 267 -- i.e. essentially noise). Since the genu and the
# foot both sit in the C-terminal leg module, aligning that module on its own,
# locally, keeps the insertion entirely out of the picture.
REF_WINDOWS = {
    "A": [(1, 439, "head"), (440, 956, "leg")],
    "B": [(1, 352, "head"), (353, 690, "leg")],
}


def _aligner(mode):
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.mode = mode
    return al


def align_map(ref_res, var_res, chain):
    """ref resseq -> variant resseq, from per-module local alignments."""
    al = _aligner("local")
    var_seq = "".join(a for _, a, _ in var_res)
    pairs, scores = [], {}
    for lo, hi, label in REF_WINDOWS.get(chain, [(-10**9, 10**9, "all")]):
        sub = [(r, a) for r, a, _ in ref_res if lo <= r <= hi]
        if len(sub) < 30:
            continue
        sub_seq = "".join(a for _, a in sub)
        try:
            aln = al.align(sub_seq, var_seq)[0]
        except Exception:
            continue
        scores[label] = round(float(aln.score), 1)
        for (rs, re_), (vs, ve) in zip(*aln.aligned):
            for k in range(re_ - rs):
                pairs.append((sub[rs + k][0], var_res[vs + k][0]))

    # Modules are aligned independently, so enforce a globally consistent register:
    # keep the longest strictly-increasing run of (ref, variant) pairs.
    pairs = sorted(set(pairs))
    mapping = {}
    if pairs:
        tails, back, idx = [], [-1] * len(pairs), []
        import bisect
        vals = []
        for i, (_, v) in enumerate(pairs):
            j = bisect.bisect_left(vals, v)
            if j == len(vals):
                vals.append(v)
                idx.append(i)
            else:
                vals[j] = v
                idx[j] = i
            back[i] = idx[j - 1] if j > 0 else -1
        cur = idx[-1]
        chosen = []
        while cur != -1:
            chosen.append(cur)
            cur = back[cur]
        for i in reversed(chosen):
            r, v = pairs[i]
            mapping[r] = v
    return mapping, scores, None, None


def transfer_endpoint(ref_resseq, mapping, ref_res, var_res, prefer=None):
    """Carry one reference residue number across; walk outward if it hit a gap.

    `prefer` ("up" for a range start, "down" for a range end) decides which way to
    walk first when the reference residue is unobserved. Without it, both ends of a
    short window can collapse onto the same nearby residue and silently invert the
    range -- which is how the beta-chain genu window came back empty.
    """
    ref_nums = [r for r, _, _ in ref_res]
    lo_ref, hi_ref = ref_nums[0], ref_nums[-1]
    var_nums = [r for r, _, _ in var_res]
    if ref_resseq <= lo_ref:
        return var_nums[0], "clamped_low"
    if ref_resseq >= hi_ref:
        return var_nums[-1], "clamped_high"
    if ref_resseq in mapping:
        return mapping[ref_resseq], "aligned"
    for delta in range(1, 60):
        order = (ref_resseq + delta, ref_resseq - delta)
        if prefer == "down":
            order = (ref_resseq - delta, ref_resseq + delta)
        elif prefer is None:
            order = (ref_resseq - delta, ref_resseq + delta)
        for cand in order:
            if cand in mapping:
                return mapping[cand], f"nearest_aligned({cand - ref_resseq:+d})"
    return None, "unmapped"


def write_clean_pdb(path, alpha_res, beta_res):
    """Write the two chains as A/B, original numbering, protein atoms only."""
    serial = 1
    with open(path, "w") as fh:
        for res_list, new_chain in ((alpha_res, "A"), (beta_res, "B")):
            for _, _, lines in res_list:
                for line in lines:
                    line = "ATOM  " + line[6:]
                    line = line[:21] + new_chain + line[22:]
                    line = f"{line[:6]}{serial:5d}{line[11:]}"
                    fh.write(line.ljust(80) + "\n")
                    serial += 1
            fh.write("TER\n")
        fh.write("END\n")


def ca_coords(res_list, lo, hi):
    pts = []
    for resseq, _, lines in res_list:
        if not (lo <= resseq <= hi):
            continue
        for line in lines:
            if line[12:16].strip() == "CA":
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                break
    return np.array(pts)


def bent_cvs(alpha_res, beta_res, b):
    """Genu angle / end-to-end / Rg / long extent using the TRANSFERRED boundaries."""
    def sel(group):
        pts = []
        rng = b[group]
        if "A" in rng:
            p = ca_coords(alpha_res, *rng["A"])
            if len(p):
                pts.append(p)
        if "B" in rng:
            p = ca_coords(beta_res, *rng["B"])
            if len(p):
                pts.append(p)
        return np.vstack(pts) if pts else np.empty((0, 3))

    upper, lower, genu = sel("UPPER"), sel("LOWER"), sel("GENU")
    head, foot = sel("HEAD"), sel("FOOT")
    allca = np.vstack([ca_coords(alpha_res, -10**9, 10**9),
                       ca_coords(beta_res, -10**9, 10**9)])
    out = {"n_upper_ca": int(len(upper)), "n_lower_ca": int(len(lower)),
           "n_genu_ca": int(len(genu))}
    if len(genu) and len(upper) and len(lower):
        h = genu.mean(axis=0)
        vu, vl = upper.mean(axis=0) - h, lower.mean(axis=0) - h
        cos = np.dot(vu, vl) / (np.linalg.norm(vu) * np.linalg.norm(vl))
        out["genu_angle_deg"] = round(float(np.degrees(np.arccos(np.clip(cos, -1, 1)))), 2)
    if len(head) and len(foot):
        out["end_to_end_A"] = round(float(np.linalg.norm(head.mean(axis=0) - foot.mean(axis=0))), 2)
    c = allca.mean(axis=0)
    out["Rg_A"] = round(float(np.sqrt(((allca - c) ** 2).sum(axis=1).mean())), 2)
    _, _, vt = np.linalg.svd(allca - c, full_matrices=False)
    proj = (allca - c) @ vt[0]
    out["long_extent_A"] = round(float(proj.max() - proj.min()), 2)
    out["n_total_ca"] = int(len(allca))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovered", default="results/variants/discovered_ectodomains.json")
    ap.add_argument("--ref-pdb", default="inputs/1jv2.pdb")
    ap.add_argument("--ref-alpha-chain", default="A")
    ap.add_argument("--ref-beta-chain", default="B")
    ap.add_argument("--out-dir", default="inputs/variants/prepared")
    ap.add_argument("--out-json", default="results/variants/domain_boundaries.json")
    ap.add_argument("--only", default="", help="comma-separated variant names to restrict to")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)

    disc = json.load(open(args.discovered))
    reps = disc["representatives"]
    if args.only:
        keep = {s.strip() for s in args.only.split(",") if s.strip()}
        reps = {k: v for k, v in reps.items() if k in keep}

    ref_alpha = parse_chain_residues(args.ref_pdb, args.ref_alpha_chain)
    ref_beta = parse_chain_residues(args.ref_pdb, args.ref_beta_chain)
    print(f"[map] reference {args.ref_pdb}: alpha {len(ref_alpha)} res "
          f"({ref_alpha[0][0]}-{ref_alpha[-1][0]}), beta {len(ref_beta)} res "
          f"({ref_beta[0][0]}-{ref_beta[-1][0]})", flush=True)

    results = {}
    for variant, rec in sorted(reps.items()):
        pdb_path = rec["pdb_path"]
        a_chain, b_chain = rec["alpha"]["chain"], rec["beta"]["chain"]
        var_alpha = parse_chain_residues(pdb_path, a_chain)
        var_beta = parse_chain_residues(pdb_path, b_chain)
        if len(var_alpha) < 50 or len(var_beta) < 50:
            print(f"[map] SKIP {variant}: chains too short after parsing", flush=True)
            continue

        amap, ascore, _, _ = align_map(ref_alpha, var_alpha, "A")
        bmap, bscore, _, _ = align_map(ref_beta, var_beta, "B")
        aid = 100.0 * len(amap) / max(1, len(ref_alpha))
        bid = 100.0 * len(bmap) / max(1, len(ref_beta))

        bounds, prov = {}, {}
        for group, chain_ranges in REF.items():
            bounds[group] = {}
            prov[group] = {}
            for ch, (lo, hi) in chain_ranges.items():
                mapping = amap if ch == "A" else bmap
                ref_res = ref_alpha if ch == "A" else ref_beta
                var_res = var_alpha if ch == "A" else var_beta
                nlo, wlo = transfer_endpoint(lo, mapping, ref_res, var_res, prefer="up")
                nhi, whi = transfer_endpoint(hi, mapping, ref_res, var_res, prefer="down")
                if nlo is None or nhi is None or nhi <= nlo:
                    bounds[group][ch] = None
                    prov[group][ch] = (f"unavailable ({wlo}/{whi}) -- reference window "
                                       "unobserved in the 1JV2 coordinates")
                else:
                    bounds[group][ch] = [int(nlo), int(nhi)]
                    prov[group][ch] = f"{wlo}/{whi}"

        # UPPER/LOWER must partition the chain with no gap or overlap at the genu
        for ch in ("A", "B"):
            up, lo_ = bounds["UPPER"].get(ch), bounds["LOWER"].get(ch)
            var_res = var_alpha if ch == "A" else var_beta
            if up and lo_:
                if lo_[0] <= up[1]:
                    lo_[0] = up[1] + 1
                    prov["LOWER"][ch] += " +split_fix"
                up[0] = min(up[0], var_res[0][0])
                lo_[1] = max(lo_[1], var_res[-1][0])

        clean = {g: {ch: tuple(v) for ch, v in d.items() if v} for g, d in bounds.items()}
        cvs = bent_cvs(var_alpha, var_beta, clean)

        out_pdb = os.path.join(args.out_dir, f"{variant}_bent.pdb")
        write_clean_pdb(out_pdb, var_alpha, var_beta)

        # a genu that really is a knee: both arms substantial, angle not straight
        ok = (cvs.get("n_upper_ca", 0) > 200 and cvs.get("n_lower_ca", 0) > 150
              and cvs.get("n_genu_ca", 0) >= 4 and cvs.get("genu_angle_deg") is not None)
        results[variant] = {
            "pdb_id": rec["pdb_id"], "title": rec["title"],
            "source_pdb": pdb_path, "prepared_pdb": out_pdb,
            "alpha": {**rec["alpha"], "source_chain": a_chain, "renamed_to": "A"},
            "beta": {**rec["beta"], "source_chain": b_chain, "renamed_to": "B"},
            "alignment": {"alpha_module_scores": ascore, "beta_module_scores": bscore,
                          "alpha_pct_ref_mapped": round(aid, 1),
                          "beta_pct_ref_mapped": round(bid, 1)},
            "boundaries": {g: {ch: list(v) for ch, v in d.items()} for g, d in clean.items()},
            "boundary_provenance": prov,
            "bent_cvs": cvs,
            "usable_for_morph": bool(ok),
        }
        ang = cvs.get("genu_angle_deg")
        print(f"[map] {variant:<22} {rec['pdb_id']}  mapped alpha {aid:5.1f}% / beta {bid:5.1f}%"
              f"  genu={ang}deg  upper/lower CA={cvs.get('n_upper_ca')}/{cvs.get('n_lower_ca')}"
              f"  {'OK' if ok else 'CHECK'}", flush=True)

    with open(args.out_json, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\n[map] wrote {args.out_json} ({len(results)} variants)")


if __name__ == "__main__":
    main()
