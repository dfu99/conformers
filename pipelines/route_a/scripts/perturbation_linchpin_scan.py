#!/usr/bin/env python3
"""Perturbation / linchpin scan: which residues govern the αVβ3 bent<->extended state?

The PI's question: which residues, if mutated, would destabilize the EXTENDED state so
it snaps back to bent? Those are the residues that LOCK the extended (open) conformation
-- i.e. residues that make favorable contacts in extended that they do NOT make in bent.

Method (structure-based, CPU-only; the definitive test is mutate->MD->snap-back on GPU):
  1. Heavy-atom residue-residue contacts (min heavy-atom distance < 4.5 A) in BOTH the
     bent crystal (1JV2, state A) and the extended model (state B), keyed by residue RANK
     so the two differently-numbered PDBs compare 1:1.
  2. EXTENDED-UNIQUE contacts (present extended, absent bent) = the locks holding the legs
     straight.  BENT-UNIQUE contacts = the head-leg clasp that re-forms on snap-back.
  3. Per-residue "extension-lock score": count extended-unique contacts, weighted by
       x3 salt bridge (Asp/Glu carboxylate <-> Lys/Arg/His cation < 4.0 A)
       x2 H-bond-capable (N/O <-> N/O < 3.5 A)
       x2 if the contact spans the genu (one partner above the knee, one below) or is
          inter-chain -- these mechanically couple the straightened legs
       +genu-proximity bonus (within +/-8 residues of the knee hinge).
  4. Rank residues; the top ones are the candidate mutations predicted to let the
     extended state relax back toward bent. Also report the bent-clasp residues (mutating
     those would instead STABILIZE extended -- the opposite lever).

Reports residues in ORIGINAL 1JV2 numbering + residue name.
"""
import argparse
import json
import numpy as np
from scipy.spatial import cKDTree

# knee split (original 1JV2 numbering, from the morph domain defs)
UPPER = {"A": (1, 592), "B": (1, 440)}
GENU_CENTER = {"A": 592, "B": 440}
CATION = {"LYS": ["NZ"], "ARG": ["NH1", "NH2", "NE"], "HIS": ["ND1", "NE2"]}
ANION = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}

# αV / β3 domain labels for reporting (original 1JV2 numbering)
def domain(chain, resid):
    if chain == "A":
        if resid <= 438: return "αV β-propeller (head)"
        if resid <= 592: return "αV thigh"
        if resid <= 746: return "αV calf-1"
        return "αV calf-2 (foot)"
    else:
        if resid <= 108: return "β3 PSI/hybrid (head)"
        if resid <= 352: return "β3 βI (head)"
        if resid <= 440: return "β3 I-EGF-1/2 (knee)"
        return "β3 I-EGF-3/4 + β-tail (leg)"


def load_residues(path):
    """Return per-chain lists of residues (sorted by resid), each a dict with atoms."""
    chains = {}
    cur = {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        nm = line[12:16].strip()
        el = line[76:78].strip() or nm[0]
        if el == "H":
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        ch = line[21]
        resid = int(line[22:26])
        rn = line[17:20].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        key = (ch, resid)
        r = cur.get(key)
        if r is None:
            r = {"chain": ch, "resid": resid, "resname": rn, "atoms": []}
            cur[key] = r
            chains.setdefault(ch, []).append(r)
        r["atoms"].append((nm, el, xyz))
    for ch in chains:
        chains[ch].sort(key=lambda r: r["resid"])
    return chains


def contact_set(residues):
    """residues: flat list. Return dict {(i,j): min_dist} for heavy-atom contacts <4.5A,
    plus per-pair charged-atom min distance and N/O min distance."""
    coords, owner = [], []
    charged, polar = [], []  # per atom flags
    for idx, r in enumerate(residues):
        for nm, el, xyz in r["atoms"]:
            coords.append(xyz); owner.append(idx)
            is_cat = r["resname"] in CATION and nm in CATION[r["resname"]]
            is_an = r["resname"] in ANION and nm in ANION[r["resname"]]
            charged.append(1 if is_cat else (-1 if is_an else 0))
            polar.append(1 if el in ("N", "O") else 0)
    coords = np.array(coords); owner = np.array(owner)
    charged = np.array(charged); polar = np.array(polar)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(4.5, output_type="ndarray")
    dmin, saltmin, hbmin = {}, {}, {}
    oi, oj = owner[pairs[:, 0]], owner[pairs[:, 1]]
    keep = oi != oj
    pairs, oi, oj = pairs[keep], oi[keep], oj[keep]
    dvec = np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1)
    ca, cb = charged[pairs[:, 0]], charged[pairs[:, 1]]
    pa, pb = polar[pairs[:, 0]], polar[pairs[:, 1]]
    is_salt = (ca * cb) < 0
    is_hb = (pa == 1) & (pb == 1)
    for k in range(len(pairs)):
        a, b = int(oi[k]), int(oj[k])
        key = (a, b) if a < b else (b, a)
        d = dvec[k]
        if key not in dmin or d < dmin[key]:
            dmin[key] = d
        if is_salt[k] and (key not in saltmin or d < saltmin[key]):
            saltmin[key] = d
        if is_hb[k] and (key not in hbmin or d < hbmin[key]):
            hbmin[key] = d
    return dmin, saltmin, hbmin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bent", default="data/reference_pdbs/1jv2.pdb")
    ap.add_argument("--extended", default="results/route_a/extended_state_b.pdb")
    ap.add_argument("--out-json", default="results/route_a/linchpin_scan.json")
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    bent_ch = load_residues(args.bent)
    ext_ch = load_residues(args.extended)

    # rank-match: same residue set/order; identity + orig numbering from the bent crystal
    meta = []          # global index -> residue metadata (orig numbering)
    bent_flat, ext_flat = [], []
    for ch in sorted(bent_ch):
        bl, el = bent_ch[ch], ext_ch[ch]
        assert len(bl) == len(el), f"chain {ch} residue count mismatch"
        for rb, re in zip(bl, el):
            gi = len(meta)
            up = rb["resid"] <= UPPER[ch][1]
            meta.append({"gidx": gi, "chain": ch, "resid": rb["resid"],
                         "resname": rb["resname"], "side": "upper" if up else "lower",
                         "domain": domain(ch, rb["resid"])})
            bent_flat.append(rb); ext_flat.append(re)

    bd, bs, bh = contact_set(bent_flat)
    ed, es, eh = contact_set(ext_flat)
    bent_pairs, ext_pairs = set(bd), set(ed)
    ext_unique = ext_pairs - bent_pairs
    bent_unique = bent_pairs - ext_pairs

    def spans_knee(i, j):
        mi, mj = meta[i], meta[j]
        return (mi["chain"] != mj["chain"]) or (mi["side"] != mj["side"])

    def near_genu(m):
        return abs(m["resid"] - GENU_CENTER[m["chain"]]) <= 8

    # ---- score residues by extended-unique contacts ----
    score = {}
    detail = {}
    for (i, j) in ext_unique:
        salt = (i, j) in es and es[(i, j)] < 4.0
        hb = (i, j) in eh and eh[(i, j)] < 3.5
        w = 3.0 if salt else (2.0 if hb else 1.0)
        if spans_knee(i, j):
            w *= 2.0
        kind = "salt-bridge" if salt else ("H-bond" if hb else "packing")
        for a, b in ((i, j), (j, i)):
            score[a] = score.get(a, 0.0) + w + (1.0 if near_genu(meta[a]) else 0.0)
            detail.setdefault(a, []).append({
                "partner": f'{meta[b]["chain"]}:{meta[b]["resname"]}{meta[b]["resid"]}',
                "partner_domain": meta[b]["domain"], "type": kind,
                "spans_knee": spans_knee(a, b)})

    ranked = sorted(score, key=lambda g: -score[g])
    linchpins = []
    for g in ranked[:args.top]:
        m = meta[g]
        ds = detail[g]
        salt = [d for d in ds if d["type"] == "salt-bridge"]
        cross = [d for d in ds if d["spans_knee"]]
        linchpins.append({
            "residue": f'{m["chain"]}:{m["resname"]}{m["resid"]}',
            "domain": m["domain"], "score": round(score[g], 1),
            "n_extended_unique_contacts": len(ds),
            "n_cross_knee": len(cross), "n_salt_bridges": len(salt),
            "near_genu": near_genu(m),
            "top_partners": [d["partner"] for d in ds[:5]],
            "salt_bridge_partners": [d["partner"] for d in salt],
        })

    # ---- bent-clasp residues (the opposite lever) ----
    clasp_score = {}
    for (i, j) in bent_unique:
        if spans_knee(i, j):
            for a in (i, j):
                clasp_score[a] = clasp_score.get(a, 0.0) + 1.0
    clasp = []
    for g in sorted(clasp_score, key=lambda g: -clasp_score[g])[:8]:
        m = meta[g]
        clasp.append({"residue": f'{m["chain"]}:{m["resname"]}{m["resid"]}',
                      "domain": m["domain"], "n_clasp_contacts": int(clasp_score[g])})

    # ---- exact bent-vs-extended distance for the top cross-knee salt-bridge locks ----
    def charged_atoms(res):
        out = []
        for nm, el, xyz in res["atoms"]:
            if res["resname"] in CATION and nm in CATION[res["resname"]]:
                out.append(("+", np.array(xyz)))
            if res["resname"] in ANION and nm in ANION[res["resname"]]:
                out.append(("-", np.array(xyz)))
        return out

    def pair_charged_dist(flat, i, j):
        a, b = charged_atoms(flat[i]), charged_atoms(flat[j])
        if not a or not b:
            return None
        return float(min(np.linalg.norm(x[1] - y[1]) for x in a for y in b
                         if x[0] != y[0]) if any(x[0] != y[0] for x in a for y in b) else 99.9)

    gi_by_label = {f'{meta[g]["chain"]}:{meta[g]["resname"]}{meta[g]["resid"]}': g
                   for g in range(len(meta))}
    key_locks = []
    seen_lock = set()
    for L in linchpins:
        gi = gi_by_label[L["residue"]]
        for partner in L["salt_bridge_partners"]:
            gj = gi_by_label.get(partner)
            if gj is None:
                continue
            key = tuple(sorted((gi, gj)))
            if key in seen_lock:
                continue
            seen_lock.add(key)
            db = pair_charged_dist(bent_flat, gi, gj)
            de = pair_charged_dist(ext_flat, gi, gj)
            if db is None or de is None:
                continue
            key_locks.append({
                "pair": f'{L["residue"]} — {partner}',
                "bent_dist_A": round(db, 2), "extended_dist_A": round(de, 2),
                "forms_on_extension": bool(db > 6.0 and de < 4.0)})
    key_locks.sort(key=lambda k: k["bent_dist_A"] - k["extended_dist_A"], reverse=True)

    # per-residue score array for the figure
    per_res = [{"chain": meta[g]["chain"], "resid": meta[g]["resid"],
                "resname": meta[g]["resname"], "score": round(score.get(g, 0.0), 2)}
               for g in range(len(meta))]

    result = {
        "analysis": "extended-state linchpin scan (contact/interaction differential)",
        "n_residues": len(meta),
        "n_contacts_bent": len(bent_pairs), "n_contacts_extended": len(ext_pairs),
        "n_extended_unique": len(ext_unique), "n_bent_unique": len(bent_unique),
        "extended_lock_linchpins": linchpins,
        "key_cross_knee_salt_bridges": key_locks,
        "bent_clasp_residues": clasp,
        "per_residue_score": per_res,
        "caveat": ("structure-based prediction from the current extended model; the "
                   "definitive test is mutate->MD and observe snap-back (GPU follow-up)."),
    }
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"residues={len(meta)}  contacts bent={len(bent_pairs)} ext={len(ext_pairs)}")
    print(f"extended-unique contacts={len(ext_unique)}  bent-unique(clasp)={len(bent_unique)}")
    print(f"\nTop {args.top} extended-lock linchpins (mutate -> predicted snap-back):")
    for L in linchpins:
        sb = f" | salt-bridge: {','.join(L['salt_bridge_partners'])}" if L["salt_bridge_partners"] else ""
        genu = " [GENU]" if L["near_genu"] else ""
        print(f"  {L['residue']:<12} {L['domain']:<28} score={L['score']:<5} "
              f"x-knee={L['n_cross_knee']}{genu}{sb}")
    print(f"\nKey cross-knee salt-bridge locks (charged-atom distance, bent -> extended):")
    for k in key_locks:
        flag = "  <== forms on extension" if k["forms_on_extension"] else ""
        print(f"  {k['pair']:<28} {k['bent_dist_A']:>6} A  ->  {k['extended_dist_A']:>5} A{flag}")
    print(f"\nBent-clasp residues (mutate -> would STABILIZE extended, opposite lever):")
    for c in clasp:
        print(f"  {c['residue']:<12} {c['domain']:<28} clasp-contacts={c['n_clasp_contacts']}")
    print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
