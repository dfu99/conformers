#!/usr/bin/env python3
"""Transfer the genu Ca2+ from the 1JV2 crystal into a morphed extended alphaV structure.

WHY. `extended_state_b.pdb` carries no metals at all, and the four "genu lock" ion pairs
obj-078..obj-083 are built on do not exist in the 1JV2 crystal OR in the pre-minimisation morph
seed -- they appear only after the ion-free, unscreened vacuum relaxation. 1JV2 puts a
well-coordinated Ca2+ (A4008) in the middle of that site:

    GLU636.OE1 2.09 A   ASP599.OD1 2.11 A   VAL601.O 2.54 A   GLY597.O 2.55 A

so E636's carboxylate -- one half of the monitored K459-E636 pair -- is a first-shell ligand of an
ion this model deletes, and Xiong 2001 predicted that this ion neutralises the thigh/calf-1 acidic
interface in an extended integrin. Putting it back is what decides whether the network is real.

HOW. Superpose 1JV2's coordinating atoms onto the same atoms of the target (chain A numbering is
identical for these residues: extended_state_b chain A resSeq 1-838 == alphaV mature 1-838), then
apply that transform to the ion. Superposing the SITE rather than the whole domain keeps the ion in
its own coordination geometry, and the fit RMSD then reports honestly whether the morph deformed
the site: a large RMSD means the crystal ion has no well-defined home here, which is itself a
finding and not something to paper over.

    python add_genu_calcium.py --target results/route_a/extended_state_b.pdb \
                               --out results/route_a/extended_state_b_ca.pdb
"""
import argparse
import numpy as np

# 1JV2 genu Ca2+ (A4008) first shell. Atom -> the residue it belongs to, in chain A.
SITE = [(597, "O"), (599, "OD1"), (599, "OD2"), (601, "O"), (636, "OE1"), (636, "OE2")]
ION_RES, ION_CHAIN = 4008, "A"


def read_atoms(path):
    """{(chain, resSeq, atomName): xyz} plus the raw lines, so we can write the file back out."""
    xyz, lines = {}, []
    for l in open(path):
        lines.append(l)
        if l.startswith(("ATOM", "HETATM")):
            try:
                xyz[(l[21], int(l[22:26]), l[12:16].strip())] = np.array(
                    [float(l[30 + 8 * i:38 + 8 * i]) for i in range(3)])
            except ValueError:
                pass
    return xyz, lines


def kabsch(P, Q):
    """Rigid transform taking P onto Q (both N x 3). Returns (R, t, rmsd)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crystal", default="data/reference_pdbs/1jv2.pdb",
                    help="1JV2, the source of the genu Ca2+ and its coordination geometry")
    ap.add_argument("--target", default="results/route_a/extended_state_b.pdb")
    ap.add_argument("--out", default="results/route_a/extended_state_b_ca.pdb")
    ap.add_argument("--max-rmsd", type=float, default=1.5,
                    help="abort if the coordination site does not superpose this well (A)")
    args = ap.parse_args()

    cx, _ = read_atoms(args.crystal)
    tx, tlines = read_atoms(args.target)

    ion = cx.get((ION_CHAIN, ION_RES, "CA"))
    if ion is None:
        raise SystemExit(f"no Ca2+ {ION_CHAIN}{ION_RES} in {args.crystal}")

    pairs = [(cx[("A", r, a)], tx[("A", r, a)]) for r, a in SITE
             if ("A", r, a) in cx and ("A", r, a) in tx]
    missing = [f"A{r}.{a}" for r, a in SITE if ("A", r, a) not in cx or ("A", r, a) not in tx]
    if len(pairs) < 4:
        raise SystemExit(f"only {len(pairs)} of {len(SITE)} site atoms present (missing {missing})")

    P = np.array([p for p, _ in pairs])
    Q = np.array([q for _, q in pairs])
    R, t, rmsd = kabsch(P, Q)
    print(f"site superposition on {len(pairs)} atoms: RMSD {rmsd:.2f} A"
          + (f"  (missing {', '.join(missing)})" if missing else ""))
    if rmsd > args.max_rmsd:
        raise SystemExit(
            f"coordination site RMSD {rmsd:.2f} A exceeds {args.max_rmsd} A — the morph deformed "
            f"this site, so the crystal ion has no well-defined position here. Refusing to invent "
            f"one; re-morph with the ion present instead.")

    pos = R @ ion + t
    print(f"placed Ca2+ at {pos[0]:8.3f} {pos[1]:8.3f} {pos[2]:8.3f}")
    print("resulting coordination in the target:")
    for r, a in SITE:
        v = tx.get(("A", r, a))
        if v is not None:
            print(f"   A{r}.{a:<4s} {np.linalg.norm(v - pos):5.2f} A")

    # Write the ion in as the last HETATM before END, using its own residue+chain so nothing
    # downstream can confuse it with a protein CA atom.
    card = (f"HETATM99999 CA    CA {ION_CHAIN}{ION_RES:4d}    "
            f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          CA\n")
    with open(args.out, "w") as fh:
        for l in tlines:
            if l.startswith("END"):
                fh.write(card)
            fh.write(l)
        if not any(l.startswith("END") for l in tlines):
            fh.write(card)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
