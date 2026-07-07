#!/usr/bin/env python3
"""Definitive snap-back test: mutate an extended-state linchpin, run MD, watch the knee.

Hypothesis (obj-078): the extended αVβ3 state is held open by a cross-knee salt-bridge
network at the genu -- αV thigh K459/D457 <-> calf-1 E598/D595 + K688/E636. Mutating
these should release the lock and let the extended state relax back toward bent (knee
angle 175° -> smaller; Rg 66 -> smaller).

This runs unbiased MD in implicit solvent (GB-OBC2, amber14/ff14SB) on the 22k-atom
ectodomain -- no explicit water (~4 GB, GPU-friendly; no solvent viscosity so large
domain motion is reachable in short MD). Tracks the knee angle, Rg, and the key
salt-bridge distances over time. Run WT (control, should stay extended) and mutants
side by side; the WT-vs-mutant DIFFERENTIAL in knee angle is the signal.

  python snapback_md.py --mutations WT              --ns 20 --platform CUDA --tag wt
  python snapback_md.py --mutations A:459:ALA       --ns 20 --platform CUDA --tag k459a
  python snapback_md.py --mutations A:459:ALA,A:598:ALA --ns 20 --platform CUDA --tag double

CPU smoke: --smoke  (tiny minimize + a few steps, to validate the build/observables).
"""
import argparse
import csv
import json
import os
import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit
from pdbfixer import PDBFixer

# knee split + genu, ORIGINAL 1JV2 numbering (== extended_state_b numbering for these residues)
UPPER = {"A": (1, 592), "B": (1, 440)}
LOWER = {"A": (593, 956), "B": (441, 690)}
GENU = {"A": (588, 596), "B": (436, 444)}
# key cross-knee salt-bridge atom pairs to monitor  (anion res/atoms, cation res/atoms)
SALT_PAIRS = [
    ("E598-K459", ("A", 598, ["OE1", "OE2"]), ("A", 459, ["NZ"])),
    ("D457-K688", ("A", 457, ["OD1", "OD2"]), ("A", 688, ["NZ"])),
    ("K459-E636", ("A", 636, ["OE1", "OE2"]), ("A", 459, ["NZ"])),
    ("D595-K688", ("A", 595, ["OD1", "OD2"]), ("A", 688, ["NZ"])),
]


def strip_h(in_pdb, out_pdb):
    with open(in_pdb) as fh, open(out_pdb, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                nm = line[12:16].strip()
                el = line[76:78].strip() or (nm[0] if nm else "")
                if el == "H" or (not el and nm.startswith("H")):
                    continue
            out.write(line)


def build(pdb, mutations, workdir, implicit=True):
    heavy = os.path.join(workdir, "heavy.pdb")
    strip_h(pdb, heavy)
    fixer = PDBFixer(filename=heavy)
    if mutations and mutations.upper() != "WT":
        by_chain = {}
        for m in mutations.split(","):
            ch, resid, new = m.split(":")
            # current resname at (ch, resid)
            by_chain.setdefault(ch, []).append((int(resid), new))
        # need old resname; read from the heavy pdb
        oldname = {}
        for line in open(heavy):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                oldname[(line[21], int(line[22:26]))] = line[17:20].strip()
        for ch, muts in by_chain.items():
            spec = [f"{oldname[(ch, rid)]}-{rid}-{new}" for rid, new in muts]
            fixer.applyMutations(spec, ch)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    if implicit:
        ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/obc2.xml")
    else:  # vacuum -- fast validation path
        ff = app.ForceField("amber14/protein.ff14SB.xml")
    system = ff.createSystem(fixer.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=2.0 * unit.nanometer, constraints=app.HBonds)
    return fixer.topology, fixer.positions, system


def ca_index(topology):
    idx = {}
    for atom in topology.atoms():
        if atom.name == "CA":
            idx[(atom.residue.chain.id, int(atom.residue.id))] = atom.index
    return idx


def in_rng(ch, r, tbl):
    return tbl[ch][0] <= r <= tbl[ch][1]


def atom_index(topology):
    idx = {}
    for atom in topology.atoms():
        idx[(atom.residue.chain.id, int(atom.residue.id), atom.name)] = atom.index
    return idx


def make_observables(topology):
    cai = ca_index(topology)
    upper = [i for (ch, r), i in cai.items() if in_rng(ch, r, UPPER)]
    lower = [i for (ch, r), i in cai.items() if in_rng(ch, r, LOWER)]
    genu = [i for (ch, r), i in cai.items() if in_rng(ch, r, GENU)]
    allca = list(cai.values())
    ai = atom_index(topology)
    salt = []
    for name, (ach, ar, aatoms), (cch, cr, catoms) in SALT_PAIRS:
        an = [ai[(ach, ar, a)] for a in aatoms if (ach, ar, a) in ai]
        cn = [ai[(cch, cr, a)] for a in catoms if (cch, cr, a) in ai]
        if an and cn:
            salt.append((name, an, cn))

    def compute(pos_nm):
        p = pos_nm
        hinge = p[genu].mean(0)
        vu = p[upper].mean(0) - hinge
        vl = p[lower].mean(0) - hinge
        cos = np.dot(vu, vl) / (np.linalg.norm(vu) * np.linalg.norm(vl))
        knee = float(np.degrees(np.arccos(np.clip(cos, -1, 1))))
        c = p[allca].mean(0)
        rg = float(np.sqrt(((p[allca] - c) ** 2).sum(1).mean()) * 10.0)  # nm->A
        sd = {}
        for name, an, cn in salt:
            d = min(np.linalg.norm(p[a] - p[c2]) for a in an for c2 in cn) * 10.0
            sd[name] = round(float(d), 2)
        return round(knee, 2), round(rg, 2), sd

    return compute


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="results/route_a/extended_state_b.pdb")
    ap.add_argument("--mutations", default="WT")
    ap.add_argument("--tag", default="wt")
    ap.add_argument("--ns", type=float, default=20.0)
    ap.add_argument("--report-ps", type=float, default=20.0)
    ap.add_argument("--platform", default="CUDA")
    ap.add_argument("--threads", default="8")
    ap.add_argument("--temperature", type=float, default=310.0)
    ap.add_argument("--out-dir", default="results/route_a/snapback")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--vacuum", action="store_true",
                    help="build in vacuum (fast; for validation). Default implicit solvent.")
    args = ap.parse_args()

    outdir = os.path.join(args.out_dir, args.tag)
    os.makedirs(outdir, exist_ok=True)

    implicit = not (args.vacuum or args.smoke)  # smoke validates in vacuum for speed
    topo, pos, system = build(args.pdb, args.mutations, outdir, implicit=implicit)
    print(f"[{args.tag}] force field: {'implicit-solvent GB-OBC2' if implicit else 'vacuum'}", flush=True)
    integ = mm.LangevinMiddleIntegrator(args.temperature * unit.kelvin,
                                        1.0 / unit.picosecond, 0.002 * unit.picosecond)
    if args.platform == "CPU":
        plat = mm.Platform.getPlatformByName("CPU")
        props = {"Threads": args.threads}
    else:
        plat = mm.Platform.getPlatformByName(args.platform)
        props = {"Precision": "mixed"} if args.platform == "CUDA" else {}
    sim = app.Simulation(topo, system, integ, plat, props)
    sim.context.setPositions(pos)

    obs = make_observables(topo)
    p0 = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.nanometer))
    knee0, rg0, sd0 = obs(p0)
    print(f"[{args.tag}] built: {system.getNumParticles()} atoms  start knee={knee0}° Rg={rg0}Å  salt={sd0}", flush=True)

    min_iter = 50 if args.smoke else 2000
    sim.minimizeEnergy(maxIterations=min_iter)
    sim.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)

    steps_total = 250 if args.smoke else int(args.ns * 1000 / 0.002)
    steps_chunk = 50 if args.smoke else int(args.report_ps * 1000 / 0.002)
    csv_path = os.path.join(outdir, "trajectory.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ps", "knee_deg", "Rg_A"] + [p[0] for p in SALT_PAIRS if p[0] in sd0])
        done = 0
        while done < steps_total:
            n = min(steps_chunk, steps_total - done)
            sim.step(n)
            done += n
            st = sim.context.getState(getPositions=True, getEnergy=True)
            p = np.array(st.getPositions().value_in_unit(unit.nanometer))
            knee, rg, sd = obs(p)
            row = [round(done * 0.002, 1), knee, rg] + [sd.get(k, "") for k in sd0]
            w.writerow(row); fh.flush()
            print(f"[{args.tag}] {done*0.002:8.1f} ps  knee={knee:6.1f}°  Rg={rg:5.1f}Å  {sd}", flush=True)
    sim.saveState(os.path.join(outdir, "final.xml"))
    with open(os.path.join(outdir, "final.pdb"), "w") as fh:
        app.PDBFile.writeFile(topo, sim.context.getState(getPositions=True).getPositions(), fh)
    summary = {"tag": args.tag, "mutations": args.mutations, "ns": args.ns,
               "start": {"knee": knee0, "Rg": rg0, "salt": sd0},
               "end": {"knee": knee, "Rg": rg, "salt": sd}}
    json.dump(summary, open(os.path.join(outdir, "summary.json"), "w"), indent=2)
    print(f"[{args.tag}] DONE  knee {knee0}->{knee}°  Rg {rg0}->{rg}Å  -> {csv_path}", flush=True)


if __name__ == "__main__":
    main()
