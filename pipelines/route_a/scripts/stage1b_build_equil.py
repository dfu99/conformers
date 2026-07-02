#!/usr/bin/env python3
"""Route-A Stage 1b: bent αVβ3 + structural Ca²⁺, solvent-only equilibration.

Refines Stage-1a (stage1_build_equil.py, solvent-only, no ions) by preserving
the 6 structural Ca²⁺ ions from 1JV2 (5 in the αV β-propeller, 1 in β3). These
rigidify the propeller/leg and are expected to reduce the cv1 β-leg drift seen
in Stage-1a (28.36 → 25.93 Å over 200 ps).

Scientific note: 1JV2 is the αVβ3 **ectodomain** (chain A ends res 956, chain B
res 690 — no transmembrane helices), so there is nothing to embed in a lipid
bilayer. The kickoff's "endothelial membrane" step is therefore inapplicable to
this construct; solvent-only is the correct model. A membrane would only be
needed for a full-length (TM-containing) integrin. This removes route A's
hardest/most error-prone build step (kickoff §2.2 Risk 2).

--build-only stops after createSystem (fast template check, no MD).
Requires the route_a conda env (openmm, pdbfixer).
"""
import argparse
import json
import os
import time

import numpy as np
import openmm as mm
from openmm import app, unit, Vec3
from pdbfixer import PDBFixer

CV_DEFS = {
    "cv0_head_calf": (("A", 1, 435), ("A", 442, 605)),
    "cv1_beta_head_tail": (("B", 1, 352), ("B", 353, 690)),
    "cv2_head_open": (("A", 1, 435), ("B", 1, 352)),
}


def parse_ca_ions(pdb_path):
    """Return list of Ca²⁺ positions (nm) from HETATM CA/CA records."""
    ions = []
    with open(pdb_path) as fh:
        for line in fh:
            if (line.startswith("HETATM") and line[12:16].strip() == "CA"
                    and line[17:20].strip() == "CA"):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                ions.append(Vec3(x / 10.0, y / 10.0, z / 10.0))
    return ions


def ca_coords_by_chain(topology, positions):
    out = {}
    pos = positions.value_in_unit(unit.nanometer)
    for atom in topology.atoms():
        if atom.name != "CA" or atom.residue.name == "CA":  # skip calcium ions
            continue
        ch = atom.residue.chain.id
        try:
            rs = int(atom.residue.id)
        except (ValueError, TypeError):
            continue
        out.setdefault(ch, []).append((rs, np.array(pos[atom.index])))
    return out


def compute_cvs(topology, positions):
    ca = ca_coords_by_chain(topology, positions)
    res = {}
    for name, (a, b) in CV_DEFS.items():
        pa = [xyz for (rs, xyz) in ca.get(a[0], []) if a[1] <= rs <= a[2]]
        pb = [xyz for (rs, xyz) in ca.get(b[0], []) if b[1] <= rs <= b[2]]
        if pa and pb:
            res[name] = round(float(np.linalg.norm(np.mean(pa, 0) - np.mean(pb, 0)) * 10.0), 2)
        else:
            res[name] = None
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--outdir", default="/workspace/route_a/stage1b")
    ap.add_argument("--nvt-ps", type=float, default=100.0)
    ap.add_argument("--npt-ps", type=float, default=100.0)
    ap.add_argument("--dt-fs", type=float, default=2.0)
    ap.add_argument("--build-only", action="store_true",
                    help="stop after createSystem (fast FF-template check)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    log = {"pdb": args.pdb, "stage": "1b_solvent_plus_Ca"}
    t0 = time.time()

    ca_ions = parse_ca_ions(args.pdb)
    log["n_ca_ions"] = len(ca_ions)
    print(f"[0] parsed {len(ca_ions)} structural Ca2+ ions", flush=True)

    print("[1] PDBFixer clean (protein)", flush=True)
    fixer = PDBFixer(filename=args.pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    modeller = app.Modeller(fixer.topology, fixer.positions)
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller.addHydrogens(ff, pH=7.0)

    # add the structural Ca2+ ions back as a dedicated chain
    ca_top = app.Topology()
    ch = ca_top.addChain("I")
    for _ in ca_ions:
        r = ca_top.addResidue("CA", ch)
        ca_top.addAtom("CA", app.element.calcium, r)
    modeller.add(ca_top, unit.Quantity(ca_ions, unit.nanometer))
    print(f"[2] added Ca2+; solvating", flush=True)

    modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer,
                        ionicStrength=0.15 * unit.molar, neutralize=True)
    log["n_atoms_solvated"] = modeller.topology.getNumAtoms()
    log["cvs_after_build"] = compute_cvs(modeller.topology, modeller.positions)
    print("    solvated atoms:", log["n_atoms_solvated"],
          "| CVs:", log["cvs_after_build"], flush=True)

    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds)
    print("[3] system created OK (FF templates matched, incl. Ca2+)", flush=True)
    log["build_seconds"] = round(time.time() - t0, 1)

    if args.build_only:
        log["mode"] = "build-only"
        with open(os.path.join(args.outdir, "stage1b_build.json"), "w") as fh:
            json.dump(log, fh, indent=2)
        print("[build-only done]", json.dumps(log, indent=2), flush=True)
        return

    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                             args.dt_fs * unit.femtoseconds)
    sim = app.Simulation(modeller.topology, system, integrator,
                         mm.Platform.getPlatformByName("CUDA"), {"Precision": "mixed"})
    sim.context.setPositions(modeller.positions)
    print("[4] minimize", flush=True)
    sim.minimizeEnergy(maxIterations=5000)
    print("[5] NVT", flush=True)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
    sim.step(int(args.nvt_ps * 1000 / args.dt_fs))
    print("[6] NPT", flush=True)
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin))
    sim.context.reinitialize(preserveState=True)
    sim.step(int(args.npt_ps * 1000 / args.dt_fs))

    state = sim.context.getState(getPositions=True, getEnergy=True)
    log["cvs_equilibrated"] = compute_cvs(modeller.topology, state.getPositions())
    log["potential_kJ_mol"] = round(
        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), 1)
    log["wall_seconds"] = round(time.time() - t0, 1)
    app.PDBFile.writeFile(modeller.topology, state.getPositions(),
                          open(os.path.join(args.outdir, "equilibrated.pdb"), "w"))
    with open(os.path.join(args.outdir, "stage1b_equil.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    print("[done]", json.dumps(log, indent=2), flush=True)


if __name__ == "__main__":
    main()
