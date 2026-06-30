#!/usr/bin/env python3
"""Route-A Stage 1 (part b): build solvated αVβ3 system + short equilibration.

Stage-1 gate (docs/route_a_kickoff.md §4): confirm the ported setup is sound
before the weeks-long Stage-2 string optimization. This first pass uses the
αVβ3 ectodomain in **explicit water + 0.15 M ions** (no membrane yet) — the
cheapest correctness test of the topology→MD→CV chain. The endothelial
membrane is the Stage-1b increment once this solvent-only run is stable.

Pipeline:
  1. PDBFixer clean (add missing heavy atoms + H at pH 7). Missing *residues*
     (crystal loop gaps) are NOT rebuilt — avoids fabricating large unphysical
     loops for a stability test. Ca2+/heteroatoms are dropped in this first
     pass (a documented simplification; structural Ca2+ is a Stage-1b add — see
     kickoff §2.4 Risk 2).
  2. amber14-all + TIP3P force field; solvate (1.0 nm padding), neutralize.
  3. Minimize → NVT (LangevinMiddle 300 K) → NPT (MonteCarloBarostat 1 bar).
  4. Recompute the route-A CVs on the equilibrated protein and compare to the
     bent baseline (stage1_cv_baseline.py) — the gate is "bent geometry holds
     at 300 K" (CVs within a few Å of baseline, no blow-up).

Requires the route_a conda env (openmm, pdbfixer). UNVALIDATED until first run
against the live env — treat numbers as provisional until confirmed.
"""
import argparse
import json
import os
import time

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer


# route-A CV residue ranges (αV numbering), from av_cvs_remapped.json
CV_DEFS = {
    "cv0_head_calf": (("A", 1, 435), ("A", 442, 605)),
    "cv1_beta_head_tail": (("B", 1, 352), ("B", 353, 690)),
    "cv2_head_open": (("A", 1, 435), ("B", 1, 352)),
}


def ca_coords_by_chain(topology, positions):
    """{chain_id: [(resseq, np.array xyz_nm)]} for CA atoms."""
    out = {}
    pos = positions.value_in_unit(unit.nanometer)
    for atom in topology.atoms():
        if atom.name != "CA":
            continue
        ch = atom.residue.chain.id
        try:
            rs = int(atom.residue.id)
        except (ValueError, TypeError):
            continue
        out.setdefault(ch, []).append((rs, np.array(pos[atom.index])))
    return out


def centroid_range(ca_map, chain, lo, hi):
    pts = [xyz for (rs, xyz) in ca_map.get(chain, []) if lo <= rs <= hi]
    return np.mean(pts, axis=0) if pts else None


def compute_cvs(topology, positions):
    ca = ca_coords_by_chain(topology, positions)
    res = {}
    for name, (a, b) in CV_DEFS.items():
        ca_a = centroid_range(ca, *a)
        ca_b = centroid_range(ca, *b)
        if ca_a is None or ca_b is None:
            res[name] = None
        else:
            res[name] = round(float(np.linalg.norm(ca_a - ca_b) * 10.0), 2)  # nm->Å
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--outdir", default="/workspace/route_a/stage1")
    ap.add_argument("--nvt-ps", type=float, default=100.0)
    ap.add_argument("--npt-ps", type=float, default=100.0)
    ap.add_argument("--dt-fs", type=float, default=2.0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    log = {"pdb": args.pdb, "stage": "1b_solvent_only"}
    t0 = time.time()

    # 1. clean
    print("[1] PDBFixer clean", flush=True)
    fixer = PDBFixer(filename=args.pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}            # do NOT rebuild crystal loop gaps
    fixer.removeHeterogens(keepWater=False)  # drop Ca2+/ligands this pass
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()  # heavy atoms only; H added once by Modeller below
    app.PDBFile.writeFile(fixer.topology, fixer.positions,
                          open(os.path.join(args.outdir, "fixed_dry.pdb"), "w"))
    log["cvs_after_fix"] = compute_cvs(fixer.topology, fixer.positions)
    print("    CVs after fix:", log["cvs_after_fix"], flush=True)

    # 2. force field + solvate
    print("[2] solvate + neutralize", flush=True)
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(ff, pH=7.0)
    modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer,
                        ionicStrength=0.15 * unit.molar, neutralize=True)
    log["n_atoms_solvated"] = modeller.topology.getNumAtoms()
    print("    solvated atoms:", log["n_atoms_solvated"], flush=True)

    # 3. system
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds)
    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin,
                                             1.0 / unit.picosecond,
                                             args.dt_fs * unit.femtoseconds)
    platform = mm.Platform.getPlatformByName("CUDA")
    sim = app.Simulation(modeller.topology, system, integrator, platform,
                         {"Precision": "mixed"})
    sim.context.setPositions(modeller.positions)

    print("[3] minimize", flush=True)
    sim.minimizeEnergy(maxIterations=5000)

    print("[4] NVT", flush=True)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
    sim.step(int(args.nvt_ps * 1000 / args.dt_fs))

    print("[5] NPT", flush=True)
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin))
    sim.context.reinitialize(preserveState=True)
    sim.step(int(args.npt_ps * 1000 / args.dt_fs))

    # 4. recompute CVs on equilibrated protein
    state = sim.context.getState(getPositions=True, getEnergy=True)
    eq_pos = state.getPositions()
    log["cvs_equilibrated"] = compute_cvs(modeller.topology, eq_pos)
    log["potential_kJ_mol"] = round(
        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), 1)
    log["wall_seconds"] = round(time.time() - t0, 1)
    app.PDBFile.writeFile(
        modeller.topology, eq_pos,
        open(os.path.join(args.outdir, "equilibrated.pdb"), "w"))
    print("    CVs equilibrated:", log["cvs_equilibrated"], flush=True)

    with open(os.path.join(args.outdir, "stage1b_equil.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    print("[done]", json.dumps(log, indent=2), flush=True)


if __name__ == "__main__":
    main()
