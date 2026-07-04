#!/usr/bin/env python3
"""Route-A Stage 2b: relax the extended-pose SEED into a usable state B.

The Stage-2 seed (extended_seed.pdb) is a rigid-body genu leg-swing, so the
genu (knee) hinge carries strained bonds/angles and steric clashes where the
rotated lower body meets the fixed headpiece. This script relaxes that strain.

CPU-ONLY BY DESIGN: runs on the OpenMM CPU platform with bounded threads so it
touches ZERO GPU — it does not contend with the GPU-bound oxDNA/FFS jobs sharing
this pod. Implicit solvent (GB-OBC2) keeps the system at ~22k atoms (protein+H)
instead of ~477k solvated, so a clash-relief minimize + short equilibration is
minutes on CPU, not hours.

Pipeline: pdbfixer (patch missing side-chain atoms + add H) -> amber14 ff14SB +
implicit OBC2 -> energy minimize (removes the hinge clashes) -> short Langevin
settle at 300 K -> report Rg / long-axis extent (real extension metrics) before
and after, and write relaxed_seed.pdb.

The relaxed seed is the endpoint for the full solvated GPU relax + string method;
implicit-solvent CPU relaxation is the non-contending first pass that makes the
seed physically valid (finite energy, no clashes) while the GPU is busy.
"""
import argparse
import json
import time

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer


def ca_coords_by_chain(topology, positions):
    """Return {chain_id: {resseq:int -> np.array xyz(Angstrom)}} for CA atoms."""
    pos = np.array(positions.value_in_unit(unit.angstrom))
    out = {}
    for atom in topology.atoms():
        if atom.name != "CA":
            continue
        cid = atom.residue.chain.id
        try:
            rs = int(atom.residue.id)
        except ValueError:
            continue
        out.setdefault(cid, {})[rs] = pos[atom.index]
    return out


def centroid(cabychain, chain, lo, hi):
    d = cabychain.get(chain, {})
    pts = [xyz for rs, xyz in d.items() if lo <= rs <= hi]
    return np.mean(pts, axis=0) if pts else None


def metrics(topology, positions):
    """Global extension metrics from CA atoms: Rg, long-axis extent, cv1."""
    cabychain = ca_coords_by_chain(topology, positions)
    allca = np.array([xyz for d in cabychain.values() for xyz in d.values()])
    com = allca.mean(axis=0)
    rg = float(np.sqrt(np.mean(np.sum((allca - com) ** 2, axis=1))))
    # long axis = first principal component; extent = spread of projections
    centred = allca - com
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    proj = centred @ vt[0]
    extent = float(proj.max() - proj.min())
    # placeholder cv1 (beta head-tail) for continuity with prior reports
    a = centroid(cabychain, "B", 1, 352)
    b = centroid(cabychain, "B", 353, 690)
    cv1 = float(np.linalg.norm(a - b)) if a is not None and b is not None else None
    return {"Rg_A": round(rg, 2), "long_axis_extent_A": round(extent, 2),
            "cv1_beta_head_tail_A": round(cv1, 2) if cv1 is not None else None,
            "n_ca": int(len(allca))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="results/route_a/extended_seed.pdb")
    ap.add_argument("--out-pdb", default="results/route_a/relaxed_seed.pdb")
    ap.add_argument("--out-json", default="results/route_a/stage2b_relax_result.json")
    ap.add_argument("--equil-ps", type=float, default=25.0)
    ap.add_argument("--threads", default="32")
    ap.add_argument("--min-iter", type=int, default=0, help="0 = minimize to convergence")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[stage2b] loading {args.seed}", flush=True)
    fixer = PDBFixer(filename=args.seed)
    fixer.findMissingResidues()
    fixer.missingResidues = {}            # do not rebuild loops, only patch atoms
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)  # drop structural Ca2+ for CPU pass
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    n_atoms = fixer.topology.getNumAtoms()
    print(f"[stage2b] patched topology: {n_atoms} atoms", flush=True)

    m_before = metrics(fixer.topology, fixer.positions)
    print(f"[stage2b] seed metrics: {m_before}", flush=True)

    ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/obc2.xml")
    system = ff.createSystem(
        fixer.topology, nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=2.0 * unit.nanometer, constraints=app.HBonds)
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)
    platform = mm.Platform.getPlatformByName("CPU")
    sim = app.Simulation(fixer.topology, system, integrator, platform,
                         {"Threads": args.threads})
    sim.context.setPositions(fixer.positions)

    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    print(f"[stage2b] E(seed) = {e0:.1f} kJ/mol -- minimizing...", flush=True)
    sim.minimizeEnergy(maxIterations=args.min_iter)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    st = sim.context.getState(getPositions=True)
    m_min = metrics(fixer.topology, st.getPositions())
    print(f"[stage2b] E(min)  = {e1:.1f} kJ/mol  metrics={m_min}", flush=True)
    with open(args.out_pdb, "w") as fh:
        app.PDBFile.writeFile(fixer.topology, st.getPositions(), fh)

    result = {
        "stage": "2b_relax_seed_cpu_implicit",
        "platform": "CPU", "threads": args.threads,
        "n_atoms": n_atoms,
        "E_seed_kJ_mol": round(e0, 1),
        "E_min_kJ_mol": round(e1, 1),
        "metrics_seed": m_before,
        "metrics_min": m_min,
        "minimize_seconds": round(time.time() - t0, 1),
    }
    # write the minimize result immediately so progress is visible before equil
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    print("[stage2b] MINIMIZE_DONE", flush=True)

    if args.equil_ps > 0:
        n_steps = int(args.equil_ps * 1000 / 2.0)  # 2 fs steps
        print(f"[stage2b] equil {args.equil_ps} ps ({n_steps} steps)...", flush=True)
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
        chunk = max(1, n_steps // 10)
        done = 0
        while done < n_steps:
            k = min(chunk, n_steps - done)
            sim.step(k)
            done += k
            e = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)
            print(f"[stage2b]   {done}/{n_steps} steps  E={e:.1f} kJ/mol", flush=True)
        st = sim.context.getState(getPositions=True, getEnergy=True)
        m_eq = metrics(fixer.topology, st.getPositions())
        e2 = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        with open(args.out_pdb, "w") as fh:
            app.PDBFile.writeFile(fixer.topology, st.getPositions(), fh)
        result.update({
            "equil_ps": args.equil_ps,
            "E_equil_kJ_mol": round(e2, 1),
            "metrics_equil": m_eq,
            "total_seconds": round(time.time() - t0, 1),
        })
        with open(args.out_json, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"[stage2b] EQUIL_DONE E={e2:.1f}  metrics={m_eq}", flush=True)

    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
