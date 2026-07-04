#!/usr/bin/env python3
"""Route-A Stage 2c: build the extended pose by INCREMENTAL leg-swing + relaxation.

The one-shot ~140 deg rigid swing (Stage 2 / 2b) leaves a genu (knee) clash so
severe the minimizer stalls: E ~ 1e13 kJ/mol, the L-BFGS line search can't get a
foothold and never descends. Instead, morph the legs open in small angular steps,
relaxing after each, so the structure never accumulates a catastrophic clash --
every small step is trivial to minimize.

CPU-ONLY (implicit solvent GB-OBC2), bounded threads -> zero GPU, no contention
with the oxDNA FFS jobs sharing this pod.

Loads the bent structure, builds one implicit-solvent context, minimizes it, then
repeats: rotate the lower-leg atoms by step_deg about the genu hinge -> minimize
-> record energy + extension. Produces the extended endpoint AND every intermediate
frame (a ready-made initial path for the string method). The swing is adaptive:
each step rotates toward the "legs opposite the head" target and stops when the
remaining angle is small.
"""
import argparse
import json
import os
import time

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

UPPER = {"A": (1, 592), "B": (1, 440)}
LOWER = {"A": (593, 956), "B": (441, 690)}
GENU = {"A": (588, 596), "B": (436, 444)}


def rot_matrix(axis, theta):
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = axis
    return np.array([
        [c + x*x*(1-c),   x*y*(1-c) - z*s, x*z*(1-c) + y*s],
        [y*x*(1-c) + z*s, c + y*y*(1-c),   y*z*(1-c) - x*s],
        [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)],
    ])


def build_index_sets(top):
    lower_idx, genu_ca, upper_ca, lower_ca, all_ca, b_head, b_tail = ([] for _ in range(7))
    for atom in top.atoms():
        cid = atom.residue.chain.id
        try:
            rs = int(atom.residue.id)
        except ValueError:
            continue
        is_ca = atom.name == "CA"
        if is_ca:
            all_ca.append(atom.index)
        lo, up, gn = LOWER.get(cid), UPPER.get(cid), GENU.get(cid)
        if lo and lo[0] <= rs <= lo[1]:
            lower_idx.append(atom.index)
            if is_ca:
                lower_ca.append(atom.index)
        if up and up[0] <= rs <= up[1] and is_ca:
            upper_ca.append(atom.index)
        if gn and gn[0] <= rs <= gn[1] and is_ca:
            genu_ca.append(atom.index)
        if cid == "B" and is_ca:
            if 1 <= rs <= 352:
                b_head.append(atom.index)
            elif 353 <= rs <= 690:
                b_tail.append(atom.index)
    return {k: np.array(v) for k, v in dict(
        lower_idx=lower_idx, genu_ca=genu_ca, upper_ca=upper_ca, lower_ca=lower_ca,
        all_ca=all_ca, b_head=b_head, b_tail=b_tail).items()}


def metrics(pos_nm, idx):
    ca = pos_nm[idx["all_ca"]] * 10.0  # Angstrom
    com = ca.mean(axis=0)
    rg = float(np.sqrt(np.mean(np.sum((ca - com) ** 2, axis=1))))
    centred = ca - com
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    proj = centred @ vt[0]
    extent = float(proj.max() - proj.min())
    cv1 = None
    if len(idx["b_head"]) and len(idx["b_tail"]):
        a = (pos_nm[idx["b_head"]] * 10).mean(axis=0)
        b = (pos_nm[idx["b_tail"]] * 10).mean(axis=0)
        cv1 = float(np.linalg.norm(a - b))
    return {"Rg_A": round(rg, 2), "long_axis_extent_A": round(extent, 2),
            "cv1_beta_head_tail_A": round(cv1, 2) if cv1 is not None else None}


def get_pos_nm(sim):
    return np.array(sim.context.getState(getPositions=True).getPositions()
                    .value_in_unit(unit.nanometer))


def energy(sim):
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="inputs/1jv2.pdb")   # bent source
    ap.add_argument("--out-pdb", default="stage2c/extended_relaxed.pdb")
    ap.add_argument("--out-json", default="stage2c/stage2c_result.json")
    ap.add_argument("--frames-dir", default="stage2c/frames")
    ap.add_argument("--step-deg", type=float, default=9.0)
    ap.add_argument("--stop-deg", type=float, default=8.0, help="stop when remaining < this")
    ap.add_argument("--max-steps", type=int, default=40)
    ap.add_argument("--min-iter", type=int, default=500)
    ap.add_argument("--threads", default="48")
    ap.add_argument("--equil-ps", type=float, default=10.0)
    args = ap.parse_args()

    os.makedirs(args.frames_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)
    t0 = time.time()

    print(f"[stage2c] loading bent {args.pdb}", flush=True)
    fixer = PDBFixer(filename=args.pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    top = fixer.topology
    idx = build_index_sets(top)
    print(f"[stage2c] atoms={top.getNumAtoms()}  lower-body atoms={len(idx['lower_idx'])}",
          flush=True)

    # VACUUM force field (no implicit solvent): GB Born-radii evaluation is the
    # expensive, poorly-threaded part on CPU (~1 s/eval for 22k atoms) and makes a
    # solvated minimize take many minutes per step. For pure clash relief during
    # the incremental swing, vacuum amber14 is ~5-10x faster per iteration and
    # relieves steric overlaps just as well -- the extended shape is held by the
    # rigid-body rotations, and the seed gets a full solvated GPU relax later.
    ff = app.ForceField("amber14/protein.ff14SB.xml")
    system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.2 * unit.nanometer, constraints=None)
    integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                        1.0 * unit.femtoseconds)
    plat = mm.Platform.getPlatformByName("CPU")
    sim = app.Simulation(top, system, integ, plat, {"Threads": args.threads})
    sim.context.setPositions(fixer.positions)

    print("[stage2c] initial minimize of bent...", flush=True)
    sim.minimizeEnergy(maxIterations=args.min_iter)
    m0 = metrics(get_pos_nm(sim), idx)
    e0 = energy(sim)
    print(f"[stage2c] bent relaxed: E={e0:.3e}  {m0}", flush=True)

    traj = [{"step": 0, "cum_deg": 0.0, "E_kJ_mol": round(e0, 1), **m0}]
    cum_deg = 0.0
    for step in range(1, args.max_steps + 1):
        pos = get_pos_nm(sim)
        hinge = pos[idx["genu_ca"]].mean(axis=0)
        v_up = pos[idx["upper_ca"]].mean(axis=0) - hinge
        v_low = pos[idx["lower_ca"]].mean(axis=0) - hinge
        u_up = v_up / np.linalg.norm(v_up)
        u_low = v_low / np.linalg.norm(v_low)
        target = -u_up
        ang_to_target = float(np.degrees(np.arccos(np.clip(np.dot(u_low, target), -1, 1))))
        if ang_to_target < args.stop_deg:
            print(f"[stage2c] reached extended (remaining {ang_to_target:.1f} deg)", flush=True)
            break
        d = min(args.step_deg, ang_to_target)
        axis = np.cross(u_low, target)
        if np.linalg.norm(axis) < 1e-6:
            axis = np.array([0.0, 0.0, 1.0])
        R = rot_matrix(axis, np.radians(d))
        newpos = pos.copy()
        newpos[idx["lower_idx"]] = hinge + (pos[idx["lower_idx"]] - hinge) @ R.T
        sim.context.setPositions(newpos * unit.nanometer)
        sim.minimizeEnergy(maxIterations=args.min_iter)
        cum_deg += d
        e = energy(sim)
        m = metrics(get_pos_nm(sim), idx)
        traj.append({"step": step, "cum_deg": round(cum_deg, 1), "E_kJ_mol": round(e, 1), **m})
        print(f"[stage2c] step {step:2d}  swung {cum_deg:5.1f} deg (remaining {ang_to_target-d:5.1f})"
              f"  E={e:.3e}  Rg={m['Rg_A']}  extent={m['long_axis_extent_A']}", flush=True)
        with open(os.path.join(args.frames_dir, f"frame_{step:02d}.pdb"), "w") as fh:
            app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)
        with open(args.out_json, "w") as fh:
            json.dump({"stage": "2c_morph_extend", "status": "morphing",
                       "trajectory": traj, "seconds": round(time.time() - t0, 1)}, fh, indent=2)

    # final structure
    with open(args.out_pdb, "w") as fh:
        app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)
    m_final = metrics(get_pos_nm(sim), idx)
    e_final = energy(sim)
    print(f"[stage2c] MORPH_DONE  E={e_final:.3e}  {m_final}", flush=True)

    result = {"stage": "2c_morph_extend", "status": "done",
              "source_pdb": args.pdb, "step_deg": args.step_deg,
              "n_steps": len(traj) - 1, "total_swing_deg": round(cum_deg, 1),
              "E_bent_kJ_mol": round(e0, 1), "E_extended_kJ_mol": round(e_final, 1),
              "metrics_bent": m0, "metrics_extended": m_final,
              "trajectory": traj, "out_pdb": args.out_pdb,
              "seconds": round(time.time() - t0, 1)}

    if args.equil_ps > 0 and e_final < 1e6:
        print(f"[stage2c] equil {args.equil_ps} ps (HBonds, 2 fs)...", flush=True)
        st = sim.context.getState(getPositions=True)
        eq_sys = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.2 * unit.nanometer, constraints=app.HBonds)
        eq_int = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                             2.0 * unit.femtoseconds)
        sim = app.Simulation(top, eq_sys, eq_int, plat, {"Threads": args.threads})
        sim.context.setPositions(st.getPositions())
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
        n = int(args.equil_ps * 1000 / 2.0)
        sim.step(n)
        with open(args.out_pdb, "w") as fh:
            app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)
        m_eq = metrics(get_pos_nm(sim), idx)
        result["metrics_equil"] = m_eq
        result["E_equil_kJ_mol"] = round(energy(sim), 1)
        print(f"[stage2c] EQUIL_DONE  {m_eq}", flush=True)

    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k != "trajectory"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
