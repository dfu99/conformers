#!/usr/bin/env python3
"""Incremental genu-hinge morph to the EXTENDED state, for any integrin variant.

Generalisation of stage2c_morph_extend.py (which hardcoded alphaVbeta3 / 1JV2
residue ranges) to arbitrary heterodimers. Domain boundaries come from
map_domain_boundaries.py, which transfers the validated alphaVbeta3 split by
sequence alignment; everything else follows the alphaVbeta3 recipe exactly:

  rotate the lower legs about the genu in small steps, vacuum-minimising after
  each, so no step ever accumulates the ~1e13 kJ/mol clash that stalls a one-shot
  ~140 deg rigid swing (obj-075).

Vacuum amber14 ff14SB, CPU only -- no GPU, so this never contends with the oxDNA
FFS jobs holding the pod's A5000. Emits the extended endpoint AND every
intermediate frame, which doubles as the string-method initial path.
"""
import argparse
import json
import os
import time

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer


def rot_matrix(axis, theta):
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = axis
    return np.array([
        [c + x*x*(1-c),   x*y*(1-c) - z*s, x*z*(1-c) + y*s],
        [y*x*(1-c) + z*s, c + y*y*(1-c),   y*z*(1-c) - x*s],
        [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)],
    ])


def build_index_sets(top, b):
    """Atom indices for the movable lower body + CA sets, from transferred bounds."""
    def rng(group, cid):
        v = b.get(group, {}).get(cid)
        return (v[0], v[1]) if v else None

    sets = {k: [] for k in ("lower_idx", "genu_ca", "upper_ca", "lower_ca",
                            "all_ca", "b_head", "b_tail")}
    for atom in top.atoms():
        cid = atom.residue.chain.id
        try:
            rs = int(atom.residue.id)
        except ValueError:
            continue
        is_ca = atom.name == "CA"
        if is_ca:
            sets["all_ca"].append(atom.index)
        lo = rng("LOWER", cid)
        if lo and lo[0] <= rs <= lo[1]:
            sets["lower_idx"].append(atom.index)
            if is_ca:
                sets["lower_ca"].append(atom.index)
        up = rng("UPPER", cid)
        if up and is_ca and up[0] <= rs <= up[1]:
            sets["upper_ca"].append(atom.index)
        gn = rng("GENU", cid)
        if gn and is_ca and gn[0] <= rs <= gn[1]:
            sets["genu_ca"].append(atom.index)
        if is_ca and cid == "B":
            bh, bt = rng("B_HEAD", "B"), rng("B_TAIL", "B")
            if bh and bh[0] <= rs <= bh[1]:
                sets["b_head"].append(atom.index)
            elif bt and bt[0] <= rs <= bt[1]:
                sets["b_tail"].append(atom.index)
    return {k: np.array(v, dtype=int) for k, v in sets.items()}


def metrics(pos_nm, idx):
    ca = pos_nm[idx["all_ca"]] * 10.0
    com = ca.mean(axis=0)
    rg = float(np.sqrt(np.mean(np.sum((ca - com) ** 2, axis=1))))
    centred = ca - com
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    proj = centred @ vt[0]
    out = {"Rg_A": round(rg, 2),
           "long_axis_extent_A": round(float(proj.max() - proj.min()), 2)}
    if len(idx["b_head"]) and len(idx["b_tail"]):
        a = (pos_nm[idx["b_head"]] * 10).mean(axis=0)
        bb = (pos_nm[idx["b_tail"]] * 10).mean(axis=0)
        out["cv1_beta_head_tail_A"] = round(float(np.linalg.norm(a - bb)), 2)
    if len(idx["genu_ca"]) and len(idx["upper_ca"]) and len(idx["lower_ca"]):
        h = (pos_nm[idx["genu_ca"]] * 10).mean(axis=0)
        vu = (pos_nm[idx["upper_ca"]] * 10).mean(axis=0) - h
        vl = (pos_nm[idx["lower_ca"]] * 10).mean(axis=0) - h
        cos = np.dot(vu, vl) / (np.linalg.norm(vu) * np.linalg.norm(vl))
        out["genu_angle_deg"] = round(float(np.degrees(np.arccos(np.clip(cos, -1, 1)))), 2)
    return out


def get_pos_nm(sim):
    return np.array(sim.context.getState(getPositions=True).getPositions()
                    .value_in_unit(unit.nanometer))


def energy(sim):
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True,
                    help="key into the boundaries JSON (e.g. alphaVbeta3__1JV2)")
    ap.add_argument("--label", default="",
                    help="output directory name; defaults to --variant")
    ap.add_argument("--boundaries", default="results/variants/boundaries_all.json")
    ap.add_argument("--out-dir", default="results/variants")
    ap.add_argument("--step-deg", type=float, default=9.0)
    ap.add_argument("--stop-deg", type=float, default=8.0)
    ap.add_argument("--max-steps", type=int, default=40)
    ap.add_argument("--min-iter", type=int, default=500)
    ap.add_argument("--threads", default="16")
    ap.add_argument("--equil-ps", type=float, default=10.0)
    args = ap.parse_args()

    spec = json.load(open(args.boundaries))[args.variant]
    b = spec["boundaries"]
    pdb_path = spec["prepared_pdb"]
    label = args.label or args.variant

    vdir = os.path.join(args.out_dir, label)
    frames_dir = os.path.join(vdir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    out_pdb = os.path.join(vdir, f"{label}_extended.pdb")
    out_json = os.path.join(vdir, f"{label}_morph.json")
    t0 = time.time()

    print(f"[{label}] loading bent {pdb_path}", flush=True)
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # keep original numbering (see tasks/lessons.md)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    top = fixer.topology
    idx = build_index_sets(top, b)
    print(f"[{label}] atoms={top.getNumAtoms()} lower-body atoms={len(idx['lower_idx'])} "
          f"CA upper/lower/genu={len(idx['upper_ca'])}/{len(idx['lower_ca'])}/{len(idx['genu_ca'])}",
          flush=True)
    if len(idx["lower_idx"]) == 0 or len(idx["genu_ca"]) == 0:
        raise SystemExit(f"[{label}] ABORT: empty lower body or genu selection")

    ff = app.ForceField("amber14/protein.ff14SB.xml")
    system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.2 * unit.nanometer, constraints=None)
    integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                        1.0 * unit.femtoseconds)
    plat = mm.Platform.getPlatformByName("CPU")
    sim = app.Simulation(top, system, integ, plat, {"Threads": args.threads})
    sim.context.setPositions(fixer.positions)

    print(f"[{label}] initial minimize of bent...", flush=True)
    sim.minimizeEnergy(maxIterations=args.min_iter)
    m0, e0 = metrics(get_pos_nm(sim), idx), energy(sim)
    print(f"[{label}] bent relaxed: E={e0:.3e} {m0}", flush=True)

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
            print(f"[{label}] reached extended (remaining {ang_to_target:.1f} deg)",
                  flush=True)
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
        traj.append({"step": step, "cum_deg": round(cum_deg, 1),
                     "E_kJ_mol": round(e, 1), **m})
        print(f"[{label}] step {step:2d} swung {cum_deg:5.1f} deg "
              f"(remaining {ang_to_target-d:5.1f}) E={e:.3e} Rg={m['Rg_A']} "
              f"extent={m['long_axis_extent_A']} genu={m.get('genu_angle_deg')}", flush=True)
        with open(os.path.join(frames_dir, f"frame_{step:02d}.pdb"), "w") as fh:
            app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)
        with open(out_json, "w") as fh:
            json.dump({"variant": label, "survey_key": args.variant, "status": "morphing", "trajectory": traj,
                       "seconds": round(time.time() - t0, 1)}, fh, indent=2)

    with open(out_pdb, "w") as fh:
        app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)
    m_final, e_final = metrics(get_pos_nm(sim), idx), energy(sim)
    print(f"[{label}] MORPH_DONE E={e_final:.3e} {m_final}", flush=True)

    result = {
        "variant": label, "survey_key": args.variant, "status": "done", "pdb_id": spec["pdb_id"],
        "source_pdb": pdb_path, "step_deg": args.step_deg,
        "n_steps": len(traj) - 1, "total_swing_deg": round(cum_deg, 1),
        "E_bent_kJ_mol": round(e0, 1), "E_extended_kJ_mol": round(e_final, 1),
        "metrics_bent": m0, "metrics_extended": m_final,
        "trajectory": traj, "out_pdb": out_pdb, "frames_dir": frames_dir,
        "seconds": round(time.time() - t0, 1),
    }

    if args.equil_ps > 0 and e_final < 1e6:
        print(f"[{label}] equil {args.equil_ps} ps...", flush=True)
        st = sim.context.getState(getPositions=True)
        eq_sys = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.2 * unit.nanometer,
                                 constraints=app.HBonds)
        eq_int = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                             2.0 * unit.femtoseconds)
        sim = app.Simulation(top, eq_sys, eq_int, plat, {"Threads": args.threads})
        sim.context.setPositions(st.getPositions())
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
        sim.step(int(args.equil_ps * 1000 / 2.0))
        with open(out_pdb, "w") as fh:
            app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)
        result["metrics_equil"] = metrics(get_pos_nm(sim), idx)
        result["E_equil_kJ_mol"] = round(energy(sim), 1)
        print(f"[{label}] EQUIL_DONE {result['metrics_equil']}", flush=True)

    result["seconds"] = round(time.time() - t0, 1)
    with open(out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"[{label}] wrote {out_json}", flush=True)


if __name__ == "__main__":
    main()
