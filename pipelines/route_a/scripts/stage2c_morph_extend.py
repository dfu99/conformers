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


ION_NAMES = {"CA", "MG", "MN", "ZN", "NA", "CL", "K"}  # residue names, not atom names


def assign_ions(top, positions_nm, lower_idx):
    """Ion atoms travel with whichever body coordinates them.

    An ion's residue id (e.g. 4008) falls outside every UPPER/LOWER residue range, so it would
    otherwise be left behind while the leg rotates away from it -- silently tearing the ion out of
    its own coordination shell. Decide by geometry instead: whichever body holds most of the atoms
    within 3.0 A of the ion carries it.
    """
    lower = set(int(i) for i in lower_idx)
    extra = []
    for res in top.residues():
        if res.name.strip().upper() not in ION_NAMES or len(list(res.atoms())) != 1:
            continue
        ion = next(res.atoms())
        d = np.linalg.norm(positions_nm - positions_nm[ion.index], axis=1) * 10.0
        shell = [i for i in np.where(d < 3.0)[0] if i != ion.index]
        if not shell:
            raise SystemExit(f"ion {res.name}{res.id} has no coordinating atom within 3 A — "
                             f"refusing to guess which body it belongs to")
        n_lower = sum(1 for i in shell if i in lower)
        side = "lower" if n_lower * 2 > len(shell) else "upper"
        print(f"[stage2c] ion {res.name}{res.id}: {len(shell)} contacts <3 A, "
              f"{n_lower} in the lower body -> rotates with the {side} body", flush=True)
        if side == "lower":
            extra.append(ion.index)
    return extra


def set_gb_screening(system, salt_M, temperature=300.0, eps=78.5):
    """Turn on Debye screening in the GB force. It is built with kappa hard-coded to 0.

    ForceField.createSystem in openmm 8.4 rejects implicitSolventSaltConc ("argument was never
    used"), and the CustomGBForce built from implicit/obc2.xml bakes `kappa=0` into every energy
    term, i.e. zero ionic strength. Rewriting the constant in place is the only route that does not
    mean hand-rebuilding the force. kappa = sqrt(2 e^2 N_A I / (eps0 eps kB T)), in nm^-1.
    """
    e, NA, eps0, kB = 1.602176634e-19, 6.02214076e23, 8.8541878128e-12, 1.380649e-23
    kappa = float(np.sqrt(2 * e * e * NA * (salt_M * 1000.0) / (eps0 * eps * kB * temperature))) / 1e9
    gb = [f for f in system.getForces() if isinstance(f, mm.CustomGBForce)]
    if not gb:
        raise SystemExit("--salt-mM given but the system has no CustomGBForce (not implicit?)")
    patched = 0
    for force in gb:
        for i in range(force.getNumEnergyTerms()):
            expr, typ = force.getEnergyTermParameters(i)
            if "kappa=0" in expr:
                force.setEnergyTermParameters(i, expr.replace("kappa=0", f"kappa={kappa}"), typ)
                patched += 1
    if not patched:
        raise SystemExit("could not find kappa=0 to patch — GB expression format changed")
    print(f"[stage2c] Debye screening on: {salt_M*1000:.0f} mM -> kappa {kappa:.4f} /nm "
          f"(Debye length {1/kappa:.2f} nm), {patched} energy terms", flush=True)


def add_ion_restraints(system, top, pos_nm, k, shell_A=3.0, tol_A=0.2):
    """Flat-bottom restraints holding each metal ion onto the ligands it has in the input.

    Without these the genu Ca2+ simply leaves: in GB implicit solvent there is no water shell, the
    ion is over-solvated relative to explicit water, and the minimiser finds it favourable to pull
    it off its carboxylates -- measured at ~10.5 A from ANY protein atom by the first morph frame,
    while the five ions in tighter sites stayed at 2.8-3.2 A.

    Flat-bottom, not harmonic: zero force while the contact is at or inside its input length, so the
    site can still relax and tighten naturally; the penalty only resists the ion LEAVING. r0 is
    taken per contact from the input geometry, so this holds the crystal coordination rather than
    imposing an idealised one.
    """
    rest = mm.CustomBondForce("0.5*k*step(r - r0)*(r - r0)^2")
    rest.addPerBondParameter("r0")
    rest.addGlobalParameter("k", k)
    ions = [a for res in top.residues() if res.name.strip().upper() in ION_NAMES
            and len(list(res.atoms())) == 1 for a in res.atoms()]
    prot = [a for a in top.atoms() if a.residue.name.strip().upper() not in ION_NAMES]
    n_total = 0
    for ion in ions:
        d = np.linalg.norm(pos_nm[[a.index for a in prot]] - pos_nm[ion.index], axis=1) * 10.0
        shell = np.where(d < shell_A)[0]
        for j in shell:
            rest.addBond(ion.index, prot[j].index, [(d[j] + tol_A) / 10.0])
        n_total += len(shell)
        if len(shell) == 0:
            raise SystemExit(f"ion {ion.residue.name}{ion.residue.id} has no ligand within "
                             f"{shell_A} A in the input — nothing to restrain it to")
    system.addForce(rest)
    print(f"[stage2c] ion restraints: {n_total} contacts over {len(ions)} ions, "
          f"k={k:g} kJ/mol/nm^2, flat-bottom +{tol_A} A", flush=True)
    return n_total


def check_ions_home(top, pos_nm, shell_A=3.5, min_ligands=2):
    """Verify every ion is still in a coordination shell. A silent escape invalidates the run."""
    prot = [a for a in top.atoms() if a.residue.name.strip().upper() not in ION_NAMES]
    P = pos_nm[[a.index for a in prot]]
    bad = []
    for res in top.residues():
        if res.name.strip().upper() not in ION_NAMES or len(list(res.atoms())) != 1:
            continue
        ion = next(res.atoms())
        d = np.linalg.norm(P - pos_nm[ion.index], axis=1) * 10.0
        n = int((d < shell_A).sum())
        print(f"[stage2c] ion {res.name}{res.id}: {n} ligands < {shell_A} A "
              f"(closest {d.min():.2f} A)", flush=True)
        if n < min_ligands:
            bad.append(f"{res.name}{res.id} ({n} ligands, closest {d.min():.2f} A)")
    if bad:
        raise SystemExit("ION(S) ESCAPED during the morph: " + "; ".join(bad) +
                         " — the seed is not usable; raise --ion-restraint-k and re-run.")


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
    ap.add_argument("--keep-ions", action="store_true",
                    help="keep structural metal ions (PDBFixer's removeHeterogens deletes them). "
                         "Sugars and waters are still dropped.")
    ap.add_argument("--implicit", action="store_true",
                    help="relax in GB-OBC2 implicit solvent instead of vacuum. Vacuum has "
                         "dielectric 1 and no counter-ions, which drives adjacent carboxylates "
                         "apart -- that is what collapsed the genu Ca site in the 2026-07 seed.")
    ap.add_argument("--salt-mM", type=float, default=0.0,
                    help="ionic strength for Debye screening in implicit solvent (e.g. 150).")
    ap.add_argument("--platform", default="CPU", help="CPU or CUDA. GB on CPU is slow.")
    ap.add_argument("--ion-restraint-k", type=float, default=0.0,
                    help="flat-bottom restraint (kJ/mol/nm^2) holding each ion on the ligands it "
                         "has in the input. Without it the genu Ca2+ leaves during minimisation.")
    args = ap.parse_args()

    os.makedirs(args.frames_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)
    t0 = time.time()

    print(f"[stage2c] loading bent {args.pdb}", flush=True)
    src = args.pdb
    if args.keep_ions:
        # removeHeterogens() is all-or-nothing, so filter by residue name up front instead and
        # skip it below: keep protein + ions, drop sugars/waters.
        src = os.path.join(os.path.dirname(args.out_json) or ".", "input_with_ions.pdb")
        kept = 0
        with open(args.pdb) as fh, open(src, "w") as out:
            for line in fh:
                if line.startswith("HETATM"):
                    if line[17:20].strip().upper() not in ION_NAMES:
                        continue
                    kept += 1
                out.write(line)
        print(f"[stage2c] kept {kept} ion atoms, dropped other heterogens -> {src}", flush=True)
    fixer = PDBFixer(filename=src)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    if not args.keep_ions:
        fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    top = fixer.topology
    idx = build_index_sets(top)
    if args.keep_ions:
        pos0_nm = np.array(fixer.positions.value_in_unit(unit.nanometer))
        extra = assign_ions(top, pos0_nm, idx["lower_idx"])
        if extra:
            idx["lower_idx"] = np.concatenate([idx["lower_idx"], np.array(extra, dtype=int)])
    print(f"[stage2c] atoms={top.getNumAtoms()}  lower-body atoms={len(idx['lower_idx'])}",
          flush=True)

    # VACUUM force field (no implicit solvent): GB Born-radii evaluation is the
    # expensive, poorly-threaded part on CPU (~1 s/eval for 22k atoms) and makes a
    # solvated minimize take many minutes per step. For pure clash relief during
    # the incremental swing, vacuum amber14 is ~5-10x faster per iteration and
    # relieves steric overlaps just as well -- the extended shape is held by the
    # rigid-body rotations, and the seed gets a full solvated GPU relax later.
    if args.implicit:
        # amber14/tip3p.xml is loaded ONLY for its ion templates -- there is no water here.
        ff = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml", "implicit/obc2.xml")
    else:
        ff = app.ForceField("amber14/protein.ff14SB.xml")
    system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.2 * unit.nanometer, constraints=None)
    if args.implicit and args.salt_mM > 0:
        set_gb_screening(system, args.salt_mM / 1000.0)
    if args.ion_restraint_k > 0:
        add_ion_restraints(system, top, np.array(fixer.positions.value_in_unit(unit.nanometer)),
                           args.ion_restraint_k)
    integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                        1.0 * unit.femtoseconds)
    plat = mm.Platform.getPlatformByName(args.platform)
    props = {"Threads": args.threads} if args.platform == "CPU" else {"Precision": "mixed"}
    sim = app.Simulation(top, system, integ, plat, props)
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
        # the equilibration system is built fresh -- every modification of the morph system has to
        # be repeated here or the last 10 ps quietly undoes it
        if args.implicit and args.salt_mM > 0:
            set_gb_screening(eq_sys, args.salt_mM / 1000.0)
        if args.ion_restraint_k > 0:
            add_ion_restraints(eq_sys, top, get_pos_nm(sim), args.ion_restraint_k)
        eq_int = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                             2.0 * unit.femtoseconds)
        sim = app.Simulation(top, eq_sys, eq_int, plat, props)
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

    if args.keep_ions:
        check_ions_home(top, get_pos_nm(sim))

    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k != "trajectory"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
