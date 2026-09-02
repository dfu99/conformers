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

FORCE MODE (obj-083). At F = 0 the above measures the wrong ensemble: integrins are
force-sensitive machines, so bent IS the ground state and extension is uphill without
tension -- obj-082 watched WT collapse alongside every mutant because there was no
barrier to cross, not because sampling was thin. --ramp-pn / --force-pn add a constant
pulling force along the physiological axis (ligand at the headpiece <-> cytoskeleton at
the C-terminal leg ends) so the knee is held against a load that can actually be varied:

  # RAMP DOWN from a holding load, after equilibrating under it -- the production protocol
  python snapback_md.py --mutations WT --ns 8 --equil-ns 2 \
                        --ramp-from-pn 60 --ramp-pn 0 --platform CUDA --tag wt_ramp
  # constant-load rung of a ladder (kappa and P(extended|F) come from these, not from a ramp)
  python snapback_md.py --mutations WT --ns 20 --equil-ns 2 --force-pn 20 --tag wt_f20

Ramp DOWN, not up. The run starts from the extended state, so an up-ramp from 0 spends its
low-force half in exactly the conditions that produced the obj-082 null: the knee collapses
on its own before the load arrives, and F½ then measures the collapse clock rather than a
force (the zero-force replicates lose half their knee drop at 3.8 +/- 2.4 ns, which a 3 pN/ns
up-ramp would report as "F½ = 12 +/- 7 pN" with no force at all). Starting held at 60 pN and
releasing measures the force at which holding stops, which is what F½ means.

Readouts (analyze_force_ramp.py): F½ from P(extended | F), and κ = k_BT/var(knee) -- κ only
from --force-pn runs, since under a ramp var(knee) is dominated by the drift the ramp drives.

CPU smoke: --smoke  (tiny minimize + a few steps, to validate the build/observables).
Force check: --check-force  (build + prove the requested pN reaches the atoms, no MD).
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

# Knee split + genu hinge, in extended_state_b.pdb's OWN numbering (chain A 1-927, chain B
# 1-539, both contiguous over observed residues) -- NOT original 1JV2 numbering. The two agree
# for chain A up to 838 and for chain B only up to 380; the ranges below stay inside that.
#
# The B entries used to be (1,440)/(441,690)/(436,444), copied from the 1JV2-numbered scripts
# (stage2_build_extended.py, stage2c_morph_extend.py, define_extension_cvs.py -- all still on
# 1JV2 numbering, correct for their own inputs). Against the renumbered file those matched real
# residues and raised nothing, but put the hinge in the wrong place: B:436-444 = beta3 mature
# 587-595 sits 61.3 A from the alphaV genu, at fraction 0.841 along the 153 A head->foot axis vs
# the real knee at 0.441, and B:381-440 (mature 532-591, lower leg, fraction 0.671) counted as
# UPPER. A 30 deg genu bend registered as 15.3 deg -- 1.8x compressed, var(knee) 3.3x too small.
#
# beta3 mature 435-531 -- which contains the beta3 genu -- is absent from this model, so there is
# no beta3 hinge to select: GENU is alphaV-only. Chain B's numbering break at 380/381 IS that
# gap, which is why it is the right place to split UPPER from LOWER.
UPPER = {"A": (1, 592), "B": (1, 380)}
LOWER = {"A": (593, 956), "B": (381, 539)}
GENU = {"A": (588, 596)}
# Expected selection sizes on extended_state_b.pdb. A range that silently matches the wrong
# residues is this repo's most repeated bug; these turn it into an abort instead of a null.
EXPECT_CA = {"UPPER": 972, "LOWER": 494, "GENU": 9}
# key cross-knee salt-bridge atom pairs to monitor  (anion res/atoms, cation res/atoms)
SALT_PAIRS = [
    ("E598-K459", ("A", 598, ["OE1", "OE2"]), ("A", 459, ["NZ"])),
    ("D457-K688", ("A", 457, ["OD1", "OD2"]), ("A", 688, ["NZ"])),
    ("K459-E636", ("A", 636, ["OE1", "OE2"]), ("A", 459, ["NZ"])),
    ("D595-K688", ("A", 595, ["OD1", "OD2"]), ("A", 688, ["NZ"])),
]

# ── constant-force pulling ──────────────────────────────────────────────────────
# pN -> kJ/(mol*nm). Identical to the fixed domain_steering.py; duplicated rather than
# imported because CLAUDE.md §6.3 keeps pipeline roots decoupled. Do NOT "correct" this
# to 6.022e-4 -- that value is low by exactly 1000x and was the bug that silently nulled
# every force run this repo has ever done.
PN_TO_KJ_PER_MOL_NM = 0.602214076
# ...and check it against openmm's own unit system rather than trusting the literal. verify_pull
# multiplies by this constant and divides by it again, so it reports the requested pN back to
# itself for ANY value -- including the historical 1000x-low one. This assert is the part that
# actually catches that bug; it costs one import-time unit conversion.
assert abs(PN_TO_KJ_PER_MOL_NM - (1.0 * unit.piconewton * unit.AVOGADRO_CONSTANT_NA).value_in_unit(
    unit.kilojoule_per_mole / unit.nanometer)) < 1e-9, "pN conversion constant is wrong"
PULL_FORCE_GROUP = 1  # own force group so its forces can be read back in isolation

# Anchors for the physiological force axis: a ligand pulls on the headpiece while the
# cytoskeleton holds the C-terminal leg ends (where the TM helices continue).
#
# extended_state_b.pdb is renumbered contiguously over OBSERVED residues, so resSeq is
# only trustworthy in the first segment of each chain (verified against P06756 / P05106):
#   chain A  resSeq 1-838   -> αV mature 1-838    (identity; αV 839-867 absent)
#            resSeq 839-927 -> αV mature 868-956  (+29)
#   chain B  resSeq 1-380   -> β3 mature 55-434   (+54; β3 435-531 absent)
#            resSeq 381-539 -> β3 mature 532-690  (+151)
# HEAD stays inside the first segment of each chain, so its resSeq ranges are safe:
# A:(1,438) = αV β-propeller, B:(1,298) = β3 βI + hybrid. Together, the ligand-binding
# headpiece. FOOT is taken POSITIONALLY as the last FOOT_NRES residues of each chain
# rather than by resSeq, because both chains are past their offset break at the C
# terminus -- "last N of the chain" is exact where a resSeq range would silently pick up
# the wrong residues. It resolves to αV mature 927-956 + β3 mature 661-690, i.e. the
# membrane-proximal ends where the TM helices would continue.
HEAD = {"A": (1, 438), "B": (1, 298)}
FOOT_NRES = 30
MIN_AXIS_A = 100.0  # head->foot separation expected in the extended state (~150 Å)


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
    return ch in tbl and tbl[ch][0] <= r <= tbl[ch][1]  # GENU is alphaV-only: chain B must miss


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
    for name, sel in (("UPPER", upper), ("LOWER", lower), ("GENU", genu)):
        if len(sel) != EXPECT_CA[name]:
            raise SystemExit(f"{name} selected {len(sel)} CA, expected {EXPECT_CA[name]} — "
                             f"the input is not numbered like extended_state_b.pdb, so knee_deg "
                             f"would be measured about the wrong hinge. Refusing to run.")
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


def make_pull_groups(topology):
    """CA indices of the two pull anchors: the headpiece and the C-terminal leg ends."""
    cai = ca_index(topology)
    head = sorted(i for (ch, r), i in cai.items()
                  if ch in HEAD and HEAD[ch][0] <= r <= HEAD[ch][1])
    by_chain = {}
    for (ch, r), i in cai.items():
        by_chain.setdefault(ch, []).append((r, i))
    foot = sorted(i for ch in by_chain for _, i in sorted(by_chain[ch])[-FOOT_NRES:])
    if not head or not foot:
        raise SystemExit(f"pull anchors empty: {len(head)} head CA, {len(foot)} foot CA")
    return head, foot


def masses_of(system):
    return np.array([system.getParticleMass(i).value_in_unit(unit.dalton)
                     for i in range(system.getNumParticles())])


def com(pos_nm, masses, idx):
    m = masses[idx]
    return (pos_nm[idx] * m[:, None]).sum(0) / m.sum()


def check_axis(pos_nm, masses, head, foot):
    """Guard against a silently wrong selection: the anchors must be at opposite ends.

    A resSeq range that matches the wrong residues (this repo's most common bug) shows up
    here as a short or nonsensical axis, whereas the extended state must span ~150 Å.
    """
    d = float(np.linalg.norm(com(pos_nm, masses, foot) - com(pos_nm, masses, head))) * 10.0
    if not np.isfinite(d) or d < MIN_AXIS_A:
        raise SystemExit(f"pull axis is {d:.1f} Å (expected >{MIN_AXIS_A:.0f} Å) — "
                         f"head/foot selection is wrong, refusing to run")
    return d


def add_pull_force(system, head, foot):
    """Constant force `f_pull` pulling the head and foot anchors apart.

    E = -f*d gives F = -dE/dd = +f along the head->foot axis regardless of separation:
    a true constant-force pull, not a moving-spring SMD. Equal and opposite on the two
    anchors, so there is no net force on the molecule and no COM runaway. f_pull is a
    global parameter in kJ/(mol*nm) so it can be ramped live via context.setParameter().
    """
    pull = mm.CustomCentroidBondForce(2, "-f_pull*distance(g1,g2)")
    pull.addGlobalParameter("f_pull", 0.0)
    pull.addGroup(head)   # g1
    pull.addGroup(foot)   # g2
    pull.addBond([0, 1], [])
    pull.setForceGroup(PULL_FORCE_GROUP)
    system.addForce(pull)
    return pull


def verify_pull(context, masses, head, foot, probe_pn=20.0, tol=0.01):
    """Prove the requested pN actually reaches the atoms. Aborts the run if it does not.

    Reads back ONLY the pull force group, sums its per-atom forces over each anchor and
    projects onto the head->foot axis. For E = -f*d the analytic answer is exactly -f on
    the head, +f on the foot, and zero net -- so any other number means the request never
    landed on the atoms. Probed at two different forces, because a single probe cannot
    tell a working force from one that is hardcoded or ignored.

    Every force run this repo has done so far was a silent null (a 1000x unit error, a
    lab-frame restraint cage, a residue-selection fallback that never fired). This is the
    gate that stops the next one.
    """
    report = {}
    for pn in (probe_pn, 2.0 * probe_pn):
        context.setParameter("f_pull", pn * PN_TO_KJ_PER_MOL_NM)
        st = context.getState(getPositions=True, getForces=True, groups={PULL_FORCE_GROUP})
        p = np.array(st.getPositions().value_in_unit(unit.nanometer))
        f = np.array(st.getForces().value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer))
        axis = com(p, masses, foot) - com(p, masses, head)
        axis /= np.linalg.norm(axis)
        to_pn = 1.0 / PN_TO_KJ_PER_MOL_NM
        got_head = float(f[head].sum(0) @ axis) * to_pn
        got_foot = float(f[foot].sum(0) @ axis) * to_pn
        net = float(np.linalg.norm(f.sum(0))) * to_pn
        report[f"{pn:g}pN"] = {"on_head_pN": round(got_head, 4),
                               "on_foot_pN": round(got_foot, 4),
                               "net_pN": round(net, 6)}
        bad = (abs(got_foot - pn) > tol * pn + 1e-3
               or abs(got_head + pn) > tol * pn + 1e-3
               or net > tol * pn + 1e-3)
        print(f"  force check @ {pn:g} pN: head={got_head:+.3f} pN  foot={got_foot:+.3f} pN  "
              f"net={net:.2e} pN  {'FAIL' if bad else 'ok'}", flush=True)
        if bad:
            raise SystemExit(
                f"FORCE CHECK FAILED: requested {pn:g} pN but the atoms feel "
                f"{got_foot:.3f} pN (net {net:.3g}). Refusing to run a silent null.")
    return report


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
    ap.add_argument("--build-only", action="store_true",
                    help="build the production implicit system + report start observables, no MD.")
    ap.add_argument("--ramp-pn", type=float, default=0.0,
                    help="force at the END of the ramp, in pN (start = --ramp-from-pn).")
    ap.add_argument("--ramp-from-pn", type=float, default=0.0,
                    help="force at the START of the ramp. Set it ABOVE --ramp-pn to ramp DOWN "
                         "from a holding load, which is the measurement that works here.")
    ap.add_argument("--equil-ns", type=float, default=0.0,
                    help="equilibrate this long at the starting force before recording anything.")
    ap.add_argument("--force-pn", type=float, default=0.0,
                    help="hold the pull at this many pN for the whole run (ladder mode).")
    ap.add_argument("--probe-pn", type=float, default=20.0,
                    help="force used by the pre-run check that the pN reaches the atoms.")
    ap.add_argument("--check-force", action="store_true",
                    help="build + verify the pull force reaches the atoms, then exit. No MD.")
    args = ap.parse_args()
    ramping = bool(args.ramp_pn or args.ramp_from_pn)
    if ramping and args.force_pn:
        ap.error("--force-pn is a constant load; it cannot be combined with a ramp")
    pulling = bool(ramping or args.force_pn or args.check_force)

    outdir = os.path.join(args.out_dir, args.tag)
    os.makedirs(outdir, exist_ok=True)

    implicit = not (args.vacuum or args.smoke)  # smoke validates in vacuum for speed
    topo, pos, system = build(args.pdb, args.mutations, outdir, implicit=implicit)
    print(f"[{args.tag}] force field: {'implicit-solvent GB-OBC2' if implicit else 'vacuum'}", flush=True)
    if args.build_only:
        obs = make_observables(topo)
        p0 = np.array(pos.value_in_unit(unit.nanometer))
        knee0, rg0, sd0 = obs(p0)
        print(f"[{args.tag}] BUILD OK: {system.getNumParticles()} atoms  "
              f"knee={knee0}° Rg={rg0}Å  salt={sd0}", flush=True)
        json.dump({"tag": args.tag, "mutations": args.mutations, "atoms": system.getNumParticles(),
                   "start": {"knee": knee0, "Rg": rg0, "salt": sd0}},
                  open(os.path.join(outdir, "build_check.json"), "w"), indent=2)
        return
    head = foot = masses = None
    if pulling:
        head, foot = make_pull_groups(topo)
        masses = masses_of(system)
        axis_a = check_axis(np.array(pos.value_in_unit(unit.nanometer)), masses, head, foot)
        add_pull_force(system, head, foot)
        print(f"[{args.tag}] pull: head {len(head)} CA <-> foot {len(foot)} CA, "
              f"axis {axis_a:.1f} Å", flush=True)

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

    check = None
    if pulling:
        check = verify_pull(sim.context, masses, head, foot, probe_pn=args.probe_pn)
        sim.context.setParameter("f_pull", 0.0)  # minimize and equilibrate unloaded
        if args.check_force:
            json.dump({"tag": args.tag, "mutations": args.mutations,
                       "head_CA": len(head), "foot_CA": len(foot),
                       "axis_A": round(axis_a, 1), "check": check},
                      open(os.path.join(outdir, "force_check.json"), "w"), indent=2)
            print(f"[{args.tag}] FORCE CHECK PASSED -> {outdir}/force_check.json", flush=True)
            return

    # 50 iterations left rmsF ~216 kJ/mol/nm on this morph and NaN'd ~half of all --smoke runs;
    # 500 gets to -6.2e4 kJ/mol (2000 reaches -6.3e4), which is enough for the pre-flight to be
    # a gate rather than a coin flip.
    min_iter = 500 if args.smoke else 2000
    sim.minimizeEnergy(maxIterations=min_iter)
    sim.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)

    steps_total = 250 if args.smoke else int(args.ns * 1000 / 0.002)
    # report_ps is already picoseconds -- do NOT apply the ns->ps *1000 factor used
    # by steps_total. With it, the chunk came out 5x longer than the whole run, so
    # min(chunk, remaining) collapsed to a single chunk and the trajectory got ONE
    # row at the final frame instead of ns/report_ps rows. The snap-back readout is
    # a time course (does the mutant knee fall while WT holds?), so a lone endpoint
    # sample cannot answer it. --smoke hardcodes its chunk, which is why the
    # pre-flight smoke test never exercised this expression.
    steps_chunk = 50 if args.smoke else int(args.report_ps / 0.002)

    # The load is piecewise-constant over a report chunk, held at the value the ramp reaches at
    # the chunk's MIDPOINT -- so the f_pN logged with a row is the mean force that acted during
    # the interval that produced it, and the ramp spans its full requested range.
    def force_pn_at(frac):
        if not ramping:
            return args.force_pn
        return args.ramp_from_pn + (args.ramp_pn - args.ramp_from_pn) * frac

    # Equilibrate UNDER the starting load before recording. Without this the run measures the
    # zero-force collapse clock and calls it a force: obj-082's 9 replicates lose half their
    # total knee drop at 3.8 +/- 2.4 ns, which on a 0 -> 60 pN / 20 ns up-ramp lands at
    # "F-half = 12 +/- 7 pN" with no force applied at all -- a believable wrong answer sitting
    # right on the literature integrin extension force.
    if pulling and args.equil_ns > 0:
        sim.context.setParameter("f_pull", force_pn_at(0.0) * PN_TO_KJ_PER_MOL_NM)
        sim.step(int(args.equil_ns * 1000 / 0.002))
        p = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.nanometer))
        print(f"[{args.tag}] equilibrated {args.equil_ns} ns at {force_pn_at(0.0):.1f} pN  "
              f"knee={obs(p)[0]:.1f}°", flush=True)

    # Frames, not just the CSV. The observables are computed from module-level residue ranges,
    # and one of those ranges was wrong for this input until 2026-08-17 -- with frames on disk a
    # bad observable is a re-analysis, without them it is 13 GPU-hours re-rented.
    sim.reporters.append(app.DCDReporter(os.path.join(outdir, "frames.dcd"), steps_chunk))

    csv_path = os.path.join(outdir, "trajectory.csv")
    f_pn = 0.0
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        cols = ["ps", "knee_deg", "Rg_A"] + (["f_pN", "head_foot_A"] if pulling else [])
        w.writerow(cols + [p[0] for p in SALT_PAIRS if p[0] in sd0])
        done = 0
        while done < steps_total:
            n = min(steps_chunk, steps_total - done)
            if pulling:
                f_pn = force_pn_at((done + n / 2) / steps_total)
                sim.context.setParameter("f_pull", f_pn * PN_TO_KJ_PER_MOL_NM)
            sim.step(n)
            done += n
            st = sim.context.getState(getPositions=True, getEnergy=True)
            p = np.array(st.getPositions().value_in_unit(unit.nanometer))
            knee, rg, sd = obs(p)
            row = [round(done * 0.002, 1), knee, rg]
            if pulling:
                hf = float(np.linalg.norm(com(p, masses, foot) - com(p, masses, head))) * 10.0
                row += [round(f_pn, 3), round(hf, 2)]
            row += [sd.get(k, "") for k in sd0]
            w.writerow(row); fh.flush()
            load = f"F={f_pn:5.1f}pN  " if pulling else ""
            print(f"[{args.tag}] {done*0.002:8.1f} ps  {load}knee={knee:6.1f}°  Rg={rg:5.1f}Å  {sd}", flush=True)
    sim.saveState(os.path.join(outdir, "final.xml"))
    with open(os.path.join(outdir, "final.pdb"), "w") as fh:
        app.PDBFile.writeFile(topo, sim.context.getState(getPositions=True).getPositions(), fh)
    summary = {"tag": args.tag, "mutations": args.mutations, "ns": args.ns,
               "start": {"knee": knee0, "Rg": rg0, "salt": sd0},
               "end": {"knee": knee, "Rg": rg, "salt": sd}}
    if pulling:
        summary["pull"] = {"mode": "ramp" if ramping else "constant",
                           "ramp_from_pn": args.ramp_from_pn, "ramp_pn": args.ramp_pn,
                           "force_pn": args.force_pn, "equil_ns": args.equil_ns,
                           "final_pn": round(f_pn, 3), "head_CA": len(head),
                           "foot_CA": len(foot), "axis_A": round(axis_a, 1),
                           "force_check": check}
    json.dump(summary, open(os.path.join(outdir, "summary.json"), "w"), indent=2)
    print(f"[{args.tag}] DONE  knee {knee0}->{knee}°  Rg {rg0}->{rg}Å  -> {csv_path}", flush=True)


if __name__ == "__main__":
    main()
