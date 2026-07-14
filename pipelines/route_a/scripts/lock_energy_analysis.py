#!/usr/bin/env python3
"""How STRONG is each genu salt-bridge lock, and how much lock energy does each mutation remove?

The snap-back MD (pending an A5000 slot) answers the *dynamical* question — does the extended
state fall back to bent when a lock is knocked out. This is its *static / thermodynamic* companion,
and it is fully CPU-only: no dynamics, no minimization, no GPU. We only need per-atom force-field
parameters, which come straight off the built System's NonbondedForce (no Context, no platform).

For the WT extended state we build the exact production system snap-back uses
(amber14/ff14SB + GB-OBC2, via snapback_md.build) and read the real ff14SB partial charges + LJ
parameters. For each of the four cross-knee salt bridges we sum the direct nonbonded interaction
energy (Coulomb + Lennard-Jones) between the two participating residues — a local, per-pair quantity
that does not depend on total atom count, so it is comparable across bridges and across mutants.

Then, using the knockout matrix already validated in obj-079 (which bridges each mutation abolishes),
the destabilization proxy for a mutation is simply the sum of the interaction energies of the bridges
it removes. That converts "K459 is the hub" (obj-080) into a quantitative, falsifiable prediction:
K459A strips two locks, E598A only one, so K459A should destabilize the extended state more — exactly
what the pending snap-back MD will test dynamically.

Caveat, stated plainly: the reported energy is the DIRECT (unscreened, eps=1) ff14SB nonbonded
interaction, an upper bound on lock strength. GB solvent screens the absolute magnitude (~2-4x); a
distance-dependent-dielectric estimate (eps=4r) is reported alongside to show the RANKING survives
screening. The ranking — driven by which and how many bridges break — is the falsifiable part.

Run under the isolated CUDA-OpenMM env (CPU here, no GPU touched):
  /home/dan/anaconda3/envs/snapback/bin/python pipelines/route_a/scripts/lock_energy_analysis.py
"""
import json
import os
import sys
import numpy as np
from openmm import NonbondedForce, unit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snapback_md as sb  # build(), strip_h()

KE = 138.935458  # kJ*nm / (mol*e^2)  -- Coulomb constant in OpenMM units

# bridge: (name, anion (chain,resid), cation (chain,resid), charged-atom pair for the min distance)
BRIDGES = [
    ("E598-K459", ("A", 598), ("A", 459), (["OE1", "OE2"], ["NZ"])),
    ("K459-E636", ("A", 636), ("A", 459), (["OE1", "OE2"], ["NZ"])),
    ("D457-K688", ("A", 457), ("A", 688), (["OD1", "OD2"], ["NZ"])),
    ("D595-K688", ("A", 595), ("A", 688), (["OD1", "OD2"], ["NZ"])),
]

# which bridges each mutation abolishes (validated in obj-079 build_check.json)
KNOCKOUT = {
    "E598A": ["E598-K459"],
    "K459A": ["E598-K459", "K459-E636"],
    "K459A+E598A": ["E598-K459", "K459-E636"],  # E598's only bridge already gone via K459A
}


def residue_atoms(topology):
    """(chain, resid) -> {atomname: index} for every residue."""
    out = {}
    for atom in topology.atoms():
        key = (atom.residue.chain.id, int(atom.residue.id))
        out.setdefault(key, {})[atom.name] = atom.index
    return out


def pair_energy(idx_a, idx_b, q, sig, eps, pos_nm, ddd=False):
    """Direct nonbonded interaction (Coulomb + LJ) between two disjoint atom sets, kJ/mol.

    ddd=False -> eps_r = 1 (direct). ddd=True -> distance-dependent dielectric eps_r = 4*r[Angstrom].
    """
    ia = np.array(idx_a); ib = np.array(idx_b)
    ra = pos_nm[ia]; rb = pos_nm[ib]
    d = np.linalg.norm(ra[:, None, :] - rb[None, :, :], axis=2)  # nm, (na,nb)
    d = np.maximum(d, 1e-4)
    qa = q[ia][:, None]; qb = q[ib][None, :]
    eps_r = (4.0 * d * 10.0) if ddd else 1.0  # 4*r in Angstrom
    coul = KE * (qa * qb) / (eps_r * d)
    sij = 0.5 * (sig[ia][:, None] + sig[ib][None, :])
    eij = np.sqrt(eps[ia][:, None] * eps[ib][None, :])
    sr6 = (sij / d) ** 6
    lj = 4.0 * eij * (sr6 * sr6 - sr6)
    return float(coul.sum()), float(lj.sum())


def min_charged_dist(ratoms_a, ratoms_b, anion_names, cation_names, pos_nm):
    an = [ratoms_a[n] for n in anion_names if n in ratoms_a]
    cn = [ratoms_b[n] for n in cation_names if n in ratoms_b]
    if not an or not cn:
        return None
    da = pos_nm[np.array(an)]; dc = pos_nm[np.array(cn)]
    return float(np.linalg.norm(da[:, None, :] - dc[None, :, :], axis=2).min() * 10.0)


def main():
    pdb = "results/route_a/extended_state_b.pdb"
    workdir = os.environ.get("LOCKE_WORKDIR", "/tmp/claude-1000/lock_energy")
    os.makedirs(workdir, exist_ok=True)

    print("building WT extended system (amber14/ff14SB + GB-OBC2) ...", flush=True)
    topo, positions, system = sb.build(pdb, "WT", workdir, implicit=True)
    pos_nm = np.array(positions.value_in_unit(unit.nanometer))

    nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    n = system.getNumParticles()
    q = np.zeros(n); sig = np.zeros(n); eps = np.zeros(n)
    for i in range(n):
        c, s, e = nb.getParticleParameters(i)
        q[i] = c.value_in_unit(unit.elementary_charge)
        sig[i] = s.value_in_unit(unit.nanometer)
        eps[i] = e.value_in_unit(unit.kilojoule_per_mole)

    ratoms = residue_atoms(topo)

    per_bridge = {}
    for name, (ach, ar), (cch, cr), (anion_names, cation_names) in BRIDGES:
        a = ratoms[(ach, ar)]; b = ratoms[(cch, cr)]
        qa = sum(q[i] for i in a.values()); qb = sum(q[i] for i in b.values())
        coul, lj = pair_energy(list(a.values()), list(b.values()), q, sig, eps, pos_nm, ddd=False)
        coul_ddd, lj_ddd = pair_energy(list(a.values()), list(b.values()), q, sig, eps, pos_nm, ddd=True)
        dmin = min_charged_dist(a, b, anion_names, cation_names, pos_nm)
        per_bridge[name] = {
            "anion": f"{ach}{ar}", "cation": f"{cch}{cr}",
            "net_charge_anion": round(qa, 3), "net_charge_cation": round(qb, 3),
            "min_charged_dist_A": round(dmin, 2) if dmin else None,
            "E_coulomb_direct_kJ": round(coul, 1),
            "E_lj_kJ": round(lj, 1),
            "E_total_direct_kJ": round(coul + lj, 1),
            "E_total_direct_kcal": round((coul + lj) / 4.184, 2),
            "E_total_ddd_kJ": round(coul_ddd + lj_ddd, 1),
            "E_total_ddd_kcal": round((coul_ddd + lj_ddd) / 4.184, 2),
        }
        print(f"  {name}: d={dmin:.2f} A  qA={qa:+.2f} qB={qb:+.2f}  "
              f"E_direct={(coul+lj)/4.184:8.2f} kcal/mol  E_ddd={(coul_ddd+lj_ddd)/4.184:7.2f}",
              flush=True)

    # per-mutation destabilization proxy = sum of interaction energies of abolished bridges
    per_mut = {}
    for mut, killed in KNOCKOUT.items():
        removed_direct = sum(per_bridge[b]["E_total_direct_kcal"] for b in killed)
        removed_ddd = sum(per_bridge[b]["E_total_ddd_kcal"] for b in killed)
        per_mut[mut] = {
            "bridges_broken": killed, "n_broken": len(killed),
            # lock energy REMOVED (structure loses this stabilization) -> positive destabilization
            "lock_energy_removed_direct_kcal": round(-removed_direct, 2),
            "lock_energy_removed_ddd_kcal": round(-removed_ddd, 2),
        }
        print(f"  {mut}: breaks {killed} -> removes "
              f"{-removed_direct:.1f} kcal/mol (direct), {-removed_ddd:.1f} (ddd)", flush=True)

    ranking = sorted(per_mut, key=lambda m: per_mut[m]["lock_energy_removed_direct_kcal"], reverse=True)
    out = {
        "model": "direct ff14SB nonbonded (Coulomb eps=1 + Lennard-Jones) between full residues; "
                 "ddd = distance-dependent dielectric eps=4r(A). WT extended, GB-OBC2 build.",
        "per_bridge": per_bridge,
        "per_mutation": per_mut,
        "destabilization_ranking": ranking,
        "note": "Static upper bound on lock strength; GB solvent screens magnitude but not ranking. "
                "The pending snap-back MD tests the ranking dynamically.",
    }
    os.makedirs("results/route_a", exist_ok=True)
    with open("results/route_a/lock_energy.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\ndestabilization ranking (most lock energy removed first):", ranking, flush=True)
    print("wrote results/route_a/lock_energy.json", flush=True)


if __name__ == "__main__":
    main()
