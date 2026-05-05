#!/usr/bin/env python3
"""Route-A day-1: build αVβ3 + endothelial-membrane topology.

Audit §12.7 + docs/route_a_kickoff.md §3 #2 deliverable. Applies
the force-field gotchas from kickoff §2 (Amber99SB-ILDN + Slipids,
endothelial lipid composition, no αIIb heavy/light split,
4 αV β-propeller Ca²⁺ sites) to a starting αVβ3 ectodomain
structure (default: 1JV2).

This is a **stub**: it documents the intended interface and prints
the steps that should run. Full implementation is on hold pending
PI sign-off and the route-A kickoff (PACE allocation).

Usage (when wired up):
  python build_av_topology.py \
      --pdb data/multi_integrin/raw_pdbs/1JV2.pdb \
      --out pipelines/route_a/topology/avb3_ready/ \
      --membrane endothelial \
      --ff amber99sb-ildn \
      --equil-ns 50
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


GOTCHAS = {
    "ff_combo": "Amber99SB-ILDN + Slipids (match Ferg-Lab; sensitivity check post-convergence)",
    "membrane": {
        "endothelial": {"POPC": 50, "POPE": 15, "POPS": 5,
                        "cholesterol": 25, "PIP2": 0},
        "platelet":    {"POPC": 30, "POPE": 25, "POPS": 15,
                        "cholesterol": 25, "PIP2": 5},
    },
    "calcium_sites": [(220, 225), (281, 286), (348, 353), (412, 417)],
    "glycosylation_sites_av_uniprot": [74, 178, 307, 401, 642, 801],
    "alpha_chain_continuity": "single peptide (no αIIb-style heavy/light split)",
    "psi_disulfide": "shared β3 (port verbatim)",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdb", type=Path, required=True,
                   help="Starting αVβ3 ectodomain PDB (e.g., 1JV2)")
    p.add_argument("--out", type=Path, required=True,
                   help="Output topology directory")
    p.add_argument("--membrane", choices=["endothelial", "platelet"],
                   default="endothelial")
    p.add_argument("--ff", default="amber99sb-ildn",
                   help="Force-field family")
    p.add_argument("--equil-ns", type=float, default=50.0,
                   help="Equilibration length in ns (Risk-2 gate: ≥ 50)")
    p.add_argument("--dry-run", action="store_true", default=True,
                   help="STUB MODE: print steps only, do not run.")
    return p.parse_args()


def stage(name: str, status: str, detail: str = "") -> None:
    marker = {"PLANNED": " ⏳", "TODO": " ⊙", "DONE": " ✓"}.get(status, "  ")
    print(f"{marker}  [{status:>7s}] {name}")
    if detail:
        print(f"            {detail}")


def main() -> int:
    args = parse_args()
    print(f"build_av_topology.py — STUB MODE\n")
    print(f"Input PDB         : {args.pdb}")
    print(f"Output directory  : {args.out}")
    print(f"Membrane          : {args.membrane}")
    print(f"Force field       : {args.ff}")
    print(f"Equilibration     : {args.equil_ns} ns\n")

    if not args.pdb.exists():
        print(f"  ⚠ PDB not found at {args.pdb}", file=sys.stderr)
        return 1

    print("Planned stages:\n")

    stage("1. Load + clean αVβ3 PDB",
          "PLANNED",
          "drop heteroatoms except Ca²⁺/Mg²⁺/Mn²⁺; "
          "use pdbtools or mdtraj.")

    stage("2. Strip αIIb heavy/light chain bond logic",
          "DONE",
          f"αV is single peptide; gotcha §2.6 — "
          f"{GOTCHAS['alpha_chain_continuity']}.")

    stage("3. Identify αV β-propeller Ca²⁺ sites",
          "PLANNED",
          f"4 sites at {GOTCHAS['calcium_sites']} per Xiong 2001.")

    stage("4. Skip αV glycosylation (Risk-4 stage 1)",
          "PLANNED",
          f"glycan sites in UniProt at residues "
          f"{GOTCHAS['glycosylation_sites_av_uniprot']} — defer.")

    stage("5. Add hydrogens via tleap",
          "PLANNED",
          f"force field: {GOTCHAS['ff_combo']}.")

    membrane_comp = GOTCHAS["membrane"][args.membrane]
    membrane_str = ", ".join(f"{k} {v}%" for k, v in membrane_comp.items() if v)
    stage("6. Build endothelial bilayer (CHARMM-GUI Membrane Builder)",
          "PLANNED",
          f"composition: {membrane_str}")

    stage("7. Solvate + ionize",
          "PLANNED",
          "TIP3P water, 150 mM NaCl, neutralize.")

    stage("8. Energy minimize",
          "PLANNED",
          "5000 steepest-descent steps.")

    stage("9. Equilibrate at 300 K",
          "PLANNED",
          f"{args.equil_ns} ns NPT, monitor "
          f"bilayer thickness (target 4.0 nm).")

    stage("10. Validate Risk-2 gate",
          "PLANNED",
          "no membrane rupture, no hexagonal phase, "
          "stable secondary structure ≥ 80% in headpiece.")

    print("\nAll stages PLANNED — full implementation activates "
          "after PI sign-off.")
    print("\nFor now this stub validates the interface. To wire up:")
    print("  1. Replace stage prints with real OpenMM / CHARMM-GUI calls.")
    print("  2. Use Ferg-Lab repo paths for stage 5 + 6 templates.")
    print("  3. Hook Risk-2 gate into return-code so SLURM can branch.")

    if args.dry_run:
        return 0
    print("\n⚠ --dry-run was disabled but no real implementation exists yet.",
          file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
