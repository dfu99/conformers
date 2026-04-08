#!/usr/bin/env python3
"""Run RoyalMD-style MD simulation with domain-preserving steering forces.

This wraps the standard OpenMM MD workflow (solvate, minimize, equilibrate,
production) but replaces the brute-force head/tail pulling with domain-aware
steering that preserves internal domain structure.

Usage:
    python run_domain_steering.py \
        --input-pdb AVB3_clean.pdb \
        --output-dir outputs/ \
        --steering-preset gentle_open \
        --production-time 1000
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-pdb", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--royalmd-root", type=Path, default=None)
    parser.add_argument("--steering-preset", default="gentle_open",
                        help="Steering preset name (gentle_open, moderate_open, "
                             "restrained_pull, cv_distance_extend)")
    parser.add_argument("--target-distances", type=str, default=None,
                        help="Custom CV target distances in nm (comma-separated). "
                             "Overrides --steering-preset with custom CV steering. "
                             "E.g. '15.0,12.0,10.0' for 3 domain pair targets.")
    parser.add_argument("--production-time", type=float, default=1000.0,
                        help="Production time in ps")
    parser.add_argument("--temperature", type=float, default=310.0)
    parser.add_argument("--timestep", type=float, default=0.002)
    parser.add_argument("--report-interval", type=int, default=1500)
    parser.add_argument("--save-frames", action="store_true", default=False,
                        help="Save PDB frames during production for AFMFold training")
    parser.add_argument("--frame-interval", type=int, default=5000,
                        help="Steps between saved PDB frames (default: 5000 → ~300 frames/ns)")
    args = parser.parse_args()

    # Add conformers scripts to path for domain_steering module
    scripts_dir = Path(__file__).resolve().parent
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))

    from domain_steering import (
        apply_steering_preset, apply_custom_cv_steering, STEERING_PRESETS,
    )

    # Validate args: either custom targets or a known preset
    use_custom_targets = args.target_distances is not None
    if not use_custom_targets and args.steering_preset not in STEERING_PRESETS:
        print(f"ERROR: Unknown preset '{args.steering_preset}'")
        print(f"Available: {', '.join(STEERING_PRESETS.keys())}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Import OpenMM
    try:
        from openmm.app import (
            PDBFile, ForceField, Modeller, PME, HBonds, Simulation,
            StateDataReporter,
        )
        from openmm import (
            LangevinMiddleIntegrator, MonteCarloBarostat, Platform,
        )
        from openmm.unit import (
            kelvin, picoseconds, atmospheres, nanometers,
            kilojoules_per_mole, molar,
        )
    except ImportError:
        from simtk.openmm.app import (
            PDBFile, ForceField, Modeller, PME, HBonds, Simulation,
            StateDataReporter,
        )
        from simtk.openmm import (
            LangevinMiddleIntegrator, MonteCarloBarostat, Platform,
        )
        from simtk.unit import (
            kelvin, picoseconds, atmospheres, nanometers,
            kilojoules_per_mole, molar,
        )

    print(f"=== Domain-Preserving Steering ===")
    print(f"Input: {args.input_pdb}")
    print(f"Preset: {args.steering_preset}")
    print(f"Production: {args.production_time} ps")

    # Load PDB
    pdb = PDBFile(str(args.input_pdb))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(pdb.topology, pdb.positions)

    # Add missing hydrogens (PDB often lacks them; Amber14 requires them)
    print("\nAdding hydrogens...")
    modeller.addHydrogens(forcefield, pH=7.0)

    # Solvate
    print("\nSolvating...")
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=1.0 * nanometers,
        ionicStrength=0.15 * molar,
    )

    # Create system
    print("Creating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometers,
        constraints=HBonds,
    )

    # Get positions as numpy array in nm
    positions = modeller.positions
    positions_nm = np.array([
        [p[0].value_in_unit(nanometers),
         p[1].value_in_unit(nanometers),
         p[2].value_in_unit(nanometers)]
        for p in positions
    ])

    # Apply domain-preserving steering forces (instead of head/tail pulling)
    if use_custom_targets:
        target_dists = [float(x) for x in args.target_distances.split(",")]
        print(f"Custom CV targets (nm): {target_dists}")
        apply_custom_cv_steering(system, modeller.topology, target_dists)
    else:
        apply_steering_preset(system, modeller.topology, positions_nm,
                              args.steering_preset)

    # Add barostat
    system.addForce(MonteCarloBarostat(
        1.0 * atmospheres, args.temperature * kelvin, 25))

    # Create integrator and simulation
    integrator = LangevinMiddleIntegrator(
        args.temperature * kelvin,
        1.0 / picoseconds,
        args.timestep * picoseconds,
    )

    # Try GPU platform
    try:
        platform = Platform.getPlatformByName("CUDA")
        simulation = Simulation(modeller.topology, system, integrator, platform)
    except Exception:
        platform = Platform.getPlatformByName("CPU")
        simulation = Simulation(modeller.topology, system, integrator, platform)

    simulation.context.setPositions(modeller.positions)

    # Minimize
    print("\nMinimizing...")
    simulation.minimizeEnergy()

    # Save minimized structure
    state = simulation.context.getState(getPositions=True)
    with open(args.output_dir / "minimized.pdb", "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Equilibrate (100 ps NVT)
    print("Equilibrating (100 ps NVT)...")
    simulation.context.setVelocitiesToTemperature(args.temperature * kelvin)

    try:
        import mdtraj.reporters
        equil_reporter = mdtraj.reporters.NetCDFReporter(
            str(args.output_dir / "equilibration.nc"), args.report_interval)
        simulation.reporters.append(equil_reporter)
    except ImportError:
        pass

    equil_steps = int(100.0 / args.timestep)
    simulation.step(equil_steps)

    # Save equilibrated
    state = simulation.context.getState(getPositions=True)
    with open(args.output_dir / "equilibrated.pdb", "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Production with steering forces active
    print(f"\nProduction ({args.production_time} ps)...")
    prod_steps = int(args.production_time / args.timestep)

    # Add production reporters
    simulation.reporters = []
    try:
        import mdtraj.reporters
        prod_reporter = mdtraj.reporters.NetCDFReporter(
            str(args.output_dir / "production.nc"), args.report_interval)
        simulation.reporters.append(prod_reporter)
    except ImportError:
        pass

    simulation.reporters.append(StateDataReporter(
        str(args.output_dir / "production.log"),
        args.report_interval,
        step=True, time=True, potentialEnergy=True,
        temperature=True, speed=True,
    ))

    t0 = time.time()

    if args.save_frames:
        # Save PDB frames at intervals for AFMFold training
        frames_dir = args.output_dir / "frames"
        frames_dir.mkdir(exist_ok=True)
        frame_idx = 0
        steps_done = 0

        # Save initial frame (equilibrated state)
        state = simulation.context.getState(getPositions=True)
        with open(frames_dir / f"frame_{frame_idx:04d}.pdb", "w") as f:
            PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        frame_idx += 1

        while steps_done < prod_steps:
            chunk = min(args.frame_interval, prod_steps - steps_done)
            simulation.step(chunk)
            steps_done += chunk
            state = simulation.context.getState(getPositions=True)
            with open(frames_dir / f"frame_{frame_idx:04d}.pdb", "w") as f:
                PDBFile.writeFile(simulation.topology, state.getPositions(), f)
            frame_idx += 1
            if frame_idx % 50 == 0:
                print(f"  Saved {frame_idx} frames ({steps_done}/{prod_steps} steps)")

        print(f"  Saved {frame_idx} frames to {frames_dir}")
    else:
        simulation.step(prod_steps)

    elapsed = time.time() - t0

    print(f"Production complete ({elapsed:.1f}s)")

    # Save final state
    state = simulation.context.getState(getPositions=True)
    with open(args.output_dir / "final.pdb", "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Write summary
    summary = {
        "steering_preset": args.steering_preset,
        "production_time_ps": args.production_time,
        "temperature_K": args.temperature,
        "elapsed_seconds": elapsed,
        "input_pdb": str(args.input_pdb),
    }
    (args.output_dir / "steering_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\nResults in {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
