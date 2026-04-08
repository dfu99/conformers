#!/usr/bin/env python3
"""Domain-preserving steering forces for integrin conformational exploration.

Three methods for steering inter-domain angles without disrupting intra-domain
folding, all implemented as OpenMM force objects:

1. Centroid angle torques (CustomCentroidBondForce)
   - Biases inter-domain angles toward target values
   - Each domain moves as a quasi-rigid body

2. Rigid-body domain restraints + inter-domain pulling
   - Adds intra-domain harmonic restraints to preserve folding
   - Pulling forces act only between domain centroids

3. Collective variable biasing
   - Defines inter-domain distances/angles as CVs
   - Applies harmonic or flat-bottom biases to steer CVs

Usage:
    from domain_steering import (
        add_centroid_angle_torque,
        add_domain_restraints_with_pulling,
        add_cv_bias,
        AVB3_DOMAINS,
    )
    # Then call any method to add forces to an OpenMM System
"""
from __future__ import annotations

import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from openmm import System
    from openmm.app import Topology

try:
    from openmm import (
        CustomCentroidBondForce,
        CustomExternalForce,
        CustomCompoundBondForce,
    )
    from openmm.unit import nanometers, kilojoules_per_mole, radians
except ImportError:
    from simtk.openmm import (
        CustomCentroidBondForce,
        CustomExternalForce,
        CustomCompoundBondForce,
    )
    from simtk.unit import nanometers, kilojoules_per_mole, radians

import numpy as np


# ── AVB3 integrin domain definitions ─────────────────────────────────────────
# (chain_id, start_resSeq, end_resSeq) — 1-indexed, inclusive
AVB3_DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 962),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}

# Inter-domain hinge angles to bias (triplets of domain names)
# angle(domain1_centroid, hinge_centroid, domain2_centroid)
AVB3_HINGE_ANGLES = [
    ("alpha_head_thigh", "alpha_calf", "alpha_tail"),   # α-leg opening
    ("beta_head",        "beta_tail",  "alpha_tail"),    # β-leg opening
    ("alpha_head_thigh", "alpha_calf", "beta_head"),     # headpiece-leg angle
]

# Inter-domain distances to bias
AVB3_HINGE_DISTANCES = [
    ("alpha_head_thigh", "alpha_tail"),   # head-tail extension
    ("beta_head",        "alpha_tail"),   # cross-chain extension
    ("alpha_head_thigh", "beta_tail"),    # head-β-tail distance
]


def _select_atoms_by_range(topology, chain_id: str, start: int, end: int) -> list[int]:
    """Select atom indices for a chain/residue range."""
    indices = []
    for atom in topology.atoms():
        if (atom.residue.chain.id == chain_id and
                start <= atom.residue.index + 1 <= end):
            indices.append(atom.index)
    # Fallback: try matching by residue.id (string PDB resSeq)
    if not indices:
        for atom in topology.atoms():
            try:
                resseq = int(atom.residue.id)
            except (ValueError, AttributeError):
                continue
            if atom.residue.chain.id == chain_id and start <= resseq <= end:
                indices.append(atom.index)
    return indices


def _get_domain_atoms(topology, domain_name: str) -> list[int]:
    """Get atom indices for a named domain."""
    chain, start, end = AVB3_DOMAINS[domain_name]
    return _select_atoms_by_range(topology, chain, start, end)


# ═══════════════════════════════════════════════════════════════════════════════
# Method 1: Centroid Angle Torques
# ═══════════════════════════════════════════════════════════════════════════════

def add_centroid_angle_torque(
    system: "System",
    topology: "Topology",
    angle_triplets: list[tuple[str, str, str]] | None = None,
    target_angles_deg: list[float] | None = None,
    force_constant: float = 100.0,  # kJ/mol/rad²
    schedule: dict[int, float] | None = None,
) -> "CustomCentroidBondForce":
    """Add centroid-based angle torques between integrin domains.

    Applies a harmonic bias on the angle formed by three domain centroids:
        E = 0.5 * k * (angle - target)²

    This torques the inter-domain hinges without applying forces to
    individual atoms. Each domain moves as a quasi-rigid body.

    Args:
        system: OpenMM System
        topology: OpenMM Topology
        angle_triplets: List of (domain1, hinge, domain2) name triplets
        target_angles_deg: Target angle per triplet (degrees). If None,
            opens each angle by 30° from initial.
        force_constant: Spring constant in kJ/mol/rad²
        schedule: Optional {step: force_constant} schedule for ramping

    Returns:
        The CustomCentroidBondForce added to the system.
    """
    if angle_triplets is None:
        angle_triplets = AVB3_HINGE_ANGLES

    # Create centroid bond force with angle bias
    # Energy: 0.5 * k * (angle(g1, g2, g3) - theta0)^2
    force = CustomCentroidBondForce(3,
        "0.5 * k * (angle(g1, g2, g3) - theta0)^2")
    force.addPerBondParameter("theta0")
    force.addGlobalParameter("k", force_constant)

    for i, (d1, d2, d3) in enumerate(angle_triplets):
        atoms1 = _get_domain_atoms(topology, d1)
        atoms2 = _get_domain_atoms(topology, d2)
        atoms3 = _get_domain_atoms(topology, d3)

        if not atoms1 or not atoms2 or not atoms3:
            print(f"  WARNING: empty domain in triplet ({d1}, {d2}, {d3})")
            continue

        # Add groups with equal weights
        g1 = force.addGroup(atoms1)
        g2 = force.addGroup(atoms2)
        g3 = force.addGroup(atoms3)

        # Default target: open by 30° from straight (180° - 30° = ~2.62 rad)
        if target_angles_deg and i < len(target_angles_deg):
            theta0 = target_angles_deg[i] * math.pi / 180.0
        else:
            theta0 = 150.0 * math.pi / 180.0  # 150° default (partially open)

        force.addBond([g1, g2, g3], [theta0])
        print(f"  Angle torque: {d1}—{d2}—{d3} → target={theta0*180/math.pi:.1f}°, "
              f"k={force_constant} kJ/mol/rad²")

    system.addForce(force)
    return force


# ═══════════════════════════════════════════════════════════════════════════════
# Method 2: Rigid-Body Domain Restraints + Inter-Domain Pulling
# ═══════════════════════════════════════════════════════════════════════════════

def add_domain_restraints_with_pulling(
    system: "System",
    topology: "Topology",
    positions_nm: np.ndarray,
    domains: dict[str, tuple[str, int, int]] | None = None,
    restraint_k: float = 1000.0,  # kJ/mol/nm² for intra-domain
    pull_pairs: list[tuple[str, str]] | None = None,
    pull_force_pn: float = 2.0,
    pull_mode: str = "centroid_distance",
) -> tuple:
    """Add intra-domain restraints and inter-domain centroid pulling.

    Within each domain: harmonic position restraints on CA atoms preserve
    the domain's internal fold. Between domains: a pulling force acts on
    centroids to change inter-domain distances.

    Args:
        system: OpenMM System
        topology: OpenMM Topology
        positions_nm: Current positions in nanometers (N×3)
        domains: Domain definitions {name: (chain, start, end)}
        restraint_k: Intra-domain restraint spring constant (kJ/mol/nm²)
        pull_pairs: List of (domain1, domain2) pairs to pull apart
        pull_force_pn: Pulling force in picoNewtons
        pull_mode: "centroid_distance" or "centroid_direction"

    Returns:
        (restraint_force, pull_force) tuple
    """
    if domains is None:
        domains = AVB3_DOMAINS
    if pull_pairs is None:
        pull_pairs = AVB3_HINGE_DISTANCES

    # Intra-domain position restraints on CA atoms
    restraint_force = CustomExternalForce(
        "0.5 * k_domain * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint_force.addGlobalParameter("k_domain", restraint_k)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    n_restrained = 0
    for dname, (chain, start, end) in domains.items():
        for atom in topology.atoms():
            if atom.name != "CA":
                continue
            try:
                resseq = int(atom.residue.id)
            except (ValueError, AttributeError):
                resseq = atom.residue.index + 1
            if atom.residue.chain.id == chain and start <= resseq <= end:
                idx = atom.index
                if idx < len(positions_nm):
                    x0, y0, z0 = positions_nm[idx]
                    restraint_force.addParticle(idx, [x0, y0, z0])
                    n_restrained += 1

    system.addForce(restraint_force)
    print(f"  Domain restraints: {n_restrained} CA atoms restrained "
          f"(k={restraint_k} kJ/mol/nm²)")

    # Inter-domain centroid pulling
    # Convert pN to kJ/(mol·nm): F_openmm = F_pN * 6.022e-4
    pull_force_value = pull_force_pn * 6.022e-4

    pull_force = CustomCentroidBondForce(2,
        "-f_pull * distance(g1, g2)")
    pull_force.addGlobalParameter("f_pull", pull_force_value)

    for d1, d2 in pull_pairs:
        atoms1 = _get_domain_atoms(topology, d1)
        atoms2 = _get_domain_atoms(topology, d2)
        if not atoms1 or not atoms2:
            continue
        g1 = pull_force.addGroup(atoms1)
        g2 = pull_force.addGroup(atoms2)
        pull_force.addBond([g1, g2], [])
        print(f"  Centroid pull: {d1} ↔ {d2}, F={pull_force_pn} pN")

    system.addForce(pull_force)
    return restraint_force, pull_force


# ═══════════════════════════════════════════════════════════════════════════════
# Method 3: Collective Variable Biasing
# ═══════════════════════════════════════════════════════════════════════════════

def add_cv_bias(
    system: "System",
    topology: "Topology",
    cv_type: str = "distance",  # "distance" or "angle"
    cv_pairs: list[tuple[str, str]] | None = None,
    cv_triplets: list[tuple[str, str, str]] | None = None,
    target_values: list[float] | None = None,
    force_constant: float = 500.0,  # kJ/mol/nm² for distance, kJ/mol/rad² for angle
    bias_type: str = "harmonic",  # "harmonic" or "flat_bottom"
    flat_bottom_width: float = 0.5,  # nm for distance, rad for angle
) -> "CustomCentroidBondForce":
    """Add collective variable biasing on inter-domain distances or angles.

    Defines inter-domain distances or angles as collective variables and
    applies harmonic or flat-bottom biases to steer them toward targets.

    Args:
        system: OpenMM System
        topology: OpenMM Topology
        cv_type: "distance" for inter-domain distances, "angle" for angles
        cv_pairs: List of (domain1, domain2) pairs (for distance CVs)
        cv_triplets: List of (d1, d2, d3) triplets (for angle CVs)
        target_values: Target CV values (nm for distance, degrees for angle)
        force_constant: Bias spring constant
        bias_type: "harmonic" or "flat_bottom"
        flat_bottom_width: Width of flat region for flat_bottom bias

    Returns:
        The CustomCentroidBondForce added to the system.
    """
    if cv_type == "distance":
        if cv_pairs is None:
            cv_pairs = AVB3_HINGE_DISTANCES

        if bias_type == "harmonic":
            energy_expr = "0.5 * k_cv * (distance(g1, g2) - d0)^2"
        else:
            energy_expr = ("0.5 * k_cv * max(0, abs(distance(g1, g2) - d0) - w)^2")

        force = CustomCentroidBondForce(2, energy_expr)
        force.addPerBondParameter("d0")
        force.addGlobalParameter("k_cv", force_constant)
        if bias_type == "flat_bottom":
            force.addGlobalParameter("w", flat_bottom_width)

        for i, (d1, d2) in enumerate(cv_pairs):
            atoms1 = _get_domain_atoms(topology, d1)
            atoms2 = _get_domain_atoms(topology, d2)
            if not atoms1 or not atoms2:
                continue

            g1 = force.addGroup(atoms1)
            g2 = force.addGroup(atoms2)

            if target_values and i < len(target_values):
                d0 = target_values[i]  # in nm
            else:
                d0 = 15.0  # default target: 15 nm (extended)

            force.addBond([g1, g2], [d0])
            print(f"  CV bias ({bias_type}): distance({d1}, {d2}) → {d0:.1f} nm, "
                  f"k={force_constant}")

    elif cv_type == "angle":
        if cv_triplets is None:
            cv_triplets = AVB3_HINGE_ANGLES

        if bias_type == "harmonic":
            energy_expr = "0.5 * k_cv * (angle(g1, g2, g3) - a0)^2"
        else:
            energy_expr = ("0.5 * k_cv * max(0, abs(angle(g1, g2, g3) - a0) - w)^2")

        force = CustomCentroidBondForce(3, energy_expr)
        force.addPerBondParameter("a0")
        force.addGlobalParameter("k_cv", force_constant)
        if bias_type == "flat_bottom":
            force.addGlobalParameter("w", flat_bottom_width)

        for i, (d1, d2, d3) in enumerate(cv_triplets):
            atoms1 = _get_domain_atoms(topology, d1)
            atoms2 = _get_domain_atoms(topology, d2)
            atoms3 = _get_domain_atoms(topology, d3)
            if not atoms1 or not atoms2 or not atoms3:
                continue

            g1 = force.addGroup(atoms1)
            g2 = force.addGroup(atoms2)
            g3 = force.addGroup(atoms3)

            if target_values and i < len(target_values):
                a0 = target_values[i] * math.pi / 180.0
            else:
                a0 = 160.0 * math.pi / 180.0  # default: nearly straight

            force.addBond([g1, g2, g3], [a0])
            print(f"  CV bias ({bias_type}): angle({d1}, {d2}, {d3}) → "
                  f"{a0*180/math.pi:.1f}°, k={force_constant}")

    else:
        raise ValueError(f"Unknown cv_type: {cv_type}")

    system.addForce(force)
    return force


# ═══════════════════════════════════════════════════════════════════════════════
# Convenience: Apply a named steering preset
# ═══════════════════════════════════════════════════════════════════════════════

STEERING_PRESETS = {
    "gentle_open": {
        "method": "centroid_angle",
        "description": "Gentle angle opening — slowly straightens inter-domain hinges",
        "force_constant": 50.0,
        "target_angles_deg": [160.0, 160.0, 150.0],
    },
    "moderate_open": {
        "method": "centroid_angle",
        "description": "Moderate opening with stronger bias",
        "force_constant": 200.0,
        "target_angles_deg": [170.0, 170.0, 160.0],
    },
    "restrained_pull": {
        "method": "domain_restrained",
        "description": "Pull centroids apart while restraining internal domain structure",
        "restraint_k": 200.0,
        "pull_force_pn": 0.5,
    },
    "cv_distance_extend": {
        "method": "cv_bias",
        "cv_type": "distance",
        "description": "Bias inter-domain distances toward extended targets",
        "force_constant": 200.0,
        "target_values": [20.0, 18.0, 16.0],  # nm
        "bias_type": "flat_bottom",
        "flat_bottom_width": 2.0,
    },
}


def apply_steering_preset(
    system: "System",
    topology: "Topology",
    positions_nm: np.ndarray,
    preset_name: str,
) -> None:
    """Apply a named steering preset to the system.

    Args:
        system: OpenMM System
        topology: OpenMM Topology
        positions_nm: Current positions in nanometers
        preset_name: One of STEERING_PRESETS keys
    """
    if preset_name not in STEERING_PRESETS:
        available = ", ".join(STEERING_PRESETS.keys())
        raise ValueError(f"Unknown preset '{preset_name}'. Available: {available}")

    preset = STEERING_PRESETS[preset_name]
    method = preset["method"]
    print(f"\n=== Steering preset: {preset_name} ===")
    print(f"  {preset['description']}")

    if method == "centroid_angle":
        add_centroid_angle_torque(
            system, topology,
            target_angles_deg=preset.get("target_angles_deg"),
            force_constant=preset.get("force_constant", 100.0),
        )
    elif method == "domain_restrained":
        add_domain_restraints_with_pulling(
            system, topology, positions_nm,
            restraint_k=preset.get("restraint_k", 500.0),
            pull_force_pn=preset.get("pull_force_pn", 1.0),
        )
    elif method == "cv_bias":
        add_cv_bias(
            system, topology,
            cv_type=preset.get("cv_type", "distance"),
            target_values=preset.get("target_values"),
            force_constant=preset.get("force_constant", 200.0),
            bias_type=preset.get("bias_type", "harmonic"),
            flat_bottom_width=preset.get("flat_bottom_width", 0.5),
        )
    else:
        raise ValueError(f"Unknown method: {method}")


def apply_custom_cv_steering(
    system: "System",
    topology: "Topology",
    target_distances_nm: list[float],
    distance_pairs: list[tuple[str, str]] | None = None,
    force_constant: float = 200.0,
    bias_type: str = "flat_bottom",
    flat_bottom_width: float = 2.0,
) -> "CustomCentroidBondForce":
    """Apply CV distance steering with arbitrary target distances.

    This is the entry point for AFMFold inference: the CNN predicts
    inter-domain distances, and this function creates the corresponding
    steering forces to generate a PDB matching those distances.

    Args:
        system: OpenMM System
        topology: OpenMM Topology
        target_distances_nm: Target inter-domain distances in nm.
            Must match the number of distance_pairs.
        distance_pairs: Domain pair names. Defaults to AVB3_HINGE_DISTANCES.
        force_constant: Bias spring constant (kJ/mol/nm²)
        bias_type: "harmonic" or "flat_bottom"
        flat_bottom_width: Width of flat region (nm)

    Returns:
        The CustomCentroidBondForce added to the system.
    """
    if distance_pairs is None:
        distance_pairs = AVB3_HINGE_DISTANCES

    if len(target_distances_nm) != len(distance_pairs):
        raise ValueError(
            f"Got {len(target_distances_nm)} target distances but "
            f"{len(distance_pairs)} domain pairs"
        )

    print(f"\n=== Custom CV Steering ===")
    print(f"  Targets (nm): {target_distances_nm}")

    return add_cv_bias(
        system, topology,
        cv_type="distance",
        cv_pairs=distance_pairs,
        target_values=target_distances_nm,
        force_constant=force_constant,
        bias_type=bias_type,
        flat_bottom_width=flat_bottom_width,
    )


if __name__ == "__main__":
    print("Domain-preserving steering methods for integrin αVβ3")
    print("\nAvailable presets:")
    for name, preset in STEERING_PRESETS.items():
        print(f"  {name}: {preset['description']}")
    print("\nAvailable domains:")
    for name, (chain, start, end) in AVB3_DOMAINS.items():
        print(f"  {name}: chain {chain}, residues {start}-{end}")
    print("\nHinge angles:")
    for d1, d2, d3 in AVB3_HINGE_ANGLES:
        print(f"  {d1} — {d2} — {d3}")
