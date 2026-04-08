#!/usr/bin/env python3
"""Analyze domain steering results — compare inter-domain geometry across presets."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Domain definitions (chain, start_res, end_res) — 1-indexed inclusive
AVB3_DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "alpha_tail":       ("A", 742, 962),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}

HINGE_ANGLES = [
    ("alpha_head_thigh", "alpha_calf", "alpha_tail"),
    ("beta_head", "beta_tail", "alpha_tail"),
    ("alpha_head_thigh", "alpha_calf", "beta_head"),
]

HINGE_DISTANCES = [
    ("alpha_head_thigh", "alpha_tail"),
    ("beta_head", "alpha_tail"),
    ("alpha_head_thigh", "beta_tail"),
]

ANGLE_LABELS = [
    "α-leg opening\n(head—calf—tail)",
    "β-leg opening\n(βhead—βtail—αtail)",
    "Headpiece-leg\n(αhead—calf—βhead)",
]

DISTANCE_LABELS = [
    "αHead ↔ αTail",
    "βHead ↔ αTail",
    "αHead ↔ βTail",
]


def parse_pdb_ca(pdb_path):
    """Parse CA atoms from PDB, return dict of {(chain, resseq): (x,y,z)}."""
    cas = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21]
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            cas[(chain, resseq)] = np.array([x, y, z])
    return cas


def domain_centroid(cas, domain_name):
    """Compute centroid of CA atoms in a domain."""
    chain, start, end = AVB3_DOMAINS[domain_name]
    coords = []
    for (c, r), xyz in cas.items():
        if c == chain and start <= r <= end:
            coords.append(xyz)
    if not coords:
        return None
    return np.mean(coords, axis=0)


def compute_angle(p1, p2, p3):
    """Angle at p2 in degrees."""
    v1 = p1 - p2
    v2 = p3 - p2
    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
    return np.degrees(np.arccos(np.clip(cos_a, -1, 1)))


def analyze_structure(pdb_path):
    """Compute inter-domain angles and distances for a PDB."""
    cas = parse_pdb_ca(pdb_path)
    if not cas:
        return None, None

    centroids = {name: domain_centroid(cas, name) for name in AVB3_DOMAINS}

    angles = []
    for d1, d2, d3 in HINGE_ANGLES:
        c1, c2, c3 = centroids[d1], centroids[d2], centroids[d3]
        if c1 is None or c2 is None or c3 is None:
            angles.append(np.nan)
        else:
            angles.append(compute_angle(c1, c2, c3))

    distances = []
    for d1, d2 in HINGE_DISTANCES:
        c1, c2 = centroids[d1], centroids[d2]
        if c1 is None or c2 is None:
            distances.append(np.nan)
        else:
            # PDB coords in Angstroms
            distances.append(np.linalg.norm(c1 - c2))

    return angles, distances


def parse_production_log(log_path):
    """Parse production.log CSV for energy and temperature time series."""
    steps, times, energies, temps, speeds = [], [], [], [], []
    with open(log_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split(",")
            if len(parts) >= 5:
                steps.append(int(parts[0]))
                times.append(float(parts[1]))
                energies.append(float(parts[2]))
                temps.append(float(parts[3]))
                speeds.append(float(parts[4]))
    return {
        "steps": np.array(steps),
        "time_ps": np.array(times),
        "energy_kjmol": np.array(energies),
        "temperature_K": np.array(temps),
        "speed_nsday": np.array(speeds),
    }


def main():
    base = Path("/home2/Documents/code/conformers/data/runs/avb3/domain_steering")
    presets = ["gentle_open", "moderate_open", "cv_distance_extend"]
    preset_labels = ["Gentle Open", "Moderate Open", "CV Distance Extend"]
    colors = ["#2196F3", "#FF9800", "#4CAF50"]

    # Analyze minimized (starting) and final structures
    results = {}
    for preset in presets:
        d = base / preset
        min_angles, min_dists = analyze_structure(d / "minimized.pdb")
        fin_angles, fin_dists = analyze_structure(d / "final.pdb")
        prod_data = parse_production_log(d / "production.log")
        results[preset] = {
            "min_angles": min_angles,
            "min_dists": min_dists,
            "fin_angles": fin_angles,
            "fin_dists": fin_dists,
            "production": prod_data,
        }

    # Also get the input structure angles (from minimized of any preset as reference)
    ref_angles = results[presets[0]]["min_angles"]
    ref_dists = results[presets[0]]["min_dists"]

    # ── Create figure ──
    fig = plt.figure(figsize=(16, 14))
    fig.suptitle("AVB3 Integrin Domain Steering — Comparison of 3 Methods",
                 fontsize=16, fontweight='bold', y=0.98)
    gs = GridSpec(3, 2, hspace=0.4, wspace=0.35, top=0.93, bottom=0.06,
                  left=0.08, right=0.95)

    # Panel 1: Inter-domain angles (grouped bar)
    ax1 = fig.add_subplot(gs[0, 0])
    x = np.arange(len(HINGE_ANGLES))
    width = 0.2
    ax1.bar(x - 1.5*width, ref_angles, width, color='#9E9E9E', label='Starting', alpha=0.7)
    for i, (preset, label, color) in enumerate(zip(presets, preset_labels, colors)):
        fin = results[preset]["fin_angles"]
        ax1.bar(x + (i-0.5)*width, fin, width, color=color, label=label, alpha=0.85)
    ax1.set_xticks(x)
    ax1.set_xticklabels(ANGLE_LABELS, fontsize=8)
    ax1.set_ylabel("Angle (degrees)")
    ax1.set_title("Inter-Domain Hinge Angles")
    ax1.legend(fontsize=7, loc='upper right')
    ax1.grid(axis='y', alpha=0.3)

    # Panel 2: Inter-domain distances (grouped bar)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.bar(x - 1.5*width, ref_dists, width, color='#9E9E9E', label='Starting', alpha=0.7)
    for i, (preset, label, color) in enumerate(zip(presets, preset_labels, colors)):
        fin = results[preset]["fin_dists"]
        ax2.bar(x + (i-0.5)*width, fin, width, color=color, label=label, alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(DISTANCE_LABELS, fontsize=8)
    ax2.set_ylabel("Distance (Å)")
    ax2.set_title("Inter-Domain Centroid Distances")
    ax2.legend(fontsize=7, loc='upper right')
    ax2.grid(axis='y', alpha=0.3)

    # Panel 3: Angle changes (delta from start)
    ax3 = fig.add_subplot(gs[1, 0])
    for i, (preset, label, color) in enumerate(zip(presets, preset_labels, colors)):
        deltas = np.array(results[preset]["fin_angles"]) - np.array(ref_angles)
        ax3.barh(x + i*0.25 - 0.25, deltas, 0.22, color=color, label=label, alpha=0.85)
    ax3.set_yticks(x)
    ax3.set_yticklabels(ANGLE_LABELS, fontsize=8)
    ax3.set_xlabel("Δ Angle (degrees)")
    ax3.set_title("Angle Change from Starting Structure")
    ax3.axvline(0, color='black', linewidth=0.5)
    ax3.legend(fontsize=7)
    ax3.grid(axis='x', alpha=0.3)

    # Panel 4: Distance changes (delta from start)
    ax4 = fig.add_subplot(gs[1, 1])
    for i, (preset, label, color) in enumerate(zip(presets, preset_labels, colors)):
        deltas = np.array(results[preset]["fin_dists"]) - np.array(ref_dists)
        ax4.barh(x + i*0.25 - 0.25, deltas, 0.22, color=color, label=label, alpha=0.85)
    ax4.set_yticks(x)
    ax4.set_yticklabels(DISTANCE_LABELS, fontsize=8)
    ax4.set_xlabel("Δ Distance (Å)")
    ax4.set_title("Distance Change from Starting Structure")
    ax4.axvline(0, color='black', linewidth=0.5)
    ax4.legend(fontsize=7)
    ax4.grid(axis='x', alpha=0.3)

    # Panel 5: Energy over time
    ax5 = fig.add_subplot(gs[2, 0])
    for preset, label, color in zip(presets, preset_labels, colors):
        d = results[preset]["production"]
        ax5.plot(d["time_ps"], d["energy_kjmol"] / 1e6, color=color, label=label,
                 linewidth=0.8, alpha=0.8)
    ax5.set_xlabel("Time (ps)")
    ax5.set_ylabel("Potential Energy (×10⁶ kJ/mol)")
    ax5.set_title("Production Energy Stability")
    ax5.legend(fontsize=7)
    ax5.grid(alpha=0.3)

    # Panel 6: Temperature over time
    ax6 = fig.add_subplot(gs[2, 1])
    for preset, label, color in zip(presets, preset_labels, colors):
        d = results[preset]["production"]
        ax6.plot(d["time_ps"], d["temperature_K"], color=color, label=label,
                 linewidth=0.8, alpha=0.8)
    ax6.axhline(310, color='red', linestyle='--', linewidth=0.5, alpha=0.5, label='Target (310K)')
    ax6.set_xlabel("Time (ps)")
    ax6.set_ylabel("Temperature (K)")
    ax6.set_title("Production Temperature Stability")
    ax6.legend(fontsize=7)
    ax6.grid(alpha=0.3)

    # Add summary table as text
    summary_lines = []
    for preset, label in zip(presets, preset_labels):
        r = results[preset]
        d_angles = np.array(r["fin_angles"]) - np.array(ref_angles)
        d_dists = np.array(r["fin_dists"]) - np.array(ref_dists)
        summary_lines.append(
            f"{label}: Δangles=[{d_angles[0]:+.1f}°, {d_angles[1]:+.1f}°, {d_angles[2]:+.1f}°]  "
            f"Δdist=[{d_dists[0]:+.1f}Å, {d_dists[1]:+.1f}Å, {d_dists[2]:+.1f}Å]"
        )

    for line in summary_lines:
        print(line)

    out_path = Path("/home2/Documents/code/conformers/figures/domain_steering_comparison.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {out_path}")
    plt.close()


if __name__ == "__main__":
    main()
