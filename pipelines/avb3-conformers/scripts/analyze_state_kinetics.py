#!/usr/bin/env python3
"""Analyze conformational state kinetics from matched CVs.

Classifies each HS-AFM frame into integrin conformational states
(BC=bent-closed, EC=extended-closed, Intermediate) based on matched
inter-domain distances from the correlation matching pipeline.

Generates:
  - State trajectory over time (color-coded strip)
  - CV trajectory with state shading
  - State percentage pie/bar chart
  - Dwell time distributions
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec


# State definitions based on inter-domain distances
# CV0: alpha_head_thigh ↔ alpha_calf (alpha leg extension)
# CV1: beta_head_hybrid_egf1 ↔ beta_tail (beta leg extension)
# CV2: alpha_head ↔ beta_head (headpiece separation, ~constant)
STATE_COLORS = {
    "BC": "#2166ac",      # blue - bent-closed
    "EC": "#b2182b",      # red - extended-closed
    "I":  "#fdae61",      # orange - intermediate
    "EO": "#1b7837",      # green - extended-open (not observed)
}

STATE_LABELS = {
    "BC": "Bent-Closed",
    "EC": "Extended-Closed",
    "I":  "Intermediate",
    "EO": "Extended-Open",
}

# Thresholds (Angstroms)
# Higher CV0/CV1 = more extended (head further from calf/tail)
BC_CV0_MAX = 65.0   # alpha leg compact
BC_CV1_MAX = 60.0   # beta leg compact
EC_CV0_MIN = 78.0   # alpha leg extended
EC_CV1_MIN = 73.0   # beta leg extended
EO_CV2_MIN = 42.0   # headpiece opening (not observed in current library)


def classify_state(cv0, cv1, cv2):
    """Classify a single frame's conformational state."""
    if cv2 > EO_CV2_MIN and cv0 > EC_CV0_MIN:
        return "EO"
    elif cv0 >= EC_CV0_MIN and cv1 >= EC_CV1_MIN:
        return "EC"
    elif cv0 <= BC_CV0_MAX and cv1 <= BC_CV1_MAX:
        return "BC"
    else:
        return "I"


def classify_trajectory(cvs):
    """Classify all frames."""
    return np.array([classify_state(c[0], c[1], c[2]) for c in cvs])


def compute_dwell_times(states):
    """Compute dwell times (consecutive frames in same state)."""
    dwells = {}
    current_state = states[0]
    current_length = 1
    for s in states[1:]:
        if s == current_state:
            current_length += 1
        else:
            dwells.setdefault(current_state, []).append(current_length)
            current_state = s
            current_length = 1
    dwells.setdefault(current_state, []).append(current_length)
    return dwells


def count_transitions(states):
    """Count state-to-state transitions."""
    transitions = {}
    for i in range(len(states) - 1):
        key = (states[i], states[i + 1])
        transitions[key] = transitions.get(key, 0) + 1
    return transitions


def plot_kinetics(cvs, states, correlations, video_name, output_dir,
                  frame_time_ms=100.0):
    """Generate composite kinetics figure."""
    n = len(states)
    time_s = np.arange(n) * frame_time_ms / 1000.0

    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(5, 2, figure=fig, height_ratios=[0.5, 1.2, 1.2, 1.0, 1.0],
                  hspace=0.35, wspace=0.3)

    # --- Row 0: State color strip ---
    ax_strip = fig.add_subplot(gs[0, :])
    for i, s in enumerate(states):
        ax_strip.axvspan(time_s[i], time_s[min(i + 1, n - 1)],
                         color=STATE_COLORS[s], alpha=0.9)
    ax_strip.set_xlim(time_s[0], time_s[-1])
    ax_strip.set_yticks([])
    ax_strip.set_xlabel("Time (s)")
    ax_strip.set_title(f"{video_name} — Conformational State Trajectory",
                       fontsize=13, fontweight="bold")
    # Legend
    patches = [Patch(facecolor=STATE_COLORS[s], label=STATE_LABELS[s])
               for s in ["BC", "I", "EC", "EO"] if s in set(states)]
    ax_strip.legend(handles=patches, loc="upper right", ncol=len(patches),
                    fontsize=9, framealpha=0.9)

    # --- Row 1: CV0 and CV1 with state shading ---
    ax_cv = fig.add_subplot(gs[1, :])
    # Background shading
    for i, s in enumerate(states):
        ax_cv.axvspan(time_s[i], time_s[min(i + 1, n - 1)],
                      color=STATE_COLORS[s], alpha=0.15)
    ax_cv.plot(time_s, cvs[:n, 0], color="#d73027", linewidth=0.8,
               label="CV0: α-leg extension", alpha=0.9)
    ax_cv.plot(time_s, cvs[:n, 1], color="#4575b4", linewidth=0.8,
               label="CV1: β-leg extension", alpha=0.9)
    ax_cv.axhline(EC_CV0_MIN, color="#d73027", linestyle="--", alpha=0.4,
                  linewidth=0.7)
    ax_cv.axhline(BC_CV0_MAX, color="#d73027", linestyle=":", alpha=0.4,
                  linewidth=0.7)
    ax_cv.axhline(EC_CV1_MIN, color="#4575b4", linestyle="--", alpha=0.4,
                  linewidth=0.7)
    ax_cv.axhline(BC_CV1_MAX, color="#4575b4", linestyle=":", alpha=0.4,
                  linewidth=0.7)
    ax_cv.set_ylabel("Distance (Å)")
    ax_cv.set_xlabel("Time (s)")
    ax_cv.legend(fontsize=9, loc="lower right")
    ax_cv.set_xlim(time_s[0], time_s[-1])
    ax_cv.set_title("Inter-domain Distances (α/β leg extension)", fontsize=11)

    # --- Row 2: Correlation over time ---
    ax_corr = fig.add_subplot(gs[2, :])
    for i, s in enumerate(states):
        ax_corr.axvspan(time_s[i], time_s[min(i + 1, n - 1)],
                        color=STATE_COLORS[s], alpha=0.15)
    ax_corr.plot(time_s, correlations[:n], color="black", linewidth=0.6,
                 alpha=0.8)
    ax_corr.set_ylabel("Correlation")
    ax_corr.set_xlabel("Time (s)")
    ax_corr.set_xlim(time_s[0], time_s[-1])
    ax_corr.set_ylim(0.85, 1.0)
    ax_corr.set_title("Fitting Correlation per Frame", fontsize=11)

    # --- Row 3 left: State percentages ---
    ax_pie = fig.add_subplot(gs[3, 0])
    unique_states = ["BC", "I", "EC"]
    counts = [np.sum(states == s) for s in unique_states]
    pcts = [100.0 * c / n for c in counts]
    # Only show states that exist
    mask = [c > 0 for c in counts]
    pie_states = [s for s, m in zip(unique_states, mask) if m]
    pie_counts = [c for c, m in zip(counts, mask) if m]
    pie_pcts = [p for p, m in zip(pcts, mask) if m]
    pie_colors = [STATE_COLORS[s] for s in pie_states]
    pie_labels = [f"{STATE_LABELS[s]}\n{p:.1f}% ({c} frames)"
                  for s, p, c in zip(pie_states, pie_pcts, pie_counts)]
    ax_pie.pie(pie_counts, labels=pie_labels, colors=pie_colors,
               autopct="", startangle=90, textprops={"fontsize": 9})
    ax_pie.set_title("State Occupancy", fontsize=11, fontweight="bold")

    # --- Row 3 right: Dwell time distributions ---
    ax_dwell = fig.add_subplot(gs[3, 1])
    dwells = compute_dwell_times(states)
    for s in ["BC", "I", "EC"]:
        if s in dwells and len(dwells[s]) > 1:
            vals = dwells[s]
            times = [v * frame_time_ms / 1000.0 for v in vals]
            ax_dwell.hist(times, bins=20, alpha=0.6, color=STATE_COLORS[s],
                          label=f"{STATE_LABELS[s]} (n={len(vals)}, "
                                f"mean={np.mean(times):.2f}s)")
    ax_dwell.set_xlabel("Dwell Time (s)")
    ax_dwell.set_ylabel("Count")
    ax_dwell.legend(fontsize=8)
    ax_dwell.set_title("Dwell Time Distribution", fontsize=11)

    # --- Row 4: Summary table ---
    ax_table = fig.add_subplot(gs[4, :])
    ax_table.axis("off")
    # Build summary text
    trans = count_transitions(states)
    lines = [f"{'State':<20} {'Frames':>8} {'%':>8} {'Mean Dwell (s)':>16} "
             f"{'Transitions In':>16} {'Transitions Out':>16}"]
    lines.append("─" * 90)
    for s in ["BC", "I", "EC"]:
        c = np.sum(states == s)
        p = 100.0 * c / n
        d = dwells.get(s, [0])
        mean_d = np.mean(d) * frame_time_ms / 1000.0 if d else 0
        t_in = sum(v for k, v in trans.items() if k[1] == s and k[0] != s)
        t_out = sum(v for k, v in trans.items() if k[0] == s and k[1] != s)
        lines.append(f"{STATE_LABELS[s]:<20} {c:>8d} {p:>7.1f}% "
                     f"{mean_d:>15.2f} {t_in:>16d} {t_out:>16d}")
    lines.append("─" * 90)
    lines.append(f"Total frames: {n}  |  "
                 f"Total transitions: {sum(1 for i in range(n-1) if states[i] != states[i+1])}  |  "
                 f"Frame rate: {1000/frame_time_ms:.0f} fps  |  "
                 f"Duration: {time_s[-1]:.1f}s")
    ax_table.text(0.02, 0.95, "\n".join(lines), transform=ax_table.transAxes,
                  fontsize=8.5, fontfamily="monospace", verticalalignment="top")

    fig.savefig(output_dir / f"state_kinetics.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_dir / 'state_kinetics.png'}")

    # Also save state assignments
    np.save(str(output_dir / "state_assignments.npy"), states)
    print(f"Saved: {output_dir / 'state_assignments.npy'}")

    return states, dwells, trans


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matched-cvs", type=Path, required=True)
    p.add_argument("--correlations", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--video-name", default="Video")
    p.add_argument("--frame-time-ms", type=float, default=100.0,
                   help="Time per frame in ms (Linz HS-AFM ~100ms)")
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    cvs = np.load(str(args.matched_cvs))
    corr = np.load(str(args.correlations))
    n = len(corr)
    cvs = cvs[:n]

    print(f"Loaded {n} frames")
    print(f"CV ranges:")
    for i, name in enumerate(["α-leg (aH↔aC)", "β-leg (bH↔bT)", "head sep (aH↔bH)"]):
        print(f"  CV{i} ({name}): [{cvs[:,i].min():.1f}, {cvs[:,i].max():.1f}] Å")

    states = classify_trajectory(cvs)
    print(f"\nState classification:")
    for s in ["BC", "I", "EC", "EO"]:
        c = np.sum(states == s)
        if c > 0:
            print(f"  {STATE_LABELS[s]}: {c} frames ({100*c/n:.1f}%)")

    plot_kinetics(cvs, states, corr, args.video_name, args.output_dir,
                  args.frame_time_ms)


if __name__ == "__main__":
    main()
