#!/usr/bin/env python3
"""Generate the 2026-05-05 audit summary figure.

4-panel project status board:
  TL: objective completion timeline (weekly cumulative)
  TR: reviewer panel — addressed/partial/open per reviewer
  BL: pipeline state — PRIMARY / supporting / drop
  BR: blocker map (rank-ordered)
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
OBJ_PATH = ROOT / "tasks" / "objectives.yaml"
OUT = ROOT / "figures" / "audit-2026-05-05.png"


def load_objectives() -> list[dict]:
    with OBJ_PATH.open() as f:
        return yaml.safe_load(f)["objectives"]


def panel_timeline(ax, objs):
    dates = []
    for o in objs:
        d = o.get("completed_date")
        if d and o.get("status") == "completed":
            dates.append(dt.date.fromisoformat(d))
    dates.sort()
    if not dates:
        return
    start = dt.date(2026, 3, 1)
    end = dt.date(2026, 5, 5)
    weeks = []
    cumulative = []
    cur = start
    while cur <= end:
        n = sum(1 for d in dates if d <= cur)
        weeks.append(cur)
        cumulative.append(n)
        cur = cur + dt.timedelta(days=7)
    ax.plot(weeks, cumulative, color="#1f77b4", linewidth=2.4)
    ax.fill_between(weeks, cumulative, color="#1f77b4", alpha=0.15)
    ax.scatter(weeks, cumulative, s=18, color="#1f77b4", zorder=5)
    today = dt.date(2026, 5, 5)
    ax.axvline(today, color="#d62728", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.text(today, max(cumulative) * 0.06, "  today",
            color="#d62728", fontsize=8, ha="left", va="bottom")
    ax.set_title("Objectives completed (cumulative)", fontsize=11)
    ax.set_ylabel("count")
    ax.set_xlabel("week")
    ax.grid(alpha=0.3)
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_ha("right")


def panel_reviewers(ax):
    reviewers = [
        ("A: Springer", 1, 1, 3),
        ("B: Jumper",   2, 1, 2),
        ("C: Ando",     1, 2, 1),
        ("D: Elber",    0, 1, 3),
        ("E: Noé",      1, 1, 3),
    ]
    labels = [r[0] for r in reviewers]
    addressed = np.array([r[1] for r in reviewers])
    partial = np.array([r[2] for r in reviewers])
    openc = np.array([r[3] for r in reviewers])
    y = np.arange(len(reviewers))
    ax.barh(y, addressed, color="#2ca02c", label="addressed")
    ax.barh(y, partial, left=addressed, color="#ff7f0e", label="partial")
    ax.barh(y, openc, left=addressed + partial, color="#d62728", label="open")
    for i, (a, p, o) in enumerate(zip(addressed, partial, openc)):
        ax.text(a + p + o + 0.1, i, f"{a}/{p}/{o}", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 7)
    ax.set_title("Reviewer concerns: addressed / partial / open", fontsize=11)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.set_xlabel("count")
    ax.grid(axis="x", alpha=0.3)


def panel_pipelines(ax):
    rows = [
        ("sim-afm-video",         "PRIMARY",    "#1b7837"),
        ("conformer-library",     "PRIMARY",    "#1b7837"),
        ("afm-overlay",           "PRIMARY",    "#1b7837"),
        ("protenix-avb3-template", "supporting", "#7570b3"),
        ("protenix-a5b1",         "supporting", "#7570b3"),
        ("boltz",                 "DROP",       "#d62728"),
        ("afcluster",             "DROP",       "#d62728"),
        ("avb3-conformers (legacy)", "DROP",   "#d62728"),
    ]
    y = np.arange(len(rows))[::-1]
    for yi, (name, label, color) in zip(y, rows):
        ax.barh(yi, 1, color=color, alpha=0.85)
        ax.text(0.04, yi, name, va="center", ha="left",
                fontsize=9.5, color="white", fontweight="bold")
        ax.text(0.96, yi, label, va="center", ha="right",
                fontsize=8, color="white")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Pipelines (3 PRIMARY • 2 supporting • 3 to drop)",
                 fontsize=11)
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_blockers(ax):
    blockers = [
        ("#1  EO state coverage (publication blocker)",
         "SMD cannot open αVβ3 headpiece in 3 ns at any k (obj-025).\n"
         "Library spans CV0=47-85 Å; EO is >85 Å. Needs metadynamics,\n"
         "REMD, or Ferg-Lab αIIbβ3 string method — all multi-week.",
         "#d62728"),
        ("#2  No ΔG(CV0) free-energy profile  ◀ NEXT SHIP",
         "Reviewer D's top concern. 1645 unbiased frames already\n"
         "exist (V1=379 + V2=1266 fitted). −kT log P(CV0) at 300 K.\n"
         "60-90 min, all CPU, no MD required.",
         "#ff7f0e"),
        ("#3  RunPod A100 GPU window",
         "Blocks αIIbβ3 steering MD launch. External, queue-managed\n"
         "via mc runpod subscribe. Not actionable on local hardware\n"
         "(12 GB local GPU is too small for solvated integrin).",
         "#fdae61"),
        ("#4  RGD docking proof-of-concept missing",
         "Reviewer E (Noé/Chodera). Conformational-selection use\n"
         "case unproven. ~1 day with smina against 5 bent + 5\n"
         "extended v7 library conformers.",
         "#fee090"),
    ]
    y = np.arange(len(blockers))[::-1]
    for yi, (title, body, color) in zip(y, blockers):
        ax.barh(yi, 1, color=color, alpha=0.85, height=0.86)
        ax.text(0.02, yi + 0.27, title, va="center", ha="left",
                fontsize=10.5, color="black", fontweight="bold")
        ax.text(0.02, yi - 0.05, body, va="top", ha="left",
                fontsize=8, color="black")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, len(blockers) - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Blockers (ranked)", fontsize=11)
    for spine in ax.spines.values():
        spine.set_visible(False)


def main() -> int:
    objs = load_objectives()
    fig, axes = plt.subplots(2, 2, figsize=(14, 9.5))
    fig.suptitle(
        "Conformers — audit 2026-05-05  •  35 objectives shipped, "
        "but reviewer panel: 5 addressed / 6 partial / 12 open",
        fontsize=13, fontweight="bold", y=0.995,
    )
    panel_timeline(axes[0, 0], objs)
    panel_reviewers(axes[0, 1])
    panel_pipelines(axes[1, 0])
    panel_blockers(axes[1, 1])
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"Saved {OUT}")

    completed = sum(1 for o in objs if o.get("status") == "completed")
    in_prog = sum(1 for o in objs if o.get("status") == "in_progress")
    print(f"objectives: {completed} completed, {in_prog} in_progress")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
