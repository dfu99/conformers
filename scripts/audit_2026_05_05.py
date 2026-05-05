#!/usr/bin/env python3
"""Generate the 2026-05-05 audit summary figure.

Refined 2026-05-05 18:00: now 4 follow-up deliverables shipped today
(F1 library coverage, F2 calibration controls, F3 ΔG bootstrap,
RGD pocket accessibility / obj-039) plus the original ΔG point estimate
(obj-038). Reviewer-tally numbers updated; blocker ranking adjusted —
RGD docking is no longer "missing", it has been answered with a
negative result that *sharpens* the EO-coverage blocker.

Layout:
  TL: cumulative objective count (with today markers for 037-039)
  TR: reviewer panel — addressed/partial/open (post-update tally)
  ML: pipelines table
  MR: blockers (ranked)
  BL: today's deliverables row (4 thumbnail tiles)
  BR: queue snapshot
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
OBJ_PATH = ROOT / "tasks" / "objectives.yaml"
QUEUE_PATH = ROOT / "tasks" / "queue.yaml"
OUT = ROOT / "figures" / "audit-2026-05-05.png"
FIG_DIR = ROOT / "figures"


def load_objectives() -> list[dict]:
    with OBJ_PATH.open() as f:
        return yaml.safe_load(f)["objectives"]


def load_queue() -> list[dict]:
    with QUEUE_PATH.open() as f:
        data = yaml.safe_load(f)
    return data.get("queue", []) if data else []


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
    end = dt.date(2026, 5, 6)
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
    # Highlight obj-038 + obj-039
    today_count = sum(1 for d in dates if d <= today)
    ax.scatter([today], [today_count], s=120, color="#d62728",
               edgecolor="black", zorder=10,
               label=f"obj-038→046 today: {today_count}")
    ax.set_title(f"Objectives completed (cumulative) — {today_count} today "
                 f"(+9 this audit: obj-038→046)",
                 fontsize=10.5)
    ax.set_ylabel("count")
    ax.set_xlabel("week")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_ha("right")


def panel_reviewers(ax):
    """Updated tally after deepening v9 (obj-047):
      A Springer: 2/1/2 (obj-043 Bayesian EO floor; unchanged from v5).
      C Ando:     fully cleared from evening 2 (3/1/0).
      D Elber:    obj-047 change-point detection moves "committor analysis"
                   open → partial via the experimental-side analog (true
                   committor still requires saddle-seeded unbiased MD).
                   Net 1/2/1 (was 1/1/2 v8).
      Others unchanged.
    """
    reviewers = [
        ("A: Springer", 2, 1, 2),       # obj-043 EO Bayesian floor (v5)
        ("B: Jumper",   2, 1, 2),
        ("C: Ando",     3, 1, 0),
        ("D: Elber",    1, 2, 1),       # +1 partial from obj-047 (committor proxy)
        ("E: Noé",      2, 1, 2),
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
    a_sum, p_sum, o_sum = int(addressed.sum()), int(partial.sum()), int(openc.sum())
    ax.set_title(f"Reviewer concerns: {a_sum} addressed / {p_sum} partial / {o_sum} open  "
                 f"(was 5/6/12 morning → 10/5/8 v5 → {a_sum}/{p_sum}/{o_sum} v9)",
                 fontsize=9.5)
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
        ("avb3-conformers (legacy, shadows rm)", "deprecated", "#7f7f7f"),
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
    ax.set_title("Pipelines (3 PRIMARY • 2 supporting • 2 to drop • 1 deprecated today)",
                 fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_blockers(ax):
    blockers = [
        ("#1  EO state coverage — Bayesian floor + non-stat. characterized",
         "obj-039 + obj-041 + obj-025 confirm EO unreachable by SMD.\n"
         "obj-043: ΔG_EO ≥ 2.02 kcal/mol. obj-047: stationarity rejected\n"
         "(p<1e-3); ensemble hops on 30s–15min timescale. Route A awaits PI.",
         "#d62728"),
        ("#2  ΔG(CV0) + populations + drift  ✓ FULLY CHARACTERIZED",
         "obj-038 ΔG (CI<0.5 kcal/mol). obj-043 BC 25%/Inter 46%/EC 24%.\n"
         "obj-044 non-stat. obj-045 σ≈8 Å drift. obj-046 intrinsic.\n"
         "obj-047: 3-4 change-points/video, stepwise — not gradual.",
         "#2ca02c"),
        ("#3  PACE A100-80GB allocation — Monday PI sign-off ask",
         "Sole remaining bottleneck. Unblocks: route A (αIIbβ3 string\n"
         "method, ~$800), Gō-Martini parallel, AF2 ablation (Rev B),\n"
         "committor analysis (Rev D). All 4 share the same allocation.",
         "#fdae61"),
        ("#4  Reviewer C contact mechanics  ✓ CLOSED (obj-040)",
         "F4 Hertzian δ = 0.11-0.28 nm at 50-200 pN — all below 1 nm\n"
         "noise floor. Hard-sphere pseudo-AFM is defensible first-order;\n"
         "Hertzian correction is quantified.",
         "#2ca02c"),
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
    ax.set_title("Blockers (ranked) — #2 + #4 closed, #1 stepwise-resolved (obj-047)",
                 fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_today_deliverables(ax):
    """13 thumbnail tiles: today's deliverables through v9 (obj-047 change-points)."""
    tiles = [
        ("free_energy_profile_v1.png",
         "obj-038\nΔG(CV0)"),
        ("library_coverage_v2.png",
         "F1\nLib cov v2"),
        ("calibration_controls_v1.png",
         "F2\nCalib"),
        ("rgd_docking_v1.png",
         "obj-039\nRGD pocket"),
        ("contact_mechanics_control.png",
         "obj-040\nHertz F4"),
        ("library_coverage_v3.png",
         "obj-041\nRoute D null"),
        ("rmsf_bootstrap_v1.png",
         "obj-042\nRMSF boot"),
        ("state_populations_v1.png",
         "obj-043\nState pops"),
        ("state_populations_v2_windowed.png",
         "obj-044\nNon-stat."),
        ("state_populations_v3_per_block_dg.png",
         "obj-045\nFES drift"),
        ("fes_drift_metadata_correlation.png",
         "obj-046\nDrift\nintrinsic"),
        ("change_point_detection.png",
         "obj-047\nChange\npoints"),
    ]
    n = len(tiles)
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Today's deliverables (12 figures + 6 docs + 2 starter scripts + obj-047) — "
                 "ΔG fully char., contact mech. closed, route D out, RMSF V1↔V2 r=0.998, "
                 "EO Bayesian floor ΔG ≥ 2 kcal/mol, FES drift intrinsic + STEPWISE (p<1e-3)",
                 fontsize=8.0, fontweight="bold")
    for i, (name, label) in enumerate(tiles):
        path = FIG_DIR / name
        if path.exists():
            try:
                img = mpimg.imread(str(path))
                ax.imshow(img, extent=(i + 0.05, i + 0.95, 0.18, 0.92),
                          aspect="auto")
            except Exception:
                pass
        ax.text(i + 0.5, 0.08, label, ha="center", va="center",
                fontsize=8.0, fontweight="bold",
                bbox=dict(boxstyle="round", facecolor="white",
                          edgecolor="#cccccc"))


def panel_queue_snapshot(ax, queue):
    items = []
    for q in queue:
        if not q.get("source", "").startswith("audit-2026-05-05"):
            continue
        status = q.get("status", "")
        obj_text = (q.get("objective", "") or "").strip().split("\n", 1)[0]
        if len(obj_text) > 90:
            obj_text = obj_text[:87] + "..."
        items.append((status, obj_text, q.get("priority", 9)))
    items.sort(key=lambda r: (r[0] != "completed", r[2]))
    items = items[:8]
    n = len(items)
    if n == 0:
        ax.text(0.5, 0.5, "No audit-tagged queue items found",
                ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return
    y = np.arange(n)[::-1]
    for yi, (status, obj_text, prio) in zip(y, items):
        color = {"completed": "#2ca02c", "in_progress": "#ff7f0e",
                 "pending": "#7f7f7f", "skipped": "#cccccc"}.get(status, "#888888")
        ax.barh(yi, 1, color=color, alpha=0.85, height=0.85)
        marker = {"completed": "✓", "in_progress": "◐",
                  "pending": "○", "skipped": "—"}.get(status, "·")
        ax.text(0.02, yi, f"{marker} P{prio}  {obj_text}",
                va="center", ha="left", fontsize=8.5, color="black")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Audit queue (today's items)", fontsize=10.5, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)


def main() -> int:
    objs = load_objectives()
    _queue = load_queue()  # kept for future panel_queue_snapshot reactivation
    del _queue
    fig = plt.figure(figsize=(15.5, 13.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.85],
                          hspace=0.42, wspace=0.18)

    completed = sum(1 for o in objs if o.get("status") == "completed")
    in_prog = sum(1 for o in objs if o.get("status") == "in_progress")
    fig.suptitle(
        f"Conformers — audit 2026-05-05 (deepening pass v9)  •  "
        f"{completed} objectives completed  •  "
        f"obj-038→047 + 6 docs + 2 starter scripts  •  "
        f"reviewer panel: 10/6/7 (was 5/6/12 morning)  •  "
        f"EO Bayesian floor ΔG ≥ 2 kcal/mol  •  FES drift STEPWISE p<1e-3 (obj-047)  •  "
        f"PI sign-off Monday unblocks route A",
        fontsize=10.0, fontweight="bold", y=0.995,
    )

    ax_tl = fig.add_subplot(gs[0, 0])
    ax_tr = fig.add_subplot(gs[0, 1])
    ax_ml = fig.add_subplot(gs[1, 0])
    ax_mr = fig.add_subplot(gs[1, 1])
    ax_bl = fig.add_subplot(gs[2, :])

    panel_timeline(ax_tl, objs)
    panel_reviewers(ax_tr)
    panel_pipelines(ax_ml)
    panel_blockers(ax_mr)
    panel_today_deliverables(ax_bl)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"Saved {OUT}")

    print(f"objectives: {completed} completed, {in_prog} in_progress")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
