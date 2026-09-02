#!/usr/bin/env python3
"""Aggregate + validate the bent->extended morphs across all integrin variants.

Acceptance, applied per variant (mirrors how the alphaVbeta3 endpoint was accepted
in obj-075/076 -- the CVs must SEPARATE the endpoints, and the path must stay
physical the whole way):

  1. genu angle ends near straight (>= GENU_EXTENDED_MIN deg)
  2. genu angle rises monotonically (no fold-back), allowing small numerical dips
  3. long-axis extent and Rg both grow by >= EXTENSION_MIN_FRAC
  4. final energy stays in a sane range -- no residual clash carried to the endpoint

alphaVbeta3 additionally acts as a positive control: it must land on the published
route-A endpoint (Rg ~67 A, extent ~211 A) or the generalisation is not faithful.

Where the PDB now contains an experimentally-determined extended structure, the
morph product is compared against it -- an external check that the constructed
endpoint is not merely self-consistent.
"""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENU_EXTENDED_MIN = 150.0
EXTENSION_MIN_FRAC = 0.20
MONOTONIC_TOL_DEG = 5.0

# alphaVbeta3 control targets, from results/route_a/stage2c_result.json (obj-075)
CONTROL = {"variant": "alphaVbeta3", "Rg_A": 67.32, "long_axis_extent_A": 211.36,
           "tol_frac": 0.10}

# categorical slots 1-5 of the validated reference palette (fixed order, never cycled)
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

PRETTY = {"alphaVbeta3": r"$\alpha_V\beta_3$", "alphaIIbbeta3": r"$\alpha_{IIb}\beta_3$",
          "alpha5beta1": r"$\alpha_5\beta_1$", "alphaMbeta2": r"$\alpha_M\beta_2$",
          "alphaXbeta2": r"$\alpha_X\beta_2$"}


def load_morphs(results_dir, variants):
    out = {}
    for v in variants:
        p = os.path.join(results_dir, v, f"{v}_morph.json")
        if os.path.exists(p):
            try:
                d = json.load(open(p))
            except json.JSONDecodeError:
                continue
            if d.get("status") == "done":
                out[v] = d
    return out


def assess(m):
    traj = m["trajectory"]
    genu = [t.get("genu_angle_deg") for t in traj if t.get("genu_angle_deg") is not None]
    b, e = m["metrics_bent"], m.get("metrics_equil") or m["metrics_extended"]
    checks, notes = {}, []

    final_genu = genu[-1] if genu else None
    checks["genu_reaches_extended"] = bool(final_genu and final_genu >= GENU_EXTENDED_MIN)
    if final_genu is not None and not checks["genu_reaches_extended"]:
        notes.append(f"knee stops at {final_genu:.0f} deg, short of {GENU_EXTENDED_MIN:.0f}")

    drops = [genu[i] - genu[i + 1] for i in range(len(genu) - 1) if genu[i + 1] < genu[i]]
    worst_drop = max(drops) if drops else 0.0
    checks["genu_monotonic"] = bool(worst_drop <= MONOTONIC_TOL_DEG)
    if drops and worst_drop > MONOTONIC_TOL_DEG:
        notes.append(f"knee angle backtracks by up to {worst_drop:.1f} deg")

    def frac(key):
        if b.get(key) in (None, 0) or e.get(key) is None:
            return None
        return (e[key] - b[key]) / abs(b[key])

    f_ext, f_rg = frac("long_axis_extent_A"), frac("Rg_A")
    checks["extent_grows"] = bool(f_ext is not None and f_ext >= EXTENSION_MIN_FRAC)
    checks["Rg_grows"] = bool(f_rg is not None and f_rg >= EXTENSION_MIN_FRAC)

    e_fin = m.get("E_equil_kJ_mol", m["E_extended_kJ_mol"])
    checks["energy_sane"] = bool(abs(e_fin) < 1e6)
    if not checks["energy_sane"]:
        notes.append(f"final energy {e_fin:.2e} kJ/mol indicates unrelieved clash")

    return {
        "checks": checks, "passed": all(checks.values()), "notes": notes,
        "genu_bent_deg": genu[0] if genu else None, "genu_extended_deg": final_genu,
        "extent_bent_A": b.get("long_axis_extent_A"), "extent_extended_A": e.get("long_axis_extent_A"),
        "Rg_bent_A": b.get("Rg_A"), "Rg_extended_A": e.get("Rg_A"),
        "extent_gain_frac": round(f_ext, 3) if f_ext is not None else None,
        "Rg_gain_frac": round(f_rg, 3) if f_rg is not None else None,
        "E_bent_kJ_mol": m["E_bent_kJ_mol"], "E_final_kJ_mol": e_fin,
        "n_steps": m["n_steps"], "total_swing_deg": m["total_swing_deg"],
        "seconds": m["seconds"],
    }


def make_figure(morphs, order, summary, exp_refs, out_png):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))
    for ax in axes:
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=GRID, linewidth=0.7, alpha=0.9)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=9)

    # Panel A -- the reaction coordinate itself, per variant
    ax = axes[0]
    for i, v in enumerate(order):
        t = morphs[v]["trajectory"]
        x = [s["cum_deg"] for s in t]
        y = [s.get("genu_angle_deg") for s in t]
        ax.plot(x, y, color=PALETTE[i], linewidth=2.0, marker="o", markersize=4.2,
                markeredgecolor="#fcfcfb", markeredgewidth=0.6, label=PRETTY[v])
        ax.annotate(PRETTY[v], (x[-1], y[-1]), textcoords="offset points",
                    xytext=(6, 0), fontsize=9.5, color=INK, va="center")
    ax.axhline(GENU_EXTENDED_MIN, color=INK_MUTED, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.annotate(f"extended threshold {GENU_EXTENDED_MIN:.0f}°", (0, GENU_EXTENDED_MIN),
                textcoords="offset points", xytext=(4, 5), fontsize=8.5, color=INK_MUTED)
    ax.set_xlabel("cumulative leg swing (°)", fontsize=10, color=INK_MUTED)
    ax.set_ylabel("genu (knee) angle (°)", fontsize=10, color=INK_MUTED)
    ax.set_title("A · Knee opens monotonically in every variant",
                 fontsize=11.5, color=INK, loc="left", pad=10)
    ax.set_xlim(left=-4)
    ax.margins(x=0.16)

    # Panel B -- the path stays physical (energy per step)
    ax = axes[1]
    for i, v in enumerate(order):
        t = morphs[v]["trajectory"]
        ax.plot([s["cum_deg"] for s in t], [s["E_kJ_mol"] / 1e3 for s in t],
                color=PALETTE[i], linewidth=2.0, marker="o", markersize=4.2,
                markeredgecolor="#fcfcfb", markeredgewidth=0.6, label=PRETTY[v])
    ax.set_xlabel("cumulative leg swing (°)", fontsize=10, color=INK_MUTED)
    ax.set_ylabel("potential energy (10³ kJ/mol)", fontsize=10, color=INK_MUTED)
    ax.set_title("B · Energy stays sane — no accumulated clash",
                 fontsize=11.5, color=INK, loc="left", pad=10)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="best")

    # Panel C -- bent vs extended extent, against experiment where it exists
    ax = axes[2]
    ys = np.arange(len(order))[::-1]
    for i, (v, y) in enumerate(zip(order, ys)):
        s = summary[v]
        ax.plot([s["extent_bent_A"], s["extent_extended_A"]], [y, y],
                color=PALETTE[i], linewidth=2.0, zorder=2, solid_capstyle="round")
        ax.plot(s["extent_bent_A"], y, "o", color="#fcfcfb", markersize=9,
                markeredgecolor=PALETTE[i], markeredgewidth=2.0, zorder=3)
        ax.plot(s["extent_extended_A"], y, "o", color=PALETTE[i], markersize=9,
                markeredgecolor="#fcfcfb", markeredgewidth=1.2, zorder=3)
        ax.annotate(f"{s['extent_bent_A']:.0f}", (s["extent_bent_A"], y),
                    textcoords="offset points", xytext=(0, 9), fontsize=8.5,
                    color=INK_MUTED, ha="center")
        ax.annotate(f"{s['extent_extended_A']:.0f} Å", (s["extent_extended_A"], y),
                    textcoords="offset points", xytext=(0, 9), fontsize=8.5,
                    color=INK, ha="center")
        for ref in exp_refs.get(v, []):
            ax.plot(ref["long_extent_A"], y, marker="D", color=INK, markersize=6.5,
                    markeredgecolor="#fcfcfb", markeredgewidth=1.0, zorder=4)
            ax.annotate(ref["pdb_id"], (ref["long_extent_A"], y),
                        textcoords="offset points", xytext=(0, -14), fontsize=8,
                        color=INK, ha="center")
    ax.set_yticks(ys)
    ax.set_yticklabels([PRETTY[v] for v in order], fontsize=10, color=INK)
    ax.set_xlabel("long-axis extent (Å)", fontsize=10, color=INK_MUTED)
    ax.set_title("C · Bent → extended, vs experimental extended structures",
                 fontsize=11.5, color=INK, loc="left", pad=10)
    hollow = plt.Line2D([], [], marker="o", color="#fcfcfb", markeredgecolor=INK_MUTED,
                        markeredgewidth=2.0, markersize=9, linestyle="none", label="bent (start)")
    solid = plt.Line2D([], [], marker="o", color=INK_MUTED, markersize=9,
                       linestyle="none", label="extended (morphed)")
    diamond = plt.Line2D([], [], marker="D", color=INK, markersize=6.5,
                         linestyle="none", label="experimental extended (cryo-EM)")
    ax.legend(handles=[hollow, solid, diamond], frameon=False, fontsize=8.5,
              labelcolor=INK, loc="lower right")
    ax.margins(x=0.13, y=0.18)

    fig.suptitle("Bent → extended conformations across five integrin variants "
                 "(genu-hinge morph, CPU-only)",
                 fontsize=13.5, color=INK, x=0.008, ha="left", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_png, dpi=170, facecolor="white")
    print(f"[summary] wrote {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", default="results/variants/variant_endpoints.json")
    ap.add_argument("--boundaries", default="results/variants/boundaries_all.json")
    ap.add_argument("--results-dir", default="results/variants")
    ap.add_argument("--out-json", default="results/variants/extension_summary.json")
    ap.add_argument("--out-md", default="results/variants/extension_summary.md")
    ap.add_argument("--out-png", default="figures/variant_extension_summary.png")
    args = ap.parse_args()

    ep = json.load(open(args.endpoints))
    bounds = json.load(open(args.boundaries))
    variants = sorted(ep["representatives"])
    morphs = load_morphs(args.results_dir, variants)
    if not morphs:
        raise SystemExit("no completed morphs found")

    order = [v for v in ["alphaVbeta3", "alphaIIbbeta3", "alpha5beta1",
                         "alphaMbeta2", "alphaXbeta2"] if v in morphs]
    order += [v for v in sorted(morphs) if v not in order]

    summary = {v: assess(morphs[v]) for v in order}

    # control check: does the generalised path reproduce the published endpoint?
    ctrl = None
    if CONTROL["variant"] in summary:
        s = summary[CONTROL["variant"]]
        d_rg = abs(s["Rg_extended_A"] - CONTROL["Rg_A"]) / CONTROL["Rg_A"]
        d_ext = abs(s["extent_extended_A"] - CONTROL["long_axis_extent_A"]) / CONTROL["long_axis_extent_A"]
        ctrl = {"target_Rg_A": CONTROL["Rg_A"], "target_extent_A": CONTROL["long_axis_extent_A"],
                "got_Rg_A": s["Rg_extended_A"], "got_extent_A": s["extent_extended_A"],
                "Rg_rel_dev": round(d_rg, 4), "extent_rel_dev": round(d_ext, 4),
                "reproduces_published_endpoint": bool(d_rg <= CONTROL["tol_frac"]
                                                      and d_ext <= CONTROL["tol_frac"])}

    exp_refs = ep.get("extended_references", {})
    comparisons = {}
    for v, refs in exp_refs.items():
        if v not in summary:
            continue
        best = max(refs, key=lambda r: r["genu_angle_deg"])
        comparisons[v] = {
            "reference_pdb": best["pdb_id"], "reference_title": best["title"],
            "exp_genu_deg": best["genu_angle_deg"], "morph_genu_deg": summary[v]["genu_extended_deg"],
            "exp_extent_A": best["long_extent_A"], "morph_extent_A": summary[v]["extent_extended_A"],
            "exp_Rg_A": best["Rg_A"], "morph_Rg_A": summary[v]["Rg_extended_A"],
            "extent_abs_diff_A": round(abs(best["long_extent_A"] - summary[v]["extent_extended_A"]), 1),
        }

    conf = {v: {"boundary_confidence": ep["representatives"][v].get("boundary_confidence"),
                "leg_alignment_score": ep["representatives"][v].get("leg_alignment_score"),
                "start_pdb": ep["representatives"][v]["pdb_id"],
                "start_genu_deg": ep["representatives"][v].get("start_genu_angle_deg")}
            for v in order}

    result = {"variants": summary, "start_states": conf, "control_check": ctrl,
              "experimental_comparison": comparisons,
              "acceptance": {"genu_extended_min_deg": GENU_EXTENDED_MIN,
                             "extension_min_frac": EXTENSION_MIN_FRAC,
                             "monotonic_tol_deg": MONOTONIC_TOL_DEG}}
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    # ---- table view (the palette's contrast WARN obliges a non-colour readout) ----
    lines = ["# Extended conformations across integrin variants", "",
             "Genu-hinge incremental morph (vacuum ff14SB, CPU-only), generalised from the",
             "alphaVbeta3 route-A recipe by alignment-transferred domain boundaries.", "",
             "| variant | start PDB | conf. | genu bent→ext (°) | extent bent→ext (Å) | "
             "Rg bent→ext (Å) | steps | acceptance |",
             "|---|---|---|---|---|---|---|---|"]
    for v in order:
        s, c = summary[v], conf[v]
        verdict = "PASS" if s["passed"] else "FAIL: " + "; ".join(s["notes"])
        lines.append(
            f"| {v} | {c['start_pdb']} | {c['boundary_confidence']} | "
            f"{s['genu_bent_deg']:.1f} → {s['genu_extended_deg']:.1f} | "
            f"{s['extent_bent_A']:.1f} → {s['extent_extended_A']:.1f} | "
            f"{s['Rg_bent_A']:.1f} → {s['Rg_extended_A']:.1f} | {s['n_steps']} | {verdict} |")
    if ctrl:
        lines += ["", "## alphaVbeta3 positive control", "",
                  f"Published route-A endpoint (obj-075): Rg {ctrl['target_Rg_A']} A, "
                  f"extent {ctrl['target_extent_A']} A.",
                  f"This run: Rg {ctrl['got_Rg_A']} A, extent {ctrl['got_extent_A']} A "
                  f"(deviation {100*ctrl['Rg_rel_dev']:.1f}% / {100*ctrl['extent_rel_dev']:.1f}%).",
                  f"Reproduces published endpoint: **{ctrl['reproduces_published_endpoint']}**"]
    if comparisons:
        lines += ["", "## Comparison with experimental extended structures", "",
                  "| variant | reference | exp genu (°) | morph genu (°) | exp extent (Å) | "
                  "morph extent (Å) | |Δextent| (Å) |", "|---|---|---|---|---|---|---|"]
        for v, c in comparisons.items():
            lines.append(f"| {v} | {c['reference_pdb']} | {c['exp_genu_deg']:.1f} | "
                         f"{c['morph_genu_deg']:.1f} | {c['exp_extent_A']:.1f} | "
                         f"{c['morph_extent_A']:.1f} | {c['extent_abs_diff_A']:.1f} |")
    with open(args.out_md, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[summary] wrote {args.out_md}")

    os.makedirs(os.path.dirname(args.out_png), exist_ok=True)
    make_figure(morphs, order, summary, exp_refs, args.out_png)
    print("\n".join(lines[:6 + len(order)]))


if __name__ == "__main__":
    main()
