#!/usr/bin/env python3
"""Force-ramp readouts (obj-083): F½, P(extended | F), and κ = k_BT/var(θ).

Reads <dir>/<tag>/trajectory.csv as written by snapback_md.py in pull mode (columns
ps, knee_deg, Rg_A, f_pN, head_foot_A, <salt bridges>). Replicates are picked up automatically:
tag `wt_ramp` pools with `wt_ramp_r2`, `wt_ramp_r3`, ... .

  F½              force at which the knee stops holding      -- magnetic tweezers
  P(extended | F)  conformational distribution vs load
  κ = kT/var(θ)    knee compliance at fixed load             -- AFM/tweezers

F½ is reported PER REPLICATE and summarised as mean ± SD. That is the whole point of running
replicates: obj-082's pooled within-genotype SD was 12.8° of knee, so a single trace per genotype
cannot rank genotypes, and a single F½ with no error bar cannot either.

κ is computed ONLY from constant-force runs (--force-pn). Under a ramp, var(knee) inside a force
bin is dominated by the drift the ramp is driving, not by thermal fluctuation, so kT/var would be
a number with no physical meaning attached. Ramp files report kappa = null on purpose.

Absolute pN from a ramp is loading-rate-dependent and ~10 orders of magnitude off any tweezers
experiment; the transferable results are the genotype ORDERING and the gaps between genotypes.

  python analyze_force_ramp.py --dir results/route_a/force_ramp
  python analyze_force_ramp.py --demo     # self-check on synthetic traces, needs no MD data
"""
import argparse
import csv
import glob
import json
import os

import numpy as np

KB_KJ = 0.00831446261815324  # kJ/(mol*K)
GENOTYPES = {"wt_ramp": ("WT", "#2166ac"), "k459a_ramp": ("K459A", "#d6604d"),
             "e598a_ramp": ("E598A", "#1b7837"), "double_ramp": ("K459A/E598A", "#762a83")}


def load(tag_dir):
    """(force_pN, knee_deg) arrays, or None if this tag has no pull-mode trajectory."""
    path = os.path.join(tag_dir, "trajectory.csv")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "f_pN" not in rows[0]:
        return None  # zero-force snap-back run -- no force column, not a ramp
    return (np.array([float(r["f_pN"]) for r in rows]),
            np.array([float(r["knee_deg"]) for r in rows]),
            {c: np.array([float(r[c]) for r in rows]) for c in rows[0] if "-" in c})


def replicate_dirs(root, tag):
    """<root>/<tag> plus its _r2/_r3/... replicates, in order."""
    return [d for d in [os.path.join(root, tag)] + sorted(glob.glob(os.path.join(root, tag + "_r*")))
            if os.path.isdir(d)]


def bin_by_force(f_pn, knee_deg, bin_pn, extended_deg, temperature, with_kappa):
    """Bin frames by applied force; per bin report P(extended), <knee>, and (optionally) κ."""
    kt = KB_KJ * temperature
    edges = np.arange(0.0, f_pn.max() + bin_pn, bin_pn)
    if edges[-1] <= f_pn.max():
        # bins are [lo, hi), so a run held at exactly a bin edge -- every constant-force rung,
        # and every F = 0 rung -- would otherwise land in no bin at all and vanish silently.
        edges = np.append(edges, edges[-1] + bin_pn)
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (f_pn >= lo) & (f_pn < hi)
        if m.sum() < 2:  # variance needs >=2 frames; a lone frame gives no spread
            continue
        var = float(np.radians(knee_deg[m]).var(ddof=1))
        bins.append({"f_pN": round(float((lo + hi) / 2), 2), "n": int(m.sum()),
                     "knee_mean_deg": round(float(knee_deg[m].mean()), 2),
                     "knee_sd_deg": round(float(knee_deg[m].std(ddof=1)), 2),
                     "P_extended": round(float((knee_deg[m] >= extended_deg).mean()), 3),
                     # kJ/(mol*rad^2); x1.66054 -> pN*nm/rad^2
                     "kappa_kJ_per_mol_rad2": (round(kt / var, 2)
                                               if with_kappa and var > 0 else None)})
    return bins


def f_half(bins, p=0.5):
    """Force where P(extended) crosses `p` going UP: the load at which the knee holds.

    Linear interpolation between the two bracketing bin centres. None means the knee never
    reached P = p over the force range explored -- report that, do not extrapolate a number.
    """
    for a, b in zip(bins[:-1], bins[1:]):
        if a["P_extended"] < p <= b["P_extended"]:
            span = b["P_extended"] - a["P_extended"]
            frac = (p - a["P_extended"]) / span if span else 0.0
            return round(a["f_pN"] + frac * (b["f_pN"] - a["f_pN"]), 2)
    if bins and bins[0]["P_extended"] >= p:
        return 0.0  # already holding at the lowest force sampled
    return None


RUPTURE_A = 5.0  # a salt bridge past this is open (clasped is 2.6-2.9 A)


def rupture_forces(runs):
    """Force at which each genu salt bridge first opens, per replicate.

    This is the readout that survives a single ramp. The knee angle is confounded -- force and
    time move together, so anything that merely decays fakes an F-half -- but bridge distances
    are direct atom-pair distances, so they compare like-for-like against the F = 0 runs, which
    carry the same columns and are unaffected by the hinge definition. None = never opened.
    """
    out = {}
    for name in runs[0][2]:
        vals = []
        for f, _, br in runs:
            open_at = np.where(br[name] > RUPTURE_A)[0]
            vals.append(round(float(f[open_at[0]]), 1) if len(open_at) else None)
        out[name] = vals
    return out


def analyse(runs, bin_pn, extended_deg, temperature):
    """runs = [(force, knee, bridges), ...] for one genotype's replicates."""
    # A run whose force never changes is a constant-load rung; only those give a real κ.
    with_kappa = all(float(f.max() - f.min()) < 1e-6 for f, _, _ in runs)
    per_rep = [f_half(bin_by_force(f, k, bin_pn, extended_deg, temperature, with_kappa))
               for f, k, _ in runs]
    pooled = bin_by_force(np.concatenate([f for f, _, _ in runs]),
                          np.concatenate([k for _, k, _ in runs]),
                          bin_pn, extended_deg, temperature, with_kappa)
    got = [x for x in per_rep if x is not None]
    return {"n_replicates": len(runs),
            "mode": "constant-force" if with_kappa else "ramp",
            "kappa_reported": with_kappa,
            "f_half_per_replicate_pN": per_rep,
            "f_half_mean_pN": round(float(np.mean(got)), 2) if got else None,
            "f_half_sd_pN": round(float(np.std(got, ddof=1)), 2) if len(got) > 1 else None,
            "f_half_pooled_pN": f_half(pooled),
            "bridge_rupture_pN": rupture_forces(runs),
            "extended_deg": extended_deg, "bin_pn": bin_pn, "temperature_K": temperature,
            "bins": pooled}


def plot(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axK, axP, axS) = plt.subplots(1, 3, figsize=(16.5, 5.0))
    any_kappa = False
    for tag, (label, colour) in GENOTYPES.items():
        r = results.get(tag)
        if not r:
            continue
        b = r["bins"]
        fs = [x["f_pN"] for x in b]
        n = r["n_replicates"]
        axK.errorbar(fs, [x["knee_mean_deg"] for x in b], yerr=[x["knee_sd_deg"] for x in b],
                     color=colour, lw=1.8, capsize=2, label=f"{label} (n={n})")
        axP.plot(fs, [x["P_extended"] for x in b], "-o", color=colour, ms=3, lw=1.8,
                 label=f"{label} (n={n})")
        for fh in [x for x in r["f_half_per_replicate_pN"] if x]:
            axP.axvline(fh, color=colour, ls=":", lw=1.0, alpha=0.7)
        k = [(x["f_pN"], x["kappa_kJ_per_mol_rad2"]) for x in b
             if x["kappa_kJ_per_mol_rad2"] is not None]
        if k:
            any_kappa = True
            axS.plot([x for x, _ in k], [y for _, y in k], "-o", color=colour, ms=3, lw=1.8,
                     label=label)

    axK.axhline(173.5, ls=":", color="gray", lw=1)
    axK.set_xlabel("applied force (pN)"); axK.set_ylabel("genu knee angle (deg)")
    axK.set_title("Knee vs load (mean ± SD per bin)", fontsize=11, weight="bold")
    axP.axhline(0.5, ls=":", color="gray", lw=1)
    axP.set_xlabel("applied force (pN)"); axP.set_ylabel("P(extended | F)")
    axP.set_title("Extension probability vs load (dotted = per-replicate F½)",
                  fontsize=11, weight="bold")
    axS.set_xlabel("applied force (pN)")
    axS.set_ylabel(r"$\kappa = k_BT/\mathrm{var}(\theta)$  (kJ/mol/rad²)")
    axS.set_title("Knee stiffness vs load", fontsize=11, weight="bold")
    if not any_kappa:
        axS.text(0.5, 0.5, "κ needs constant-force runs (--force-pn):\nunder a ramp var(θ) is "
                           "ramp drift, not\nthermal fluctuation",
                 transform=axS.transAxes, ha="center", va="center", fontsize=10, color="dimgray")
    for ax in (axK, axP, axS):
        ax.grid(alpha=0.25)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=9)
    fig.suptitle("αVβ3 genu lock under load: F½, P(extended | F), compliance",
                 fontsize=12.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=125)
    print(f"wrote {out}")


def demo():
    """Self-check: synthesise runs with a KNOWN transition and confirm F½ recovers it."""
    rng = np.random.default_rng(0)
    for truth in (15.0, 30.0):
        # ramp DOWN, like the production protocol: force is descending in time, F½ must not care
        f = np.linspace(60, 0, 600)
        knee = 120 + 55 / (1 + np.exp(-(f - truth) / 1.5)) + rng.normal(0, 4, f.size)
        r = analyse([(f, knee, {"X-Y": np.full(f.size, 2.8)})], bin_pn=2.0,
                    extended_deg=150.0, temperature=310.0)
        assert r["bridge_rupture_pN"]["X-Y"] == [None], "a bridge that never opens must be None"
        got = r["f_half_pooled_pN"]
        assert got is not None and abs(got - truth) < 3.0, f"F½ {got} != {truth}"
        assert r["kappa_reported"] is False, "a ramp must not report kappa"
    # replicates: mean/SD across runs, and a constant-force run DOES report kappa
    reps = [(np.linspace(60, 0, 600),
             120 + 55 / (1 + np.exp(-(np.linspace(60, 0, 600) - t) / 1.5)) + rng.normal(0, 4, 600),
             {"X-Y": np.where(np.linspace(60, 0, 600) > t, 2.8, 9.0)})
            for t in (14.0, 16.0, 18.0)]
    r = analyse(reps, bin_pn=2.0, extended_deg=150.0, temperature=310.0)
    assert r["n_replicates"] == 3 and abs(r["f_half_mean_pN"] - 16.0) < 2.0, r
    assert r["f_half_sd_pN"] > 0, "replicates must produce an error bar"
    # the bridge opens as the ramp passes each replicate's threshold, descending
    got = r["bridge_rupture_pN"]["X-Y"]
    assert all(abs(g - t) < 1.0 for g, t in zip(got, (14.0, 16.0, 18.0))), got
    const = (np.full(400, 20.0), 160 + rng.normal(0, 5, 400), {"X-Y": np.full(400, 2.8)})
    rc = analyse([const], bin_pn=5.0, extended_deg=150.0, temperature=310.0)
    assert rc["kappa_reported"] and rc["bins"][0]["kappa_kJ_per_mol_rad2"] > 0, rc
    # a knee that never holds reports None, not a fabricated number
    f = np.linspace(60, 0, 600)
    assert analyse([(f, np.full(600, 120.0) + rng.normal(0, 4, 600),
                     {"X-Y": np.full(600, 9.0)})], 5.0, 150.0, 310.0)["f_half_pooled_pN"] is None
    print("demo OK: F½ within 3 pN on down-ramps, error bar across replicates, "
          "κ only for constant force, never-holds -> None, bridge rupture force recovered")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/route_a/force_ramp")
    ap.add_argument("--out-json", default="results/route_a/force_ramp/readouts.json")
    ap.add_argument("--out-fig", default="figures/route_a_force_ramp.png")
    ap.add_argument("--bin-pn", type=float, default=5.0)
    ap.add_argument("--extended-deg", type=float, default=150.0,
                    help="knee angle at or above which a frame counts as extended")
    ap.add_argument("--temperature", type=float, default=310.0)
    ap.add_argument("--demo", action="store_true", help="run the self-check and exit")
    args = ap.parse_args()
    if args.demo:
        return demo()

    results = {}
    for tag, (label, _) in GENOTYPES.items():
        runs = [d for d in (load(p) for p in replicate_dirs(args.dir, tag)) if d]
        if not runs:
            continue
        results[tag] = analyse(runs, args.bin_pn, args.extended_deg, args.temperature)
        r = results[tag]
        mean, sd = r["f_half_mean_pN"], r["f_half_sd_pN"]
        val = "never holds" if mean is None else (f"{mean:.1f} ± {sd:.1f} pN" if sd is not None
                                                  else f"{mean:.1f} pN")
        print(f"{label:<12} F½ = {val:<18} (n={r['n_replicates']}, {r['mode']}, "
              f"{len(r['bins'])} force bins)")
    if not results:
        raise SystemExit(f"no pull-mode trajectories under {args.dir} "
                         f"(expected {'/'.join(GENOTYPES)}[_r2...]/trajectory.csv with an f_pN column)")
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    json.dump(results, open(args.out_json, "w"), indent=2)
    print(f"wrote {args.out_json}")
    plot(results, args.out_fig)


if __name__ == "__main__":
    main()
