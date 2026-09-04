#!/usr/bin/env python3
"""Step 0 of the calcium objection: do our four genu ion pairs exist OUTSIDE our own morph?

obj-078..083 all rest on four αV genu ion pairs that first appeared during OUR incremental
leg-swing morph, run in vacuum with the 6 structural Ca2+ deleted and no ionic screening.
Rebuilding that morph correctly (Ca2+ kept, GB-OBC2 + 150 mM) opened three of the four —
including both cross-knee pairs that carried the entire force-protection result.

Before funding a corrected GPU campaign, ask the two free questions:

  (a) do the pairs appear in DEPOSITED extended αVβ3 (8XEN, cryo-EM 3.2 A)?
  (b) do they appear in an INDEPENDENT published explicit-solvent force-clamp trajectory of
      full-length αVβ3 with metals and 150 mM salt (Kolasangiani 2025, Structure;
      github.com/tamarabidone/alphaV_vs_alphaIIB)?

If they appear there, the calcium objection collapses and obj-078..083 survives a referee.
If they do not, the thread is a null and should be reported as one.

Pair definitions are imported from linchpin_engagement_order so they cannot drift from
obj-080/081/083. Distance = min(anion carboxylate O ... Lys NZ), the same definition used
throughout the thread. Thresholds: 4 A = engaged (obj-080), 5 A = ruptured (obj-083).

Residue numbering is NOT assumed. Every structure is identity-checked (D457/K459/D595/
E598/E636/K688) across all chains and a window of offsets; a structure that cannot be
matched unambiguously is reported as such rather than silently mismeasured. This is the
lesson from the GENU['B'] mis-selection that made a 30 deg bend read as 15.3 deg.
"""
import argparse
import json
import os
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from linchpin_engagement_order import BRIDGES  # noqa: E402  single source of truth

# BRIDGES uses en-dashes; obj-082/083 tables use ASCII hyphens. Normalise for comparability.
PAIRS = [(n.replace("–", "-"), a, c) for n, a, c in BRIDGES]
PAIR_NAMES = [p[0] for p in PAIRS]

ENGAGED_A = 4.0   # obj-080 engagement criterion
RUPTURED_A = 5.0  # obj-083 rupture criterion

# Identity fingerprint: our mature-αV numbering -> expected residue. All six are chain-A αV.
FINGERPRINT = {457: "ASP", 459: "LYS", 595: "ASP", 598: "GLU", 636: "GLU", 688: "LYS"}

# 1JV2 bent crystal, measured in the novelty audit. The self-test target: if the measurement
# path is right, these come back. Recorded in memory/route-a-genu-lock-novelty.md.
SELFTEST_1JV2 = {"E598-K459": 22.95, "K459-E636": 25.42, "D457-K688": 19.36, "D595-K688": 7.39}

RCSB = "https://files.rcsb.org/download/{}.pdb"
KOLA_RAW = "https://raw.githubusercontent.com/tamarabidone/alphaV_vs_alphaIIB/main/{}"
KOLA_API = "https://api.github.com/repos/tamarabidone/alphaV_vs_alphaIIB/contents/{}"
KOLA_TOP = "proteinalphaVbeta3.pdb"
KOLA_TRAJ_DIRS = ["AlphaVBeta3/Extension_Trajectories", "AlphaVBeta3/Bending_Trajectories"]


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def fetch(url, dest):
    """Download once. Returns dest, or None if the resource is unavailable."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=120) as r, open(dest, "wb") as f:
            f.write(r.read())
        return dest
    except Exception as e:  # a missing reference must not kill the whole run
        print(f"  ! could not fetch {url}: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return None


def read_pdb(path):
    """-> (coords {(chain,resid,atom): xyz}, resnames {(chain,resid): resname})."""
    coords, resnames = {}, {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            ch = line[21]
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            name = line[12:16].strip()
            resnames[(ch, resid)] = line[17:20].strip()
            coords[(ch, resid, name)] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return coords, resnames


def resolve_numbering(resnames, offsets=range(-60, 210)):
    """Find every (chain, offset) where all six fingerprint residues match.

    Returns (chain, offset, candidates). offset is added to OUR numbering to get the
    file's numbering, so file_resid = our_resid + offset.
    """
    chains = sorted({c for c, _ in resnames})
    hits = []
    for ch in chains:
        for off in offsets:
            if all(resnames.get((ch, r + off)) == aa for r, aa in FINGERPRINT.items()):
                hits.append((ch, off))
    if len(hits) == 1:
        return hits[0][0], hits[0][1], hits
    return None, None, hits


def pair_distances(coords, chain, offset):
    """-> {pair_name: min carboxylate-O..NZ distance in A, or None if atoms absent}."""
    out = {}
    for name, (ar, aatoms), (cr, catoms) in PAIRS:
        apts = [coords[(chain, ar + offset, a)] for a in aatoms
                if (chain, ar + offset, a) in coords]
        cpts = [coords[(chain, cr + offset, a)] for a in catoms
                if (chain, cr + offset, a) in coords]
        out[name] = (float(min(np.linalg.norm(x - y) for x in apts for y in cpts))
                     if apts and cpts else None)
    return out


def measure_structure(path, label):
    coords, resnames = read_pdb(path)
    chain, offset, hits = resolve_numbering(resnames)
    if chain is None:
        return {"label": label, "path": path, "error":
                f"numbering unresolved ({len(hits)} candidate chain/offset matches: {hits[:6]})"}
    d = pair_distances(coords, chain, offset)
    return {"label": label, "path": path, "chain": chain, "offset": offset,
            "distances_A": d,
            "engaged": {k: (None if v is None else v < ENGAGED_A) for k, v in d.items()}}


def measure_trajectory(top, traj, label, stride=1):
    """Per-frame pair distances over an xtc. Streams frames; one frame in memory at a time."""
    import MDAnalysis as mda

    u = mda.Universe(top, traj)
    # Identity-check against the topology, using MDAnalysis' own tables so a PDB that the
    # plain parser and MDA disagree about cannot slip through.
    ch = getattr(u.atoms, "chainIDs", None)
    resnames = {}
    for i, at in enumerate(u.atoms):
        key = ((ch[i] if ch is not None else " "), int(at.resid))
        resnames[key] = at.resname
    chain, offset, hits = resolve_numbering(resnames)
    if chain is None:
        return {"label": label, "topology": top, "trajectory": traj, "error":
                f"numbering unresolved ({len(hits)} candidates: {hits[:6]})"}

    idx = {}
    for name, (ar, aatoms), (cr, catoms) in PAIRS:
        def sel(resid, names):
            return [i for i, at in enumerate(u.atoms)
                    if int(at.resid) == resid + offset and at.name in names
                    and (ch is None or ch[i] == chain)]
        idx[name] = (sel(ar, aatoms), sel(cr, catoms))

    series = {n: [] for n in PAIR_NAMES}
    nframes = 0
    for ts in u.trajectory[::stride]:
        pos = ts.positions
        for name, (ai, ci) in idx.items():
            if not ai or not ci:
                series[name].append(np.nan)
                continue
            series[name].append(float(np.min(np.linalg.norm(
                pos[ai][:, None, :] - pos[ci][None, :, :], axis=-1))))
        nframes += 1

    stats = {}
    for n in PAIR_NAMES:
        a = np.array(series[n], dtype=float)
        ok = a[~np.isnan(a)]
        stats[n] = ({"n": 0} if ok.size == 0 else {
            "n": int(ok.size),
            "median_A": round(float(np.median(ok)), 2),
            "min_A": round(float(ok.min()), 2),
            "p05_A": round(float(np.percentile(ok, 5)), 2),
            "frac_engaged_lt4A": round(float((ok < ENGAGED_A).mean()), 3),
            "frac_intact_lt5A": round(float((ok < RUPTURED_A).mean()), 3)})
    return {"label": label, "topology": top, "trajectory": traj, "chain": chain,
            "offset": offset, "frames": nframes, "stats": stats,
            "series": {n: [None if np.isnan(x) else round(x, 2) for x in series[n]]
                       for n in PAIR_NAMES}}


def list_github_dir(path):
    try:
        with urllib.request.urlopen(KOLA_API.format(path), timeout=60) as r:
            return [e["name"] for e in json.load(r) if e["type"] == "file"]
    except Exception as e:
        print(f"  ! could not list {path}: {e}")
        return []


def selftest(cache):
    """One runnable check: the measurement path must reproduce the recorded 1JV2 values."""
    p = fetch(RCSB.format("1JV2"), os.path.join(cache, "1JV2.pdb"))
    assert p, "self-test needs 1JV2 from RCSB"
    res = measure_structure(p, "1JV2 bent crystal")
    assert "error" not in res, res
    assert res["offset"] == 0, f"1JV2 should be identity-numbered, got offset {res['offset']}"
    for name, want in SELFTEST_1JV2.items():
        got = res["distances_A"][name]
        assert got is not None, f"{name} not measurable in 1JV2"
        assert abs(got - want) < 0.5, f"{name}: got {got:.2f} A, audit recorded {want} A"
    d = res["distances_A"]
    print("selftest OK — 1JV2 reproduces all four recorded distances: "
          + ", ".join(f"{k} {d[k]:.2f} Å" for k in SELFTEST_1JV2))
    return res


def plot(report, out):
    """Left: every static reference, one group per pair. Right: the published
    trajectories pooled into Extension vs Bending per pair — 8 boxes, not 28.
    Byte-identical inputs are dropped so a duplicated file cannot widen a box."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    statics = [r for r in report["structures"] if "error" not in r]
    seen, trajs = set(), []
    for r in report["trajectories"]:
        if "stats" not in r or r.get("sha256_16") in seen:
            continue
        seen.add(r.get("sha256_16"))
        trajs.append(r)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8), gridspec_kw={"width_ratios": [1.25, 1]})

    ax = axes[0]
    x = np.arange(len(PAIR_NAMES))
    w = 0.8 / max(len(statics), 1)
    for i, r in enumerate(statics):
        vals = [r["distances_A"][n] or np.nan for n in PAIR_NAMES]
        bars = ax.bar(x + i * w - 0.4 + w / 2, vals, w, label=r["label"])
        for b, v in zip(bars, vals):
            if v == v and v < ENGAGED_A:  # annotate only what is actually engaged
                ax.text(b.get_x() + b.get_width() / 2, v + 0.4, f"{v:.1f}",
                        ha="center", fontsize=7, weight="bold")
    ax.axhline(ENGAGED_A, ls="--", color="k", lw=1)
    ax.axhline(RUPTURED_A, ls=":", color="crimson", lw=1)
    ax.text(-0.46, ENGAGED_A - 1.3, "engaged (4 Å)", fontsize=8)
    ax.text(-0.46, RUPTURED_A + 0.5, "ruptured (5 Å)", fontsize=8, color="crimson")
    ax.set_xticks(x); ax.set_xticklabels(PAIR_NAMES, fontsize=10)
    ax.set_ylabel("min carboxylate-O ··· Lys-NZ distance (Å)")
    ax.set_title("Only our own calcium-free morph makes these contacts",
                 fontsize=11, weight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.25, axis="y")

    ax = axes[1]
    groups = {"Extension": "#2166ac", "Bending": "#d6604d"}
    data, pos, colors = [], [], []
    for j, n in enumerate(PAIR_NAMES):
        for k, (g, c) in enumerate(groups.items()):
            pooled = [v for r in trajs if r["label"].startswith(g)
                      for v in r["series"][n] if v is not None]
            if pooled:
                data.append(pooled); pos.append(j + (k - 0.5) * 0.34); colors.append(c)
    if data:
        bp = ax.boxplot(data, positions=pos, widths=0.3, showfliers=False, patch_artist=True)
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.6)
        for med in bp["medians"]:
            med.set_color("k")
    ax.axhline(ENGAGED_A, ls="--", color="k", lw=1)
    ax.axhline(RUPTURED_A, ls=":", color="crimson", lw=1)
    ax.set_xticks(np.arange(len(PAIR_NAMES)))
    ax.set_xticklabels(PAIR_NAMES, fontsize=10)
    ax.set_xlim(-0.6, len(PAIR_NAMES) - 0.4)
    ax.set_ylabel("distance (Å)")
    ax.set_ylim(0, None)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.6) for c in groups.values()]
    ax.legend(handles, [f"{g} replicas" for g in groups], fontsize=8, loc="lower right")
    ax.set_title(f"Independent explicit-solvent αVβ3 under load\n"
                 f"(Kolasangiani 2025, metals + 150 mM, {len(trajs)} distinct trajectories)",
                 fontsize=11, weight="bold")
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle("Step 0 — the four αV genu ion pairs do not exist outside our own "
                 "calcium-free morph", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default="data/external/step0")
    ap.add_argument("--out-json", default="results/route_a/step0_external_pairs.json")
    ap.add_argument("--out-fig", default="figures/route_a_step0_external_pairs.png")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--selftest", action="store_true", help="only run the 1JV2 check")
    ap.add_argument("--replot", action="store_true",
                    help="rebuild the figure from --out-json; no network, no trajectories")
    a = ap.parse_args()

    if a.replot:
        plot(json.load(open(a.out_json)), a.out_fig)
        return

    os.makedirs(a.cache, exist_ok=True)
    if a.selftest:
        selftest(a.cache)
        return

    report = {"pairs": PAIR_NAMES, "engaged_A": ENGAGED_A, "ruptured_A": RUPTURED_A,
              "structures": [], "trajectories": []}

    print("== self-test ==")
    report["structures"].append(selftest(a.cache))

    print("== deposited references ==")
    for pdb, label in [("8XEN", "8XEN extended αVβ3 cryo-EM 3.2 Å"),
                       ("6DJP", "6DJP αVβ8 cryo-EM 4.8 Å (weak: side chains unresolved)")]:
        p = fetch(RCSB.format(pdb), os.path.join(a.cache, f"{pdb}.pdb"))
        if p:
            r = measure_structure(p, label)
            report["structures"].append(r)
            print(f"  {label}: {r.get('distances_A', r.get('error'))}")

    print("== our own structures ==")
    for path, label in [("results/route_a/extended_state_b.pdb", "ours: vacuum morph, no Ca²⁺"),
                        ("results/route_a/extended_state_ca.pdb", "ours: Ca²⁺ + 150 mM rebuild"),
                        ("results/route_a/extended_seed.pdb", "ours: pre-relax seed")]:
        if os.path.exists(path):
            r = measure_structure(path, label)
            report["structures"].append(r)
            print(f"  {label}: {r.get('distances_A', r.get('error'))}")

    print("== Kolasangiani 2025 published trajectories ==")
    top = fetch(KOLA_RAW.format(KOLA_TOP), os.path.join(a.cache, KOLA_TOP))
    if not top:
        print("  ! topology unavailable; skipping trajectory arm")
    else:
        for d in KOLA_TRAJ_DIRS:
            kind = d.split("/")[-1].replace("_Trajectories", "")
            for fn in list_github_dir(d):
                if not fn.endswith(".xtc"):
                    continue
                # Cache under the SOURCE DIRECTORY. Both trajectory dirs contain a
                # "Replica1.xtc"; keying the cache on the basename alone silently
                # served the Extension file for the Bending arm and inflated n.
                p = fetch(KOLA_RAW.format(f"{d}/{fn}"),
                          os.path.join(a.cache, d.replace("/", "_"), fn))
                if not p:
                    continue
                try:
                    r = measure_trajectory(top, p, f"{kind} {fn}", stride=a.stride)
                except Exception as e:
                    r = {"label": f"{kind} {fn}", "error": f"{type(e).__name__}: {e}"}
                r["sha256_16"] = sha256(p)
                report["trajectories"].append(r)
                print(f"  {kind} {fn}: {r.get('stats', r.get('error'))}")

    os.makedirs(os.path.dirname(a.out_json) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(a.out_fig) or ".", exist_ok=True)
    plot(report, a.out_fig)

    # Verdict, stated plainly so the log answers the question without re-reading the JSON.
    print("\n== VERDICT ==")
    ext = next((r for r in report["structures"] if r.get("label", "").startswith("8XEN")), None)
    if ext and "error" not in ext:
        eng = [n for n, v in ext["engaged"].items() if v]
        print(f"8XEN (deposited extended αVβ3): {len(eng)}/4 pairs engaged <4 Å -> {eng or 'none'}")
    tstats = [r for r in report["trajectories"] if "stats" in r]
    seen, uniq = set(), []
    for r in tstats:  # identical files must not be counted as independent replicas
        if r.get("sha256_16") in seen:
            print(f"  (duplicate input skipped: {r['label']} == "
                  f"{next(x['label'] for x in uniq if x['sha256_16'] == r['sha256_16'])})")
            continue
        seen.add(r.get("sha256_16"))
        uniq.append(r)
    if uniq:
        for n in PAIR_NAMES:
            fr = [r["stats"][n].get("frac_intact_lt5A") for r in uniq
                  if r["stats"][n].get("n")]
            fr = [x for x in fr if x is not None]
            if fr:
                print(f"  {n}: intact <5 Å in {np.mean(fr)*100:.1f}% of frames "
                      f"across {len(fr)} DISTINCT published trajectories")
    report["distinct_trajectories"] = len(uniq)

    with open(a.out_json, "w") as f:
        json.dump(report, f, indent=2)
    print(f"wrote {a.out_json}")


if __name__ == "__main__":
    main()
