#!/usr/bin/env python3
"""Cryptic-binding-pocket scan over the v7 615-frame αVβ3 library.

Reviewer-E follow-up (audit-2026-05-05 §13.5). The earlier queue entry
was auto-marked completed via obj-046 fallout, but no actual analysis
exists on disk. This script delivers the genuine scan: residue-level
SASA differential between bent and extended library frames, with spatial
clustering into candidate cryptic-pocket sites + a bounded null result
when none is found.

Approach:
  1. Reuse the CV0 cache from obj-039 to pick the N most-bent + N
     most-extended frames (default N=20 per state).
  2. For each picked PDB compute Shrake-Rupley residue-level SASA
     (1654, ) plus the residue Cα coordinate.
  3. Per-residue: mean SASA in bent vs extended, ΔSASA, two-sample
     Welch t-statistic. Flag "cryptic-opening" residues with |ΔSASA|
     > τ (default τ=20 Å² to be sensitive without overflagging).
  4. Spatial clustering: agglomerative-style merge of opening
     residues whose Cα separation is < 10 Å into candidate pockets.
     Filter to pockets with ≥ 3 residues and aggregate ΔSASA > 100 Å²
     (the success-criterion threshold from queue.yaml).
  5. Sensitivity bound: report what aggregate ΔSASA the design
     could resolve at α=0.05 given the bent/extended noise σ.
  6. Save pocket_summary.json + 4-panel figure.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LIB_DIR = ROOT / "data" / "runs" / "avb3" / "conformers" / "all_frames_bent_extended"
CV0_CACHE = ROOT / "results" / "afm_pipeline" / "rgd_docking" / "library_cv0_cache.npy.npz"
OUT_DIR = ROOT / "results" / "cryptic_pockets"
FIG = ROOT / "figures" / "cryptic_pockets_v1.png"

N_PER_STATE = 20
DSASA_RESIDUE_THRESH = 20.0
POCKET_DSASA_THRESH = 100.0
POCKET_MIN_RES = 3
POCKET_MAX_RES = 25
SPATIAL_CUTOFF_A = 7.0
SPATIAL_RESPLIT_FLOOR_A = 4.0

DOMAINS = [
    ("alpha_propeller_thigh", "A", 1, 435),
    ("alpha_calf", "A", 436, 741),
    ("alpha_tail", "A", 742, 956),
    ("beta_psi_hybrid_betaA", "B", 1, 352),
    ("beta_egf_betatail", "B", 353, 692),
]
HYDROPHOBIC = set("VILMFYWA")
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    if cid:
        return cid
    return chr(ord("A") + atom.residue.chain.index)


def residue_ca_coords(traj: md.Trajectory) -> tuple[np.ndarray, list[tuple[str, int, str]]]:
    topology = traj.topology
    assert topology is not None
    coords = np.full((traj.n_residues, 3), np.nan, dtype=np.float64)
    labels: list[tuple[str, int, str]] = []
    res_atom_lookup: dict[int, int] = {}
    for atom in topology.atoms:
        if atom.name == "CA":
            res_atom_lookup.setdefault(atom.residue.index, atom.index)
    xyz = traj.xyz
    assert xyz is not None
    for res in topology.residues:
        ai = res_atom_lookup.get(res.index)
        if ai is not None:
            coords[res.index] = xyz[0, ai] * 10.0
        first_atom = next(iter(res.atoms))
        labels.append((chain_id(first_atom), int(res.resSeq), res.name))
    return coords, labels


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG.parent.mkdir(parents=True, exist_ok=True)

    cache = np.load(CV0_CACHE, allow_pickle=True)
    names = list(cache["names"])
    cv0 = np.array(cache["cv0"], dtype=np.float64)
    sort_idx = np.argsort(cv0)
    bent_idx = sort_idx[:N_PER_STATE]
    ext_idx = sort_idx[-N_PER_STATE:]
    print(f"library frames {len(names)}; bent CV0 ≤ {cv0[bent_idx][-1]:.2f} Å; "
          f"extended CV0 ≥ {cv0[ext_idx][0]:.2f} Å")

    cache_file = OUT_DIR / "sasa_cache.npz"
    sasa_bent = sasa_ext = None
    cv0_bent = cv0_ext = None
    ref_coords = None
    residue_labels: list[tuple[str, int, str]] = []
    n_res = 0
    if cache_file.exists():
        c = np.load(cache_file, allow_pickle=True)
        if (int(c["n_per_state"]) == N_PER_STATE
                and np.array_equal(c["bent_idx"], bent_idx)
                and np.array_equal(c["ext_idx"], ext_idx)):
            sasa_bent = c["sasa_bent"]
            sasa_ext = c["sasa_ext"]
            cv0_bent = c["cv0_bent"]
            cv0_ext = c["cv0_ext"]
            ref_coords = c["ref_coords"]
            residue_labels = list(c["residue_labels"])
            n_res = sasa_bent.shape[1]
            print(f"  using cached SASA ({2*N_PER_STATE} frames × {n_res} residues)")
        else:
            cache_file.unlink()

    if sasa_bent is None:
        sample = md.load(str(LIB_DIR / names[bent_idx[0]]))
        n_res = sample.n_residues
        ref_coords, residue_labels = residue_ca_coords(sample)
        print(f"residues per frame: {n_res}")

        sasa_bent = np.zeros((N_PER_STATE, n_res))
        sasa_ext = np.zeros((N_PER_STATE, n_res))
        cv0_bent = np.zeros(N_PER_STATE)
        cv0_ext = np.zeros(N_PER_STATE)

        t0 = time.time()
        for k, idx in enumerate(bent_idx):
            traj = md.load(str(LIB_DIR / names[idx]))
            sasa_bent[k] = md.shrake_rupley(traj, mode="residue")[0] * 100.0
            cv0_bent[k] = cv0[idx]
            if (k + 1) % 5 == 0:
                print(f"  bent {k+1}/{N_PER_STATE} ({time.time()-t0:.1f}s)")
        for k, idx in enumerate(ext_idx):
            traj = md.load(str(LIB_DIR / names[idx]))
            sasa_ext[k] = md.shrake_rupley(traj, mode="residue")[0] * 100.0
            cv0_ext[k] = cv0[idx]
            if (k + 1) % 5 == 0:
                print(f"  extended {k+1}/{N_PER_STATE} ({time.time()-t0:.1f}s)")
        print(f"SASA scan: {time.time()-t0:.1f}s for {2*N_PER_STATE} frames × {n_res} residues")
        np.savez(cache_file,
                 n_per_state=N_PER_STATE,
                 bent_idx=bent_idx, ext_idx=ext_idx,
                 sasa_bent=sasa_bent, sasa_ext=sasa_ext,
                 cv0_bent=cv0_bent, cv0_ext=cv0_ext,
                 ref_coords=ref_coords,
                 residue_labels=np.array(residue_labels, dtype=object))

    assert sasa_bent is not None and sasa_ext is not None
    assert cv0_bent is not None and cv0_ext is not None
    assert ref_coords is not None

    mean_bent = sasa_bent.mean(axis=0)
    mean_ext = sasa_ext.mean(axis=0)
    delta = mean_ext - mean_bent
    sd_bent = sasa_bent.std(axis=0, ddof=1)
    sd_ext = sasa_ext.std(axis=0, ddof=1)
    pooled_se = np.sqrt(sd_bent ** 2 / N_PER_STATE + sd_ext ** 2 / N_PER_STATE)
    welch_t = np.divide(delta, pooled_se,
                        out=np.zeros_like(delta), where=pooled_se > 1e-9)

    sensitivity = {
        "median_per_residue_se_A2": float(np.median(pooled_se)),
        "p95_per_residue_se_A2": float(np.percentile(pooled_se, 95)),
        "min_detectable_single_residue_dsasa_A2_alpha05": 1.96 * float(np.median(pooled_se)),
        "min_detectable_pocket_dsasa_A2_3res": float(np.sqrt(3.0) * 1.96 * np.median(pooled_se)),
        "min_detectable_pocket_dsasa_A2_5res": float(np.sqrt(5.0) * 1.96 * np.median(pooled_se)),
    }

    opening_mask = delta > DSASA_RESIDUE_THRESH
    closing_mask = delta < -DSASA_RESIDUE_THRESH
    open_idx = np.where(opening_mask)[0]
    close_idx = np.where(closing_mask)[0]
    print(f"\nresidues with ΔSASA > {DSASA_RESIDUE_THRESH:+.0f} Å² "
          f"(opening on extension): {opening_mask.sum()}")
    print(f"residues with ΔSASA < -{DSASA_RESIDUE_THRESH:.0f} Å² "
          f"(closing on extension): {closing_mask.sum()}")

    pockets, opening_megas = cluster_pockets(
        open_idx, ref_coords, residue_labels, delta,
        spatial_cutoff_A=SPATIAL_CUTOFF_A,
        min_residues=POCKET_MIN_RES,
        aggregate_thresh=POCKET_DSASA_THRESH)
    print(f"\ncandidate opening pockets (3-{POCKET_MAX_RES} residues, "
          f"aggregate ΔSASA > {POCKET_DSASA_THRESH:.0f} Å²): {len(pockets)}")
    intra_count = sum(1 for p in pockets if p["intra_domain"])
    drug_count = sum(1 for p in pockets if p["druggable_candidate"])
    print(f"  of which intra-domain: {intra_count}; druggable candidates "
          f"(intra-domain, ≥40% hydrophobic, ≥5 residues): {drug_count}")
    for p in pockets[:12]:
        flag = ("DRUGGABLE" if p["druggable_candidate"]
                else "intra" if p["intra_domain"]
                else "interface")
        print(f"  pocket {p['pocket_id']} [{flag}]: {p['n_residues']} res "
              f"({p['primary_domain']}), agg ΔSASA = {p['aggregate_delta_sasa_A2']:.1f} Å² "
              f"hydroph={p['hydrophobic_fraction']:.0%}; top: {p['top_residues'][0]}")
    if opening_megas:
        print(f"  opening megaclusters (>{POCKET_MAX_RES} residues, bulk interface): "
              f"{len(opening_megas)}")
        for m in opening_megas[:3]:
            print(f"    mega: {m['n_residues']} res, agg ΔSASA = "
                  f"{m['aggregate_delta_sasa_A2']:.0f} Å² (head-leg interface)")

    closing_pockets, closing_megas = cluster_pockets(
        close_idx, ref_coords, residue_labels, -delta,
        spatial_cutoff_A=SPATIAL_CUTOFF_A,
        min_residues=POCKET_MIN_RES,
        aggregate_thresh=POCKET_DSASA_THRESH)
    print(f"\ncandidate closing pockets (3-{POCKET_MAX_RES} residues, |aggregate| > "
          f"{POCKET_DSASA_THRESH:.0f}): {len(closing_pockets)}")
    for p in closing_pockets[:5]:
        print(f"  pocket {p['pocket_id']}: {p['n_residues']} res, "
              f"agg |ΔSASA| = {p['aggregate_delta_sasa_A2']:.1f} Å² "
              f"({p['chain_a_count']} αV, {p['chain_b_count']} β3); "
              f"top: {p['top_residues'][0]}")

    summary = {
        "library_dir": str(LIB_DIR.relative_to(ROOT)),
        "n_frames_bent": N_PER_STATE,
        "n_frames_extended": N_PER_STATE,
        "mean_cv0_bent_A": float(cv0_bent.mean()),
        "mean_cv0_extended_A": float(cv0_ext.mean()),
        "n_residues_total": int(n_res),
        "n_residues_opening_dsasa_gt_20": int(opening_mask.sum()),
        "n_residues_closing_dsasa_lt_minus_20": int(closing_mask.sum()),
        "n_candidate_opening_pockets": len(pockets),
        "n_candidate_closing_pockets": len(closing_pockets),
        "n_intra_domain_opening_pockets": sum(1 for p in pockets if p["intra_domain"]),
        "n_druggable_opening_pockets": sum(1 for p in pockets if p["druggable_candidate"]),
        "opening_pockets": pockets,
        "closing_pockets": closing_pockets,
        "opening_megaclusters": opening_megas,
        "closing_megaclusters": closing_megas,
        "sensitivity_bound": sensitivity,
        "thresholds": {
            "dsasa_residue_A2": DSASA_RESIDUE_THRESH,
            "pocket_dsasa_A2": POCKET_DSASA_THRESH,
            "pocket_min_residues": POCKET_MIN_RES,
            "pocket_max_residues": POCKET_MAX_RES,
            "spatial_cutoff_A": SPATIAL_CUTOFF_A,
            "spatial_resplit_floor_A": SPATIAL_RESPLIT_FLOOR_A,
        },
        "interpretation": (
            f"With N={N_PER_STATE}+{N_PER_STATE} frames, the design resolves "
            f"single-residue ΔSASA ≥ {sensitivity['min_detectable_single_residue_dsasa_A2_alpha05']:.1f} Å² "
            f"and 5-residue aggregate ΔSASA ≥ {sensitivity['min_detectable_pocket_dsasa_A2_5res']:.1f} Å² "
            f"at α=0.05. "
            f"Found {len(pockets)} pocket(s) opening on extension; of these, "
            f"{sum(1 for p in pockets if p['intra_domain'])} are intra-domain "
            f"(rule out rigid-body inter-domain separation as the cause) and "
            f"{sum(1 for p in pockets if p['druggable_candidate'])} are druggable "
            f"candidates (intra-domain × ≥40% hydrophobic × ≥5 residues). "
            f"Bulk-interface megaclusters (which open trivially during BC→EC "
            f"head-leg separation) are reported separately and excluded from "
            f"the cryptic-binding-site count."
        ),
    }
    summary_path = OUT_DIR / "pocket_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWrote {summary_path}")

    np.savez(OUT_DIR / "residue_dsasa.npz",
             names=np.array(names),
             bent_idx=bent_idx,
             ext_idx=ext_idx,
             mean_bent=mean_bent,
             mean_ext=mean_ext,
             delta=delta,
             pooled_se=pooled_se,
             welch_t=welch_t,
             ref_coords=ref_coords,
             residue_labels=np.array(residue_labels, dtype=object))

    plot(summary, mean_bent, mean_ext, delta, welch_t,
         residue_labels, cv0_bent, cv0_ext)
    print(f"Wrote {FIG}")
    return 0


def _single_link_groups(member_arr: np.ndarray, coords: np.ndarray,
                        cutoff: float) -> list[list[int]]:
    """Single-linkage clustering: members within `cutoff` Å are connected."""
    n = len(member_arr)
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    pos = coords[member_arr]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(pos[i] - pos[j]) < cutoff:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[ra] = rb
    g: dict[int, list[int]] = {}
    for i in range(n):
        g.setdefault(find(i), []).append(int(member_arr[i]))
    return list(g.values())


def cluster_pockets(idx: np.ndarray, coords: np.ndarray,
                    labels: list[tuple[str, int, str]], delta: np.ndarray,
                    spatial_cutoff_A: float, min_residues: int,
                    aggregate_thresh: float,
                    max_residues: int = POCKET_MAX_RES,
                    resplit_floor_A: float = SPATIAL_RESPLIT_FLOOR_A
                    ) -> tuple[list[dict], list[dict]]:
    """Spatially cluster residues into pockets. Returns (pockets, megaclusters).

    Megaclusters (n_residues > max_residues) are reported separately as
    bulk-interface gains rather than discrete cryptic pockets — the head-leg
    interface that opens during BC→EC extension dominates by raw SASA gain
    but is not a binding-pocket candidate.
    """
    if len(idx) == 0:
        return [], []
    valid = ~np.isnan(coords[idx, 0])
    idx = idx[valid]
    if len(idx) == 0:
        return [], []

    queue = [(np.asarray(idx), spatial_cutoff_A)]
    final_clusters: list[list[int]] = []
    megaclusters: list[list[int]] = []
    while queue:
        member_arr, cutoff = queue.pop(0)
        groups = _single_link_groups(member_arr, coords, cutoff)
        for g in groups:
            if len(g) <= max_residues:
                final_clusters.append(g)
            elif cutoff > resplit_floor_A:
                queue.append((np.asarray(g), max(resplit_floor_A, cutoff - 1.5)))
            else:
                megaclusters.append(g)

    def _format(members: list[int], k: int) -> dict:
        chain_a = sum(1 for r in members if labels[r][0] == "A")
        chain_b = sum(1 for r in members if labels[r][0] == "B")
        sorted_m = sorted(members, key=lambda r: -float(np.abs(delta[r])))
        top = [
            f"{labels[r][0]}:{labels[r][2]}{labels[r][1]}  ΔSASA={delta[r]:+.1f} Å²"
            for r in sorted_m[:5]
        ]
        member_domains = []
        for r in members:
            ch, num, _ = labels[r]
            d = next((dn for dn, dch, lo, hi in DOMAINS
                      if dch == ch and lo <= num <= hi), "unknown")
            member_domains.append(d)
        domain_set = set(member_domains)
        intra_domain = len(domain_set) == 1
        domains_dict: dict[str, int] = {}
        for d in member_domains:
            domains_dict[d] = domains_dict.get(d, 0) + 1
        hydroph_count = 0
        aromatic_count = 0
        for r in members:
            three = labels[r][2]
            one = THREE_TO_ONE.get(three, "X")
            if one in HYDROPHOBIC:
                hydroph_count += 1
            if one in {"F", "Y", "W"}:
                aromatic_count += 1
        return {
            "pocket_id": k,
            "n_residues": len(members),
            "chain_a_count": chain_a,
            "chain_b_count": chain_b,
            "aggregate_delta_sasa_A2": float(np.abs(delta[members]).sum()),
            "centroid_A": [float(x) for x in coords[members].mean(axis=0).tolist()],
            "top_residues": top,
            "all_member_residues": [
                f"{labels[r][0]}:{labels[r][2]}{labels[r][1]}" for r in sorted_m
            ],
            "domain_counts": domains_dict,
            "intra_domain": intra_domain,
            "primary_domain": max(domains_dict, key=lambda d: domains_dict[d]),
            "hydrophobic_fraction": hydroph_count / len(members),
            "aromatic_count": aromatic_count,
            "druggable_candidate": (intra_domain
                                    and (hydroph_count / len(members) >= 0.4)
                                    and len(members) >= 5),
        }

    pockets = []
    for k, members in enumerate(sorted(final_clusters,
                                       key=lambda m: -float(np.abs(delta[m]).sum()))):
        if len(members) < min_residues:
            continue
        if float(np.abs(delta[members]).sum()) < aggregate_thresh:
            continue
        pockets.append(_format(members, k))

    megalist = []
    for k, members in enumerate(sorted(megaclusters,
                                       key=lambda m: -float(np.abs(delta[m]).sum()))):
        megalist.append(_format(members, k))

    return pockets, megalist


def plot(summary: dict, mean_bent: np.ndarray, mean_ext: np.ndarray,
         delta: np.ndarray, welch_t: np.ndarray,
         labels: list[tuple[str, int, str]], cv0_bent: np.ndarray,
         cv0_ext: np.ndarray) -> None:
    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95], hspace=0.36, wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(mean_bent, mean_ext, s=4, alpha=0.35, color="#7570b3")
    lim = max(mean_bent.max(), mean_ext.max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="#888", linewidth=0.8)
    open_mask = delta > DSASA_RESIDUE_THRESH
    close_mask = delta < -DSASA_RESIDUE_THRESH
    ax.scatter(mean_bent[open_mask], mean_ext[open_mask],
               s=12, color="#d62728", label=f"opening (n={int(open_mask.sum())})")
    ax.scatter(mean_bent[close_mask], mean_ext[close_mask],
               s=12, color="#1f77b4", label=f"closing (n={int(close_mask.sum())})")
    ax.set_xlabel("mean SASA bent (Å²)")
    ax.set_ylabel("mean SASA extended (Å²)")
    ax.set_title(f"Per-residue SASA: bent vs extended (CV0 {cv0_bent.mean():.1f} → {cv0_ext.mean():.1f} Å)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    bins = np.linspace(-150, 150, 61).tolist()
    ax.hist(delta, bins=bins, color="#cccccc", edgecolor="#666",
            label=f"all residues (n={len(delta)})")
    ax.axvline(DSASA_RESIDUE_THRESH, color="#d62728", linestyle="--", linewidth=1.0)
    ax.axvline(-DSASA_RESIDUE_THRESH, color="#1f77b4", linestyle="--", linewidth=1.0)
    ax.axvline(0, color="black", linewidth=0.6, alpha=0.5)
    ax.set_yscale("log")
    ax.set_xlabel("ΔSASA = SASA_extended − SASA_bent  (Å²)")
    ax.set_ylabel("residues  (log)")
    ax.set_title("Distribution of per-residue ΔSASA")
    sd = float(np.std(delta))
    ax.text(0.02, 0.95,
            f"mean Δ = {delta.mean():+.2f} Å²\nstd Δ  = {sd:.2f} Å²\n"
            f"|Δ|>{DSASA_RESIDUE_THRESH:.0f}: {(np.abs(delta) > DSASA_RESIDUE_THRESH).sum()} of {len(delta)}",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 2])
    res_idx = np.arange(len(delta))
    ax.scatter(res_idx, delta, s=2, alpha=0.5, color="#888")
    open_idx = res_idx[open_mask]
    close_idx = res_idx[close_mask]
    ax.scatter(open_idx, delta[open_idx], s=10, color="#d62728")
    ax.scatter(close_idx, delta[close_idx], s=10, color="#1f77b4")
    ax.axhline(DSASA_RESIDUE_THRESH, color="#d62728", linestyle="--", linewidth=0.8)
    ax.axhline(-DSASA_RESIDUE_THRESH, color="#1f77b4", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="black", linewidth=0.6, alpha=0.4)
    chain_b_start = next((i for i, l in enumerate(labels) if l[0] == "B"), len(labels))
    ax.axvline(chain_b_start, color="#666", linestyle=":", linewidth=1.0)
    ax.text(chain_b_start * 0.5, ax.get_ylim()[1] * 0.85, "αV chain", ha="center", fontsize=9)
    ax.text((chain_b_start + len(labels)) * 0.5, ax.get_ylim()[1] * 0.85, "β3 chain", ha="center", fontsize=9)
    ax.set_xlabel("residue index")
    ax.set_ylabel("ΔSASA (Å²)")
    ax.set_title("Per-residue ΔSASA along the chain")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 0])
    pockets = summary["opening_pockets"]
    if pockets:
        agg = np.array([p["aggregate_delta_sasa_A2"] for p in pockets])
        hyd = np.array([p["hydrophobic_fraction"] for p in pockets])
        intra = np.array([p["intra_domain"] for p in pockets])
        drug = np.array([p["druggable_candidate"] for p in pockets])
        colors = []
        for i in range(len(pockets)):
            if drug[i]:
                colors.append("#2ca02c")
            elif intra[i]:
                colors.append("#ff7f0e")
            else:
                colors.append("#7f7f7f")
        ax.scatter(agg, hyd * 100, c=colors, s=80, edgecolor="black", linewidth=0.6)
        for i, p in enumerate(pockets[:18]):
            ax.annotate(f"P{i}\n{p['n_residues']}r", (agg[i], hyd[i] * 100),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=7, color="#333")
        ax.axhline(40, color="#888", linestyle=":", linewidth=0.8,
                   label="40% hydroph (drug threshold)")
        ax.axvline(POCKET_DSASA_THRESH, color="#888", linestyle=":", linewidth=0.8)
        ax.set_xlabel("aggregate ΔSASA (Å²)")
        ax.set_ylabel("hydrophobic-residue fraction (%)")
        n_drug = int(drug.sum())
        n_intra = int(intra.sum())
        ax.set_title(f"Pockets: green=druggable ({n_drug}), orange=intra-domain non-drug "
                     f"({n_intra - n_drug}), grey=interface ({len(pockets) - n_intra})",
                     fontsize=9)
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5,
                "NULL RESULT: no candidate pocket meets\n"
                f"≥{POCKET_MIN_RES} residues × aggregate >{POCKET_DSASA_THRESH:.0f} Å²",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=12,
                bbox=dict(boxstyle="round", facecolor="#fff3cd", edgecolor="#856404"))
        ax.axis("off")

    ax = fig.add_subplot(gs[1, 1])
    sens = summary["sensitivity_bound"]
    bins = np.linspace(0, max(15, np.percentile(welch_t, 99.5)), 60).tolist()
    ax.hist(np.abs(welch_t), bins=bins, color="#7570b3", edgecolor="#444")
    ax.axvline(1.96, color="#d62728", linestyle="--", label="|t|=1.96 (α=0.05)")
    ax.axvline(2.58, color="#fd8d3c", linestyle="--", label="|t|=2.58 (α=0.01)")
    ax.set_xlabel("Welch |t| per residue")
    ax.set_ylabel("residues")
    ax.set_yscale("log")
    n_sig = int((np.abs(welch_t) > 1.96).sum())
    expected = len(welch_t) * 0.05
    ax.set_title(f"Welch-t distribution — {n_sig} sig at α=0.05 "
                 f"(expected by chance: {expected:.0f})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    rows = [
        ["frames per state", f"{N_PER_STATE}"],
        ["bent CV0 mean", f"{summary['mean_cv0_bent_A']:.1f} Å"],
        ["extended CV0 mean", f"{summary['mean_cv0_extended_A']:.1f} Å"],
        ["residues opening", f"{summary['n_residues_opening_dsasa_gt_20']}"],
        ["residues closing", f"{summary['n_residues_closing_dsasa_lt_minus_20']}"],
        ["opening pockets (all)", f"{summary['n_candidate_opening_pockets']}"],
        ["intra-domain opening", f"{summary['n_intra_domain_opening_pockets']}"],
        ["DRUGGABLE candidates", f"{summary['n_druggable_opening_pockets']}"],
        ["closing pockets", f"{summary['n_candidate_closing_pockets']}"],
        ["min detect. single-res", f"{sens['min_detectable_single_residue_dsasa_A2_alpha05']:.1f} Å²"],
        ["min detect. 5-res pocket", f"{sens['min_detectable_pocket_dsasa_A2_5res']:.1f} Å²"],
    ]
    cell_text = [[k, v] for k, v in rows]
    table = ax.table(cellText=cell_text, colLabels=["metric", "value"],
                     loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)
    ax.set_title("Sensitivity + summary stats")

    fig.suptitle(
        "Cryptic-pocket scan v7 αVβ3 library — "
        f"{summary['n_candidate_opening_pockets']} opening pocket(s) → "
        f"{summary['n_intra_domain_opening_pockets']} intra-domain → "
        f"{summary['n_druggable_opening_pockets']} druggable cryptic candidate(s) "
        f"({2*N_PER_STATE} frames × {summary['n_residues_total']} residues)",
        fontsize=11, fontweight="bold", y=0.995,
    )
    fig.savefig(FIG, dpi=140, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
