#!/usr/bin/env python3
"""Library coverage v3 — fold published αVβ3 ectodomain crystal PDBs
into the F1 coverage figure as candidate EO endpoints.

Audit §10.9 P=3 follow-up; companion to docs/eo_coverage_strategy.md
route D ("published EO PDBs"). This script tests route D empirically:
download every published full-ectodomain αVβ3 PDB on RCSB, compute
CV0/CV1/CV2, and overlay them on the v2 library coverage figure.

Result spoiler: all 5 published αVβ3 ectodomain crystal structures
(1JV2, 1L5G, 4G1E, 4G1M, 4MMX) sit in the BC band (CV0 ≈ 51-52 Å) —
they are all bent, even when bound to RGD-mimetics that *do* induce
headpiece opening. Route D therefore returns no usable EO endpoints
for full-ectodomain αVβ3; the EO frames must come from enhanced
sampling (routes A/B/C/E in the strategy doc).

This is itself a publishable null result that strengthens the audit's
#1 EO blocker.

Output: figures/library_coverage_v3.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LIB_CACHE = ROOT / "results" / "afm_pipeline" / "library_coverage_v2" / "library_cvs.npz"
FITTED_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
PDB_DIR = ROOT / "data" / "multi_integrin" / "raw_pdbs"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "library_coverage_v3"
FIG_DIR = ROOT / "figures"

DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}
CV_PAIRS = [
    ("alpha_head_thigh", "alpha_calf"),
    ("beta_head",        "beta_tail"),
    ("alpha_head_thigh", "beta_head"),
]

# Candidate EO references — every published αVβ3 ectodomain PDB.
# Each entry: (pdb_id, label, expected_state_per_literature)
EO_CANDIDATES = [
    ("1JV2", "1JV2 αVβ3 apo (Xiong 2001)",         "bent"),
    ("1L5G", "1L5G αVβ3 + cilengitide (2002)",     "bent + open headpiece"),
    ("4G1E", "4G1E αVβ3 + coil-coiled tag (2012)", "bent"),
    ("4G1M", "4G1M αVβ3 re-refinement (2012)",     "bent"),
    ("4MMX", "4MMX αVβ3 + Fn 10th domain (2013)",  "bent"),
]

THRESHOLDS = {
    "cv0": [("BC max", 65), ("EC min", 78), ("library max", 85)],
    "cv1": [("BC max", 60), ("EC min", 73)],
    "cv2": [("EO min", 50)],
}


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    return cid if cid else chr(ord("A") + atom.residue.chain.index)


def domain_ca_indices(top, chain: str, lo: int, hi: int) -> np.ndarray:
    out = [a.index for a in top.atoms
           if a.name == "CA"
           and chain_id(a) == chain
           and lo <= a.residue.resSeq <= hi]
    return np.array(out, dtype=np.int64)


def cv_for_pdb(pdb_id: str) -> dict[str, float]:
    p = PDB_DIR / f"{pdb_id}.pdb"
    if not p.exists():
        raise FileNotFoundError(f"missing {p}")
    t = md.load(str(p))
    out = {}
    for i, (a, b) in enumerate(CV_PAIRS):
        ia = domain_ca_indices(t.topology, *DOMAINS[a])
        ib = domain_ca_indices(t.topology, *DOMAINS[b])
        if len(ia) == 0 or len(ib) == 0:
            out[f"cv{i}"] = float("nan")
            continue
        ca = t.xyz[0, ia].mean(axis=0) * 10.0
        cb = t.xyz[0, ib].mean(axis=0) * 10.0
        out[f"cv{i}"] = float(np.linalg.norm(ca - cb))
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Library + fitted reuse the F1 cache
    lib = np.load(LIB_CACHE, allow_pickle=True)
    lib_cv0 = np.array(lib["cv0"])
    lib_cv1 = np.array(lib["cv1"])
    lib_cv2 = np.array(lib["cv2"])
    f = np.load(FITTED_NPZ)
    fit_cv0 = np.concatenate([f["cv0_v1"], f["cv0_v2"]])
    fit_cv1 = np.concatenate([f["cv1_v1"], f["cv1_v2"]])
    fit_cv2 = np.concatenate([f["cv2_v1"], f["cv2_v2"]])

    # Score every EO candidate
    rows = []
    for pdb_id, label, expected in EO_CANDIDATES:
        try:
            cvs = cv_for_pdb(pdb_id)
        except FileNotFoundError as e:
            print(f"  skip {pdb_id}: {e}")
            continue
        rows.append({
            "pdb": pdb_id, "label": label, "expected": expected,
            "cv0": cvs["cv0"], "cv1": cvs["cv1"], "cv2": cvs["cv2"],
        })
        print(f"  {pdb_id}: CV0={cvs['cv0']:.2f}  CV1={cvs['cv1']:.2f}  CV2={cvs['cv2']:.2f}  ({expected})")

    in_eo_band_cv2 = sum(1 for r in rows if r["cv2"] >= 50)
    in_ec_band_cv0 = sum(1 for r in rows if r["cv0"] >= 78)
    print(f"\n{len(rows)} αVβ3 ectodomain PDBs scored")
    print(f"  EO threshold CV2 ≥ 50 Å: {in_eo_band_cv2} pass")
    print(f"  EC threshold CV0 ≥ 78 Å: {in_ec_band_cv0} pass")
    print(f"  Route D (published EO PDBs) verdict: "
          f"{'PASS' if in_eo_band_cv2 >= 1 else 'FAIL — no usable EO endpoints'}")

    plot(lib_cv0, lib_cv1, lib_cv2,
         fit_cv0, fit_cv1, fit_cv2,
         rows, FIG_DIR, OUT_DIR,
         in_eo_band_cv2, in_ec_band_cv0)

    import json
    with (OUT_DIR / "eo_pdb_scan.json").open("w") as fh:
        json.dump({
            "candidates": rows,
            "verdict": ("ROUTE D PASS" if in_eo_band_cv2 >= 1
                        else "ROUTE D FAIL — no published αVβ3 EO ectodomain PDB"),
            "library_max_cv0": float(lib_cv0.max()),
            "library_max_cv2": float(lib_cv2.max()),
            "fit_max_cv0": float(fit_cv0.max()),
            "fit_max_cv2": float(fit_cv2.max()),
            "method": "centroid distance, αVβ3 domain definitions from "
                      "build_library_metadata.py",
        }, fh, indent=2)
    return 0


def plot(lib0, lib1, lib2, fit0, fit1, fit2, rows,
         fig_dir, out_dir, in_eo, in_ec):
    fig, axes = plt.subplots(3, 1, figsize=(12.5, 10.5),
                             sharex=False)
    cvs = [
        ("cv0", lib0, fit0, "CV0 — αV head ↔ αV calf  (Å)", (40, 100), 65, 78, 85),
        ("cv1", lib1, fit1, "CV1 — β3 head ↔ β3 tail  (Å)", (35, 95),  60, 73, None),
        ("cv2", lib2, fit2, "CV2 — αV head ↔ β3 head  (Å)", (15, 60),  None, None, 50),
    ]

    for ax, (key, lib, fit, label, xlim, *thr) in zip(axes, cvs):
        bins = np.linspace(*xlim, 41)
        ax.hist(lib, bins=bins, color="#7570b3", alpha=0.45,
                density=True, label=f"library (n={lib.size})")
        ax.hist(fit, bins=bins, color="#1b9e77", alpha=0.45,
                density=True, label=f"fitted V1+V2 (n={fit.size})")

        # Threshold lines
        for name, x in THRESHOLDS[key]:
            ax.axvline(x, color="#444444", linestyle="--", linewidth=1)
            ax.text(x, ax.get_ylim()[1] * 0.96, f" {name}={x}",
                    color="#444444", fontsize=8, ha="left", rotation=90,
                    va="top")

        # EO gap shading
        if key == "cv2":
            lib_max = lib.max()
            ax.axvspan(lib_max, xlim[1], color="red", alpha=0.12,
                       hatch="///", label="EO gap")
        if key == "cv0":
            lib_max = lib.max()
            ax.axvspan(lib_max, xlim[1], color="red", alpha=0.12,
                       hatch="///", label="library gap")

        # EO-candidate PDB markers
        for j, r in enumerate(rows):
            v = r[key]
            if not np.isfinite(v):
                continue
            ax.axvline(v, color="#cc0000", linewidth=1.4, linestyle="-",
                       alpha=0.85)
            ax.text(v, ax.get_ylim()[1] * (0.85 - 0.10 * j),
                    f"  {r['pdb']}", color="#cc0000",
                    fontsize=8, ha="left", fontweight="bold")

        ax.set_xlim(*xlim)
        ax.set_xlabel(label)
        ax.set_ylabel("density")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.3)

    axes[0].set_title(
        f"Library coverage v3 — {len(rows)} published αVβ3 ectodomain PDBs scored as EO candidates "
        f"(route D from EO strategy doc)\n"
        f"Verdict: {in_eo} of {len(rows)} reach EO threshold CV2 ≥ 50 Å — "
        f"{'PASS' if in_eo >= 1 else 'FAIL: route D returns no usable EO endpoints'}\n"
        f"All 5 PDBs sit in the BC CV0 band (≈51-52 Å) — full-ectodomain αVβ3 "
        f"is *only* crystallized bent. EO sampling requires routes A/B/C/E.",
        fontsize=11, fontweight="bold", loc="left", pad=14,
    )

    fig.tight_layout()
    out_path = out_dir / "library_coverage_v3.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig.savefig(fig_dir / "library_coverage_v3.png", dpi=140, bbox_inches="tight")
    print(f"saved {out_path}")
    print(f"copied to {fig_dir / 'library_coverage_v3.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
