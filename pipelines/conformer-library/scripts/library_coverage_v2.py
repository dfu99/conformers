#!/usr/bin/env python3
"""Library coverage v2 — CV0/CV1/CV2 distributions across the 615-frame
v7 library and the V1+V2 fitted trajectories, side by side.

Audit follow-up F1 (audit-2026-05-05 §8). Makes the EO coverage gap
visually unambiguous and shows that the fitted trajectory CVs sit
*inside* the library convex hull on CV0 + CV1 but escape it on CV2.

Output: figures/library_coverage_v2.png plus per-CV array dumps in
results/afm_pipeline/library_coverage_v2/.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LIB_DIR = ROOT / "data" / "runs" / "avb3" / "conformers" / "all_frames_bent_extended"
FITTED_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "cv_series.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "library_coverage_v2"
FIG_DIR = ROOT / "figures"

DOMAINS = {
    "alpha_head_thigh": ("A", 1, 435),
    "alpha_calf":       ("A", 436, 741),
    "beta_head":        ("B", 1, 352),
    "beta_tail":        ("B", 353, 692),
}

CV_PAIRS = [
    ("alpha_head_thigh", "alpha_calf"),     # CV0
    ("beta_head",        "beta_tail"),      # CV1
    ("alpha_head_thigh", "beta_head"),      # CV2
]

CV_LABELS = {
    "cv0": "CV0 — αV head ↔ αV calf  (Å)",
    "cv1": "CV1 — β3 head ↔ β3 tail  (Å)",
    "cv2": "CV2 — αV head ↔ β3 head  (Å)",
}

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


def cv_for_pdb(traj: md.Trajectory) -> dict[str, float]:
    out = {}
    for i, (a, b) in enumerate(CV_PAIRS):
        ia = domain_ca_indices(traj.topology, *DOMAINS[a])
        ib = domain_ca_indices(traj.topology, *DOMAINS[b])
        ca = traj.xyz[0, ia].mean(axis=0) * 10.0
        cb = traj.xyz[0, ib].mean(axis=0) * 10.0
        out[f"cv{i}"] = float(np.linalg.norm(ca - cb))
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    cache = OUT_DIR / "library_cvs.npz"
    if cache.exists():
        data = np.load(cache, allow_pickle=True)
        names = list(data["names"])
        lib_cv0 = np.array(data["cv0"])
        lib_cv1 = np.array(data["cv1"])
        lib_cv2 = np.array(data["cv2"])
        print(f"Loaded cached library CVs ({len(names)} frames)")
    else:
        pdbs = sorted(LIB_DIR.glob("*.pdb"))
        names = [p.name for p in pdbs]
        lib_cv0, lib_cv1, lib_cv2 = [], [], []
        for k, p in enumerate(pdbs):
            t = md.load(str(p))
            cvs = cv_for_pdb(t)
            lib_cv0.append(cvs["cv0"])
            lib_cv1.append(cvs["cv1"])
            lib_cv2.append(cvs["cv2"])
            if (k + 1) % 100 == 0:
                print(f"  scored {k+1}/{len(pdbs)}")
        lib_cv0 = np.array(lib_cv0)
        lib_cv1 = np.array(lib_cv1)
        lib_cv2 = np.array(lib_cv2)
        np.savez(cache, names=np.array(names),
                 cv0=lib_cv0, cv1=lib_cv1, cv2=lib_cv2)
        print(f"Cached to {cache}")

    print(f"Library: n={len(names)}  CV0 [{lib_cv0.min():.1f},{lib_cv0.max():.1f}]")
    print(f"          CV1 [{lib_cv1.min():.1f},{lib_cv1.max():.1f}]")
    print(f"          CV2 [{lib_cv2.min():.1f},{lib_cv2.max():.1f}]")

    f = np.load(FITTED_NPZ)
    fit_cv0 = np.concatenate([f["cv0_v1"], f["cv0_v2"]])
    fit_cv1 = np.concatenate([f["cv1_v1"], f["cv1_v2"]])
    fit_cv2 = np.concatenate([f["cv2_v1"], f["cv2_v2"]])
    n_fitted = fit_cv0.size
    print(f"Fitted: n={n_fitted}")
    print(f"  CV0 [{fit_cv0.min():.1f},{fit_cv0.max():.1f}]")
    print(f"  CV1 [{fit_cv1.min():.1f},{fit_cv1.max():.1f}]")
    print(f"  CV2 [{fit_cv2.min():.1f},{fit_cv2.max():.1f}]")

    plot(lib_cv0, lib_cv1, lib_cv2,
         fit_cv0, fit_cv1, fit_cv2,
         f["cv0_v1"], f["cv0_v2"],
         FIG_DIR, OUT_DIR)
    return 0


def plot(lib0, lib1, lib2, fit0, fit1, fit2, fit0_v1, fit0_v2,
         fig_dir: Path, out_dir: Path) -> None:
    fig = plt.figure(figsize=(13.5, 9.5))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.0, 1.0, 1.2],
                          hspace=0.55, wspace=0.40)

    # Row 0 = CV0, row 1 = CV1, row 2 = CV2
    cvs = [
        ("cv0", lib0, fit0, "CV0 — αV head ↔ αV calf  (Å)", (40, 100)),
        ("cv1", lib1, fit1, "CV1 — β3 head ↔ β3 tail  (Å)", (35, 95)),
        ("cv2", lib2, fit2, "CV2 — αV head ↔ β3 head  (Å)", (15, 60)),
    ]
    for r, (key, lib, fit, label, xlim) in enumerate(cvs):
        # Library hist (left)
        ax = fig.add_subplot(gs[r, 0])
        ax.hist(lib, bins=40, color="#7570b3", alpha=0.85,
                label=f"library (n={lib.size})")
        for name, x in THRESHOLDS[key]:
            ax.axvline(x, color="#444444", linestyle="--", linewidth=1)
            ax.text(x, ax.get_ylim()[1] * 0.92,
                    f" {name}={x}", color="#444444", fontsize=7,
                    ha="left", rotation=90, va="top")
        ax.set_xlim(*xlim)
        ax.set_xlabel(label)
        ax.set_ylabel("library count")
        if r == 0:
            ax.set_title("v7 library distribution\n(steering MD samples)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")

        # Fitted hist (middle)
        ax2 = fig.add_subplot(gs[r, 1])
        ax2.hist(fit, bins=40, color="#1b9e77", alpha=0.85,
                 label=f"fitted V1+V2 (n={fit.size})")
        for name, x in THRESHOLDS[key]:
            ax2.axvline(x, color="#444444", linestyle="--", linewidth=1)
        ax2.set_xlim(*xlim)
        ax2.set_xlabel(label)
        ax2.set_ylabel("fitted count")
        if r == 0:
            ax2.set_title("V1+V2 fitted-trajectory distribution\n(unbiased HS-AFM samples)")
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=8, loc="upper right")

        # Overlay (right) — both normalized
        ax3 = fig.add_subplot(gs[r, 2])
        bins = np.linspace(*xlim, 41)
        ax3.hist(lib, bins=bins, color="#7570b3", alpha=0.45,
                 density=True, label="library")
        ax3.hist(fit, bins=bins, color="#1b9e77", alpha=0.45,
                 density=True, label="fitted V1+V2")
        # Highlight gap regions
        if key == "cv0":
            lib_max = lib.max()
            ax3.axvspan(lib_max, xlim[1], color="red", alpha=0.12, hatch="///",
                        label="library gap")
        if key == "cv2":
            lib_max = lib.max()
            eo_min = 50
            if lib_max < eo_min:
                ax3.axvspan(lib_max, xlim[1], color="red", alpha=0.12,
                            hatch="///", label="EO gap")
        for name, x in THRESHOLDS[key]:
            ax3.axvline(x, color="#444444", linestyle="--", linewidth=1)
        ax3.set_xlim(*xlim)
        ax3.set_xlabel(label)
        ax3.set_ylabel("density")
        if r == 0:
            ax3.set_title("Overlay (density) — library vs fitted")
        ax3.grid(alpha=0.3)
        ax3.legend(fontsize=8, loc="upper right")

    # Footer summary text
    n_v1 = fit0_v1.size
    n_v2 = fit0_v2.size
    eo_in_fit = int((fit2 > 50).sum())
    eo_in_lib = int((lib2 > 50).sum())
    above_lib_max = int((fit0 > lib0.max()).sum())

    fig.suptitle("Library coverage v2 — CV0/CV1/CV2 distributions across 615 v7 library frames "
                 "and V1+V2 fitted trajectories\n"
                 f"V1={n_v1}  V2={n_v2}  fitted total={fit0.size}  |  "
                 f"library CV0 max={lib0.max():.1f} Å, fitted frames > library max: {above_lib_max}  |  "
                 f"library CV2 max={lib2.max():.1f} Å (EO threshold 50 Å) — EO frames in library: {eo_in_lib}, in fitted: {eo_in_fit}",
                 fontsize=11, fontweight="bold", y=0.995)

    out_path = out_dir / "library_coverage_v2.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig.savefig(fig_dir / "library_coverage_v2.png", dpi=140, bbox_inches="tight")
    print(f"  saved {out_path}")
    print(f"  copied to {fig_dir / 'library_coverage_v2.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
