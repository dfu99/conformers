#!/usr/bin/env python3
"""F4 Reviewer-C contact-mechanics control (audit §10.8).

Compares the hard-sphere geometric pseudo-AFM dilation (currently used
in `simulate_afm_video.py`) against a Hertzian contact-mechanics model
on a 2 nm probe sphere imaged with a 1.5 nm Si tip at HS-AFM-typical
imaging forces (50, 100, 200 pN).

Physics
-------
For two spheres of radii R_p (protein, soft) and R_t (tip, hard):
  R_eff = R_p R_t / (R_p + R_t)
  1/E*  = (1 - ν_p²)/E_p + (1 - ν_t²)/E_t
  δ     = (9 F² / (16 E*² R_eff))^(1/3)             (indentation)
  a     = (3 F R_eff / (4 E*))^(1/3)                (contact radius)

Material constants (from Müller/Ando HS-AFM literature on proteins):
  E_protein ≈ 1 GPa  (range 0.3-3 GPa across proteins)
  E_tip     ≈ 130 GPa (Si)
  ν         ≈ 0.30   for both

For HS-AFM imaging on the αVβ3 head domain at 100 pN, the predicted
indentation δ ≈ 0.18 nm — quantifiable and below the typical 1 nm
vertical noise floor.

Output: figures/contact_mechanics_control.png +
results/afm_pipeline/contact_mechanics/{summary.json, summary.npz}.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "results" / "afm_pipeline" / "contact_mechanics"
FIG_DIR = ROOT / "figures"

# Material + geometry
R_PROTEIN_NM = 2.0
R_TIP_NM = 1.5
E_PROTEIN_PA = 1.0e9        # 1 GPa  (Müller/Ando; 0.3-3 typical)
E_TIP_PA = 130.0e9          # Silicon tip
NU = 0.30                   # Poisson, both materials

# Imaging
PIXEL_PITCH_NM = 0.98       # standard sim-AFM
CANVAS_NM = 16.0            # focused view
FORCES_PN = np.array([50.0, 100.0, 200.0])   # HS-AFM imaging force range


def hertzian_indentation_nm(force_n: np.ndarray, r_eff_m: float,
                            e_star_pa: float) -> np.ndarray:
    """Hertz indentation δ in metres, returned in nm."""
    delta_m = (9.0 * force_n**2 / (16.0 * e_star_pa**2 * r_eff_m)) ** (1.0 / 3.0)
    return delta_m * 1e9


def hertzian_contact_radius_nm(force_n: np.ndarray, r_eff_m: float,
                               e_star_pa: float) -> np.ndarray:
    """Hertz contact radius a in m, returned in nm."""
    a_m = (3.0 * force_n * r_eff_m / (4.0 * e_star_pa)) ** (1.0 / 3.0)
    return a_m * 1e9


def render_hard_sphere(r_part_nm: float, r_tip_nm: float,
                       pixel_nm: float, canvas_nm: float) -> tuple[np.ndarray, np.ndarray]:
    """Hard-sphere paraboloid dilation (matches calibration_controls.py).
    h(x,y) = sqrt((R+r)² - x² - y²) - r   for r² < (R+r)²."""
    n = int(round(canvas_nm / pixel_nm))
    xs = (np.arange(n) - n / 2.0) * pixel_nm
    ys = (np.arange(n) - n / 2.0) * pixel_nm
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    R = np.sqrt(X**2 + Y**2)
    rc = r_part_nm + r_tip_nm
    h = np.zeros_like(R)
    inside = R < rc
    h[inside] = np.sqrt(rc**2 - R[inside]**2) - r_tip_nm
    return np.clip(h, 0, None), xs


def render_hertzian(r_part_nm: float, r_tip_nm: float,
                    pixel_nm: float, canvas_nm: float,
                    force_pn: float, e_star_pa: float,
                    r_eff_m: float) -> np.ndarray:
    """Apparent height under Hertzian compression at constant imaging
    force. The whole height map is depressed by δ(F) at the contact
    point; off-contact (where local curvature does not put the tip in
    contact at all) reverts to hard-sphere geometric.

    Implementation: for each pixel where the hard-sphere height is
    above the half-indentation threshold, subtract δ(F). Below that
    threshold the tip is "off" the molecule and δ → 0 smoothly.
    """
    h_hard, _ = render_hard_sphere(r_part_nm, r_tip_nm, pixel_nm, canvas_nm)
    delta = hertzian_indentation_nm(np.array([force_pn * 1e-12]),
                                    r_eff_m, e_star_pa)[0]
    # Smooth transition: full indentation in contact zone, taper at periphery
    a_nm = hertzian_contact_radius_nm(np.array([force_pn * 1e-12]),
                                      r_eff_m, e_star_pa)[0]
    n = h_hard.shape[0]
    xs = (np.arange(n) - n / 2.0) * pixel_nm
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    radial = np.sqrt(X**2 + Y**2)
    # Contact zone: indentation = δ at center, decays smoothly with
    # tip footprint (diameter ≈ 2a). Use a cosine-tapered window.
    contact_mask = radial < (r_part_nm + r_tip_nm)
    taper = np.where(radial < a_nm, 1.0,
                     np.where(radial < (r_part_nm + r_tip_nm),
                              0.5 * (1 + np.cos(np.pi
                                                * (radial - a_nm)
                                                / max(r_part_nm + r_tip_nm - a_nm, 1e-6))),
                              0.0))
    indentation = delta * taper * contact_mask
    h = np.clip(h_hard - indentation, 0, None)
    return h


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    r_eff_m = (R_PROTEIN_NM * R_TIP_NM
               / (R_PROTEIN_NM + R_TIP_NM)) * 1e-9
    e_star_pa = 1.0 / ((1 - NU**2) / E_PROTEIN_PA
                       + (1 - NU**2) / E_TIP_PA)
    print(f"R_eff = {r_eff_m * 1e9:.3f} nm   E* = {e_star_pa / 1e9:.2f} GPa")

    # Indentation table
    forces_n = FORCES_PN * 1e-12
    deltas = hertzian_indentation_nm(forces_n, r_eff_m, e_star_pa)
    contact_radii = hertzian_contact_radius_nm(forces_n, r_eff_m, e_star_pa)
    pmax_pa = (6.0 * forces_n * e_star_pa**2 / (np.pi**3 * r_eff_m**2)) ** (1.0 / 3.0)
    print("Force(pN)  δ(nm)   a(nm)   p_max(MPa)")
    for f, d, a, p in zip(FORCES_PN, deltas, contact_radii, pmax_pa):
        print(f"  {f:6.1f}    {d:.3f}   {a:.3f}   {p / 1e6:.1f}")

    # Render fields
    h_hard, xs = render_hard_sphere(R_PROTEIN_NM, R_TIP_NM,
                                    PIXEL_PITCH_NM, CANVAS_NM)
    h_hertz = {f: render_hertzian(R_PROTEIN_NM, R_TIP_NM,
                                  PIXEL_PITCH_NM, CANVAS_NM,
                                  f, e_star_pa, r_eff_m)
               for f in FORCES_PN}

    # Hard-sphere FWHM + apparent height for sanity
    hard_max = float(h_hard.max())
    hertz_max = {f: float(h.max()) for f, h in h_hertz.items()}

    plot(h_hard, h_hertz, xs, deltas, contact_radii, pmax_pa,
         hard_max, hertz_max, e_star_pa, r_eff_m)

    summary = {
        "params": {
            "R_protein_nm": R_PROTEIN_NM,
            "R_tip_nm": R_TIP_NM,
            "E_protein_GPa": E_PROTEIN_PA / 1e9,
            "E_tip_GPa": E_TIP_PA / 1e9,
            "nu": NU,
            "R_eff_nm": r_eff_m * 1e9,
            "E_star_GPa": e_star_pa / 1e9,
            "pixel_pitch_nm": PIXEL_PITCH_NM,
            "canvas_nm": CANVAS_NM,
        },
        "force_pN": FORCES_PN.tolist(),
        "indentation_nm": deltas.tolist(),
        "contact_radius_nm": contact_radii.tolist(),
        "p_max_MPa": (pmax_pa / 1e6).tolist(),
        "hard_sphere_apparent_height_nm": hard_max,
        "hertzian_apparent_height_nm": hertz_max,
        "systematic_offset_nm": {f: hard_max - hertz_max[f]
                                 for f in FORCES_PN},
    }
    import json
    with (OUT_DIR / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2,
                  default=lambda o: float(o) if hasattr(o, "item") else str(o))
    np.savez(OUT_DIR / "summary.npz",
             h_hard=h_hard,
             h_hertz_50=h_hertz[50.0],
             h_hertz_100=h_hertz[100.0],
             h_hertz_200=h_hertz[200.0],
             xs=xs, forces_pn=FORCES_PN,
             deltas_nm=deltas, contact_radii_nm=contact_radii,
             pmax_mpa=pmax_pa / 1e6,
             E_star_pa=e_star_pa, R_eff_m=r_eff_m)
    print(f"saved {OUT_DIR / 'summary.json'}")
    return 0


def plot(h_hard, h_hertz, xs, deltas, contact_radii, pmax_pa,
         hard_max, hertz_max, e_star_pa, r_eff_m):
    fig = plt.figure(figsize=(13.5, 9.0))
    gs = fig.add_gridspec(3, 4,
                          height_ratios=[1.0, 0.8, 0.6],
                          hspace=0.45, wspace=0.34)

    # --- Top row: 4 height maps (hard, F=50, 100, 200)
    panels = [("Hard-sphere\n(no contact mechanics)", h_hard, None),
              ("Hertzian, F = 50 pN",  h_hertz[50.0],  50.0),
              ("Hertzian, F = 100 pN", h_hertz[100.0], 100.0),
              ("Hertzian, F = 200 pN", h_hertz[200.0], 200.0)]
    extent = (xs[0], xs[-1], xs[0], xs[-1])
    vmax = h_hard.max() * 1.02
    for i, (title, h, force) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i])
        im = ax.imshow(h, cmap="copper", origin="lower",
                       extent=extent, vmin=0, vmax=vmax, aspect="equal")
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel("x (nm)")
        if i == 0:
            ax.set_ylabel("y (nm)")
        if force is not None:
            d_idx = list(panels).index((title, h, force)) - 1
            d = deltas[d_idx]
            ax.text(0.04, 0.94, f"δ={d:.2f} nm",
                    transform=ax.transAxes, color="white", fontsize=9,
                    fontweight="bold", va="top",
                    bbox=dict(boxstyle="round,pad=0.3",
                              facecolor="black", alpha=0.55))
    cbar = fig.colorbar(im, ax=fig.axes[:4], shrink=0.62, pad=0.02,
                        location="right")
    cbar.set_label("Apparent height (nm)")

    # --- Mid row L: cross-section line plots (y=0 across panels)
    ax = fig.add_subplot(gs[1, :2])
    n = h_hard.shape[0]
    mid = n // 2
    ax.plot(xs, h_hard[mid], color="black", linewidth=2.2,
            label="hard-sphere")
    colors = ["#1b9e77", "#d95f02", "#7570b3"]
    for color, (force, h) in zip(colors, h_hertz.items()):
        ax.plot(xs, h[mid], linewidth=1.8, linestyle="--",
                color=color, label=f"Hertzian F={int(force)} pN")
    ax.axhline(R_PROTEIN_NM, color="#888888", linestyle=":", linewidth=1)
    ax.text(xs[-1] * 0.98, R_PROTEIN_NM,
            f" R_part={R_PROTEIN_NM} nm", ha="right", fontsize=8, va="bottom")
    ax.set_xlabel("x (nm)  (cross-section through y=0)")
    ax.set_ylabel("Apparent height (nm)")
    ax.set_title("Cross-section: indentation reduces apparent height",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)

    # --- Mid row R: indentation vs force curve
    ax = fig.add_subplot(gs[1, 2:])
    f_grid = np.linspace(10, 400, 100)
    delta_grid = hertzian_indentation_nm(f_grid * 1e-12, r_eff_m, e_star_pa)
    ax.plot(f_grid, delta_grid, color="black", linewidth=2,
            label=r"Hertz $\delta = (9F^2/(16 E^{*2} R_{eff}))^{1/3}$")
    ax.scatter([50, 100, 200], deltas, s=80, color="#d95f02",
               zorder=5, edgecolor="black", linewidth=0.6,
               label="HS-AFM imaging force range")
    ax.axhline(1.0, color="red", linestyle="--", linewidth=1,
               label="HS-AFM noise floor ≈ 1 nm")
    ax.set_xlabel("Imaging force F (pN)")
    ax.set_ylabel(r"Indentation $\delta$ (nm)")
    ax.set_title("Hertz indentation vs imaging force",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)

    # --- Bottom row: systematic-offset summary table + take-aways
    ax = fig.add_subplot(gs[2, :2])
    ax.axis("off")
    rows = [
        ["F (pN)", "δ (nm)", "contact radius a (nm)", "p_max (MPa)",
         "hard − Hertz Δh (nm)"],
    ]
    for f, d, a, p in zip([50.0, 100.0, 200.0], deltas, contact_radii, pmax_pa):
        rows.append([f"{f:.0f}", f"{d:.3f}", f"{a:.3f}",
                     f"{p / 1e6:.0f}", f"{hard_max - hertz_max[f]:.3f}"])
    table = ax.table(cellText=rows, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 1.55)
    for k in range(len(rows[0])):
        table[(0, k)].set_text_props(weight="bold", color="white")
        table[(0, k)].set_facecolor("#444444")
    ax.set_title("Hertzian summary at HS-AFM imaging forces",
                 fontsize=10, fontweight="bold")

    ax = fig.add_subplot(gs[2, 2:])
    ax.axis("off")
    bullets = (
        "Reviewer-C contact-mechanics control:\n\n"
        f"  • R_eff = {r_eff_m * 1e9:.2f} nm,  E* = {e_star_pa / 1e9:.2f} GPa\n"
        "    (E_protein=1 GPa, E_tip=130 GPa, ν=0.3)\n\n"
        "  • At HS-AFM imaging force 50-200 pN, indentation\n"
        "    δ = 0.11 - 0.28 nm  →  apparent height drops\n"
        "    by 5-14 % for a 2 nm probe sphere.\n\n"
        "  • All values sit BELOW the HS-AFM noise floor\n"
        "    (1 nm vertical, per Ando 2013 review).\n\n"
        "  • Conclusion: hard-sphere pseudo-AFM is a defensible\n"
        "    first-order model for fitting; Hertzian correction\n"
        "    is a known systematic in the noise-bounded regime.\n\n"
        "  • Recommended pipeline change: optional --hertz-force\n"
        "    argument in `simulate_afm_video.py` for downstream\n"
        "    correction; document offset in figure captions."
    )
    ax.text(0.0, 1.0, bullets, transform=ax.transAxes, fontsize=9.5,
            ha="left", va="top", family="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#fafafa",
                      edgecolor="#cccccc"))

    fig.suptitle("F4 — Pseudo-AFM contact-mechanics control "
                 "(Reviewer C, audit-2026-05-05 §10.8)\n"
                 "Hertzian indentation on a 2 nm probe sphere "
                 "at HS-AFM imaging forces",
                 fontsize=12, fontweight="bold", y=0.995)
    out_path = OUT_DIR / "contact_mechanics_control.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig.savefig(FIG_DIR / "contact_mechanics_control.png",
                dpi=140, bbox_inches="tight")
    print(f"saved {out_path}")
    print(f"copied to {FIG_DIR / 'contact_mechanics_control.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
