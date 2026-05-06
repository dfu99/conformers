#!/usr/bin/env python3
"""Per-state RMSF decomposition — audit §32 v23 / obj-063.

Stratifies the rotation-corrected per-residue RMSF (obj-027 / obj-042
methodology) by HMM Viterbi state assignment (obj-048: BC / Inter / EC).
Tests whether the state-resolved flexibility structure differs from
the V1+V2 pooled structure of obj-042.

Method:
  1. Load V7 fitted_coords_median_reanchored.npy for V1 (379 frames)
     and V2 (1266 frames). Both share the same topology.
  2. Kabsch-align all frames to V1 frame 0 using 790 headpiece CAs
     (αV 1-440 + β3 1-350). Same alignment as obj-042.
  3. Concatenate to combined (1645 frames) and tag each frame with
     its HMM Viterbi state (0=BC, 1=Inter, 2=EC) from
     results/afm_pipeline/free_energy_profile/hmm_states.npz.
  4. Per state s ∈ {BC, Inter, EC}: compute per-residue RMSF using
     only the frames in s (n_BC=658, n_Inter=587, n_EC=400 combined).
  5. Bootstrap (n=500) per state on the state-conditional frame
     subset → 95 % CI per residue per state.
  6. State-difference test: per-residue Δ(EC-BC) RMSF and Δ(Inter-BC)
     RMSF + bootstrap CI on Δ. Residues whose Δ-CI excludes 0 are
     state-differential.

Output:
  - figures/rmsf_per_state_v1.png
  - results/afm_pipeline/rmsf_per_state/{rmsf_per_state.npz,
      state_differential.json, summary.json}
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
V7 = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final"
HMM_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "hmm_states.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "rmsf_per_state"
FIG_DIR = ROOT / "figures"

HEAD_A = (1, 440)
HEAD_B = (1, 350)

STATE_NAMES = ["BC", "Inter", "EC"]
STATE_COLORS = ["#4a90e2", "#f5a623", "#d0021b"]


def chain_id(atom) -> str:
    cid = atom.residue.chain.chain_id
    return cid if cid else chr(ord("A") + atom.residue.chain.index)


def collect_ca_meta(top):
    ca_idx, ca_chain, ca_resseq = [], [], []
    for atom in top.atoms:
        if atom.name != "CA":
            continue
        ca_idx.append(atom.index)
        ca_chain.append(0 if chain_id(atom) == "A" else 1)
        ca_resseq.append(atom.residue.resSeq)
    return np.array(ca_idx), np.array(ca_chain), np.array(ca_resseq)


def head_mask(ca_chain, ca_resseq):
    return (
        ((ca_chain == 0) & (ca_resseq >= HEAD_A[0]) & (ca_resseq <= HEAD_A[1]))
        | ((ca_chain == 1) & (ca_resseq >= HEAD_B[0]) & (ca_resseq <= HEAD_B[1]))
    )


def kabsch_R(P, Q):
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return U @ D @ Vt


def head_align_to_ref(coords_ca_nm, head_mask_ca, ref_head_centered, ref_centroid):
    """Kabsch-align all frames to a fixed reference using only headpiece CAs."""
    aligned = np.empty_like(coords_ca_nm)
    for t in range(coords_ca_nm.shape[0]):
        head = coords_ca_nm[t, head_mask_ca]
        head_centroid = head.mean(axis=0)
        moving_centered = head - head_centroid
        R = kabsch_R(moving_centered, ref_head_centered)
        translated = coords_ca_nm[t] - head_centroid
        aligned[t] = translated @ R + ref_centroid
    return aligned


def rmsf_per_residue(coords_ca_aligned_nm, idx=None):
    """RMSF in Å per residue."""
    if idx is not None:
        coords = coords_ca_aligned_nm[idx]
    else:
        coords = coords_ca_aligned_nm
    mean = coords.mean(axis=0)
    delta = coords - mean
    return np.sqrt((delta ** 2).sum(axis=-1).mean(axis=0)) * 10.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-boot", type=int, default=500)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--top-k", type=int, default=10)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    top1 = md.load(str(V7 / "video1" / "topology.pdb")).topology
    ca_idx, ca_chain, ca_resseq = collect_ca_meta(top1)
    n_ca = len(ca_idx)
    head_mask_ca = head_mask(ca_chain, ca_resseq)
    print(f"{n_ca} CAs total, {int(head_mask_ca.sum())} headpiece CAs")

    print("Loading V1 + V2 ...")
    v1 = np.load(V7 / "video1" / "fitted_coords_median_reanchored.npy")
    v2 = np.load(V7 / "video2" / "fitted_coords_median_reanchored.npy")
    v1_ca = v1[:, ca_idx, :]
    v2_ca = v2[:, ca_idx, :]

    # Common reference: V1 frame 0 headpiece centroid + centered headpiece.
    # All V1 + V2 frames align to the same reference so per-state subgroups
    # share a coordinate frame.
    ref_head = v1_ca[0, head_mask_ca]
    ref_centroid = ref_head.mean(axis=0)
    ref_centered = ref_head - ref_centroid

    print("Aligning V1 ...")
    v1_aligned = head_align_to_ref(v1_ca, head_mask_ca, ref_centered, ref_centroid)
    print("Aligning V2 ...")
    v2_aligned = head_align_to_ref(v2_ca, head_mask_ca, ref_centered, ref_centroid)
    combined = np.concatenate([v1_aligned, v2_aligned], axis=0)
    print(f"Combined: {combined.shape[0]} frames")

    print("Loading HMM Viterbi states ...")
    hmm = np.load(HMM_NPZ)
    states = np.concatenate([hmm["viterbi_v1"], hmm["viterbi_v2"]]).astype(int)
    assert states.shape[0] == combined.shape[0], "frame/state count mismatch"

    state_idx = {s: np.where(states == s)[0] for s in (0, 1, 2)}
    for s in (0, 1, 2):
        print(f"  state {STATE_NAMES[s]} (id={s}): n_frames = {state_idx[s].size}")

    # Pooled (control) RMSF for sanity check
    rmsf_pooled = rmsf_per_residue(combined)
    print(f"Pooled RMSF: mean={rmsf_pooled.mean():.2f} Å, max={rmsf_pooled.max():.2f}")

    # Per-state point estimates
    rmsf_state = np.zeros((3, n_ca))
    for s in (0, 1, 2):
        rmsf_state[s] = rmsf_per_residue(combined[state_idx[s]])
        print(f"  RMSF[{STATE_NAMES[s]}] mean={rmsf_state[s].mean():.2f}, "
              f"max={rmsf_state[s].max():.2f}")

    # Bootstrap per state (resample frames within each state)
    print(f"\nBootstrap {args.n_boot} per state ...")
    rng = np.random.default_rng(args.seed)
    boot = np.zeros((3, args.n_boot, n_ca))
    for s in (0, 1, 2):
        n_s = state_idx[s].size
        for k in range(args.n_boot):
            sub_idx = rng.integers(0, n_s, size=n_s)
            boot[s, k] = rmsf_per_residue(combined[state_idx[s][sub_idx]])
        print(f"  state {STATE_NAMES[s]} bootstrap done")

    rmsf_lo = np.quantile(boot, 0.025, axis=1)
    rmsf_hi = np.quantile(boot, 0.975, axis=1)
    rmsf_med = np.quantile(boot, 0.50, axis=1)

    # State-difference Δ_RMSF per residue + bootstrap CI on Δ
    # Use paired difference within bootstrap iterations because per-iter
    # frame draws are independent across states; the across-iteration
    # quantiles bound the per-residue Δ sampling distribution.
    print("\nState-difference contrasts (EC-BC) and (Inter-BC) ...")
    delta_ec_bc_boot = boot[2] - boot[0]            # [n_boot, n_ca]
    delta_int_bc_boot = boot[1] - boot[0]
    delta_ec_int_boot = boot[2] - boot[1]
    delta_ec_bc = delta_ec_bc_boot.mean(axis=0)
    delta_int_bc = delta_int_bc_boot.mean(axis=0)
    delta_ec_int = delta_ec_int_boot.mean(axis=0)

    delta_ec_bc_lo = np.quantile(delta_ec_bc_boot, 0.025, axis=0)
    delta_ec_bc_hi = np.quantile(delta_ec_bc_boot, 0.975, axis=0)
    delta_int_bc_lo = np.quantile(delta_int_bc_boot, 0.025, axis=0)
    delta_int_bc_hi = np.quantile(delta_int_bc_boot, 0.975, axis=0)
    delta_ec_int_lo = np.quantile(delta_ec_int_boot, 0.025, axis=0)
    delta_ec_int_hi = np.quantile(delta_ec_int_boot, 0.975, axis=0)

    sig_ec_bc = (delta_ec_bc_lo > 0) | (delta_ec_bc_hi < 0)
    sig_int_bc = (delta_int_bc_lo > 0) | (delta_int_bc_hi < 0)
    sig_ec_int = (delta_ec_int_lo > 0) | (delta_ec_int_hi < 0)

    n_sig_ec_bc = int(sig_ec_bc.sum())
    n_sig_int_bc = int(sig_int_bc.sum())
    n_sig_ec_int = int(sig_ec_int.sum())
    print(f"  Δ(EC-BC) sig: {n_sig_ec_bc} / {n_ca} residues")
    print(f"  Δ(Inter-BC) sig: {n_sig_int_bc} / {n_ca} residues")
    print(f"  Δ(EC-Inter) sig: {n_sig_ec_int} / {n_ca} residues")

    # Top-K differential residues by |Δ EC-BC|
    top_diff_idx = np.argsort(-np.abs(delta_ec_bc))[: args.top_k]
    top_diff = []
    for k, i in enumerate(top_diff_idx):
        top_diff.append({
            "rank": k + 1,
            "ca_idx": int(i),
            "chain": "A" if ca_chain[i] == 0 else "B",
            "resSeq": int(ca_resseq[i]),
            "rmsf_BC": float(rmsf_state[0, i]),
            "rmsf_Inter": float(rmsf_state[1, i]),
            "rmsf_EC": float(rmsf_state[2, i]),
            "delta_EC_BC": float(delta_ec_bc[i]),
            "delta_EC_BC_ci": [float(delta_ec_bc_lo[i]),
                                float(delta_ec_bc_hi[i])],
            "sig_EC_BC": bool(sig_ec_bc[i]),
        })

    np.savez(
        OUT_DIR / "rmsf_per_state.npz",
        rmsf_state=rmsf_state,
        rmsf_lo=rmsf_lo,
        rmsf_hi=rmsf_hi,
        rmsf_med=rmsf_med,
        rmsf_pooled=rmsf_pooled,
        delta_ec_bc=delta_ec_bc,
        delta_int_bc=delta_int_bc,
        delta_ec_int=delta_ec_int,
        delta_ec_bc_lo=delta_ec_bc_lo,
        delta_ec_bc_hi=delta_ec_bc_hi,
        ca_chain=ca_chain,
        ca_resseq=ca_resseq,
        n_BC=state_idx[0].size,
        n_Inter=state_idx[1].size,
        n_EC=state_idx[2].size,
        n_boot=args.n_boot,
    )
    print(f"saved {OUT_DIR / 'rmsf_per_state.npz'}")

    summary = {
        "n_frames_combined": int(combined.shape[0]),
        "n_BC": int(state_idx[0].size),
        "n_Inter": int(state_idx[1].size),
        "n_EC": int(state_idx[2].size),
        "n_ca_total": int(n_ca),
        "n_boot": int(args.n_boot),
        "rmsf_pooled_mean_A": float(rmsf_pooled.mean()),
        "rmsf_pooled_max_A": float(rmsf_pooled.max()),
        "rmsf_BC_mean_A": float(rmsf_state[0].mean()),
        "rmsf_Inter_mean_A": float(rmsf_state[1].mean()),
        "rmsf_EC_mean_A": float(rmsf_state[2].mean()),
        "rmsf_BC_max_A": float(rmsf_state[0].max()),
        "rmsf_Inter_max_A": float(rmsf_state[1].max()),
        "rmsf_EC_max_A": float(rmsf_state[2].max()),
        "n_sig_ec_bc": n_sig_ec_bc,
        "n_sig_int_bc": n_sig_int_bc,
        "n_sig_ec_int": n_sig_ec_int,
        "frac_sig_ec_bc": n_sig_ec_bc / n_ca,
        "max_abs_delta_ec_bc_A": float(np.abs(delta_ec_bc).max()),
        "median_delta_ec_bc_A": float(np.median(delta_ec_bc)),
        "ranking_pearson_BC_EC": float(
            np.corrcoef(rmsf_state[0], rmsf_state[2])[0, 1]),
        "ranking_pearson_BC_Inter": float(
            np.corrcoef(rmsf_state[0], rmsf_state[1])[0, 1]),
        "ranking_pearson_Inter_EC": float(
            np.corrcoef(rmsf_state[1], rmsf_state[2])[0, 1]),
    }
    with (OUT_DIR / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
    with (OUT_DIR / "state_differential.json").open("w") as f:
        json.dump({"top_diff_EC_BC": top_diff,
                   "n_sig_ec_bc": n_sig_ec_bc,
                   "n_sig_int_bc": n_sig_int_bc,
                   "n_sig_ec_int": n_sig_ec_int,
                   "n_ca_total": int(n_ca),
                   "method": "Bootstrap n=500 per state on V1+V2 "
                             "Kabsch-aligned headpiece coords; Δ-CI "
                             "excludes 0 → state-differential"},
                  f, indent=2)
    print(f"saved summary + state_differential JSON")

    print(f"\nTop-{args.top_k} differential residues (|Δ EC-BC|):")
    print(f"{'rank':>4s} {'chain':>5s} {'resSeq':>6s} "
          f"{'BC':>7s} {'Inter':>7s} {'EC':>7s} "
          f"{'Δ(EC-BC)':>10s} {'95% CI':>16s} {'sig':>4s}")
    for r in top_diff:
        print(f"  {r['rank']:>2d}  {r['chain']:>4s} {r['resSeq']:>5d}  "
              f"{r['rmsf_BC']:>6.1f}  {r['rmsf_Inter']:>6.1f}  {r['rmsf_EC']:>6.1f}  "
              f"{r['delta_EC_BC']:>+9.2f}  "
              f"[{r['delta_EC_BC_ci'][0]:+.2f}, {r['delta_EC_BC_ci'][1]:+.2f}]  "
              f"{'✓' if r['sig_EC_BC'] else '·':>4s}")

    plot(rmsf_state, rmsf_lo, rmsf_hi, rmsf_pooled,
         delta_ec_bc, delta_ec_bc_lo, delta_ec_bc_hi, sig_ec_bc,
         ca_chain, top_diff,
         summary, n_boot=args.n_boot)
    return 0


def plot(rmsf_state, rmsf_lo, rmsf_hi, rmsf_pooled,
         delta_ec_bc, delta_ec_bc_lo, delta_ec_bc_hi, sig_ec_bc,
         ca_chain, top_diff,
         summary, n_boot):
    fig = plt.figure(figsize=(15.5, 12))
    gs = fig.add_gridspec(3, 2,
                          height_ratios=[1.6, 1.0, 1.2],
                          hspace=0.55, wspace=0.18)

    n_ca = rmsf_state.shape[1]
    x = np.arange(n_ca)
    chain_a = ca_chain == 0

    # Top: per-residue RMSF curves per state, both chains overlaid with
    # subtle separator
    ax = fig.add_subplot(gs[0, :])
    chain_a_max_idx = int(np.where(chain_a)[0].max()) if chain_a.any() else 0
    for s in range(3):
        ax.fill_between(x, rmsf_lo[s], rmsf_hi[s],
                        color=STATE_COLORS[s], alpha=0.18)
        ax.plot(x, rmsf_state[s], color=STATE_COLORS[s], linewidth=0.95,
                label=(f"{STATE_NAMES[s]}  "
                       f"(n_frames={summary[f'n_{STATE_NAMES[s]}']:>4d}, "
                       f"mean RMSF "
                       f"{summary[f'rmsf_{STATE_NAMES[s]}_mean_A']:.1f} Å)"))
    # Faint pooled overlay for comparison
    ax.plot(x, rmsf_pooled, color="black", linewidth=0.6, linestyle=":",
            alpha=0.55, label=f"Pooled V1+V2 (n=1645, mean "
                              f"{summary['rmsf_pooled_mean_A']:.1f} Å)")
    ax.axvline(chain_a_max_idx + 0.5, color="#888", linewidth=0.7,
               linestyle="--", alpha=0.6)
    ax.text(chain_a_max_idx * 0.5, ax.get_ylim()[1] * 0.94, "αV (chain A)",
            fontsize=9, ha="center", color="#1f4970")
    ax.text((chain_a_max_idx + n_ca) * 0.5, ax.get_ylim()[1] * 0.94,
            "β3 (chain B)",
            fontsize=9, ha="center", color="#7a1410")
    ax.set_xlabel("CA index (chain A → chain B)")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title(
        f"Per-state RMSF decomposition (obj-063) — "
        f"{n_boot} bootstrap per state, V1+V2 fitted (n=1645)\n"
        f"State-resolved per-residue flexibility under HMM Viterbi "
        f"labels (obj-048): does the BC↔EC ranking depend on state?",
        fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=9.0, framealpha=0.92)
    ax.grid(alpha=0.3)

    # Mid-left: Δ(EC-BC) RMSF with significance markers
    ax = fig.add_subplot(gs[1, 0])
    ax.fill_between(x, delta_ec_bc_lo, delta_ec_bc_hi,
                    color="#666666", alpha=0.30)
    ax.plot(x, delta_ec_bc, color="#222222", linewidth=0.85,
            label="Δ(EC − BC) point")
    ax.scatter(x[sig_ec_bc], delta_ec_bc[sig_ec_bc],
               s=8, color="#d62728", zorder=10, alpha=0.8,
               label=f"Δ-CI excludes 0 ({int(sig_ec_bc.sum())} residues)")
    ax.axhline(0, color="#888", linewidth=0.6)
    ax.set_xlabel("CA index")
    ax.set_ylabel("Δ RMSF (Å)")
    ax.set_title(
        f"State-differential flexibility: Δ(EC − BC) per residue\n"
        f"{int(sig_ec_bc.sum())} / {n_ca} residues with 95 %-CI "
        f"excluding 0 ({100 * sig_ec_bc.mean():.1f} %)",
        fontsize=10, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)

    # Mid-right: BC RMSF vs EC RMSF scatter (state-correlation diagnosis)
    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(rmsf_state[0], rmsf_state[2], s=4, alpha=0.45,
               color="#7570b3")
    lim = max(rmsf_state[0].max(), rmsf_state[2].max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="#444444", linewidth=1,
            label="y = x")
    r = summary["ranking_pearson_BC_EC"]
    ax.text(0.05, 0.95,
            f"Pearson r(BC, EC) = {r:.3f}\n"
            f"r(BC, Inter) = {summary['ranking_pearson_BC_Inter']:.3f}\n"
            f"r(Inter, EC) = {summary['ranking_pearson_Inter_EC']:.3f}",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round", facecolor="white",
                      edgecolor="#cccccc"))
    ax.set_xlabel("BC RMSF (Å)")
    ax.set_ylabel("EC RMSF (Å)")
    ax.set_title("BC vs EC per-residue RMSF — does the ranking change?",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)

    # Bottom: top-10 differential residues table
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    rows = [["rank", "chain:resSeq", "BC RMSF (Å)", "Inter RMSF (Å)",
             "EC RMSF (Å)", "Δ(EC−BC)", "95 % CI", "sig"]]
    for r in top_diff:
        rows.append([
            f"{r['rank']}",
            f"{r['chain']}:{r['resSeq']}",
            f"{r['rmsf_BC']:.1f}",
            f"{r['rmsf_Inter']:.1f}",
            f"{r['rmsf_EC']:.1f}",
            f"{r['delta_EC_BC']:+.2f}",
            f"[{r['delta_EC_BC_ci'][0]:+.2f}, {r['delta_EC_BC_ci'][1]:+.2f}]",
            "✓" if r["sig_EC_BC"] else "·",
        ])
    table = ax.table(cellText=rows, loc="center", cellLoc="center",
                     colWidths=[0.05, 0.13, 0.10, 0.10, 0.10, 0.10,
                                 0.18, 0.05])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.30)
    for k in range(len(rows[0])):
        table[(0, k)].set_text_props(weight="bold", color="white")
        table[(0, k)].set_facecolor("#444444")
    for i, r in enumerate(top_diff):
        if r["sig_EC_BC"]:
            for k in range(len(rows[0])):
                table[(i + 1, k)].set_facecolor("#e8f4ea")
    ax.set_title(
        f"Top-10 most state-differential residues by |Δ(EC−BC)| RMSF — "
        f"green rows are 95 %-CI-significant\n"
        f"({summary['n_sig_ec_bc']} residues sig. EC-BC, "
        f"{summary['n_sig_int_bc']} sig. Inter-BC, "
        f"{summary['n_sig_ec_int']} sig. EC-Inter; max |Δ| = "
        f"{summary['max_abs_delta_ec_bc_A']:.1f} Å)",
        fontsize=10, fontweight="bold")

    out_png = OUT_DIR / "rmsf_per_state.png"
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    fig.savefig(FIG_DIR / "rmsf_per_state_v1.png",
                dpi=140, bbox_inches="tight")
    print(f"saved {out_png}")
    print(f"copied to {FIG_DIR / 'rmsf_per_state_v1.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
