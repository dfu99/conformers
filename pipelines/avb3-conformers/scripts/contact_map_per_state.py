#!/usr/bin/env python3
"""Per-state Cα-Cα contact-map differential — audit §34 v25 / obj-065.

Computes per-frame Cα-Cα contact maps over the 790 headpiece CAs
(αV 1-440 + β3 1-350), averages per HMM state, and reports the
differential ΔC(EC - BC), ΔC(Inter - BC), ΔC(EC - Inter) plus the
top-K most-disrupted contact pairs.

Method:
  1. Load V7 fitted_coords_median_reanchored.npy for V1 (379) + V2
     (1266), Kabsch-aligned to V1 frame 0 using 790 headpiece CAs.
  2. Subset to headpiece CAs (790 residues).
  3. Per frame: Cα-Cα distance matrix; threshold at 8 Å → binary
     contact matrix.
  4. Average per HMM state (obj-048 Viterbi labels) → C_s[i,j] =
     P(contact | state s).
  5. Differential maps + top-K most-disrupted contacts (|ΔC| > 0.3).
  6. Per-residue contact-disruption score Σ_j |ΔC_{ij}| highlights
     residues whose contact pattern changes most across states.

Output:
  - figures/contact_map_per_state_v1.png
  - results/afm_pipeline/contact_map_per_state/{contact_maps.npz,
      top_contacts.json, summary.json}
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
V7 = ROOT / "results" / "afm_pipeline" / "v7_smoothed_final"
HMM_NPZ = ROOT / "results" / "afm_pipeline" / "free_energy_profile" / "hmm_states.npz"
OUT_DIR = ROOT / "results" / "afm_pipeline" / "contact_map_per_state"
FIG_DIR = ROOT / "figures"

HEAD_A = (1, 440)
HEAD_B = (1, 350)
CONTACT_CUTOFF_NM = 0.8  # 8 Å Cα cutoff
DIFF_THRESHOLD = 0.30  # contact-prob change > 0.3 = "disrupted"
TOP_K = 20

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


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    top = md.load(str(V7 / "video1" / "topology.pdb")).topology
    ca_idx, ca_chain, ca_resseq = collect_ca_meta(top)
    head_idx = np.where(head_mask(ca_chain, ca_resseq))[0]
    n_head = head_idx.size
    print(f"{ca_idx.size} CAs total, {n_head} headpiece CAs")

    print("Loading V1 + V2 ...")
    v1 = np.load(V7 / "video1" / "fitted_coords_median_reanchored.npy")
    v2 = np.load(V7 / "video2" / "fitted_coords_median_reanchored.npy")
    v1_head = v1[:, ca_idx[head_idx], :]
    v2_head = v2[:, ca_idx[head_idx], :]
    head_coords = np.concatenate([v1_head, v2_head], axis=0)
    n_frames = head_coords.shape[0]
    print(f"Combined: {n_frames} frames × {n_head} CAs × 3")

    print("Loading HMM Viterbi labels ...")
    hmm = np.load(HMM_NPZ)
    states = np.concatenate([hmm["viterbi_v1"], hmm["viterbi_v2"]]).astype(int)
    assert states.size == n_frames

    state_n = {s: int((states == s).sum()) for s in (0, 1, 2)}
    print(f"States: BC={state_n[0]}, Inter={state_n[1]}, EC={state_n[2]}")

    # Compute per-state contact maps incrementally to avoid storing
    # 1645 × 790 × 790 contact tensor (~3 GB)
    print("\nComputing per-state contact-prob maps ...")
    C_state = np.zeros((3, n_head, n_head), dtype=np.float32)

    for t in range(n_frames):
        coords_nm = head_coords[t]
        diff = coords_nm[:, None, :] - coords_nm[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        contact = (dist < CONTACT_CUTOFF_NM).astype(np.float32)
        np.fill_diagonal(contact, 0.0)  # exclude self
        C_state[states[t]] += contact
        if (t + 1) % 200 == 0:
            print(f"  {t+1}/{n_frames} frames")

    for s in (0, 1, 2):
        C_state[s] /= max(state_n[s], 1)

    # Differential maps
    delta_ec_bc = C_state[2] - C_state[0]
    delta_int_bc = C_state[1] - C_state[0]
    delta_ec_int = C_state[2] - C_state[1]

    # Top-K most-disrupted contacts (rank by |Δ EC-BC|)
    iu = np.triu_indices(n_head, k=1)
    delta_pairs = delta_ec_bc[iu]
    abs_delta = np.abs(delta_pairs)
    order = np.argsort(-abs_delta)[: TOP_K]
    top_pairs = []
    for k, p in enumerate(order):
        i, j = iu[0][p], iu[1][p]
        i_ca = head_idx[i]
        j_ca = head_idx[j]
        top_pairs.append({
            "rank": k + 1,
            "ca_idx_i": int(i_ca),
            "ca_idx_j": int(j_ca),
            "chain_i": "A" if ca_chain[i_ca] == 0 else "B",
            "chain_j": "A" if ca_chain[j_ca] == 0 else "B",
            "resSeq_i": int(ca_resseq[i_ca]),
            "resSeq_j": int(ca_resseq[j_ca]),
            "C_BC": float(C_state[0, i, j]),
            "C_Inter": float(C_state[1, i, j]),
            "C_EC": float(C_state[2, i, j]),
            "delta_EC_BC": float(delta_ec_bc[i, j]),
            "abs_delta_EC_BC": float(abs_delta[p]),
            "direction": "form" if delta_ec_bc[i, j] > 0 else "break",
        })

    # Per-residue contact-disruption scores
    disrupt_score = np.abs(delta_ec_bc).sum(axis=1)

    n_break = int((delta_ec_bc < -DIFF_THRESHOLD).sum() / 2)  # symmetric
    n_form = int((delta_ec_bc > DIFF_THRESHOLD).sum() / 2)

    print(f"\nContacts disrupted ({DIFF_THRESHOLD} threshold):")
    print(f"  Break (BC-only):   {n_break}")
    print(f"  Form (EC-only):    {n_form}")
    print(f"  Net change EC-BC:  {n_form - n_break}")
    print(f"  Sum |ΔC|:          {np.abs(delta_ec_bc).sum() / 2:.1f}")

    summary = {
        "n_frames": n_frames,
        "n_BC": state_n[0],
        "n_Inter": state_n[1],
        "n_EC": state_n[2],
        "n_head_ca": n_head,
        "contact_cutoff_A": CONTACT_CUTOFF_NM * 10,
        "diff_threshold": DIFF_THRESHOLD,
        "n_contacts_break_EC_BC": n_break,
        "n_contacts_form_EC_BC": n_form,
        "sum_abs_delta_EC_BC": float(np.abs(delta_ec_bc).sum() / 2),
        "max_abs_delta_EC_BC": float(np.abs(delta_ec_bc).max()),
        "BC_total_contacts": float(C_state[0].sum() / 2),
        "Inter_total_contacts": float(C_state[1].sum() / 2),
        "EC_total_contacts": float(C_state[2].sum() / 2),
        "top_residue_disruption_idx": int(np.argmax(disrupt_score)),
        "top_residue_disruption_score": float(disrupt_score.max()),
    }
    with (OUT_DIR / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
    with (OUT_DIR / "top_contacts.json").open("w") as f:
        json.dump({"top_K": TOP_K,
                   "diff_threshold": DIFF_THRESHOLD,
                   "ranked_by": "|ΔC EC-BC|",
                   "contacts": top_pairs}, f, indent=2)
    np.savez(OUT_DIR / "contact_maps.npz",
             C_state=C_state,
             delta_ec_bc=delta_ec_bc,
             delta_int_bc=delta_int_bc,
             delta_ec_int=delta_ec_int,
             disrupt_score=disrupt_score,
             head_ca_idx=head_idx,
             head_ca_chain=ca_chain[head_idx],
             head_ca_resseq=ca_resseq[head_idx])
    print(f"\nsaved {OUT_DIR}/{{contact_maps.npz, top_contacts.json, summary.json}}")

    print(f"\nTop-{TOP_K} most-disrupted contacts (|ΔC EC-BC|):")
    print(f"{'rank':>4s} {'pair':>20s} {'C_BC':>6s} {'C_Inter':>8s} "
          f"{'C_EC':>6s} {'ΔC':>7s} {'dir':>5s}")
    for r in top_pairs:
        pair = f"{r['chain_i']}:{r['resSeq_i']}↔{r['chain_j']}:{r['resSeq_j']}"
        print(f"  {r['rank']:>2d}  {pair:>17s}  "
              f"{r['C_BC']:>5.2f}  {r['C_Inter']:>7.2f}  "
              f"{r['C_EC']:>5.2f}  {r['delta_EC_BC']:>+6.2f}  "
              f"{r['direction']:>5s}")

    plot(C_state, delta_ec_bc, delta_int_bc, delta_ec_int,
         disrupt_score, ca_chain[head_idx], ca_resseq[head_idx],
         top_pairs, summary)
    return 0


def plot(C_state, delta_ec_bc, delta_int_bc, delta_ec_int,
         disrupt_score, head_chain, head_resseq, top_pairs, summary):
    fig = plt.figure(figsize=(15.5, 13))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.6, 1.2, 1.0],
                          hspace=0.45, wspace=0.30)

    n = C_state.shape[1]
    chain_b_start = int(np.where(head_chain == 1)[0].min())

    # Row 1: per-state contact-prob heatmaps
    for s in range(3):
        ax = fig.add_subplot(gs[0, s])
        im = ax.imshow(C_state[s], cmap="Blues", vmin=0.0, vmax=1.0,
                       origin="lower", aspect="equal")
        ax.axhline(chain_b_start, color="red", linewidth=0.6,
                   alpha=0.6)
        ax.axvline(chain_b_start, color="red", linewidth=0.6,
                   alpha=0.6)
        ax.set_xlabel("CA index (αV → β3)")
        ax.set_ylabel("CA index")
        ax.set_title(f"({chr(ord('a') + s)}) {STATE_NAMES[s]} state — "
                     f"P(contact, < 8 Å)\n"
                     f"n_frames = {summary[f'n_{STATE_NAMES[s]}']}, "
                     f"total contacts = "
                     f"{summary[f'{STATE_NAMES[s]}_total_contacts']:.0f}",
                     fontsize=9.5, fontweight="bold")
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)

    # Row 2: differential maps (ΔC EC-BC, Inter-BC, EC-Inter)
    for j, (delta, label) in enumerate([
        (delta_int_bc, "Inter − BC"),
        (delta_ec_bc, "EC − BC"),
        (delta_ec_int, "EC − Inter"),
    ]):
        ax = fig.add_subplot(gs[1, j])
        # Diverging colormap centered at 0
        vmax = max(0.5, float(np.abs(delta).max()))
        im = ax.imshow(delta, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       origin="lower", aspect="equal")
        ax.axhline(chain_b_start, color="black", linewidth=0.6,
                   alpha=0.5)
        ax.axvline(chain_b_start, color="black", linewidth=0.6,
                   alpha=0.5)
        ax.set_xlabel("CA index")
        ax.set_ylabel("CA index")
        n_break = int((delta < -DIFF_THRESHOLD).sum() / 2)
        n_form = int((delta > DIFF_THRESHOLD).sum() / 2)
        ax.set_title(f"({chr(ord('d') + j)}) ΔC ({label})\n"
                     f"break: {n_break}, form: {n_form} contacts "
                     f"(|Δ| > {DIFF_THRESHOLD})",
                     fontsize=9.5, fontweight="bold")
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)

    # Row 3: per-residue disruption score + top-K table
    ax_disrupt = fig.add_subplot(gs[2, 0])
    ax_disrupt.fill_between(np.arange(n), 0, disrupt_score,
                             color="#7570b3", alpha=0.65)
    ax_disrupt.axvline(chain_b_start, color="red", linewidth=0.6,
                        alpha=0.6, label="α/β boundary")
    ax_disrupt.set_xlabel("CA index")
    ax_disrupt.set_ylabel(r"$\Sigma_j |\Delta C_{ij}|$ (EC-BC)")
    top_disrupt_i = int(np.argmax(disrupt_score))
    chain_label = "A" if head_chain[top_disrupt_i] == 0 else "B"
    ax_disrupt.set_title(
        f"Per-residue contact-disruption score\n"
        f"(top: {chain_label}:{head_resseq[top_disrupt_i]}, "
        f"score = {disrupt_score.max():.2f})",
        fontsize=9.5, fontweight="bold")
    ax_disrupt.legend(fontsize=8)
    ax_disrupt.grid(alpha=0.3)

    # Top-K table spanning two columns
    ax_table = fig.add_subplot(gs[2, 1:])
    ax_table.axis("off")
    rows = [["#", "pair", "C_BC", "C_Inter", "C_EC", "ΔC", "dir"]]
    for r in top_pairs[:15]:
        rows.append([
            f"{r['rank']}",
            f"{r['chain_i']}:{r['resSeq_i']}↔{r['chain_j']}:{r['resSeq_j']}",
            f"{r['C_BC']:.2f}",
            f"{r['C_Inter']:.2f}",
            f"{r['C_EC']:.2f}",
            f"{r['delta_EC_BC']:+.2f}",
            r["direction"],
        ])
    table = ax_table.table(cellText=rows, loc="center", cellLoc="center",
                            colWidths=[0.04, 0.18, 0.07, 0.08, 0.07,
                                       0.07, 0.07])
    table.auto_set_font_size(False)
    table.set_fontsize(8.8)
    table.scale(1.0, 1.20)
    for k in range(len(rows[0])):
        table[(0, k)].set_text_props(weight="bold", color="white")
        table[(0, k)].set_facecolor("#444444")
    for i, r in enumerate(top_pairs[:15]):
        if r["direction"] == "break":
            for k in range(len(rows[0])):
                table[(i + 1, k)].set_facecolor("#fff0e8")
        else:
            for k in range(len(rows[0])):
                table[(i + 1, k)].set_facecolor("#e8f4ea")
    ax_table.set_title(
        f"Top-15 most-disrupted contacts by |ΔC(EC-BC)|  "
        f"(orange = break on extension, green = form)\n"
        f"Total break/form counts at |ΔC| > {DIFF_THRESHOLD}: "
        f"{summary['n_contacts_break_EC_BC']} break / "
        f"{summary['n_contacts_form_EC_BC']} form  •  "
        f"max |ΔC| = {summary['max_abs_delta_EC_BC']:.2f}",
        fontsize=9.5, fontweight="bold")

    fig.suptitle(
        f"Per-state Cα-Cα contact-map differential (obj-065) — "
        f"V1+V2 fitted, n=1645, {summary['n_head_ca']} headpiece CAs at "
        f"{summary['contact_cutoff_A']:.0f} Å cutoff\n"
        f"BC contacts {summary['BC_total_contacts']:.0f}, "
        f"Inter {summary['Inter_total_contacts']:.0f}, "
        f"EC {summary['EC_total_contacts']:.0f}  •  "
        f"top disruption residue at "
        f"{('A' if head_chain[summary['top_residue_disruption_idx']] == 0 else 'B')}:"
        f"{int(head_resseq[summary['top_residue_disruption_idx']])}, "
        f"score {summary['top_residue_disruption_score']:.2f}",
        fontsize=10.5, fontweight="bold", y=0.995,
    )

    out_png = OUT_DIR / "contact_map_per_state.png"
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    fig.savefig(FIG_DIR / "contact_map_per_state_v1.png",
                dpi=140, bbox_inches="tight")
    print(f"saved {out_png}")
    print(f"copied to {FIG_DIR / 'contact_map_per_state_v1.png'}")


if __name__ == "__main__":
    raise SystemExit(main())
