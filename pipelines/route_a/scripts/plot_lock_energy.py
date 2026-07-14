#!/usr/bin/env python3
"""Figure for lock_energy_analysis.py: per-bridge lock strengths + energy each mutation removes.

Reads results/route_a/lock_energy.json (produced under the snapback env) and renders under base
python3 (matplotlib). Two panels: (a) the four genu salt-bridge lock strengths in the WT extended
state; (b) the lock energy each mutation removes = the static snap-back prediction.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BRIDGE_ORDER = ["E598-K459", "K459-E636", "D457-K688", "D595-K688"]
COLORS = {"E598-K459": "#b2182b", "K459-E636": "#d6604d",
          "D457-K688": "#2166ac", "D595-K688": "#1b7837"}
MUT_COLORS = {"E598A": "#d6604d", "K459A": "#b2182b", "K459A+E598A": "#762a83"}

d = json.load(open("results/route_a/lock_energy.json"))
pb = d["per_bridge"]; pm = d["per_mutation"]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 5.4))

# ---- (a) per-bridge lock strength (magnitude of direct interaction energy) ----
names = sorted(BRIDGE_ORDER, key=lambda b: abs(pb[b]["E_total_direct_kcal"]))
y = range(len(names))
mags = [abs(pb[b]["E_total_direct_kcal"]) for b in names]
axL.barh(list(y), mags, color=[COLORS[b] for b in names], edgecolor="black", lw=0.6)
for i, b in enumerate(names):
    axL.text(mags[i] + 1.5, i, f'd={pb[b]["min_charged_dist_A"]} A   (screened {abs(pb[b]["E_total_ddd_kcal"]):.1f})',
             va="center", fontsize=8.5, color="#333")
    if "K459" in b:  # tag the two bridges K459A removes together
        axL.text(mags[i] - 2, i, "K459 hub", va="center", ha="right", fontsize=8.5,
                 color="white", weight="bold")
axL.set_yticks(list(y)); axL.set_yticklabels(names, fontsize=10)
axL.set_xlabel("lock strength  |E_interaction|  (kcal/mol, direct ff14SB)", fontsize=10)
axL.set_title("(a) Genu salt-bridge lock strengths — WT extended", fontsize=11, weight="bold")
axL.set_xlim(0, max(mags) * 1.28)
axL.grid(axis="x", alpha=0.25)

# ---- (b) lock energy removed per mutation = snap-back prediction ----
muts = d["destabilization_ranking"]  # most-removed first
vals = [pm[m]["lock_energy_removed_direct_kcal"] for m in muts]
bars = axR.bar(range(len(muts)), vals, color=[MUT_COLORS[m] for m in muts],
               edgecolor="black", lw=0.6, width=0.62)
for i, m in enumerate(muts):
    nb = pm[m]["n_broken"]; brk = ", ".join(pm[m]["bridges_broken"])
    axR.text(i, vals[i] + 4, f'{vals[i]:.0f}', ha="center", fontsize=10, weight="bold")
    axR.text(i, vals[i] / 2, f'{nb} lock{"s" if nb > 1 else ""}\nbroken', ha="center",
             va="center", fontsize=8.5, color="white", weight="bold")
    axR.text(i, -12, brk.replace(", ", "\n"), ha="center", va="top", fontsize=7.4, color="#444")
axR.set_xticks(range(len(muts)))
axR.set_xticklabels([m.replace("+", "\n+") for m in muts], fontsize=9.5)
axR.set_ylabel("genu lock energy removed  (kcal/mol, direct)", fontsize=10)
axR.set_title("(b) Destabilization by mutation → snap-back prediction", fontsize=11, weight="bold")
axR.set_ylim(-34, max(vals) * 1.2)
axR.axhline(0, color="black", lw=0.8)
axR.grid(axis="y", alpha=0.25)

fig.suptitle("Which residues energetically lock αVβ3 open — K459 is the hub (removes 2 locks, ~1.8× E598A)",
             fontsize=12.5, weight="bold", y=0.99)
fig.text(0.5, 0.008,
         "Direct ff14SB Coulomb+LJ interaction between full residues (WT extended, GB-OBC2 build); "
         "screened = distance-dependent dielectric eps=4r. Static upper bound — the pending A5000 snap-back MD tests the ranking dynamically.",
         ha="center", fontsize=7.6, color="#666")
fig.tight_layout(rect=(0, 0.03, 1, 0.96))
fig.savefig("figures/route_a_lock_energy.png", dpi=130)
print("wrote figures/route_a_lock_energy.png")
