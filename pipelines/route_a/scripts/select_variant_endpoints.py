#!/usr/bin/env python3
"""Choose, per integrin variant, the BENT structure to morph from -- and any
experimentally-determined EXTENDED structure available to validate against.

Two separate jobs:

1. Bent start state. The alphaVbeta3 recipe morphs *from a bent crystal*, so the
   start must be the most compact full ectodomain available, not merely the most
   complete one. Ranking rule: smallest long-axis extent, subject to the
   completeness floors already applied during discovery. One documented exception:
   alphaVbeta3 is pinned to 1JV2 so this run reproduces the established route-A
   endpoint (obj-075: Rg 39 -> 67 A, extent 130 -> 211 A) and thereby acts as a
   positive control on the generalised code path.

2. Extended validation references. Discovery turned up alphaVbeta3 cryo-EM
   entries at 196-214 A long-axis extent (8XEZ / 8XEN, "conformation 5 / 4") and
   comparably elongated alphaXbeta2 and alphaIIbbeta3 entries. tasks/planning.md
   still records "no extended alphaVbeta3 crystal exists (obj-041)", which was
   true for crystallography but is no longer true for cryo-EM. Any structure whose
   extent exceeds EXTENDED_EXTENT_A is therefore flagged as a candidate
   experimental extended state, to be compared against the morph product rather
   than used as its input.
"""
import argparse
import json
import os

# a bent ectodomain spans ~110-140 A; anything much longer has already unbent
EXTENDED_EXTENT_A = 150.0

# on the validated primary coordinate, bent alphaVbeta3 sits at ~41 deg and the
# morphed extended endpoint at ~174 deg; treat a wide knee as already-extended
EXTENDED_GENU_DEG = 100.0

# completeness floors, relative to the best structure available for that variant
MIN_ALPHA_FRAC = 0.90
MIN_BETA_FRAC = 0.80

# leg-module local-alignment score below which the transferred genu is only a
# homology guess: alphaVbeta3 scores ~2500/1200 (alpha/beta), alphaIIbbeta3 663/1210,
# alpha5beta1 888/353, but the beta2 integrins bottom out near 140/244
MIN_LEG_SCORE_HIGH = 300.0

# pinned starts, with reasons (kept explicit so the choice is auditable)
PINNED_BENT = {
    "alphaVbeta3": ("1JV2", "established route-A state A; reproduces obj-072..081 "
                            "so the generalised pipeline has a positive control"),
}


def emit_all(disc, out_json):
    """Mode 1: every accepted structure, keyed variant__PDBID, for a genu-angle survey.

    Long-axis extent is only a proxy for bentness -- it also grows when a leg
    straightens without the knee opening. The genu angle is the validated primary
    coordinate (obj-076), but it can only be measured once boundaries have been
    transferred. So emit all candidates here, run them through
    map_domain_boundaries.py, and let choose_from_boundaries() pick on the real CV.
    """
    reps = {f"{r['variant']}__{r['pdb_id']}": r for r in disc["accepted"]}
    with open(out_json, "w") as fh:
        json.dump({"_purpose": "all accepted structures, for genu-angle survey",
                   "representatives": reps}, fh, indent=2)
    print(f"wrote {out_json} ({len(reps)} structures for survey)")


def choose_from_boundaries(disc, bounds_path, out_json):
    """Mode 2: pick each variant's bent start by measured genu angle."""
    bounds = json.load(open(bounds_path))

    # Completeness guard. The genu angle is a centroid construction, so a chain
    # missing a couple of hundred residues reports a knee angle for a molecule that
    # is partly absent -- alphaXbeta2 3K72 (884 of 1084 alpha CA) comes back at 3 deg
    # while carrying the LARGEST long-axis extent of its variant, which cannot both
    # be true. Require each candidate to be nearly as complete as the best structure
    # available for the same variant before its angle is allowed to decide anything.
    max_ca = {}
    for key, spec in bounds.items():
        v = key.split("__")[0]
        a, bta = spec["alpha"]["observed_ca"], spec["beta"]["observed_ca"]
        cur = max_ca.setdefault(v, [0, 0])
        cur[0], cur[1] = max(cur[0], a), max(cur[1], bta)

    surveyed = {}
    for key, spec in bounds.items():
        variant = key.split("__")[0]
        ang = spec.get("bent_cvs", {}).get("genu_angle_deg")
        if ang is None or not spec.get("usable_for_morph"):
            continue
        a_ca, b_ca = spec["alpha"]["observed_ca"], spec["beta"]["observed_ca"]
        a_frac = a_ca / max(1, max_ca[variant][0])
        b_frac = b_ca / max(1, max_ca[variant][1])
        aln = spec.get("alignment", {})
        a_mod = aln.get("alpha_module_scores", {}) or {}
        b_mod = aln.get("beta_module_scores", {}) or {}
        leg_score = min(a_mod.get("leg", 0), b_mod.get("leg", 0))
        surveyed.setdefault(variant, []).append({
            "pdb_id": spec["pdb_id"], "title": spec["title"],
            "genu_angle_deg": ang,
            "long_extent_A": spec["bent_cvs"].get("long_extent_A"),
            "Rg_A": spec["bent_cvs"].get("Rg_A"),
            "end_to_end_A": spec["bent_cvs"].get("end_to_end_A"),
            "alpha_completeness": round(a_frac, 3),
            "beta_completeness": round(b_frac, 3),
            "leg_alignment_score": leg_score,
            "boundary_confidence": ("high" if leg_score >= MIN_LEG_SCORE_HIGH
                                    else "low"),
            "complete_enough": bool(a_frac >= MIN_ALPHA_FRAC and b_frac >= MIN_BETA_FRAC),
            "survey_key": key,
        })

    by_variant = {}
    for r in disc["accepted"]:
        by_variant.setdefault(r["variant"], []).append(r)

    starts, notes = {}, {}
    for variant, rows in sorted(by_variant.items()):
        allc = sorted(surveyed.get(variant, []), key=lambda d: d["genu_angle_deg"])
        cands = [c for c in allc if c["complete_enough"]] or allc
        if not cands:
            continue
        dropped = [c["pdb_id"] for c in allc if not c["complete_enough"]]
        if variant in PINNED_BENT:
            pdb_id, reason = PINNED_BENT[variant]
            pick = next((c for c in allc if c["pdb_id"] == pdb_id), cands[0])
            notes[variant] = (f"pinned to {pick['pdb_id']}: {reason} "
                              f"(genu {pick['genu_angle_deg']} deg; most-bent complete "
                              f"candidate was {cands[0]['pdb_id']} at "
                              f"{cands[0]['genu_angle_deg']} deg)")
        else:
            pick = cands[0]
            notes[variant] = (f"smallest genu angle ({pick['genu_angle_deg']} deg) "
                              f"of {len(cands)} sufficiently complete candidates")
        if dropped:
            notes[variant] += f"; excluded as incomplete: {','.join(dropped)}"
        if pick["boundary_confidence"] == "low":
            notes[variant] += ("; BOUNDARY CONFIDENCE LOW -- leg-module alignment to "
                               f"alphaVbeta3 scores only {pick['leg_alignment_score']}, "
                               "so the transferred genu is a homology guess and the "
                               "morph must be judged on its own trajectory")
        rec = next(r for r in rows if r["pdb_id"] == pick["pdb_id"])
        starts[variant] = {**rec,
                           "start_genu_angle_deg": pick["genu_angle_deg"],
                           "boundary_confidence": pick["boundary_confidence"],
                           "leg_alignment_score": pick["leg_alignment_score"],
                           "survey_key": pick["survey_key"]}

    extended_refs = {}
    for variant, cands in surveyed.items():
        wide = sorted([c for c in cands if c["genu_angle_deg"] >= EXTENDED_GENU_DEG],
                      key=lambda c: -c["genu_angle_deg"])
        if wide:
            extended_refs[variant] = wide

    result = {
        "_purpose": "bent morph inputs chosen on measured genu angle + experimental "
                    "extended references",
        "extended_genu_threshold_deg": EXTENDED_GENU_DEG,
        "representatives": starts,
        "selection_notes": notes,
        "genu_survey": {k: sorted(v, key=lambda d: d["genu_angle_deg"])
                        for k, v in surveyed.items()},
        "extended_references": extended_refs,
    }
    with open(out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"{'variant':<16}{'start':<7}{'genu':>7}{'extent':>8}{'conf':>7}   selection")
    for v, r in sorted(starts.items()):
        print(f"{v:<16}{r['pdb_id']:<7}{r['start_genu_angle_deg']:>7}"
              f"{r['long_axis_extent_A']:>8}{r['boundary_confidence']:>7}   {notes[v]}")
    print("\nGenu-angle survey (all usable structures):")
    for v, lst in sorted(result["genu_survey"].items()):
        line = "  ".join(f"{c['pdb_id']}:{c['genu_angle_deg']:.0f}" for c in lst)
        print(f"  {v:<16}{line}")
    print("\nExperimental extended references (validation targets, NOT inputs):")
    for v, lst in sorted(extended_refs.items()):
        for c in lst:
            print(f"  {v:<16}{c['pdb_id']:<7}genu={c['genu_angle_deg']:>6} deg  "
                  f"extent={c['long_extent_A']} A  {c['title'][:52]}")
    print(f"\nwrote {out_json}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovered", default="results/variants/discovered_ectodomains.json")
    ap.add_argument("--out-json", default="results/variants/variant_endpoints.json")
    ap.add_argument("--all", action="store_true",
                    help="emit every accepted structure for the genu-angle survey")
    ap.add_argument("--from-boundaries", default="",
                    help="survey boundaries JSON; pick starts by measured genu angle")
    args = ap.parse_args()

    disc = json.load(open(args.discovered))
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)

    if args.all:
        emit_all(disc, args.out_json)
        return
    if args.from_boundaries:
        choose_from_boundaries(disc, args.from_boundaries, args.out_json)
        return

    accepted = disc["accepted"]

    by_variant = {}
    for r in accepted:
        by_variant.setdefault(r["variant"], []).append(r)

    starts, extended_refs, notes = {}, {}, {}
    for variant, rows in sorted(by_variant.items()):
        bent_pool = [r for r in rows if r["long_axis_extent_A"] < EXTENDED_EXTENT_A]
        pool = bent_pool or rows
        if variant in PINNED_BENT:
            pdb_id, reason = PINNED_BENT[variant]
            match = [r for r in rows if r["pdb_id"] == pdb_id]
            if match:
                chosen = match[0]
                notes[variant] = f"pinned to {pdb_id}: {reason}"
            else:
                chosen = min(pool, key=lambda r: r["long_axis_extent_A"])
                notes[variant] = (f"pinned {pdb_id} not among accepted structures; "
                                  f"fell back to most-bent {chosen['pdb_id']}")
        else:
            chosen = min(pool, key=lambda r: r["long_axis_extent_A"])
            notes[variant] = (f"most compact full ectodomain "
                              f"({chosen['long_axis_extent_A']} A) of {len(rows)} candidates")
            if not bent_pool:
                notes[variant] += (" -- WARNING: every candidate exceeds "
                                   f"{EXTENDED_EXTENT_A} A, so this start is not fully bent")
        starts[variant] = chosen

        ext = sorted([r for r in rows if r["long_axis_extent_A"] >= EXTENDED_EXTENT_A],
                     key=lambda r: -r["long_axis_extent_A"])
        if ext:
            extended_refs[variant] = [
                {"pdb_id": r["pdb_id"], "title": r["title"],
                 "long_axis_extent_A": r["long_axis_extent_A"],
                 "pdb_path": r["pdb_path"],
                 "alpha_chain": r["alpha"]["chain"], "beta_chain": r["beta"]["chain"]}
                for r in ext]

    result = {
        "_purpose": "bent morph inputs + experimental extended references per variant",
        "extended_extent_threshold_A": EXTENDED_EXTENT_A,
        "representatives": starts,          # consumed by map_domain_boundaries.py
        "selection_notes": notes,
        "extended_references": extended_refs,
    }
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"{'variant':<16}{'start':<7}{'extent':>8}   selection")
    for v, r in sorted(starts.items()):
        print(f"{v:<16}{r['pdb_id']:<7}{r['long_axis_extent_A']:>8}   {notes[v]}")
    print("\nExperimental extended references (validation targets, NOT inputs):")
    if not extended_refs:
        print("  none")
    for v, lst in sorted(extended_refs.items()):
        for r in lst:
            print(f"  {v:<16}{r['pdb_id']:<7}{r['long_axis_extent_A']:>8} A  {r['title'][:64]}")
    print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
