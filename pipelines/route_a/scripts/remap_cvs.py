#!/usr/bin/env python3
"""Route-A day-1: remap αIIb-numbered CV definitions to αV numbering.

Audit §12.7 + docs/route_a_kickoff.md Risk-1 mitigation. Consumes
results/route_a/av_aiib_alignment.json (the BLOSUM62 αV(1JV2) ↔
αIIb(3FCS) per-position lookup) and rewrites any CV definition
file using αIIb residues into αV residues.

CV definitions are read as JSON of the form:

  {
    "cvs": [
      {"name": "cv0", "type": "centroid_distance",
       "domain_a": {"chain": "A", "residues": "1-452"},
       "domain_b": {"chain": "A", "residues": "455-620"}}
    ]
  }

Chain B definitions pass through unchanged because both αVβ3 and
αIIbβ3 share the β3 chain.

Usage:
  python remap_cvs.py \
      --alignment results/route_a/av_aiib_alignment.json \
      --in   ferg_lab_cvs.json \
      --out  pipelines/route_a/cvs_av.json
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_residue_range(spec: str) -> tuple[int, int]:
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", spec)
    if m:
        return int(m.group(1)), int(m.group(2))
    raise ValueError(f"Cannot parse residue range: {spec!r}")


def load_alignment(path: Path) -> dict[int, int]:
    """Return αIIb→αV resSeq dict.

    The alignment JSON maps αV positions to αIIb in the `landmarks`
    list, so we invert it across the full alignment by walking the
    pairwise lookup; for landmarks we have direct entries, for
    unlandmarked residues we interpolate using neighbouring deltas.
    """
    with path.open() as f:
        data = json.load(f)

    # Build αV → αIIb from the landmarks (sparse, exact at landmarks)
    landmarks = sorted(
        ((lm["av_resSeq"], lm["aiib_resSeq"])
         for lm in data["landmarks"]
         if lm["aiib_resSeq"] is not None),
        key=lambda p: p[0],
    )
    if not landmarks:
        raise RuntimeError("No usable landmarks in alignment JSON")

    # For positions between landmarks, linear-interpolate the offset
    av_to_aiib = {}
    for i, (av_a, aiib_a) in enumerate(landmarks):
        av_to_aiib[av_a] = aiib_a
        if i + 1 < len(landmarks):
            av_b, aiib_b = landmarks[i + 1]
            for av_x in range(av_a + 1, av_b):
                # Linear interpolation between offset_a and offset_b
                f = (av_x - av_a) / (av_b - av_a)
                aiib_x = int(round(aiib_a + f * (aiib_b - aiib_a)))
                av_to_aiib[av_x] = aiib_x

    # Inverse: αIIb → αV (single-valued per αIIb for our purposes)
    aiib_to_av = {}
    for av_pos, aiib_pos in av_to_aiib.items():
        if aiib_pos is None:
            continue
        # Multiple αV positions may map to same αIIb → keep the closest match
        if aiib_pos not in aiib_to_av:
            aiib_to_av[aiib_pos] = av_pos
        else:
            # If a duplicate, keep the first (the gap creates one-to-many;
            # we resolve by smallest αV resSeq for determinism).
            if av_pos < aiib_to_av[aiib_pos]:
                aiib_to_av[aiib_pos] = av_pos
    return aiib_to_av


def remap_residue_range(spec: str, aiib_to_av: dict[int, int]) -> str:
    """Convert e.g. '1-452' (αIIb) to the αV equivalent."""
    lo, hi = parse_residue_range(spec)
    av_lo = aiib_to_av.get(lo)
    av_hi = aiib_to_av.get(hi)
    if av_lo is None:
        # Step inward until we find a mapped αIIb pos
        while lo <= hi and av_lo is None:
            lo += 1
            av_lo = aiib_to_av.get(lo)
    if av_hi is None:
        while hi >= lo and av_hi is None:
            hi -= 1
            av_hi = aiib_to_av.get(hi)
    if av_lo is None or av_hi is None:
        raise RuntimeError(f"Cannot remap residue range {spec}: no αIIb→αV "
                           f"mapping in [{lo}, {hi}]")
    return f"{min(av_lo, av_hi)}-{max(av_lo, av_hi)}"


def remap_cv(cv: dict, aiib_to_av: dict[int, int]) -> dict:
    out = dict(cv)
    for key in ("domain_a", "domain_b"):
        d = cv.get(key)
        if not d:
            continue
        new = dict(d)
        chain = d.get("chain", "A")
        if chain == "A":
            new["residues"] = remap_residue_range(d["residues"], aiib_to_av)
            new["_remapped_from_aiib"] = d["residues"]
        # Chain B (β3) shared between αVβ3 and αIIbβ3 → no remap
        out[key] = new
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--alignment", type=Path,
                   default=Path("results/route_a/av_aiib_alignment.json"))
    p.add_argument("--in", dest="input", type=Path, required=True,
                   help="αIIb-numbered CV definition JSON")
    p.add_argument("--out", type=Path, required=True,
                   help="Output αV-numbered CV definition JSON")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    aiib_to_av = load_alignment(args.alignment)
    print(f"Loaded {len(aiib_to_av)} αIIb→αV residue mappings")

    with args.input.open() as f:
        data = json.load(f)
    cvs = data.get("cvs", [])
    if not cvs:
        raise RuntimeError(f"No 'cvs' list found in {args.input}")
    print(f"Remapping {len(cvs)} CV definitions ...")

    remapped = [remap_cv(cv, aiib_to_av) for cv in cvs]
    out = {**data, "cvs": remapped,
           "_alignment_source": str(args.alignment),
           "_remap_method": "BLOSUM62 αV-αIIb pairwise + linear interpolation between landmarks"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        json.dump(out, f, indent=2)
    print(f"Wrote {args.out}")

    # Echo the residue swaps for the kickoff log
    for old, new in zip(cvs, remapped):
        for key in ("domain_a", "domain_b"):
            if old.get(key, {}).get("chain") == "A":
                old_r = old[key].get("residues")
                new_r = new[key].get("residues")
                print(f"  {old.get('name', '?')}.{key}: αIIb {old_r}  →  αV {new_r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
