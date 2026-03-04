#!/usr/bin/env python3
"""Rank Boltz outputs by an extension proxy (max CA span)."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable, List, Tuple


def ca_positions(struct_path: Path) -> List[Tuple[str, float, float, float]]:
    try:
        import gemmi  # type: ignore
    except Exception as exc:
        raise RuntimeError("gemmi is required. Install with: pip install gemmi") from exc

    structure = gemmi.read_structure(str(struct_path))
    model = structure[0]
    pts: List[Tuple[str, float, float, float]] = []
    for chain in model:
        for res in chain:
            atom = res.find_atom("CA", "\0")
            if atom is None:
                continue
            pos = atom.pos
            pts.append((chain.name, float(pos.x), float(pos.y), float(pos.z)))
    return pts


def dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def max_span(points: Iterable[Tuple[str, float, float, float]]) -> float:
    xyz = [(x, y, z) for _, x, y, z in points]
    if len(xyz) < 2:
        return 0.0
    best = 0.0
    for i in range(len(xyz)):
        ai = xyz[i]
        for j in range(i + 1, len(xyz)):
            d = dist(ai, xyz[j])
            if d > best:
                best = d
    return best


def find_confidence_json(struct_path: Path) -> Path | None:
    m = re.search(r"sample_(\d+)", struct_path.name)
    if not m:
        return None
    sample_idx = m.group(1)
    candidates = list(struct_path.parent.glob(f"*confidence*sample_{sample_idx}.json"))
    if not candidates:
        return None
    return candidates[0]


def read_confidence(path: Path | None) -> Tuple[str, str]:
    if path is None or not path.exists():
        return ("", "")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return ("", "")
    iptm = data.get("iptm", data.get("i_ptm", ""))
    ptm = data.get("ptm", data.get("p_tm", ""))
    return (str(iptm), str(ptm))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pred-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    struct_files = sorted(list(args.pred_root.rglob("*.cif")) + list(args.pred_root.rglob("*.pdb")))
    if not struct_files:
        raise RuntimeError(f"No structure files found under {args.pred_root}")

    rows = []
    for path in struct_files:
        points = ca_positions(path)
        span = max_span(points)
        conf = find_confidence_json(path)
        iptm, ptm = read_confidence(conf)
        rows.append((span, path, len(points), iptm, ptm))
    rows.sort(key=lambda x: x[0], reverse=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        handle.write("rank\tmax_ca_span\tca_atoms\tiptm\tptm\tstructure_path\n")
        for i, (span, path, n_ca, iptm, ptm) in enumerate(rows, start=1):
            handle.write(f"{i}\t{span:.3f}\t{n_ca}\t{iptm}\t{ptm}\t{path}\n")
    print(f"[rank] wrote {len(rows)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
