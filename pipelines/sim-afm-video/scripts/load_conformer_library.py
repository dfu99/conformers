#!/usr/bin/env python3
"""Load a conformer library into a single (T, N, 3) coordinate array.

Input contract:
- A directory containing PDB files for each conformer.
- A `library.json` (or `library.csv`) metadata file in that directory:

  {
    "conformers": [
      {"file": "frame_0000.pdb", "cv0_A": 52.3, "cv1_A": 42.1, "state": "BC"},
      {"file": "frame_0001.pdb", "cv0_A": 53.1, "cv1_A": 42.5, "state": "BC"},
      ...
    ],
    "trajectory_order": [0, 1, 2, ...],
    "adjacency": [[0, 1], [1, 2], ...],
    "topology_pdb": "frame_0000.pdb"
  }

`trajectory_order` is a sequence of indices defining playback order
(consecutive entries are biophysically adjacent). `adjacency` is an
optional more-general edge list. Either one is sufficient; if both are
present, `trajectory_order` takes precedence.

Output:
- coords.npy (T, N, 3) — atom xyz in nm for each frame in order
- topology.pdb — copied from metadata.topology_pdb
- order.npy — the playback order (same length as T)
- meta.json — passthrough of conformers list, in playback order

Usage:
    python load_conformer_library.py \\
        --library data/runs/avb3/conformers/all_frames_bent_extended_v7 \\
        --output results/sim_afm_pipeline/v7_input
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import mdtraj as md


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--library", type=Path, required=True,
                   help="Directory with PDBs + library.json")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--max-frames", type=int, default=0,
                   help="0 = all")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    meta_path = args.library / "library.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"No library.json in {args.library}")
    with meta_path.open() as f:
        meta = json.load(f)

    conformers = meta["conformers"]
    n = len(conformers)
    if "trajectory_order" in meta and meta["trajectory_order"]:
        order = meta["trajectory_order"]
    elif "adjacency" in meta and meta["adjacency"]:
        # Build a path through the graph (simple greedy).
        adj = {i: set() for i in range(n)}
        for a, b in meta["adjacency"]:
            adj[a].add(b)
            adj[b].add(a)
        # Start from node 0 and walk via DFS, breaking at dead ends.
        order = []
        seen = set()
        stack = [0]
        while stack:
            node = stack[-1]
            if node not in seen:
                seen.add(node)
                order.append(node)
                next_neighbors = sorted(adj[node] - seen)
                if next_neighbors:
                    stack.append(next_neighbors[0])
                else:
                    stack.pop()
            else:
                stack.pop()
        if len(order) < n:
            unseen = sorted(set(range(n)) - seen)
            order.extend(unseen)
    else:
        order = list(range(n))

    if args.max_frames > 0:
        order = order[: args.max_frames]
    print(f"Loading {len(order)} frames in playback order")

    topo_pdb_name = meta.get("topology_pdb", conformers[0]["file"])
    topo_path = args.library / topo_pdb_name
    if not topo_path.exists():
        raise FileNotFoundError(f"Topology PDB {topo_path} missing")
    topo_traj = md.load(str(topo_path))
    n_atoms = topo_traj.n_atoms
    print(f"Topology has {n_atoms} atoms")
    shutil.copy(str(topo_path), str(args.output / "topology.pdb"))

    coords = np.zeros((len(order), n_atoms, 3), dtype=np.float32)
    ordered_meta = []
    for k, idx in enumerate(order):
        entry = conformers[idx]
        pdb_path = args.library / entry["file"]
        if not pdb_path.exists():
            raise FileNotFoundError(f"Missing {pdb_path}")
        t = md.load(str(pdb_path))
        if t.n_atoms != n_atoms:
            raise ValueError(
                f"Atom count mismatch: {pdb_path} has {t.n_atoms}, "
                f"topology has {n_atoms}"
            )
        coords[k] = t.xyz[0]
        ordered_meta.append(entry)
        if (k + 1) % 100 == 0:
            print(f"  loaded {k+1}/{len(order)}")

    np.save(str(args.output / "coords.npy"), coords)
    np.save(str(args.output / "order.npy"), np.array(order, dtype=np.int64))
    with (args.output / "meta.json").open("w") as f:
        json.dump({"playback": ordered_meta, "n_frames": len(order)}, f, indent=2)
    print(f"Wrote {args.output}/coords.npy, order.npy, meta.json, topology.pdb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
