#!/usr/bin/env python3
"""Extract one-letter chain sequences from a PDB into .seq files."""

from __future__ import annotations

import argparse
from pathlib import Path


AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}


def parse_pdb(path: Path) -> dict[str, str]:
    chains: dict[str, list[str]] = {}
    seen: dict[str, set[tuple[str, str]]] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain = (line[21].strip() or "_")
        resname = line[17:20].strip().upper()
        resseq = line[22:26].strip()
        icode = line[26].strip()
        key = (resseq, icode)
        chains.setdefault(chain, [])
        seen.setdefault(chain, set())
        if key in seen[chain]:
            continue
        seen[chain].add(key)
        chains[chain].append(AA3_TO_AA1.get(resname, "X"))
    return {k: "".join(v) for k, v in chains.items()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--chain-a", default="A")
    parser.add_argument("--chain-b", default="B")
    args = parser.parse_args()

    seqs = parse_pdb(args.pdb)
    for chain in (args.chain_a, args.chain_b):
        if chain not in seqs:
            raise ValueError(f"Chain {chain} not found in {args.pdb}. Available: {sorted(seqs)}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "chain_A.seq").write_text(seqs[args.chain_a], encoding="utf-8")
    (args.outdir / "chain_B.seq").write_text(seqs[args.chain_b], encoding="utf-8")
    print(f"Wrote {args.outdir / 'chain_A.seq'} ({len(seqs[args.chain_a])} aa)")
    print(f"Wrote {args.outdir / 'chain_B.seq'} ({len(seqs[args.chain_b])} aa)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
