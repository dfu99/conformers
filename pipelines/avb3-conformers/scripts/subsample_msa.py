#!/usr/bin/env python3
"""Subsample an A3M MSA file at different depths for conformer validation.

Creates multiple A3M files with progressively fewer sequences, preserving
the query sequence (first entry) in all subsamples.

Usage:
    python subsample_msa.py --input full.a3m --output-dir msa_subsamples/ \
        --fractions 1.0,0.5,0.25,0.1,0.05,0.01
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path


def parse_a3m(path: Path) -> list[tuple[str, str]]:
    """Parse A3M file into list of (header, sequence) tuples."""
    records: list[tuple[str, str]] = []
    header = None
    seq_parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line
            seq_parts = []
        elif header is not None:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def write_a3m(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f"{header}\n{seq}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Full MSA in A3M format.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory for subsampled A3M files.")
    parser.add_argument("--fractions", default="1.0,0.5,0.25,0.1,0.05,0.01",
                        help="Comma-separated fractions of MSA to keep.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility.")
    args = parser.parse_args()

    records = parse_a3m(args.input)
    if not records:
        print(f"ERROR: No sequences found in {args.input}")
        return 1

    query = records[0]
    others = records[1:]
    total = len(others)
    print(f"Full MSA: {total + 1} sequences (1 query + {total} hits)")

    rng = random.Random(args.seed)
    fractions = [float(f) for f in args.fractions.split(",")]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []

    for frac in sorted(fractions, reverse=True):
        n_keep = max(1, int(total * frac)) if frac < 1.0 else total
        label = f"depth_{frac:.2f}"
        out_path = args.output_dir / f"{label}.a3m"

        sampled = rng.sample(others, min(n_keep, total))
        write_a3m(out_path, [query] + sampled)

        manifest.append({"fraction": frac, "n_sequences": len(sampled) + 1,
                         "path": str(out_path)})
        print(f"  {label}: {len(sampled) + 1} sequences → {out_path}")

    # Write manifest
    import json
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nManifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
