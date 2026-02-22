#!/usr/bin/env python3
"""Extract integrin alpha5/beta1 sequences from FASTA-like file."""

from __future__ import annotations

import argparse
import re
import unicodedata
from pathlib import Path


def parse_fasta_like(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    seq_parts: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or set(line) == {"="}:
            continue
        if line.startswith(">"):
            if header and seq_parts:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
            continue
        cleaned = re.sub(r"[^A-Za-z]", "", line).upper()
        if cleaned:
            seq_parts.append(cleaned)
    if header and seq_parts:
        records.append((header, "".join(seq_parts)))
    if not records:
        raise ValueError(f"No sequence records found in {path}")
    return records


def canonical(text: str) -> str:
    text = (
        text.replace("α", "alpha")
        .replace("β", "beta")
        .replace("Α", "alpha")
        .replace("Β", "beta")
    )
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    return re.sub(r"[^A-Za-z0-9]+", "", text).lower()


def pick(records: list[tuple[str, str]], target: str) -> tuple[str, str]:
    t = canonical(target)
    for name, seq in records:
        if canonical(name) == t:
            return name, seq
    names = ", ".join(n for n, _ in records)
    raise ValueError(f"Could not find target {target}. Available: {names}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-file", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--name-a", default="Integrin alpha5-Avi")
    parser.add_argument("--name-b", default="Integrin beta1-spycatcher")
    args = parser.parse_args()

    records = parse_fasta_like(args.sequence_file)
    name_a, seq_a = pick(records, args.name_a)
    name_b, seq_b = pick(records, args.name_b)

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "chain_A.seq").write_text(seq_a, encoding="utf-8")
    (args.outdir / "chain_B.seq").write_text(seq_b, encoding="utf-8")
    (args.outdir / "chain_A.name.txt").write_text(name_a + "\n", encoding="utf-8")
    (args.outdir / "chain_B.name.txt").write_text(name_b + "\n", encoding="utf-8")
    print(f"Wrote chain A ({len(seq_a)} aa): {name_a}")
    print(f"Wrote chain B ({len(seq_b)} aa): {name_b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
