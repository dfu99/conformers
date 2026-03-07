#!/usr/bin/env python3
"""Build staged Protenix input JSON files for A5B1 + tag attachment."""

from __future__ import annotations

import argparse
import json
import re
import unicodedata
from pathlib import Path
from typing import Dict, List, Tuple


def parse_fasta_like(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = None
    seq_parts: List[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or set(line) == {"="}:
            continue
        if line.startswith(">"):
            if header is not None and seq_parts:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
            continue
        cleaned = re.sub(r"[^A-Za-z]", "", line).upper()
        if cleaned:
            seq_parts.append(cleaned)
    if header is not None and seq_parts:
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


def pick_sequence(records: List[Tuple[str, str]], target_name: str) -> Tuple[str, str]:
    t = canonical(target_name)
    for name, seq in records:
        if canonical(name) == t:
            return name, seq
    available = ", ".join(name for name, _ in records)
    raise ValueError(f"Could not find component '{target_name}'. Available: {available}")


def make_payload(job_name: str, chains: List[Tuple[str, str]]) -> List[dict]:
    return [
        {
            "name": job_name,
            "covalent_bonds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "count": 1,
                        "sequence": seq,
                        "modifications": [],
                    }
                }
                for _, seq in chains
            ],
        }
    ]


def write_json(path: Path, payload: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sequence-file",
        type=Path,
        default=Path("data/a5b1/sequences/sequences_updated"),
        help="FASTA-like sequence file with all components.",
    )
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path("data/runs/a5b1/staged_attachment"),
        help="Workflow root for staged attachment artifacts.",
    )
    parser.add_argument("--alpha-name", default="Integrin alpha5-Avi")
    parser.add_argument("--beta-name", default="Integrin beta1-spycatcher")
    parser.add_argument("--spytag-name", default="Spytag")
    parser.add_argument("--streptavidin-name", default="Streptavidin")
    parser.add_argument("--stage1-job-name", default="a5b1_stage1_spytag")
    parser.add_argument("--stage2-job-name", default="a5b1_stage2_streptavidin")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    records = parse_fasta_like(args.sequence_file)
    alpha_name, alpha_seq = pick_sequence(records, args.alpha_name)
    beta_name, beta_seq = pick_sequence(records, args.beta_name)
    spytag_name, spytag_seq = pick_sequence(records, args.spytag_name)
    streptavidin_name, streptavidin_seq = pick_sequence(records, args.streptavidin_name)

    out_dir = args.workflow_dir.resolve() / "inputs" / "protenix"
    out_dir.mkdir(parents=True, exist_ok=True)

    stage1_chains = [
        (alpha_name, alpha_seq),
        (beta_name, beta_seq),
        (spytag_name, spytag_seq),
    ]
    stage2_chains = [
        (alpha_name, alpha_seq),
        (beta_name, beta_seq),
        (streptavidin_name, streptavidin_seq),
    ]

    stage1_json = out_dir / "stage1_spytag_input.json"
    stage2_json = out_dir / "stage2_streptavidin_input.json"

    write_json(stage1_json, make_payload(args.stage1_job_name, stage1_chains))
    write_json(stage2_json, make_payload(args.stage2_job_name, stage2_chains))

    chain_roles = {
        "stage1": {
            "job_name": args.stage1_job_name,
            "chain_roles": {"A": "alpha", "B": "beta", "C": "spytag"},
        },
        "stage2": {
            "job_name": args.stage2_job_name,
            "chain_roles": {"A": "alpha", "B": "beta", "C": "streptavidin"},
        },
    }
    write_json(out_dir / "chain_roles.json", chain_roles)

    print("Wrote staged Protenix inputs:")
    print(f"  {stage1_json}")
    print(f"  {stage2_json}")
    print(f"  {out_dir / 'chain_roles.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
