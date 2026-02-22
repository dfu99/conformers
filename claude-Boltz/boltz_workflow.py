#!/usr/bin/env python3
"""Boltz-1 conformation sampling for integrin alpha/beta heterodimer.

Boltz-1 (MIT, 2024) is a fully independent diffusion-based structure
prediction model — separate training data, separate architecture, separate
inference engine from AF2/Protenix. Running many seeds × diffusion samples
gives broad conformational coverage without template or training-bias
toward the bent integrin state.

This script is self-contained: it reads sequences from a PDB, writes a
Boltz YAML input, and drives `boltz predict` directly.

Install:
    pip install boltz

Docs / source:
    https://github.com/jwohlwend/boltz
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

try:
    import yaml  # PyYAML, pulled in by boltz
except ImportError:
    yaml = None  # type: ignore  — handled at runtime


# ---------------------------------------------------------------------------
# PDB sequence parser
# ---------------------------------------------------------------------------

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K", "ASX": "X", "GLX": "X", "UNK": "X",
}


def parse_pdb_sequences(pdb_path: Path) -> Dict[str, str]:
    chains: Dict[str, List[str]] = {}
    seen: Dict[str, set] = {}
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM  ") or line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A", "1"):
                continue
            chain_id = line[21].strip() or "_"
            res_key = (line[22:26].strip(), line[26].strip())
            chains.setdefault(chain_id, [])
            seen.setdefault(chain_id, set())
            if res_key not in seen[chain_id]:
                seen[chain_id].add(res_key)
                chains[chain_id].append(
                    AA3_TO_AA1.get(line[17:20].strip().upper(), "X")
                )
    if not chains:
        raise ValueError(f"No CA ATOM records in {pdb_path}")
    return {c: "".join(s) for c, s in chains.items()}


# ---------------------------------------------------------------------------
# Boltz YAML writer
# ---------------------------------------------------------------------------

def write_boltz_yaml(
    path: Path,
    chain_order: List[str],
    chain_to_seq: Dict[str, str],
) -> None:
    """Write a minimal Boltz-1 YAML input (sequences only, no MSA override).

    Boltz runs its own internal MSA search (via ColabFold databases or
    bundled MMseqs2) — we don't need to pre-supply an MSA.
    """
    if yaml is None:
        raise ImportError("PyYAML not installed. Run: pip install boltz")

    sequences = [
        {"protein": {"id": cid, "sequence": chain_to_seq[cid]}}
        for cid in chain_order
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        yaml.dump(
            {"version": 1, "sequences": sequences},
            fh,
            default_flow_style=False,
            sort_keys=False,
        )


# ---------------------------------------------------------------------------
# Boltz runner
# ---------------------------------------------------------------------------

def run_boltz_seed(
    input_yaml: Path,
    out_dir: Path,
    seed: int,
    diffusion_samples: int,
    recycling_steps: int,
    sampling_steps: int,
    devices: int,
    accelerator: str,
    use_msa_server: bool,
) -> None:
    """Run `boltz predict` for a single seed."""
    cmd = [
        "boltz", "predict", str(input_yaml),
        "--out_dir", str(out_dir),
        "--seed", str(seed),
        "--diffusion_samples", str(diffusion_samples),
        "--recycling_steps", str(recycling_steps),
        "--sampling_steps", str(sampling_steps),
        "--devices", str(devices),
        "--accelerator", accelerator,
        "--output_format", "mmcif",
    ]
    if use_msa_server:
        cmd.append("--use_msa_server")

    print("\n[run]", " ".join(cmd))
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Boltz-1 independent conformation sampling for integrin."
    )
    p.add_argument("--input_pdb", type=Path, required=True,
                   help="PDB to extract chain sequences from.")
    p.add_argument("--chain_order", type=str, default="A,B",
                   help="Comma-separated chain IDs to include.")
    p.add_argument("--workflow_dir", type=Path, required=True,
                   help="Root for inputs and outputs.")
    p.add_argument("--job_name", type=str, default="integrin_boltz")
    p.add_argument("--seeds", type=str, default="101,202,303,404,505,606,707,808,909",
                   help="Comma-separated seeds. One independent boltz run per seed.")
    p.add_argument("--diffusion_samples", type=int, default=5,
                   help="Diffusion samples per seed.")
    p.add_argument("--recycling_steps", type=int, default=3)
    p.add_argument("--sampling_steps", type=int, default=200)
    p.add_argument("--devices", type=int, default=1)
    p.add_argument("--accelerator", type=str, default="gpu",
                   choices=("gpu", "cpu"))
    p.add_argument("--use_msa_server", action="store_true",
                   help="Use ColabFold MSA server for on-the-fly MSA search "
                        "(requires internet; omit on cluster).")
    p.add_argument("--run", action="store_true",
                   help="Execute boltz predict (else just write YAML).")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    chain_order = [c.strip() for c in args.chain_order.split(",") if c.strip()]
    input_pdb = args.input_pdb.resolve()
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

    workflow_dir = args.workflow_dir.resolve()
    inputs_dir = workflow_dir / "inputs"
    outputs_dir = workflow_dir / "outputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    chain_to_seq = parse_pdb_sequences(input_pdb)
    print("Chains from PDB:")
    for cid in chain_order:
        if cid not in chain_to_seq:
            raise ValueError(f"Chain '{cid}' missing. Available: {sorted(chain_to_seq)}")
        print(f"  {cid}: len={len(chain_to_seq[cid])}")

    input_yaml = inputs_dir / f"{args.job_name}.yaml"
    write_boltz_yaml(input_yaml, chain_order, chain_to_seq)
    print(f"\nWrote Boltz input YAML: {input_yaml}")

    seeds = [s.strip() for s in args.seeds.split(",") if s.strip()]
    total = len(seeds) * args.diffusion_samples
    print(f"{len(seeds)} seeds × {args.diffusion_samples} samples = {total} conformations")

    if not args.run:
        print("Dry run — pass --run to execute.")
        return 0

    import shutil
    if not shutil.which("boltz"):
        raise FileNotFoundError(
            "boltz not found in PATH.\n"
            "Install with: pip install boltz"
        )

    for seed in seeds:
        run_boltz_seed(
            input_yaml=input_yaml,
            out_dir=outputs_dir / f"seed_{seed}",
            seed=int(seed),
            diffusion_samples=args.diffusion_samples,
            recycling_steps=args.recycling_steps,
            sampling_steps=args.sampling_steps,
            devices=args.devices,
            accelerator=args.accelerator,
            use_msa_server=args.use_msa_server,
        )

    cif_count = sum(1 for _ in outputs_dir.rglob("*.cif"))
    print(f"\nDone. {cif_count} .cif files under {outputs_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
