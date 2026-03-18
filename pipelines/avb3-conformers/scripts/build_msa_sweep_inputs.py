#!/usr/bin/env python3
"""Build Protenix input JSONs for MSA-subsampled conformer validation.

For each MSA depth level, creates a Protenix input JSON that uses the
subsampled MSA for structure prediction. The resulting predictions are
then compared against pulled conformer frames to validate realism.

Usage:
    python build_msa_sweep_inputs.py \
        --msa-manifest msa_subsamples/manifest.json \
        --sequence-file data/avb3/sequences/avb3_sequence.fasta \
        --output-dir data/runs/avb3/msa_validation/inputs
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def extract_sequence_from_pdb(pdb_path: Path) -> dict[str, str]:
    """Extract sequences per chain from a PDB file using CA atoms."""
    chains: dict[str, list[tuple[int, str]]] = {}
    aa3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
    }
    seen = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        chain_id = line[21]
        resseq = int(line[22:26])
        resname = line[17:20].strip()
        key = (chain_id, resseq)
        if key in seen:
            continue
        seen.add(key)
        aa = aa3to1.get(resname, "X")
        chains.setdefault(chain_id, []).append((resseq, aa))

    result = {}
    for chain_id, residues in sorted(chains.items()):
        residues.sort()
        result[chain_id] = "".join(aa for _, aa in residues)
    return result


def build_protenix_input(job_name: str, chain_seqs: dict[str, str],
                         msa_paths: dict[str, str] | None = None) -> list[dict]:
    """Build a Protenix input JSON payload."""
    sequences = []
    for chain_id, seq in sorted(chain_seqs.items()):
        chain_spec: dict = {
            "proteinChain": {
                "count": 1,
                "sequence": seq,
                "modifications": [],
            }
        }
        # If MSA paths provided, add them
        if msa_paths and chain_id in msa_paths:
            chain_spec["proteinChain"]["msa"] = {
                "precomputed_msa_dir": "",
                "pairedMsaPath": "",
                "unpairedMsaPath": msa_paths[chain_id],
            }
        sequences.append(chain_spec)

    return [{
        "name": job_name,
        "covalent_bonds": [],
        "sequences": sequences,
    }]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-pdb", type=Path, required=True,
                        help="Reference AVB3 PDB (frame 0) for sequence extraction.")
    parser.add_argument("--msa-dir-a", type=Path, default=None,
                        help="Directory with subsampled A3M files for chain A.")
    parser.add_argument("--msa-dir-b", type=Path, default=None,
                        help="Directory with subsampled A3M files for chain B.")
    parser.add_argument("--msa-manifest", type=Path, default=None,
                        help="Manifest JSON from subsample_msa.py (for single-chain or shared MSA).")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for Protenix input JSONs.")
    parser.add_argument("--depths", default="1.0,0.5,0.25,0.1,0.05,0.01",
                        help="MSA depth fractions to generate inputs for.")
    parser.add_argument("--no-msa-baseline", action="store_true",
                        help="Also generate a no-MSA baseline input.")
    args = parser.parse_args()

    # Extract sequences from reference PDB
    chain_seqs = extract_sequence_from_pdb(args.reference_pdb)
    print(f"Extracted sequences from {args.reference_pdb}:")
    for chain_id, seq in chain_seqs.items():
        print(f"  Chain {chain_id}: {len(seq)} residues")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fractions = [float(f) for f in args.depths.split(",")]
    configs = []

    # Generate input for each MSA depth
    for frac in fractions:
        label = f"depth_{frac:.2f}"
        job_name = f"avb3_msa_{label}"

        # Build MSA path mapping if MSA dirs provided
        msa_paths = {}
        if args.msa_dir_a:
            msa_a = args.msa_dir_a / f"{label}.a3m"
            if msa_a.exists():
                msa_paths["A"] = str(msa_a.resolve())
        if args.msa_dir_b:
            msa_b = args.msa_dir_b / f"{label}.a3m"
            if msa_b.exists():
                msa_paths["B"] = str(msa_b.resolve())

        payload = build_protenix_input(job_name, chain_seqs,
                                       msa_paths if msa_paths else None)
        out_path = args.output_dir / f"{label}_input.json"
        out_path.write_text(json.dumps(payload, indent=2))
        configs.append({
            "fraction": frac, "label": label,
            "input_json": str(out_path), "job_name": job_name,
            "use_msa": bool(msa_paths),
        })
        print(f"  {label} → {out_path}")

    # No-MSA baseline
    if args.no_msa_baseline:
        payload = build_protenix_input("avb3_no_msa_baseline", chain_seqs)
        out_path = args.output_dir / "no_msa_input.json"
        out_path.write_text(json.dumps(payload, indent=2))
        configs.append({
            "fraction": 0.0, "label": "no_msa",
            "input_json": str(out_path), "job_name": "avb3_no_msa_baseline",
            "use_msa": False,
        })
        print(f"  no_msa → {out_path}")

    # Write sweep config
    sweep_path = args.output_dir / "sweep_config.json"
    sweep_path.write_text(json.dumps(configs, indent=2))
    print(f"\nSweep config: {sweep_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
