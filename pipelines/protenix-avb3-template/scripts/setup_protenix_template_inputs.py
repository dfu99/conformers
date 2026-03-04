#!/usr/bin/env python3
"""Generate a minimal Protenix template inference input from a seed PDB.

This script does setup only:
1) Parse chain sequences from input PDB.
2) Build per-chain template .a3m files.
3) Convert PDB -> mmCIF for template structure.
4) Register template mmCIF entry in Protenix mmcif DB layout.
5) Write a single template-enabled Protenix input JSON.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Sequence


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
    "SEC": "C",
    "PYL": "K",
    "ASX": "X",
    "GLX": "X",
    "UNK": "X",
}

ENTRY_ID_RE = re.compile(r"^[A-Za-z0-9_]{3,}$")


def parse_chain_order(value: str) -> List[str]:
    chains = [c.strip() for c in value.split(",") if c.strip()]
    if not chains:
        raise ValueError("--chain_order must include at least one chain.")
    return chains


def parse_pdb_sequences(pdb_path: Path) -> Dict[str, str]:
    chain_to_res: Dict[str, List[str]] = {}
    seen: Dict[str, set] = {}

    with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            if line[12:16].strip() != "CA":
                continue
            altloc = line[16]
            if altloc not in (" ", "A", "1"):
                continue

            chain_id = line[21].strip() or "_"
            resname = line[17:20].strip().upper()
            resseq = line[22:26].strip()
            icode = line[26].strip()
            residue_key = (resseq, icode)

            if chain_id not in chain_to_res:
                chain_to_res[chain_id] = []
                seen[chain_id] = set()
            if residue_key in seen[chain_id]:
                continue
            seen[chain_id].add(residue_key)
            chain_to_res[chain_id].append(AA3_TO_AA1.get(resname, "X"))

    if not chain_to_res:
        raise ValueError(f"No protein CA records found in {pdb_path}")
    return {chain: "".join(res) for chain, res in chain_to_res.items()}


def resolve_mkdssp(explicit: str) -> str | None:
    if explicit:
        path = Path(explicit).resolve()
        if not path.exists():
            raise FileNotFoundError(f"--mkdssp_binary_path does not exist: {path}")
        return str(path)
    found = shutil.which("mkdssp")
    return found


def convert_pdb_to_mmcif(
    pdb_path: Path,
    cif_path: Path,
    converter: str,
    mkdssp_binary_path: str,
) -> str:
    cif_path.parent.mkdir(parents=True, exist_ok=True)
    if converter not in {"auto", "mkdssp", "gemmi"}:
        raise ValueError(f"Unsupported converter: {converter}")

    if converter in {"auto", "mkdssp"}:
        mkdssp = resolve_mkdssp(mkdssp_binary_path)
        if mkdssp:
            cmd = [mkdssp, "--output-format", "mmcif", str(pdb_path), str(cif_path)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0 and cif_path.exists() and cif_path.stat().st_size > 0:
                return "mkdssp"
            if converter == "mkdssp":
                stderr = result.stderr.strip() or "(no stderr)"
                raise RuntimeError(
                    "mkdssp conversion failed.\n"
                    f"Command: {' '.join(cmd)}\n"
                    f"stderr: {stderr}"
                )
        elif converter == "mkdssp":
            raise RuntimeError(
                "--template_converter mkdssp was requested, but mkdssp was not found."
            )

    try:
        import gemmi  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "PDB->mmCIF conversion failed: mkdssp unavailable and gemmi not installed.\n"
            "Install gemmi (`pip install gemmi`) or provide mkdssp."
        ) from exc

    structure = gemmi.read_structure(str(pdb_path))
    structure.make_mmcif_document().write_file(str(cif_path))
    return "gemmi"


def resolve_template_mmcif_dir(explicit_dir: Path | None) -> Path:
    if explicit_dir is not None:
        return explicit_dir.resolve()
    root = os.environ.get("PROTENIX_ROOT_DIR", "").strip()
    if not root:
        raise RuntimeError(
            "Template mmcif DB root is unknown.\n"
            "Set --template_mmcif_dir or export PROTENIX_ROOT_DIR."
        )
    return (Path(root) / "mmcif").resolve()


def register_template_mmcif_entry(
    template_cif: Path,
    template_entry_id: str,
    mmcif_dir: Path,
) -> Path:
    if not template_cif.exists():
        raise FileNotFoundError(f"Template CIF not found: {template_cif}")
    if len(template_entry_id) < 3:
        raise ValueError("template_entry_id must be at least 3 characters long.")

    subdir = mmcif_dir / template_entry_id[1:3]
    subdir.mkdir(parents=True, exist_ok=True)
    target = subdir / f"{template_entry_id}.cif.gz"
    with template_cif.open("rb") as src, gzip.open(target, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return target.resolve()


def write_template_a3m(
    chain_id: str,
    seq: str,
    template_entry_id: str,
    out_dir: Path,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{template_entry_id}_{chain_id}.a3m"
    # Keep a single hit header form to avoid parser ambiguity.
    text = f">query\n{seq}\n>{template_entry_id}_{chain_id}\n{seq}\n"
    path.write_text(text, encoding="utf-8")
    return path.resolve()


def write_input_json(
    path: Path,
    job_name: str,
    chain_order: Sequence[str],
    chain_to_seq: Dict[str, str],
    msa_root: Path,
    template_a3m_by_chain: Dict[str, Path],
) -> None:
    payload = [
        {
            "name": job_name,
            "sequences": [],
            "covalent_bonds": [],
        }
    ]

    for idx, chain_id in enumerate(chain_order):
        seq = chain_to_seq[chain_id]
        paired = (msa_root / str(idx) / "pairing.a3m").resolve()
        unpaired = (msa_root / str(idx) / "non_pairing.a3m").resolve()
        if not paired.exists():
            raise FileNotFoundError(f"Missing MSA file: {paired}")
        if not unpaired.exists():
            raise FileNotFoundError(f"Missing MSA file: {unpaired}")

        payload[0]["sequences"].append(
            {
                "proteinChain": {
                    "sequence": seq,
                    "count": 1,
                    "pairedMsaPath": str(paired),
                    "unpairedMsaPath": str(unpaired),
                    "templatesPath": str(template_a3m_by_chain[chain_id]),
                }
            }
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create minimal template-enabled Protenix input from a seed PDB."
    )
    parser.add_argument(
        "--input_pdb",
        type=Path,
        default=Path("data/avb3/template_example/seed_090_frame_000.pdb"),
        help="Input seed PDB used for sequence extraction.",
    )
    parser.add_argument(
        "--chain_order",
        type=str,
        default="A,B",
        help="Comma-separated chain order for JSON + MSA mapping.",
    )
    parser.add_argument(
        "--msa_root",
        type=Path,
        default=Path("data/avb3/template_example/msa"),
        help="MSA root with per-chain folders: 0/, 1/, ...",
    )
    parser.add_argument(
        "--output_root",
        type=Path,
        default=Path("data/runs/avb3/protenix_template/simple_pipeline"),
        help="Output root for generated inputs.",
    )
    parser.add_argument(
        "--job_name",
        type=str,
        default="seed_090_frame_000_template",
        help="Job name written to Protenix input JSON.",
    )
    parser.add_argument(
        "--template_entry_id",
        type=str,
        default="s090",
        help="Template entry ID used in a3m headers and mmcif DB registration.",
    )
    parser.add_argument(
        "--template_cif",
        type=Path,
        default=None,
        help="Template CIF output path. Defaults to <output_root>/templates/<entry>.cif",
    )
    parser.add_argument(
        "--template_converter",
        type=str,
        default="auto",
        choices=("auto", "mkdssp", "gemmi"),
        help="PDB->mmCIF converter. auto prefers mkdssp and falls back to gemmi.",
    )
    parser.add_argument(
        "--mkdssp_binary_path",
        type=str,
        default="",
        help="Optional explicit path to mkdssp binary.",
    )
    parser.add_argument(
        "--template_mmcif_dir",
        type=Path,
        default=None,
        help="Protenix mmcif DB root. Defaults to $PROTENIX_ROOT_DIR/mmcif.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    input_pdb = args.input_pdb.resolve()
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

    chain_order = parse_chain_order(args.chain_order)
    template_entry_id = args.template_entry_id.strip().lower()
    if not ENTRY_ID_RE.fullmatch(template_entry_id):
        raise ValueError("--template_entry_id must match [A-Za-z0-9_]{3,}.")

    output_root = args.output_root.resolve()
    inputs_dir = output_root / "inputs"
    template_dir = output_root / "templates"
    template_a3m_dir = inputs_dir / "template_search_results"
    input_json = inputs_dir / f"{args.job_name}.json"
    template_cif = (
        args.template_cif.resolve()
        if args.template_cif is not None
        else (template_dir / f"{template_entry_id}.cif").resolve()
    )

    chain_to_seq = parse_pdb_sequences(input_pdb)
    missing_chains = [chain for chain in chain_order if chain not in chain_to_seq]
    if missing_chains:
        available = ",".join(sorted(chain_to_seq))
        raise ValueError(
            "Missing expected chain(s): "
            + ",".join(missing_chains)
            + f". Available chains: {available}"
        )

    print("Parsed chains:")
    for chain_id in chain_order:
        print(f"  chain {chain_id}: length={len(chain_to_seq[chain_id])}")

    template_a3m_by_chain: Dict[str, Path] = {}
    for chain_id in chain_order:
        template_a3m_by_chain[chain_id] = write_template_a3m(
            chain_id=chain_id,
            seq=chain_to_seq[chain_id],
            template_entry_id=template_entry_id,
            out_dir=template_a3m_dir,
        )

    converter_used = convert_pdb_to_mmcif(
        pdb_path=input_pdb,
        cif_path=template_cif,
        converter=args.template_converter,
        mkdssp_binary_path=args.mkdssp_binary_path,
    )

    mmcif_dir = resolve_template_mmcif_dir(args.template_mmcif_dir)
    registered_mmcif = register_template_mmcif_entry(
        template_cif=template_cif,
        template_entry_id=template_entry_id,
        mmcif_dir=mmcif_dir,
    )

    write_input_json(
        path=input_json,
        job_name=args.job_name,
        chain_order=chain_order,
        chain_to_seq=chain_to_seq,
        msa_root=args.msa_root.resolve(),
        template_a3m_by_chain=template_a3m_by_chain,
    )

    print("\nGenerated files:")
    print(f"  input_json:       {input_json}")
    print(f"  template_cif:     {template_cif}")
    print(f"  mmcif_registered: {registered_mmcif}")
    print(f"  converter_used:   {converter_used}")
    for chain_id in chain_order:
        print(f"  template_a3m[{chain_id}]: {template_a3m_by_chain[chain_id]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
