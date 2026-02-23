#!/usr/bin/env python3
"""Protenix workflow for integrin conformation sampling from a self-template PDB.

This script creates two inference inputs from a starting PDB:
1) MSA-only input
2) MSA + self-template input (inline template mapping with explicit chain IDs)

It can also launch multi-seed Protenix inference for both inputs so you can
compare conformational diversity with and without template guidance.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Sequence


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
    # Common substitutions; Protenix expects 20 aa + X.
    "MSE": "M",
    "SEC": "C",
    "PYL": "K",
    "ASX": "X",
    "GLX": "X",
    "UNK": "X",
}

TEMPLATE_ENTRY_ID_PATTERN = re.compile(r"^[A-Za-z0-9_]+$")


def parse_chain_order(value: str) -> List[str]:
    chains = [x.strip() for x in value.split(",") if x.strip()]
    if not chains:
        raise ValueError("--chain_order must include at least one chain.")
    return chains


def parse_pdb_sequences(pdb_path: Path) -> Dict[str, str]:
    """Extract protein sequence per chain from ATOM CA records."""
    chains: Dict[str, List[str]] = {}
    seen: Dict[str, set] = {}

    with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            altloc = line[16]
            if altloc not in (" ", "A", "1"):
                continue

            resname = line[17:20].strip().upper()
            chain_id = line[21].strip() or "_"
            resseq = line[22:26].strip()
            icode = line[26].strip()
            residue_key = (resseq, icode)

            if chain_id not in chains:
                chains[chain_id] = []
                seen[chain_id] = set()
            if residue_key in seen[chain_id]:
                continue
            seen[chain_id].add(residue_key)

            aa1 = AA3_TO_AA1.get(resname, "X")
            chains[chain_id].append(aa1)

    if not chains:
        raise ValueError(f"No protein CA records found in {pdb_path}")

    return {chain: "".join(seq) for chain, seq in chains.items()}


def resolve_mkdssp_binary(explicit_path: str) -> str | None:
    if explicit_path:
        path = Path(explicit_path)
        if not path.exists():
            raise FileNotFoundError(
                f"--mkdssp_binary_path was provided but does not exist: {explicit_path}"
            )
        return str(path.resolve())
    discovered = shutil.which("mkdssp")
    return discovered


def convert_pdb_to_mmcif(
    pdb_path: Path,
    cif_path: Path,
    converter: str,
    mkdssp_binary_path: str,
) -> str:
    """Convert PDB to mmCIF. Prefer mkdssp for Protenix template compatibility."""
    cif_path.parent.mkdir(parents=True, exist_ok=True)

    if converter not in {"auto", "mkdssp", "gemmi"}:
        raise ValueError(f"Unsupported converter: {converter}")

    if converter in {"auto", "mkdssp"}:
        mkdssp = resolve_mkdssp_binary(mkdssp_binary_path)
        if mkdssp is not None:
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
                "--template_converter mkdssp was requested but mkdssp was not found in PATH.\n"
                "Provide --mkdssp_binary_path /path/to/mkdssp."
            )

    try:
        import gemmi  # type: ignore
    except Exception as exc:  # pragma: no cover - runtime dependency
        raise RuntimeError(
            "PDB->mmCIF conversion failed: mkdssp unavailable/failed and gemmi not installed.\n"
            "Install one of:\n"
            "  1) mkdssp (preferred)\n"
            "  2) gemmi (`pip install gemmi`)"
        ) from exc

    structure = gemmi.read_structure(str(pdb_path))
    structure.make_mmcif_document().write_file(str(cif_path))
    return "gemmi"


def normalize_template_entry_id(value: str) -> str:
    """Normalize and validate entry ID used in template .a3m headers."""
    entry_id = value.strip()
    if not entry_id:
        raise ValueError("--template_entry_id cannot be empty.")
    if not TEMPLATE_ENTRY_ID_PATTERN.fullmatch(entry_id):
        raise ValueError(
            "--template_entry_id must match [A-Za-z0-9_]+ so Protenix can parse it."
        )
    if len(entry_id) < 3:
        raise ValueError("--template_entry_id must be at least 3 characters long.")
    return entry_id.lower()


def write_self_template_a3m_files(
    chain_ids: Sequence[str],
    chain_to_seq: Dict[str, str],
    template_entry_id: str,
    output_dir: Path,
) -> Dict[str, Path]:
    """Write Protenix-compatible hmmsearch .a3m files per chain."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Path] = {}

    for chain_id in chain_ids:
        if chain_id not in chain_to_seq:
            available = ", ".join(sorted(chain_to_seq))
            raise ValueError(
                f"Chain '{chain_id}' not found in parsed PDB chains. Available: {available}"
            )
        seq = chain_to_seq[chain_id]
        # Some Protenix builds parse ENTRY_CHAIN while others parse ENTRY_CHAIN/start-end.
        # Emit both to maximize compatibility with the installed template parser.
        hit_header_plain = f">{template_entry_id}_{chain_id}"
        hit_header_with_range = f">{template_entry_id}_{chain_id}/1-{len(seq)}"
        path = output_dir / f"{template_entry_id}_{chain_id}.a3m"
        path.write_text(
            (
                f">query\n{seq}\n"
                f"{hit_header_plain}\n{seq}\n"
                f"{hit_header_with_range}\n{seq}\n"
            ),
            encoding="utf-8",
        )
        out[chain_id] = path.resolve()

    return out


def resolve_template_mmcif_dir(explicit_dir: Path | None) -> Path:
    """Resolve template mmcif DB root used by Protenix template retrieval."""
    if explicit_dir is not None:
        return explicit_dir.resolve()

    protenix_root_dir = os.environ.get("PROTENIX_ROOT_DIR", "").strip()
    if not protenix_root_dir:
        raise RuntimeError(
            "Cannot resolve template mmcif DB directory.\n"
            "Provide --template_mmcif_dir or set PROTENIX_ROOT_DIR."
        )
    return (Path(protenix_root_dir) / "mmcif").resolve()


def register_template_mmcif_entry(
    template_cif: Path,
    template_entry_id: str,
    mmcif_dir: Path,
) -> Path:
    """Register a custom template mmCIF under Protenix mmcif/<id[1:3]>/<id>.cif.gz."""
    if not template_cif.exists():
        raise FileNotFoundError(f"Template CIF not found: {template_cif}")
    if len(template_entry_id) < 3:
        raise ValueError("template_entry_id must be at least 3 characters long.")

    target_subdir = mmcif_dir / template_entry_id[1:3]
    target_subdir.mkdir(parents=True, exist_ok=True)
    target = target_subdir / f"{template_entry_id}.cif.gz"

    with template_cif.open("rb") as src, gzip.open(target, "wb") as dst:
        shutil.copyfileobj(src, dst)

    return target.resolve()


def build_input_payload(
    job_name: str,
    chain_ids: Sequence[str],
    chain_to_seq: Dict[str, str],
    msa_root: Path | None,
    template_cif: Path | None,
    template_paths_by_chain: Mapping[str, Path] | None = None,
) -> List[dict]:
    sequences: List[dict] = []

    for idx, chain_id in enumerate(chain_ids):
        if chain_id not in chain_to_seq:
            available = ", ".join(sorted(chain_to_seq))
            raise ValueError(
                f"Chain '{chain_id}' not found in parsed PDB chains. Available: {available}"
            )

        chain_obj = {
            "sequence": chain_to_seq[chain_id],
            "count": 1,
        }
        if msa_root is not None:
            chain_obj["pairedMsaPath"] = str((msa_root / str(idx) / "pairing.a3m").resolve())
            chain_obj["unpairedMsaPath"] = str(
                (msa_root / str(idx) / "non_pairing.a3m").resolve()
            )
        if template_paths_by_chain is not None:
            template_path = template_paths_by_chain.get(chain_id)
            if template_path is None:
                raise ValueError(f"Missing templatesPath for chain '{chain_id}'.")
            chain_obj["templatesPath"] = str(template_path.resolve())
        elif template_cif is not None:
            chain_obj["templates"] = [
                {
                    "mmcif_path": str(template_cif.resolve()),
                    "chain_id": chain_id,
                }
            ]

        sequences.append({"proteinChain": chain_obj})

    return [
        {
            "name": job_name,
            "sequences": sequences,
            "covalent_bonds": [],
        }
    ]


def write_json(path: Path, payload: List[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_protenix_pred(
    input_json: Path,
    out_dir: Path,
    model_name: str,
    seeds: str,
    use_msa: bool,
    use_template: bool,
    dtype: str,
    samples: int,
    triatt_kernel: str,
    trimul_kernel: str,
    kalign_binary_path: str | None = None,
) -> None:
    cmd = [
        "protenix",
        "pred",
        "-i",
        str(input_json),
        "-o",
        str(out_dir),
        "-n",
        model_name,
        "-s",
        seeds,
        "-e",
        str(samples),
        "--dtype",
        dtype,
        "--use_msa",
        str(use_msa).lower(),
        "--use_template",
        str(use_template).lower(),
        "--triatt_kernel",
        triatt_kernel,
        "--trimul_kernel",
        trimul_kernel,
    ]
    if kalign_binary_path:
        cmd.extend(["--kalign_binary_path", kalign_binary_path])
    print("\n[run]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create and optionally run Protenix self-template conformation sampling "
            "for integrin-like multichain proteins."
        )
    )
    parser.add_argument(
        "--input_pdb",
        type=Path,
        default=Path("data/template_example/seed_090_frame_000.pdb"),
        help="Starting PDB used for sequence extraction and self-template construction.",
    )
    parser.add_argument(
        "--chain_order",
        type=str,
        default="A,B",
        help="Comma-separated chain order to include in inference JSON and MSA index mapping.",
    )
    parser.add_argument(
        "--job_name",
        type=str,
        default="seed_090_frame_000",
        help="Base job name used in generated JSON files.",
    )
    parser.add_argument(
        "--workflow_dir",
        type=Path,
        default=Path("data/template_example/workflow_outputs"),
        help="Directory where generated inputs and outputs are stored.",
    )
    parser.add_argument(
        "--msa_root",
        type=Path,
        default=Path(
            "/content/drive/MyDrive/colab_cache/afmfold-data/AVB3/seed_090_frame_000/msa"
        ),
        help=(
            "Root folder containing per-chain MSA subfolders "
            "(e.g. 0/pairing.a3m, 0/non_pairing.a3m, 1/...)."
        ),
    )
    parser.add_argument(
        "--skip_msa_paths",
        action="store_true",
        help="Do not write pairedMsaPath/unpairedMsaPath fields into generated JSON.",
    )
    parser.add_argument(
        "--template_cif",
        type=Path,
        default=Path("data/template_example/seed_090_frame_000.cif"),
        help="Path for generated mmCIF self-template.",
    )
    parser.add_argument(
        "--template_json_mode",
        type=str,
        default="templatesPath",
        choices=("templatesPath", "legacy_templates"),
        help=(
            "Template JSON schema mode. `templatesPath` uses Protenix v1 format "
            "(.a3m/.hhr path); `legacy_templates` emits inline mmcif_path blocks."
        ),
    )
    parser.add_argument(
        "--template_entry_id",
        type=str,
        default="s090",
        help=(
            "Template entry ID used in generated .a3m hit headers (ENTRY_CHAIN/start-end). "
            "Example: s090 -> mmcif/09/s090.cif.gz"
        ),
    )
    parser.add_argument(
        "--template_a3m_dir",
        type=Path,
        default=None,
        help=(
            "Directory to write generated per-chain template .a3m files. "
            "Default: <workflow_dir>/inputs/template_search_results"
        ),
    )
    parser.add_argument(
        "--register_template_mmcif",
        action="store_true",
        help=(
            "Register generated template mmCIF as "
            "mmcif/<id[1:3]>/<id>.cif.gz for templatesPath-based retrieval."
        ),
    )
    parser.add_argument(
        "--template_mmcif_dir",
        type=Path,
        default=None,
        help=(
            "mmcif DB root to install template entry into. "
            "Defaults to $PROTENIX_ROOT_DIR/mmcif."
        ),
    )
    parser.add_argument(
        "--no_convert_template",
        action="store_true",
        help="Skip PDB->mmCIF conversion and use existing --template_cif file.",
    )
    parser.add_argument(
        "--template_converter",
        type=str,
        default="auto",
        choices=("auto", "mkdssp", "gemmi"),
        help=(
            "PDB->mmCIF converter. `auto` prefers mkdssp (recommended by Protenix maintainers) "
            "and falls back to gemmi."
        ),
    )
    parser.add_argument(
        "--mkdssp_binary_path",
        type=str,
        default="",
        help=(
            "Optional explicit path to mkdssp binary used when --template_converter "
            "is auto/mkdssp."
        ),
    )
    parser.add_argument(
        "--model_name",
        type=str,
        default="protenix_base_default_v1.0.0",
        help="Protenix model name.",
    )
    parser.add_argument(
        "--seeds",
        type=str,
        default="101,202,303,404,505",
        help="Comma-separated seeds for conformational sampling.",
    )
    parser.add_argument(
        "--samples_per_seed",
        type=int,
        default=5,
        help="Number of diffusion samples per seed (-e / --sample).",
    )
    parser.add_argument(
        "--dtype",
        type=str,
        default="bf16",
        choices=("bf16", "fp32"),
        help="Inference dtype.",
    )
    parser.add_argument(
        "--triatt_kernel",
        type=str,
        default="torch",
        help="Triangle attention kernel for inference.",
    )
    parser.add_argument(
        "--trimul_kernel",
        type=str,
        default="torch",
        help="Triangle multiplicative kernel for inference.",
    )
    parser.add_argument(
        "--kalign_binary_path",
        type=str,
        default="",
        help=(
            "Optional explicit path to kalign binary. If omitted, script uses "
            "the first `kalign` found in PATH."
        ),
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run `protenix pred` for generated inputs.",
    )
    parser.add_argument(
        "--run_msa_only",
        action="store_true",
        help="Run only the MSA-only job when --run is enabled.",
    )
    parser.add_argument(
        "--run_template_only",
        action="store_true",
        help="Run only the self-template job when --run is enabled.",
    )
    return parser.parse_args()


def resolve_kalign_binary(use_template: bool, explicit_path: str) -> str | None:
    if not use_template:
        return None

    if explicit_path:
        path = Path(explicit_path)
        if not path.exists():
            raise FileNotFoundError(
                f"--kalign_binary_path was provided but does not exist: {explicit_path}"
            )
        return str(path.resolve())

    discovered = shutil.which("kalign") or shutil.which("kalign3")
    if discovered:
        return discovered

    # In Colab, try to self-heal by installing kalign automatically.
    in_colab = "COLAB_GPU" in os.environ or Path("/content").exists()
    if in_colab:
        print("kalign not found; attempting installation via apt-get ...")
        subprocess.run(["apt-get", "update", "-y"], check=False)
        subprocess.run(["apt-get", "install", "-y", "kalign"], check=False)
        discovered = shutil.which("kalign") or shutil.which("kalign3")
        if discovered:
            return discovered

    raise RuntimeError(
        "Template-enabled inference requires kalign, but no kalign binary was found.\n"
        "Install it first, then rerun. For Colab:\n"
        "  apt-get update -y && apt-get install -y kalign\n"
        "Or provide --kalign_binary_path /path/to/kalign."
    )


def ensure_template_mmcif_dir(use_template: bool) -> None:
    if not use_template:
        return

    protenix_root_dir = os.environ.get("PROTENIX_ROOT_DIR", "").strip()
    if not protenix_root_dir:
        raise RuntimeError(
            "Template-enabled inference requires environment variable PROTENIX_ROOT_DIR.\n"
            "Set it to your Protenix runtime root that contains the mmcif database.\n"
            "Expected location: $PROTENIX_ROOT_DIR/mmcif"
        )

    mmcif_dir = Path(protenix_root_dir) / "mmcif"
    if not mmcif_dir.exists():
        raise RuntimeError(
            "Template-enabled inference requires mmcif directory at:\n"
            f"  {mmcif_dir}\n"
            "Create/symlink this folder or set PROTENIX_ROOT_DIR accordingly."
        )


def main() -> int:
    args = parse_args()

    chain_order = parse_chain_order(args.chain_order)
    input_pdb = args.input_pdb.resolve()
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

    workflow_dir = args.workflow_dir.resolve()
    inputs_dir = workflow_dir / "inputs"
    outputs_dir = workflow_dir / "outputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    chain_to_seq = parse_pdb_sequences(input_pdb)
    print("Parsed chains from PDB:")
    for chain_id in chain_order:
        seq = chain_to_seq.get(chain_id)
        if seq is None:
            available = ", ".join(sorted(chain_to_seq))
            raise ValueError(f"Missing chain '{chain_id}'. Available chains: {available}")
        print(f"  chain {chain_id}: length={len(seq)}")

    template_cif = args.template_cif.resolve()
    needs_template_cif = (
        args.template_json_mode == "legacy_templates" or args.register_template_mmcif
    )
    template_converter_used: str | None = None
    if needs_template_cif:
        if not args.no_convert_template:
            print(f"\nConverting template PDB -> mmCIF: {input_pdb} -> {template_cif}")
            template_converter_used = convert_pdb_to_mmcif(
                pdb_path=input_pdb,
                cif_path=template_cif,
                converter=args.template_converter,
                mkdssp_binary_path=args.mkdssp_binary_path,
            )
        elif not template_cif.exists():
            raise FileNotFoundError(
                f"--no_convert_template was set, but template CIF does not exist: {template_cif}"
            )
    else:
        print("\nSkipping template mmCIF conversion (not required for current options).")

    msa_root = None if args.skip_msa_paths else args.msa_root.resolve()
    if msa_root is not None:
        for i in range(len(chain_order)):
            for filename in ("pairing.a3m", "non_pairing.a3m"):
                candidate = msa_root / str(i) / filename
                if not candidate.exists():
                    raise FileNotFoundError(
                        f"Missing expected MSA file for chain index {i}: {candidate}"
                    )

    msa_only_json = inputs_dir / f"{args.job_name}_msa_only.json"
    self_template_json = inputs_dir / f"{args.job_name}_msa_self_template.json"

    msa_only_payload = build_input_payload(
        job_name=f"{args.job_name}_msa_only",
        chain_ids=chain_order,
        chain_to_seq=chain_to_seq,
        msa_root=msa_root,
        template_cif=None,
    )
    write_json(msa_only_json, msa_only_payload)

    template_paths_by_chain: Dict[str, Path] | None = None
    registered_template_path: Path | None = None
    if args.template_json_mode == "templatesPath":
        template_entry_id = normalize_template_entry_id(args.template_entry_id)
        template_a3m_dir = (
            args.template_a3m_dir.resolve()
            if args.template_a3m_dir is not None
            else (inputs_dir / "template_search_results").resolve()
        )
        template_paths_by_chain = write_self_template_a3m_files(
            chain_ids=chain_order,
            chain_to_seq=chain_to_seq,
            template_entry_id=template_entry_id,
            output_dir=template_a3m_dir,
        )
        if args.register_template_mmcif:
            mmcif_db_dir = resolve_template_mmcif_dir(args.template_mmcif_dir)
            registered_template_path = register_template_mmcif_entry(
                template_cif=template_cif,
                template_entry_id=template_entry_id,
                mmcif_dir=mmcif_db_dir,
            )

    self_template_payload = build_input_payload(
        job_name=f"{args.job_name}_msa_self_template",
        chain_ids=chain_order,
        chain_to_seq=chain_to_seq,
        msa_root=msa_root,
        template_cif=template_cif if args.template_json_mode == "legacy_templates" else None,
        template_paths_by_chain=template_paths_by_chain,
    )
    write_json(self_template_json, self_template_payload)

    print("\nGenerated inputs:")
    print(f"  MSA-only:      {msa_only_json}")
    print(f"  Self-template: {self_template_json}")
    print(f"  Template CIF:  {template_cif if needs_template_cif else '(not required)'}")
    if template_converter_used is not None:
        print(f"  Template converter: {template_converter_used}")
    print(f"  Template mode: {args.template_json_mode}")
    if template_paths_by_chain:
        print("  Template .a3m:")
        for chain_id in chain_order:
            print(f"    chain {chain_id}: {template_paths_by_chain[chain_id]}")
    if registered_template_path is not None:
        print(f"  Registered mmcif entry: {registered_template_path}")
    elif args.template_json_mode == "templatesPath":
        print(
            "  NOTE: template .a3m files were generated, but mmcif entry was not "
            "registered. Ensure the entry exists in runtime mmcif DB."
        )

    if not args.run:
        print("\nRun commands:")
        print(
            "  protenix pred"
            f" -i {msa_only_json}"
            f" -o {outputs_dir / 'msa_only'}"
            f" -n {args.model_name} -s {args.seeds} -e {args.samples_per_seed}"
            f" --dtype {args.dtype} --use_msa true --use_template false"
        )
        print(
            "  protenix pred"
            f" -i {self_template_json}"
            f" -o {outputs_dir / 'msa_self_template'}"
            f" -n {args.model_name} -s {args.seeds} -e {args.samples_per_seed}"
            f" --dtype {args.dtype} --use_msa true --use_template true"
        )
        return 0

    if args.run_msa_only and args.run_template_only:
        raise ValueError(
            "Choose at most one of --run_msa_only and --run_template_only, or neither for both."
        )

    run_msa_only = args.run_msa_only or not args.run_template_only
    run_template = args.run_template_only or not args.run_msa_only
    kalign_binary_path = resolve_kalign_binary(
        use_template=run_template, explicit_path=args.kalign_binary_path
    )
    ensure_template_mmcif_dir(use_template=run_template)

    if run_msa_only:
        run_protenix_pred(
            input_json=msa_only_json,
            out_dir=outputs_dir / "msa_only",
            model_name=args.model_name,
            seeds=args.seeds,
            use_msa=True,
            use_template=False,
            dtype=args.dtype,
            samples=args.samples_per_seed,
            triatt_kernel=args.triatt_kernel,
            trimul_kernel=args.trimul_kernel,
        )
    if run_template:
        run_protenix_pred(
            input_json=self_template_json,
            out_dir=outputs_dir / "msa_self_template",
            model_name=args.model_name,
            seeds=args.seeds,
            use_msa=True,
            use_template=True,
            dtype=args.dtype,
            samples=args.samples_per_seed,
            triatt_kernel=args.triatt_kernel,
            trimul_kernel=args.trimul_kernel,
            kalign_binary_path=kalign_binary_path,
        )

    print("\nDone.")
    print(f"Outputs root: {outputs_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
