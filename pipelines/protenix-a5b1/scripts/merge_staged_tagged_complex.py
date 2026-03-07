#!/usr/bin/env python3
"""Merge staged A5B1 tag-attachment predictions into one final complex.

Inputs:
- base heterodimer CIF (accepted integrin A/B structure)
- stage-1 predictions dir (A/B + SpyTag)
- stage-2 predictions dir (A/B + Streptavidin)

Output:
- merged complex with receptor chains from base plus transformed stage ligands.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


def select_best_prediction(predictions_dir: Path) -> Path:
    summaries = sorted(predictions_dir.glob("*_summary_confidence_sample_*.json"))
    if not summaries:
        raise FileNotFoundError(
            f"No summary confidence json files found in predictions dir: {predictions_dir}"
        )

    best: Tuple[float, Path] | None = None
    for summary in summaries:
        data = json.loads(summary.read_text(encoding="utf-8"))
        score = float(data.get("ranking_score", float("-inf")))
        cif = summary.with_name(summary.name.replace("_summary_confidence", "")).with_suffix(".cif")
        if not cif.exists():
            raise FileNotFoundError(f"Expected prediction cif not found for {summary}: {cif}")
        if best is None or score > best[0]:
            best = (score, cif)

    if best is None:
        raise RuntimeError(f"Failed to select best prediction from {predictions_dir}")
    return best[1]


def parse_chain_list(value: str) -> List[str]:
    chains = [x.strip() for x in value.split(",") if x.strip()]
    if not chains:
        raise ValueError("At least one chain ID is required.")
    return chains


def load_structure(path: Path):
    import gemmi  # type: ignore

    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        raise ValueError(f"No model found in {path}")
    return st, st[0]


def residue_key(chain_id: str, residue) -> Tuple[str, int, str]:
    icode = str(getattr(residue.seqid, "icode", "") or "").strip()
    if icode in {".", "?"}:
        icode = ""
    return chain_id, int(residue.seqid.num), icode


def collect_ca_positions(model, chain_ids: Sequence[str]) -> Dict[Tuple[str, int, str], np.ndarray]:
    out: Dict[Tuple[str, int, str], np.ndarray] = {}
    for chain_id in chain_ids:
        chain = model.find_chain(chain_id)
        if chain is None:
            available = [c.name for c in model]
            raise ValueError(
                f"Chain '{chain_id}' not found in model. Available chains: {available}"
            )
        for residue in chain:
            for atom in residue:
                if atom.name.strip() == "CA":
                    out[residue_key(chain_id, residue)] = np.array(
                        [atom.pos.x, atom.pos.y, atom.pos.z], dtype=float
                    )
                    break
    return out


def kabsch_transform(moving: np.ndarray, reference: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    if moving.shape != reference.shape:
        raise ValueError("moving/reference shape mismatch")
    if moving.ndim != 2 or moving.shape[1] != 3:
        raise ValueError("Expected Nx3 coordinates")
    if moving.shape[0] < 3:
        raise ValueError("Need at least 3 points for rigid alignment")

    cm = moving.mean(axis=0)
    cr = reference.mean(axis=0)

    xm = moving - cm
    xr = reference - cr

    h = xm.T @ xr
    u, _s, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1.0
        r = vt.T @ u.T
    t = cr - r @ cm

    moved = (r @ moving.T).T + t
    rmsd = float(np.sqrt(np.mean(np.sum((moved - reference) ** 2, axis=1))))
    return r, t, rmsd


def compute_alignment(base_model, stage_model, receptor_chains: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, float, int]:
    base_map = collect_ca_positions(base_model, receptor_chains)
    stage_map = collect_ca_positions(stage_model, receptor_chains)

    common_keys = sorted(set(base_map) & set(stage_map))
    if len(common_keys) < 3:
        raise ValueError(
            f"Insufficient receptor overlap for alignment: found {len(common_keys)} common CA residues"
        )

    moving = np.stack([stage_map[k] for k in common_keys], axis=0)
    reference = np.stack([base_map[k] for k in common_keys], axis=0)
    r, t, rmsd = kabsch_transform(moving, reference)
    return r, t, rmsd, len(common_keys)


def transform_chain(chain, rotation: np.ndarray, translation: np.ndarray) -> None:
    import gemmi  # type: ignore

    for residue in chain:
        for atom in residue:
            xyz = np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
            new_xyz = rotation @ xyz + translation
            atom.pos = gemmi.Position(float(new_xyz[0]), float(new_xyz[1]), float(new_xyz[2]))


def ensure_chain_free(model, chain_id: str) -> None:
    existing = [c.name for c in model]
    if chain_id in existing:
        raise ValueError(
            f"Target output chain ID '{chain_id}' already exists in base model chains: {existing}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-cif", type=Path, required=True)
    parser.add_argument("--stage1-predictions-dir", type=Path, required=True)
    parser.add_argument("--stage2-predictions-dir", type=Path, required=True)
    parser.add_argument(
        "--receptor-chains",
        default="A,B",
        help="Comma-separated receptor chains used for alignment.",
    )
    parser.add_argument("--stage1-ligand-chain", default="C")
    parser.add_argument("--stage2-ligand-chain", default="C")
    parser.add_argument("--out-stage1-chain", default="C")
    parser.add_argument("--out-stage2-chain", default="D")
    parser.add_argument(
        "--out-cif",
        type=Path,
        default=Path("data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.cif"),
    )
    parser.add_argument(
        "--out-pdb",
        type=Path,
        default=Path("data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.pdb"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receptor_chains = parse_chain_list(args.receptor_chains)

    stage1_cif = select_best_prediction(args.stage1_predictions_dir.resolve())
    stage2_cif = select_best_prediction(args.stage2_predictions_dir.resolve())

    base_st, base_model = load_structure(args.base_cif.resolve())
    stage1_st, stage1_model = load_structure(stage1_cif)
    stage2_st, stage2_model = load_structure(stage2_cif)

    r1, t1, rmsd1, n1 = compute_alignment(base_model, stage1_model, receptor_chains)
    r2, t2, rmsd2, n2 = compute_alignment(base_model, stage2_model, receptor_chains)

    stage1_lig = stage1_model.find_chain(args.stage1_ligand_chain)
    if stage1_lig is None:
        available = [c.name for c in stage1_model]
        raise ValueError(
            f"Stage1 ligand chain '{args.stage1_ligand_chain}' not found. Available: {available}"
        )

    stage2_lig = stage2_model.find_chain(args.stage2_ligand_chain)
    if stage2_lig is None:
        available = [c.name for c in stage2_model]
        raise ValueError(
            f"Stage2 ligand chain '{args.stage2_ligand_chain}' not found. Available: {available}"
        )

    transform_chain(stage1_lig, r1, t1)
    transform_chain(stage2_lig, r2, t2)

    stage1_lig.name = args.out_stage1_chain
    stage2_lig.name = args.out_stage2_chain

    ensure_chain_free(base_model, stage1_lig.name)
    ensure_chain_free(base_model, stage2_lig.name)

    # gemmi model.add_chain() copies the chain object.
    base_model.add_chain(stage1_lig)
    base_model.add_chain(stage2_lig)

    args.out_cif.parent.mkdir(parents=True, exist_ok=True)
    args.out_pdb.parent.mkdir(parents=True, exist_ok=True)

    base_st.make_mmcif_document().write_file(str(args.out_cif.resolve()))
    base_st.write_pdb(str(args.out_pdb.resolve()))

    summary = {
        "base_cif": str(args.base_cif.resolve()),
        "stage1_best_cif": str(stage1_cif),
        "stage2_best_cif": str(stage2_cif),
        "receptor_chains": receptor_chains,
        "stage1_alignment": {
            "common_ca_count": n1,
            "rmsd": rmsd1,
            "input_ligand_chain": args.stage1_ligand_chain,
            "output_ligand_chain": args.out_stage1_chain,
        },
        "stage2_alignment": {
            "common_ca_count": n2,
            "rmsd": rmsd2,
            "input_ligand_chain": args.stage2_ligand_chain,
            "output_ligand_chain": args.out_stage2_chain,
        },
        "out_cif": str(args.out_cif.resolve()),
        "out_pdb": str(args.out_pdb.resolve()),
    }

    summary_path = args.out_cif.with_suffix(".merge_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Merged final tagged complex:")
    print(f"  out_cif: {args.out_cif.resolve()}")
    print(f"  out_pdb: {args.out_pdb.resolve()}")
    print(f"  merge_summary: {summary_path.resolve()}")
    print(f"  stage1_rmsd={rmsd1:.3f} over {n1} CA")
    print(f"  stage2_rmsd={rmsd2:.3f} over {n2} CA")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
