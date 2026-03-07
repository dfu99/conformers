#!/usr/bin/env python3
"""Prepare staged A5B1 partner-attachment workflow inputs.

Stages generated:
1) Heterodimer + SpyTag (beta-chain tail anchor)
2) Heterodimer + Streptavidin (alpha-chain tail anchor)
3) Merge plan for combining stage 1/2 docked outputs by receptor overlap
"""

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


def select_best_prediction(predictions_dir: Path) -> Path:
    summaries = sorted(predictions_dir.glob("*_summary_confidence_sample_*.json"))
    if not summaries:
        raise FileNotFoundError(
            f"No summary confidence json files found in predictions dir: {predictions_dir}"
        )

    best = None
    for summary in summaries:
        data = json.loads(summary.read_text(encoding="utf-8"))
        score = float(data.get("ranking_score", float("-inf")))
        cif = summary.with_name(summary.name.replace("_summary_confidence", ""))
        cif = cif.with_suffix(".cif")
        if not cif.exists():
            raise FileNotFoundError(f"Expected prediction cif not found for {summary}: {cif}")
        if best is None or score > best[0]:
            best = (score, cif)

    assert best is not None
    return best[1]


def parse_tail_anchors_from_cif(cif_path: Path, chain_ids: List[str]) -> Dict[str, Dict[str, int]]:
    try:
        import gemmi  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "gemmi is required for chain/tail parsing from cif. Install with: pip install gemmi"
        ) from exc

    structure = gemmi.read_structure(str(cif_path))
    if len(structure) == 0:
        raise ValueError(f"No model found in {cif_path}")
    model = structure[0]

    out: Dict[str, Dict[str, int]] = {}
    for chain_id in chain_ids:
        chain = model.find_chain(chain_id)
        if chain is None:
            available = [c.name for c in model]
            raise ValueError(f"Chain '{chain_id}' not found in {cif_path}. Available: {available}")

        seen = set()
        seq_nums = []
        for residue in chain:
            has_ca = False
            for atom in residue:
                if atom.name.strip() == "CA":
                    has_ca = True
                    break
            if not has_ca:
                continue
            key = (residue.seqid.num, residue.seqid.icode)
            if key in seen:
                continue
            seen.add(key)
            seq_nums.append(int(residue.seqid.num))

        if not seq_nums:
            raise ValueError(f"No protein CA residues found in chain '{chain_id}' of {cif_path}")

        out[chain_id] = {
            "length": len(seq_nums),
            "tail_residue": max(seq_nums),
        }

    return out


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare staged A5B1 SpyTag/Streptavidin attachment workflow inputs."
    )
    parser.add_argument(
        "--sequence-file",
        type=Path,
        default=Path("data/a5b1/sequences/sequences_updated"),
        help="FASTA-like component sequence file.",
    )
    parser.add_argument(
        "--predictions-dir",
        type=Path,
        default=Path(
            "data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/"
            "integrin_alpha5_beta1/seed_101/predictions"
        ),
        help="Directory containing heterodimer prediction cifs + summary confidence json files.",
    )
    parser.add_argument(
        "--heterodimer-cif",
        type=Path,
        default=None,
        help="Optional explicit heterodimer cif. If omitted, best ranking_score in predictions-dir is used.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("data/runs/a5b1/staged_attachment"),
        help="Output workflow directory.",
    )
    parser.add_argument("--alpha-name", default="Integrin alpha5-Avi")
    parser.add_argument("--beta-name", default="Integrin beta1-spycatcher")
    parser.add_argument("--spytag-name", default="Spytag")
    parser.add_argument("--streptavidin-name", default="Streptavidin")
    parser.add_argument("--alpha-chain", default="A")
    parser.add_argument("--beta-chain", default="B")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.sequence_file.exists():
        raise FileNotFoundError(f"Sequence file not found: {args.sequence_file}")
    if not args.predictions_dir.exists() and args.heterodimer_cif is None:
        raise FileNotFoundError(f"Predictions dir not found: {args.predictions_dir}")

    records = parse_fasta_like(args.sequence_file)
    alpha_name, alpha_seq = pick_sequence(records, args.alpha_name)
    beta_name, beta_seq = pick_sequence(records, args.beta_name)
    spytag_name, spytag_seq = pick_sequence(records, args.spytag_name)
    streptavidin_name, streptavidin_seq = pick_sequence(records, args.streptavidin_name)

    heterodimer_cif = (
        args.heterodimer_cif.resolve()
        if args.heterodimer_cif is not None
        else select_best_prediction(args.predictions_dir.resolve())
    )
    if not heterodimer_cif.exists():
        raise FileNotFoundError(f"Heterodimer cif not found: {heterodimer_cif}")

    chain_ids = [args.alpha_chain, args.beta_chain]
    tail_anchors = parse_tail_anchors_from_cif(heterodimer_cif, chain_ids)

    outdir = args.outdir.resolve()
    inputs_dir = outdir / "inputs"
    seq_dir = inputs_dir / "sequences"
    constraints_dir = inputs_dir / "constraints"
    manifests_dir = inputs_dir / "manifests"

    for d in (seq_dir, constraints_dir, manifests_dir):
        d.mkdir(parents=True, exist_ok=True)

    (seq_dir / "integrin_alpha5_avi.seq").write_text(alpha_seq + "\n", encoding="utf-8")
    (seq_dir / "integrin_beta1_spycatcher.seq").write_text(beta_seq + "\n", encoding="utf-8")
    (seq_dir / "spytag.seq").write_text(spytag_seq + "\n", encoding="utf-8")
    (seq_dir / "streptavidin.seq").write_text(streptavidin_seq + "\n", encoding="utf-8")

    stage1_constraint = {
        "type": "tail_proximity",
        "receptor_chain": args.beta_chain,
        "receptor_tail_residue": tail_anchors[args.beta_chain]["tail_residue"],
        "ligand_name": "spytag",
        "ligand_residue": 1,
        "target_distance_angstrom": 6.0,
        "distance_tolerance_angstrom": 3.0,
        "note": "Bias SpyTag near beta1-SpyCatcher tail.",
    }
    stage2_constraint = {
        "type": "tail_proximity",
        "receptor_chain": args.alpha_chain,
        "receptor_tail_residue": tail_anchors[args.alpha_chain]["tail_residue"],
        "ligand_name": "streptavidin",
        "ligand_residue": 1,
        "target_distance_angstrom": 8.0,
        "distance_tolerance_angstrom": 4.0,
        "note": "Bias Streptavidin near alpha5-Avi tail.",
    }

    write_json(constraints_dir / "stage1_spytag_constraint.json", stage1_constraint)
    write_json(constraints_dir / "stage2_streptavidin_constraint.json", stage2_constraint)

    stage1_manifest = {
        "stage": "stage1_spytag",
        "objective": "Dock/refine SpyTag attachment to beta1-SpyCatcher tail",
        "receptor": {
            "heterodimer_cif": str(heterodimer_cif),
            "chains": [args.alpha_chain, args.beta_chain],
        },
        "ligand": {
            "name": spytag_name,
            "sequence_file": str((seq_dir / "spytag.seq").resolve()),
        },
        "constraint_file": str((constraints_dir / "stage1_spytag_constraint.json").resolve()),
        "output_dir": str((outdir / "outputs" / "stage1_spytag").resolve()),
    }

    stage2_manifest = {
        "stage": "stage2_streptavidin",
        "objective": "Dock/refine Streptavidin attachment to alpha5-Avi tail",
        "receptor": {
            "heterodimer_cif": str(heterodimer_cif),
            "chains": [args.alpha_chain, args.beta_chain],
        },
        "ligand": {
            "name": streptavidin_name,
            "sequence_file": str((seq_dir / "streptavidin.seq").resolve()),
        },
        "constraint_file": str((constraints_dir / "stage2_streptavidin_constraint.json").resolve()),
        "output_dir": str((outdir / "outputs" / "stage2_streptavidin").resolve()),
    }

    stage3_manifest = {
        "stage": "stage3_merge",
        "objective": "Merge stage1 and stage2 docked complexes by overlapping A5B1 receptor residues",
        "reference_heterodimer_cif": str(heterodimer_cif),
        "stage1_expected_complex": str(
            (outdir / "outputs" / "stage1_spytag" / "best_complex.cif").resolve()
        ),
        "stage2_expected_complex": str(
            (outdir / "outputs" / "stage2_streptavidin" / "best_complex.cif").resolve()
        ),
        "receptor_chains_for_alignment": [args.alpha_chain, args.beta_chain],
        "merged_output": str((outdir / "outputs" / "stage3_merged" / "a5b1_tagged_merged.cif").resolve()),
        "note": "Implement structural merge after stage1/2 docking outputs are available.",
    }

    write_json(manifests_dir / "stage1_spytag_dock.json", stage1_manifest)
    write_json(manifests_dir / "stage2_streptavidin_dock.json", stage2_manifest)
    write_json(manifests_dir / "stage3_merge_plan.json", stage3_manifest)

    components = {
        "sequence_file": str(args.sequence_file.resolve()),
        "selected_components": {
            "alpha": {"name": alpha_name, "length": len(alpha_seq)},
            "beta": {"name": beta_name, "length": len(beta_seq)},
            "spytag": {"name": spytag_name, "length": len(spytag_seq)},
            "streptavidin": {"name": streptavidin_name, "length": len(streptavidin_seq)},
        },
        "heterodimer_cif": str(heterodimer_cif),
        "tail_anchors": tail_anchors,
    }
    write_json(inputs_dir / "components.json", components)

    runbook = f"""# A5B1 Staged Attachment Runbook

Generated by: setup_staged_attachment_workflow.py

## Selected receptor
- Heterodimer CIF: {heterodimer_cif}
- Alpha chain: {args.alpha_chain} (tail residue {tail_anchors[args.alpha_chain]['tail_residue']})
- Beta chain: {args.beta_chain} (tail residue {tail_anchors[args.beta_chain]['tail_residue']})

## Stage files
- Stage1 manifest: {manifests_dir / 'stage1_spytag_dock.json'}
- Stage2 manifest: {manifests_dir / 'stage2_streptavidin_dock.json'}
- Stage3 merge plan: {manifests_dir / 'stage3_merge_plan.json'}

## Suggested next actions
1. Run constrained docking/refinement for Stage1 (SpyTag) using stage1 manifest + constraint json.
2. Run constrained docking/refinement for Stage2 (Streptavidin) using stage2 manifest + constraint json.
3. Merge Stage1/Stage2 complexes using overlap of receptor chains {args.alpha_chain},{args.beta_chain}.

Note: This script prepares structured inputs and constraints only; it does not run docking.
"""
    (outdir / "RUNBOOK.md").write_text(runbook, encoding="utf-8")

    print("Prepared staged attachment workflow inputs:")
    print(f"  outdir: {outdir}")
    print(f"  heterodimer_cif: {heterodimer_cif}")
    print(f"  stage1: {manifests_dir / 'stage1_spytag_dock.json'}")
    print(f"  stage2: {manifests_dir / 'stage2_streptavidin_dock.json'}")
    print(f"  stage3: {manifests_dir / 'stage3_merge_plan.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
