#!/usr/bin/env python3
"""Generate BoltzGen-compatible YAML specs for heterodimer sampling."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List

import yaml


def read_seq(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"Empty sequence file: {path}")
    return text


def parse_thresholds(raw: str) -> List[float]:
    if not raw.strip():
        return []
    out = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        out.append(float(token))
    return out


def build_templates(template_cif: Path, chain_id: str, force: bool, threshold: float | None) -> list[dict]:
    # Legacy helper retained only to keep import-level compatibility if referenced elsewhere.
    # Current BoltzGen schema does not accept nested protein.templates blocks.
    return []


def write_yaml(path: Path, doc: dict) -> None:
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")


def build_protein_doc(seq_a: str, seq_b: str, msa_a: str, msa_b: str) -> dict:
    return {
        "entities": [
            {"protein": {"id": "A", "sequence": seq_a, "msa": msa_a}},
            {"protein": {"id": "B", "sequence": seq_b, "msa": msa_b}},
        ]
    }


def build_file_context_doc(template_cif: Path, msa_a: str, msa_b: str) -> dict:
    return {
        "entities": [
            {
                "file": {
                    "path": str(template_cif.resolve()),
                    "include": [
                        {"chain": {"id": "A", "msa": msa_a}},
                        {"chain": {"id": "B", "msa": msa_b}},
                    ],
                    "structure_groups": "all",
                }
            }
        ]
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain-a-seq-file", required=True, type=Path)
    parser.add_argument("--chain-b-seq-file", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--job-prefix", default="integrin_ab_boltz")
    parser.add_argument("--chain-a-msa", default="")
    parser.add_argument("--chain-b-msa", default="")
    parser.add_argument("--template-cif", default="")
    parser.add_argument(
        "--force-thresholds",
        default="0.60,0.75,0.90",
        help=(
            "Legacy template-force labels. Current BoltzGen schema does not expose "
            "protein.templates force thresholds, so template-guided generation emits "
            "a single file-context job when --template-cif is provided."
        ),
    )
    parser.add_argument(
        "--extended-only",
        action="store_true",
        help=(
            "Generate only template-force jobs (skip baseline, empty-MSA control, and soft template) "
            "to focus purely on extended-state search."
        ),
    )
    args = parser.parse_args()

    seq_a = read_seq(args.chain_a_seq_file)
    seq_b = read_seq(args.chain_b_seq_file)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    msa_a = args.chain_a_msa.strip() or "empty"
    msa_b = args.chain_b_msa.strip() or "empty"
    has_template = bool(args.template_cif.strip())
    template_cif = Path(args.template_cif) if has_template else None
    if template_cif is not None and not template_cif.exists():
        raise FileNotFoundError(f"Template CIF not found: {template_cif}")
    if args.extended_only and template_cif is None:
        raise ValueError("--extended-only requires --template-cif.")

    thresholds = parse_thresholds(args.force_thresholds)
    if args.extended_only and not thresholds:
        raise ValueError("--extended-only requires at least one --force-thresholds value.")

    manifest_path = outdir / "jobs_manifest.tsv"
    job_count = 0

    with manifest_path.open("w", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["job_name", "mode", "yaml_path"])

        if not args.extended_only:
            # Baseline with user-provided MSA paths (or empty).
            baseline_name = f"{args.job_prefix}_msa_baseline"
            baseline_doc = build_protein_doc(seq_a, seq_b, msa_a, msa_b)
            baseline_yaml = outdir / f"{baseline_name}.yaml"
            write_yaml(baseline_yaml, baseline_doc)
            writer.writerow([baseline_name, "msa_baseline", baseline_yaml])
            job_count += 1

            # No-MSA control.
            no_msa_name = f"{args.job_prefix}_empty_msa_control"
            no_msa_doc = build_protein_doc(seq_a, seq_b, "empty", "empty")
            no_msa_yaml = outdir / f"{no_msa_name}.yaml"
            write_yaml(no_msa_yaml, no_msa_doc)
            writer.writerow([no_msa_name, "empty_msa_control", no_msa_yaml])
            job_count += 1

        if template_cif is not None:
            if len(thresholds) > 1:
                print(
                    "[sweep] Note: template force thresholds are not represented in the "
                    "installed BoltzGen schema; emitting one template_context job."
                )
            context_name = f"{args.job_prefix}_template_context"
            context_doc = build_file_context_doc(template_cif, msa_a, msa_b)
            context_yaml = outdir / f"{context_name}.yaml"
            write_yaml(context_yaml, context_doc)
            writer.writerow([context_name, "template_context", context_yaml])
            job_count += 1

    print(f"[sweep] wrote {job_count} jobs to {outdir}")
    print(f"[sweep] manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
