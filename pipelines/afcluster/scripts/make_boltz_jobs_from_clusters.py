#!/usr/bin/env python3
"""Generate BoltzGen YAML jobs from AF-Cluster outputs for chain A/B."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List

import yaml


def read_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"Empty sequence file: {path}")
    return text


def is_outlier_cluster(path: Path) -> bool:
    return path.name == "cluster_outliers.a3m"


def load_ranked_clusters(cluster_dir: Path, top_n: int) -> List[Path]:
    summary = cluster_dir / "clusters.tsv"
    paths: List[Path] = []
    if summary.exists():
        with summary.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                path = cluster_dir / row["filename"]
                if path.exists():
                    paths.append(path)
            real_clusters = [path for path in paths if not is_outlier_cluster(path)]
            outlier_clusters = [path for path in paths if is_outlier_cluster(path)]
            paths = real_clusters + outlier_clusters
            if top_n > 0:
                paths = paths[:top_n]
    else:
        paths = sorted(
            cluster_dir.glob("*.a3m"),
            key=lambda p: (is_outlier_cluster(p), -p.stat().st_size, p.name),
        )
        if top_n > 0:
            paths = paths[:top_n]
    if not paths:
        raise RuntimeError(f"No cluster .a3m files found in {cluster_dir}")
    return paths


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
    parser = argparse.ArgumentParser(description="Build BoltzGen YAML jobs from clustered MSAs.")
    parser.add_argument("--chain-a-seq-file", type=Path, required=True)
    parser.add_argument("--chain-b-seq-file", type=Path, required=True)
    parser.add_argument("--chain-a-cluster-dir", type=Path, required=True)
    parser.add_argument("--chain-b-cluster-dir", type=Path, required=True)
    parser.add_argument("--top-a", type=int, default=8)
    parser.add_argument("--top-b", type=int, default=8)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--job-prefix", type=str, default="integrin_ab")
    parser.add_argument("--template-cif", type=Path, default=None)
    parser.add_argument(
        "--template-force-threshold",
        type=float,
        default=0.80,
        help=(
            "Legacy compatibility flag. The installed BoltzGen schema does not expose "
            "per-chain force-template thresholds, so template-guided jobs emit one "
            "file-context spec whenever --template-cif is provided."
        ),
    )
    parser.add_argument(
        "--include-empty-msa-control",
        action="store_true",
        help="Also write one control YAML with msa=empty for both chains.",
    )
    args = parser.parse_args()

    seq_a = read_text(args.chain_a_seq_file)
    seq_b = read_text(args.chain_b_seq_file)
    args.outdir.mkdir(parents=True, exist_ok=True)

    clusters_a = load_ranked_clusters(args.chain_a_cluster_dir, args.top_a)
    clusters_b = load_ranked_clusters(args.chain_b_cluster_dir, args.top_b)

    if args.template_cif is not None and not args.template_cif.exists():
        raise FileNotFoundError(f"Template CIF not found: {args.template_cif}")

    manifest_path = args.outdir / "jobs_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8") as manifest:
        manifest.write("job_name\tchain_a_msa\tchain_b_msa\tyaml_path\n")

        job_count = 0
        for i, msa_a in enumerate(clusters_a, start=1):
            for j, msa_b in enumerate(clusters_b, start=1):
                job_name = f"{args.job_prefix}_a{i:02d}_b{j:02d}"
                yaml_path = args.outdir / f"{job_name}.yaml"
                msa_a_path = str(msa_a.resolve())
                msa_b_path = str(msa_b.resolve())
                if args.template_cif is None:
                    doc = build_protein_doc(seq_a, seq_b, msa_a_path, msa_b_path)
                else:
                    doc = build_file_context_doc(args.template_cif, msa_a_path, msa_b_path)

                write_yaml(yaml_path, doc)
                manifest.write(
                    f"{job_name}\t{msa_a}\t{msa_b}\t{yaml_path}\n"
                )
                job_count += 1

        if args.include_empty_msa_control:
            job_name = f"{args.job_prefix}_control_empty_msa"
            yaml_path = args.outdir / f"{job_name}.yaml"
            if args.template_cif is None:
                doc = build_protein_doc(seq_a, seq_b, "empty", "empty")
            else:
                doc = build_file_context_doc(args.template_cif, "empty", "empty")
            write_yaml(yaml_path, doc)
            manifest.write(f"{job_name}\tempty\tempty\t{yaml_path}\n")
            job_count += 1

    if args.template_cif is not None:
        print(
            "[jobs] Note: template_force_threshold is retained for CLI compatibility "
            "but is not represented in the installed BoltzGen schema."
        )
    print(f"[jobs] wrote {job_count} BoltzGen YAML files in {args.outdir}")
    print(f"[jobs] manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
