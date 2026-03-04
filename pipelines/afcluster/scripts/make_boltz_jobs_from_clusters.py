#!/usr/bin/env python3
"""Generate Boltz YAML jobs from AF-Cluster outputs for chain A/B."""

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
                if top_n > 0 and len(paths) >= top_n:
                    break
    else:
        paths = sorted(cluster_dir.glob("*.a3m"), key=lambda p: p.stat().st_size, reverse=True)
        if top_n > 0:
            paths = paths[:top_n]
    if not paths:
        raise RuntimeError(f"No cluster .a3m files found in {cluster_dir}")
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Boltz YAML jobs from clustered MSAs.")
    parser.add_argument("--chain-a-seq-file", type=Path, required=True)
    parser.add_argument("--chain-b-seq-file", type=Path, required=True)
    parser.add_argument("--chain-a-cluster-dir", type=Path, required=True)
    parser.add_argument("--chain-b-cluster-dir", type=Path, required=True)
    parser.add_argument("--top-a", type=int, default=8)
    parser.add_argument("--top-b", type=int, default=8)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--job-prefix", type=str, default="integrin_ab")
    parser.add_argument("--template-cif", type=Path, default=None)
    parser.add_argument("--template-force-threshold", type=float, default=0.80)
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

    template_block_a = []
    template_block_b = []
    if args.template_cif is not None:
        template_path = str(args.template_cif.resolve())
        template_block_a = [
            {
                "cif": template_path,
                "chain_id": "A",
                "template_inclusion_prob": 1.0,
                "template_filter": True,
                "force_template": True,
                "force_template_threshold": float(args.template_force_threshold),
            }
        ]
        template_block_b = [
            {
                "cif": template_path,
                "chain_id": "B",
                "template_inclusion_prob": 1.0,
                "template_filter": True,
                "force_template": True,
                "force_template_threshold": float(args.template_force_threshold),
            }
        ]

    manifest_path = args.outdir / "jobs_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8") as manifest:
        manifest.write("job_name\tchain_a_msa\tchain_b_msa\tyaml_path\n")

        job_count = 0
        for i, msa_a in enumerate(clusters_a, start=1):
            for j, msa_b in enumerate(clusters_b, start=1):
                job_name = f"{args.job_prefix}_a{i:02d}_b{j:02d}"
                yaml_path = args.outdir / f"{job_name}.yaml"

                doc = {
                    "version": 1,
                    "sequences": [
                        {
                            "protein": {
                                "id": "A",
                                "sequence": seq_a,
                                "msa": str(msa_a.resolve()),
                            }
                        },
                        {
                            "protein": {
                                "id": "B",
                                "sequence": seq_b,
                                "msa": str(msa_b.resolve()),
                            }
                        },
                    ],
                }
                if template_block_a:
                    doc["sequences"][0]["protein"]["templates"] = template_block_a
                if template_block_b:
                    doc["sequences"][1]["protein"]["templates"] = template_block_b

                yaml_path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
                manifest.write(
                    f"{job_name}\t{msa_a}\t{msa_b}\t{yaml_path}\n"
                )
                job_count += 1

        if args.include_empty_msa_control:
            job_name = f"{args.job_prefix}_control_empty_msa"
            yaml_path = args.outdir / f"{job_name}.yaml"
            doc = {
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": seq_a, "msa": "empty"}},
                    {"protein": {"id": "B", "sequence": seq_b, "msa": "empty"}},
                ],
            }
            if template_block_a:
                doc["sequences"][0]["protein"]["templates"] = template_block_a
            if template_block_b:
                doc["sequences"][1]["protein"]["templates"] = template_block_b
            yaml_path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
            manifest.write(f"{job_name}\tempty\tempty\t{yaml_path}\n")
            job_count += 1

    print(f"[jobs] wrote {job_count} Boltz YAML files in {args.outdir}")
    print(f"[jobs] manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
