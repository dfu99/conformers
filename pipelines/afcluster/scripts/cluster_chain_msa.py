#!/usr/bin/env python3
"""Cluster one chain MSA with AF-Cluster and write per-cluster A3Ms.

Compatible with multiple AFCluster API variants.
"""

from __future__ import annotations

import argparse
import inspect
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cluster a chain MSA with AF-Cluster.")
    parser.add_argument("--input-a3m", required=True, type=Path, help="Input .a3m file")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument(
        "--max-dist",
        type=float,
        default=0.9,
        help="Compatibility parameter: used only if AFCluster API supports max_dist.",
    )
    parser.add_argument("--min-samples", type=int, default=2, help="DBSCAN min_samples")
    parser.add_argument(
        "--n-controls",
        type=int,
        default=0,
        help="AFCluster n_controls (0 disables random controls)",
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=-1.0,
        help="DBSCAN eps; if <0, script tries AFCluster gridsearch_eps",
    )
    parser.add_argument(
        "--keep-top",
        type=int,
        default=12,
        help="Keep only top N cluster a3m files by size (<=0 keeps all)",
    )
    return parser.parse_args()


def numeric_hint(path: Path) -> int:
    m = re.search(r"(\d+)", path.stem)
    if not m:
        return 10**9
    return int(m.group(1))


def cluster_sort_key(path: Path) -> tuple[int, int, int, str]:
    is_outlier = 1 if path.name == "cluster_outliers.a3m" else 0
    return (is_outlier, -path.stat().st_size, numeric_hint(path), path.name)


def filter_supported_kwargs(callable_obj, kwargs: dict) -> dict:
    """Keep only kwargs supported by callable_obj signature.

    If signature cannot be inspected, return input kwargs as-is.
    """
    try:
        sig = inspect.signature(callable_obj)
    except (TypeError, ValueError):
        return kwargs

    params = sig.parameters
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params.values()):
        return kwargs
    return {k: v for k, v in kwargs.items() if k in params}


def main() -> int:
    args = parse_args()

    if not args.input_a3m.exists():
        raise FileNotFoundError(f"Input MSA not found: {args.input_a3m}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    try:
        from afcluster import AFCluster, read_a3m
    except Exception as exc:
        raise RuntimeError(
            "Could not import afcluster. Install with: pip install afcluster"
        ) from exc

    try:
        msa = read_a3m(str(args.input_a3m))
    except TypeError:
        msa = read_a3m(args.input_a3m)

    ctor_kwargs = filter_supported_kwargs(
        AFCluster.__init__,
        {
            "max_dist": args.max_dist,
            "min_samples": args.min_samples,
            "n_controls": args.n_controls,
        },
    )

    # AFCluster API changed across versions; pass only supported constructor kwargs.
    clusterer = AFCluster(**ctor_kwargs)

    eps = args.eps
    if eps < 0:
        if hasattr(clusterer, "gridsearch_eps"):
            try:
                eps = float(clusterer.gridsearch_eps(msa))
            except TypeError:
                eps = float(clusterer.gridsearch_eps(str(args.input_a3m)))
            print(f"[afcluster] gridsearch eps={eps:.4f}")
        else:
            eps = 0.6
            print(f"[afcluster] gridsearch_eps unavailable; using eps={eps:.4f}")

    cluster_kwargs = filter_supported_kwargs(
        clusterer.cluster,
        {
            "eps": eps,
            "min_samples": args.min_samples,
            "max_dist": args.max_dist,
        },
    )

    try:
        cluster_result = clusterer.cluster(msa, **cluster_kwargs)
    except TypeError:
        try:
            cluster_result = clusterer.cluster(msa)
        except TypeError:
            cluster_result = clusterer.cluster(str(args.input_a3m), **cluster_kwargs)

    if hasattr(cluster_result, "to_csv"):
        cluster_result.to_csv(args.outdir / "cluster_result.tsv", sep="\t", index=False)

    try:
        clusterer.write_a3m(str(args.outdir))
    except TypeError:
        clusterer.write_a3m(args.outdir)

    a3m_files = sorted(args.outdir.glob("*.a3m"), key=lambda p: (numeric_hint(p), p.name))
    if not a3m_files:
        raise RuntimeError(
            f"No cluster .a3m files were written to {args.outdir}. "
            "This often means an AFCluster API/data format mismatch."
        )

    # Keep bona fide clusters ahead of the outlier bucket, then break ties by size.
    ranked = sorted(a3m_files, key=cluster_sort_key)
    if args.keep_top > 0:
        keep = set(ranked[: args.keep_top])
        for path in ranked[args.keep_top :]:
            path.unlink(missing_ok=True)
        ranked = [p for p in ranked if p in keep]

    summary = args.outdir / "clusters.tsv"
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("rank\tfilename\tbytes\n")
        for i, path in enumerate(ranked, start=1):
            handle.write(f"{i}\t{path.name}\t{path.stat().st_size}\n")

    print(f"[afcluster] wrote {len(ranked)} cluster MSAs to: {args.outdir}")
    print(f"[afcluster] summary: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
