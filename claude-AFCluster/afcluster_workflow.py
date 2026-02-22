#!/usr/bin/env python3
"""AF-Cluster conformation sampling for integrin via LocalColabFold (AF2).

Clusters the MSA by gap-pattern, then runs ColabFold on each cluster subset.
AF2's recycling through the MSA representation module is sensitive to MSA
composition — different evolutionary subsets encode different structural
contexts and can unlock alternative conformations (e.g. extended integrin)
that the full MSA suppresses.

Reference: Del Alamo et al. (2022) eLife. https://doi.org/10.7554/eLife.75751

Requirements:
  - LocalColabFold installed and `colabfold_batch` in PATH
    https://github.com/YoshitakaMo/localcolabfold
  - scikit-learn  (pip install scikit-learn)
  - Per-chain a3m files from a prior MSA search (ColabFold or Protenix output)
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# A3M I/O
# ---------------------------------------------------------------------------

def parse_a3m(path: Path) -> List[Tuple[str, str]]:
    """Parse a3m → list of (header, ungapped_sequence).

    Lowercase letters are insertion columns in a3m format; strip them so
    all sequences share the same column width for distance computation.
    """
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    parts: List[str] = []

    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.rstrip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if header is not None and parts:
                    records.append((header, "".join(parts)))
                header = line[1:].strip() or f"seq_{len(records)}"
                parts = []
            else:
                parts.append(re.sub(r"[a-z]", "", line))   # strip inserts

    if header is not None and parts:
        records.append((header, "".join(parts)))

    return records


def write_a3m(path: Path, records: List[Tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# MSA clustering
# ---------------------------------------------------------------------------

def cluster_msa(
    records: List[Tuple[str, str]],
    n_clusters: int,
    random_state: int = 42,
) -> List[int]:
    """K-means on binary gap/residue encoding (1=residue, 0=gap).

    Returns integer cluster labels (one per record).
    Falls back to round-robin if sklearn is missing or N < n_clusters.
    """
    if len(records) <= n_clusters:
        print(f"  [warn] {len(records)} seqs ≤ n_clusters={n_clusters}; "
              "assigning one-per-cluster.")
        return list(range(len(records)))

    try:
        from sklearn.cluster import MiniBatchKMeans  # type: ignore
    except ImportError:
        print("  [warn] scikit-learn not found; using round-robin assignment.")
        return [i % n_clusters for i in range(len(records))]

    seqs = [s for _, s in records]
    L = len(seqs[0])
    X = np.zeros((len(seqs), L), dtype=np.float32)
    for i, s in enumerate(seqs):
        for j, c in enumerate(s):
            if c not in ("-", "X", "."):
                X[i, j] = 1.0

    km = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=random_state,
        n_init=10,
        batch_size=min(1024, len(seqs)),
    )
    return km.fit_predict(X).tolist()


# ---------------------------------------------------------------------------
# ColabFold input (paired a3m, FASTA, or query FASTA)
# ---------------------------------------------------------------------------

def build_colabfold_a3m(
    query_records: List[Tuple[str, str]],     # one per chain: (header, seq)
    subset_records: List[List[Tuple[str, str]]],  # per-chain MSA subsets
) -> str:
    """Build a single combined a3m string for ColabFold paired-chain input.

    ColabFold multi-chain format: chains separated by ':' in the query header,
    then individual per-chain MSA blocks separated by '#'.
    Format spec: https://github.com/sokrypton/ColabFold#faq
    """
    lines: List[str] = []

    # Query header: chain seqs joined by ':'
    query_seq = ":".join(seq for _, seq in query_records)
    query_header = ":".join(hdr for hdr, _ in query_records)
    lines.append(f">{query_header}")
    lines.append(query_seq)

    # Per-chain MSA blocks delimited by '#'
    for chain_idx, chain_recs in enumerate(subset_records):
        lines.append("#")
        for hdr, seq in chain_recs:
            lines.append(f">{hdr}")
            lines.append(seq)

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# ColabFold runner
# ---------------------------------------------------------------------------

def run_colabfold(
    input_a3m: Path,
    out_dir: Path,
    num_seeds: int,
    num_recycles: int,
    num_models: int,
    num_relax: int,
    amber: bool,
    colabfold_bin: str,
) -> None:
    """Run colabfold_batch on a pre-built a3m (no MSA search step)."""
    cmd = [
        colabfold_bin,
        str(input_a3m),
        str(out_dir),
        "--num-seeds", str(num_seeds),
        "--num-recycle", str(num_recycles),
        "--num-models", str(num_models),
        "--model-type", "alphafold2_multimer_v3",
        "--msa-mode", "custom",          # skip MMseqs2 search; use provided a3m
    ]
    if num_relax > 0:
        cmd += ["--num-relax", str(num_relax)]
    if not amber:
        cmd += ["--amber"]               # no-op placeholder; amber is default off

    print("\n[run]", " ".join(cmd))
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)


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
            resname = line[17:20].strip().upper()
            chain_id = line[21].strip() or "_"
            res_key = (line[22:26].strip(), line[26].strip())
            chains.setdefault(chain_id, [])
            seen.setdefault(chain_id, set())
            if res_key not in seen[chain_id]:
                seen[chain_id].add(res_key)
                chains[chain_id].append(AA3_TO_AA1.get(resname, "X"))
    if not chains:
        raise ValueError(f"No CA ATOM records in {pdb_path}")
    return {c: "".join(s) for c, s in chains.items()}


# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="AF-Cluster: cluster MSA subsets → LocalColabFold (AF2) inference."
    )
    p.add_argument("--input_pdb", type=Path,
                   default=Path("../data/template_example/seed_090_frame_000.pdb"))
    p.add_argument("--chain_order", type=str, default="A,B",
                   help="Comma-separated chain IDs in PDB order.")
    p.add_argument("--msa_root", type=Path, required=True,
                   help="MSA root with per-chain non_pairing.a3m files "
                        "(subdirs 0/, 1/, ...). Compatible with Protenix or "
                        "ColabFold --save-all outputs.")
    p.add_argument("--workflow_dir", type=Path,
                   default=Path("../data/afcluster_outputs"))
    p.add_argument("--n_clusters", type=int, default=8)
    p.add_argument("--min_cluster_size", type=int, default=5)
    p.add_argument("--cluster_on_chain", type=int, default=0,
                   help="Chain index (0-based) whose MSA drives clustering.")
    p.add_argument("--colabfold_bin", type=str, default="colabfold_batch",
                   help="Path to colabfold_batch executable.")
    p.add_argument("--num_seeds", type=int, default=3,
                   help="ColabFold --num-seeds (random seeds per model).")
    p.add_argument("--num_recycles", type=int, default=3)
    p.add_argument("--num_models", type=int, default=5,
                   help="Number of AF2 models to run (1–5).")
    p.add_argument("--num_relax", type=int, default=0,
                   help="Amber relaxation steps (0 = skip, slow).")
    p.add_argument("--run", action="store_true",
                   help="Execute colabfold_batch (else generate inputs only).")
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

    msa_root = args.msa_root.resolve()
    workflow_dir = args.workflow_dir.resolve()
    inputs_dir = workflow_dir / "inputs"
    outputs_dir = workflow_dir / "outputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    chain_to_seq = parse_pdb_sequences(input_pdb)
    print("Chains from PDB:")
    for cid in chain_order:
        if cid not in chain_to_seq:
            raise ValueError(f"Chain '{cid}' not found. Available: {sorted(chain_to_seq)}")
        print(f"  {cid}: len={len(chain_to_seq[cid])}")

    # Load per-chain unpaired MSAs
    unpaired: List[List[Tuple[str, str]]] = []
    for i in range(len(chain_order)):
        for candidate in [
            msa_root / str(i) / "non_pairing.a3m",
            msa_root / f"{i}.a3m",
        ]:
            if candidate.exists():
                recs = parse_a3m(candidate)
                print(f"  Chain {chain_order[i]} MSA: {len(recs)} seqs from {candidate}")
                unpaired.append(recs)
                break
        else:
            raise FileNotFoundError(
                f"No MSA found for chain index {i} under {msa_root}"
            )

    # Cluster on pivot chain
    pivot = args.cluster_on_chain
    print(f"\nClustering chain {chain_order[pivot]} MSA ({len(unpaired[pivot])} seqs) "
          f"into {args.n_clusters} clusters ...")
    labels = cluster_msa(unpaired[pivot], n_clusters=args.n_clusters)

    cluster_ids = sorted(set(labels))
    sizes = {c: labels.count(c) for c in cluster_ids}
    print("Cluster sizes:", sizes)

    # Build one ColabFold a3m per cluster
    query_records = [(f"chain_{cid}", chain_to_seq[cid]) for cid in chain_order]
    generated: List[Tuple[int, Path]] = []

    for cid in cluster_ids:
        indices = [i for i, lbl in enumerate(labels) if lbl == cid]
        if len(indices) < args.min_cluster_size:
            print(f"  Cluster {cid}: {len(indices)} seqs < min={args.min_cluster_size}, skip.")
            continue

        # Build per-chain subset (reindex safely for chains with fewer seqs)
        subset_per_chain: List[List[Tuple[str, str]]] = []
        for chain_recs in unpaired:
            safe = [chain_recs[i] for i in indices if i < len(chain_recs)]
            subset_per_chain.append(safe if safe else [chain_recs[0]])

        a3m_content = build_colabfold_a3m(query_records, subset_per_chain)
        a3m_path = inputs_dir / f"cluster{cid:02d}.a3m"
        a3m_path.write_text(a3m_content, encoding="utf-8")
        generated.append((cid, a3m_path))
        print(f"  Cluster {cid}: wrote {a3m_path.name} "
              f"({len(indices)} seqs per chain)")

    print(f"\nGenerated {len(generated)} cluster a3m inputs under {inputs_dir}")

    if not args.run:
        print("Dry run — pass --run to execute colabfold_batch per cluster.")
        return 0

    import shutil
    colabfold_bin = shutil.which(args.colabfold_bin) or args.colabfold_bin
    if not Path(colabfold_bin).exists() and not shutil.which(colabfold_bin):
        raise FileNotFoundError(
            f"colabfold_batch not found: {colabfold_bin}\n"
            "Install LocalColabFold: https://github.com/YoshitakaMo/localcolabfold"
        )

    for cid, a3m_path in generated:
        run_colabfold(
            input_a3m=a3m_path,
            out_dir=outputs_dir / f"cluster{cid:02d}",
            num_seeds=args.num_seeds,
            num_recycles=args.num_recycles,
            num_models=args.num_models,
            num_relax=args.num_relax,
            amber=False,
            colabfold_bin=colabfold_bin,
        )

    pdb_count = sum(1 for _ in outputs_dir.rglob("*.pdb"))
    cif_count = sum(1 for _ in outputs_dir.rglob("*.cif"))
    print(f"\nDone. {pdb_count} PDB / {cif_count} CIF files under {outputs_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
