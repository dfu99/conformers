#!/usr/bin/env python3
"""Run ProteinTTT (test-time training) on AVB3 integrin chains.

Extracts AVB3 chain sequences from a reference PDB, runs ESMFold TTT
per-chain (since full heterodimer exceeds single-sequence ESMFold limits),
and outputs predicted PDB structures with pLDDT scores.

Usage:
    python run_proteinttt.py \
        --conformers-root ~/scratch/conformers \
        --proteinttt-root ~/scratch/ProteinTTT \
        --work-dir ~/scratch/conformers/data/runs/avb3/proteinttt
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


def extract_sequences_from_pdb(pdb_path: Path) -> dict[str, str]:
    """Extract amino acid sequences per chain from PDB using CA atoms."""
    aa3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
    }
    chains: dict[str, list[tuple[int, str]]] = {}
    seen = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain_id = line[21]
        resseq = int(line[22:26])
        resname = line[17:20].strip()
        key = (chain_id, resseq)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(chain_id, []).append((resseq, aa3to1.get(resname, "X")))

    return {c: "".join(aa for _, aa in sorted(res)) for c, res in sorted(chains.items())}


def run_ttt_on_chain(chain_id: str, sequence: str, work_dir: Path,
                     ttt_steps: int, ttt_lr: float, ttt_lora_rank: int,
                     ttt_batch_size: int, ttt_crop_size: int,
                     msa_path: str | None, reference_pdb: str | None) -> dict:
    """Run ProteinTTT on a single chain."""
    import torch
    from proteinttt.base import TTTConfig
    from proteinttt.models.esmfold import ESMFoldTTT

    chain_dir = work_dir / f"chain_{chain_id}"
    chain_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Chain {chain_id}: {len(sequence)} residues")
    print(f"TTT: steps={ttt_steps}, lr={ttt_lr}, lora_rank={ttt_lora_rank}")
    print(f"{'='*60}")

    cfg = TTTConfig(
        lr=ttt_lr,
        steps=ttt_steps,
        lora_rank=ttt_lora_rank,
        lora_alpha=32.0,
        batch_size=ttt_batch_size,
        crop_size=ttt_crop_size,
        optimizer="sgd",
        eval_each_step=True,
    )

    # Load ESMFold with TTT
    print("Loading ESMFold model...")
    model = ESMFoldTTT.from_pretrained(ttt_cfg=cfg)
    model = model.cuda()
    model.eval()

    # Baseline prediction (before TTT)
    print("Running baseline prediction (no TTT)...")
    t0 = time.time()
    with torch.no_grad():
        baseline_output = model.infer(sequence)
    baseline_plddt = baseline_output["mean_plddt"].item()
    baseline_pdb = model.output_to_pdb(baseline_output)
    baseline_time = time.time() - t0

    baseline_pdb_str = baseline_pdb[0] if isinstance(baseline_pdb, list) else baseline_pdb
    (chain_dir / f"baseline_chain_{chain_id}.pdb").write_text(baseline_pdb_str)
    print(f"  Baseline pLDDT: {baseline_plddt:.2f} ({baseline_time:.1f}s)")

    # Run TTT
    print(f"Running TTT ({ttt_steps} steps)...")
    t0 = time.time()
    results = model.ttt(
        seq=sequence,
        msa_pth=Path(msa_path) if msa_path else None,
        correct_pdb_path=Path(reference_pdb) if reference_pdb else None,
    )
    ttt_time = time.time() - t0

    # Post-TTT prediction
    print("Running post-TTT prediction...")
    with torch.no_grad():
        ttt_output = model.infer(sequence)
    ttt_plddt = ttt_output["mean_plddt"].item()
    ttt_pdb = model.output_to_pdb(ttt_output)

    ttt_pdb_str = ttt_pdb[0] if isinstance(ttt_pdb, list) else ttt_pdb
    (chain_dir / f"ttt_chain_{chain_id}.pdb").write_text(ttt_pdb_str)
    print(f"  Post-TTT pLDDT: {ttt_plddt:.2f} ({ttt_time:.1f}s)")
    print(f"  Improvement: {ttt_plddt - baseline_plddt:+.2f}")

    # Save per-residue pLDDT
    if "plddt" in ttt_output:
        plddt_per_res = ttt_output["plddt"].cpu().numpy().flatten().tolist()
    else:
        plddt_per_res = []

    # Reset model for next chain
    model.ttt_reset()

    result = {
        "chain": chain_id,
        "sequence_length": len(sequence),
        "baseline_plddt": baseline_plddt,
        "ttt_plddt": ttt_plddt,
        "plddt_improvement": ttt_plddt - baseline_plddt,
        "ttt_time_seconds": ttt_time,
        "baseline_pdb": str(chain_dir / f"baseline_chain_{chain_id}.pdb"),
        "ttt_pdb": str(chain_dir / f"ttt_chain_{chain_id}.pdb"),
        "plddt_per_residue": plddt_per_res,
    }

    (chain_dir / "result.json").write_text(json.dumps(result, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conformers-root", type=Path, required=True)
    parser.add_argument("--proteinttt-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--ttt-steps", type=int, default=30)
    parser.add_argument("--ttt-lr", type=float, default=4e-4)
    parser.add_argument("--ttt-lora-rank", type=int, default=8)
    parser.add_argument("--ttt-batch-size", type=int, default=4)
    parser.add_argument("--ttt-crop-size", type=int, default=1024)
    parser.add_argument("--msa-path", default=None)
    parser.add_argument("--reference-pdb", default=None)
    args = parser.parse_args()

    # Add ProteinTTT to path
    sys.path.insert(0, str(args.proteinttt_root))

    args.work_dir.mkdir(parents=True, exist_ok=True)

    # Find reference PDB for sequence extraction
    ref_pdb = None
    for candidate in [
        args.conformers_root / "data/avb3/template_example/seed_090_frame_000.pdb",
        Path("/home2/Documents/code/RoyalMD/test_systems/AVB3_clean.pdb"),
        Path("/home2/Documents/code/RoyalMD/results-avb3_04/splits/production_trajectory_frame_000.pdb"),
    ]:
        if candidate.exists():
            ref_pdb = candidate
            break

    # Also check PACE paths
    if ref_pdb is None:
        for candidate in [
            args.conformers_root / "data/avb3/template_example/seed_090_frame_000.pdb",
            Path.home() / "scratch/RoyalMD/test_systems/AVB3_clean.pdb",
        ]:
            if candidate.exists():
                ref_pdb = candidate
                break

    if ref_pdb is None:
        print("ERROR: No reference AVB3 PDB found for sequence extraction.")
        return 1

    print(f"Reference PDB: {ref_pdb}")
    chain_seqs = extract_sequences_from_pdb(ref_pdb)
    for c, s in chain_seqs.items():
        print(f"  Chain {c}: {len(s)} residues")

    # Run TTT per chain (full heterodimer is too large for single ESMFold pass)
    all_results = []
    for chain_id, seq in chain_seqs.items():
        result = run_ttt_on_chain(
            chain_id, seq, args.work_dir,
            args.ttt_steps, args.ttt_lr, args.ttt_lora_rank,
            args.ttt_batch_size, args.ttt_crop_size,
            args.msa_path, args.reference_pdb,
        )
        all_results.append(result)

    # Save combined results
    summary = {
        "reference_pdb": str(ref_pdb),
        "chains": all_results,
        "total_plddt_improvement": sum(r["plddt_improvement"] for r in all_results),
    }
    summary_path = args.work_dir / "proteinttt_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"\nSummary saved to {summary_path}")

    for r in all_results:
        print(f"  Chain {r['chain']}: baseline={r['baseline_plddt']:.2f} → "
              f"TTT={r['ttt_plddt']:.2f} ({r['plddt_improvement']:+.2f})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
