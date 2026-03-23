#!/usr/bin/env python3
"""AlphaFold2 RandomWalk — weight perturbation for conformational diversity.

Implements the core idea from Taneja et al. (JCIM 2026): add Gaussian noise
to AF2 model parameters before each inference pass, producing diverse predicted
conformations instead of the single dominant state.

This script wraps the standard AF2 installation by:
1. Loading AF2 model parameters from checkpoint
2. Adding scaled Gaussian noise to model weights
3. Saving perturbed parameters as a new checkpoint
4. Running AF2 inference with the perturbed checkpoint
5. Repeating N times with different noise seeds

Usage (on PACE):
    module load alphafold/2.3.2
    python af2_random_walk.py \
        --fasta avb3_multimer.fasta \
        --output-dir af2_rw_outputs/ \
        --n-perturbations 10 \
        --noise-scale 0.01

Note: This modifies AF2 model weights in-place (using copies). It requires
write access to a working directory for perturbed checkpoints.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


def find_af2_params(data_dir: Path) -> list[Path]:
    """Find AlphaFold2 model parameter files."""
    params_dir = data_dir / "params"
    if not params_dir.exists():
        return []
    return sorted(params_dir.glob("params_model_*_multimer_v3.npz"))


def perturb_params(src_path: Path, dst_path: Path, noise_scale: float,
                   seed: int) -> dict:
    """Load AF2 params, add Gaussian noise, save perturbed version."""
    import numpy as np

    rng = np.random.RandomState(seed)
    params = dict(np.load(src_path, allow_pickle=True))

    n_params = 0
    n_perturbed = 0
    for key in params:
        arr = params[key]
        if not isinstance(arr, np.ndarray):
            continue
        n_params += arr.size
        if arr.dtype in (np.float32, np.float64) and arr.ndim >= 1:
            noise = rng.normal(0, noise_scale * np.std(arr), arr.shape).astype(arr.dtype)
            params[key] = arr + noise
            n_perturbed += arr.size

    np.savez(dst_path, **params)
    return {"n_params": n_params, "n_perturbed": n_perturbed,
            "noise_scale": noise_scale, "seed": seed}


def run_af2_with_params(fasta_path: Path, output_dir: Path,
                        data_dir: Path, params_dir: Path,
                        model_preset: str = "multimer",
                        db_preset: str = "reduced_dbs",
                        max_template_date: str = "2020-05-14",
                        n_cpus: int = 8) -> int:
    """Run AF2 inference using a specific params directory."""
    cmd = [
        "alphafold",
        f"--fasta_paths={fasta_path}",
        f"--max_template_date={max_template_date}",
        f"--model_preset={model_preset}",
        f"--db_preset={db_preset}",
        f"--output_dir={output_dir}",
        f"--data_dir={data_dir}",
        f"--uniref90_database_path={data_dir}/uniref90/uniref90.fasta",
        f"--mgnify_database_path={data_dir}/mgnify/mgy_clusters_2022_05.fa",
        f"--template_mmcif_dir={data_dir}/pdb_mmcif/mmcif_files",
        f"--obsolete_pdbs_path={data_dir}/pdb_mmcif/obsolete.dat",
        f"--uniprot_database_path={data_dir}/uniprot/uniprot.fasta",
        f"--pdb_seqres_database_path={data_dir}/pdb_seqres/pdb_seqres.txt",
        f"--small_bfd_database_path={data_dir}/small_bfd/bfd-first_non_consensus_sequences.fasta",
        "--use_gpu_relax=false",  # Skip relaxation for speed
        f"--num_multimer_predictions_per_model=1",  # 1 pred per model (5 total)
    ]

    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"

    result = subprocess.run(cmd, env=env, timeout=7200)  # 2h timeout per run
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", type=Path, required=True,
                        help="Input FASTA file (multimer format)")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for all perturbation runs")
    parser.add_argument("--data-dir", type=Path,
                        default=Path("/storage/coda1/d-pace_community/0/alphafold/alphafold_2.3.2_data"),
                        help="AF2 data directory on PACE")
    parser.add_argument("--n-perturbations", type=int, default=10,
                        help="Number of perturbed inference runs (default: 10)")
    parser.add_argument("--noise-scale", type=float, default=0.01,
                        help="Noise scale relative to param std (default: 0.01)")
    parser.add_argument("--base-seed", type=int, default=42,
                        help="Base random seed (default: 42)")
    parser.add_argument("--model-preset", default="multimer")
    parser.add_argument("--db-preset", default="reduced_dbs")
    parser.add_argument("--skip-unperturbed", action="store_true",
                        help="Skip the baseline (unperturbed) run")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Find AF2 model params
    param_files = find_af2_params(args.data_dir)
    if not param_files:
        print(f"ERROR: No AF2 params found in {args.data_dir}/params/")
        return 1
    print(f"Found {len(param_files)} AF2 model param files")

    # Step 0: Run baseline (unperturbed) if not skipping
    if not args.skip_unperturbed:
        baseline_dir = args.output_dir / "baseline"
        if not list(baseline_dir.glob("*/ranked_0.pdb")):
            print("\n=== Baseline (unperturbed) run ===")
            rc = run_af2_with_params(args.fasta, baseline_dir, args.data_dir,
                                     args.data_dir / "params",
                                     args.model_preset, args.db_preset)
            if rc != 0:
                print(f"  Baseline failed (exit {rc})")
        else:
            print("  Baseline already exists, skipping")

    # Step 1-N: Perturbed runs
    manifest = []
    for i in range(args.n_perturbations):
        seed = args.base_seed + i
        run_dir = args.output_dir / f"perturb_{i:03d}_seed{seed}"
        perturbed_params_dir = run_dir / "params"

        if list(run_dir.glob("*/ranked_0.pdb")):
            print(f"\n=== Perturbation {i} (seed={seed}) — already done ===")
            manifest.append({"index": i, "seed": seed, "dir": str(run_dir), "status": "cached"})
            continue

        print(f"\n=== Perturbation {i}/{args.n_perturbations} (seed={seed}, "
              f"noise={args.noise_scale}) ===")

        # Create perturbed params
        perturbed_params_dir.mkdir(parents=True, exist_ok=True)
        for pf in param_files:
            dst = perturbed_params_dir / pf.name
            if not dst.exists():
                info = perturb_params(pf, dst, args.noise_scale, seed)
                print(f"  Perturbed {pf.name}: {info['n_perturbed']:,} params")

        # Copy non-model files (e.g., params metadata)
        for f in (args.data_dir / "params").iterdir():
            if f.suffix != ".npz":
                dst = perturbed_params_dir / f.name
                if not dst.exists():
                    shutil.copy2(f, dst)

        # Create a modified data_dir symlink structure
        perturbed_data = run_dir / "data"
        perturbed_data.mkdir(parents=True, exist_ok=True)
        # Symlink everything except params
        for item in args.data_dir.iterdir():
            link = perturbed_data / item.name
            if not link.exists():
                if item.name == "params":
                    link.symlink_to(perturbed_params_dir)
                else:
                    link.symlink_to(item)

        # Run AF2 with perturbed params
        rc = run_af2_with_params(args.fasta, run_dir, perturbed_data,
                                 perturbed_params_dir,
                                 args.model_preset, args.db_preset)

        status = "completed" if rc == 0 else f"failed_exit_{rc}"
        manifest.append({"index": i, "seed": seed, "dir": str(run_dir),
                         "noise_scale": args.noise_scale, "status": status})
        print(f"  Status: {status}")

    # Write manifest
    manifest_path = args.output_dir / "af2_rw_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nManifest: {manifest_path}")
    print(f"Completed: {sum(1 for m in manifest if m['status'] in ('completed', 'cached'))}"
          f"/{len(manifest)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
