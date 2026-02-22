# codex-AFCluster

AF-Cluster-driven conformational sampling pipeline for integrin alpha/beta.

This branch clusters chain-level MSAs with AF-Cluster, then feeds cluster-specific MSAs into Boltz `predict` to search for extended states.

## What is implemented

1. `scripts/cluster_chain_msa.py`
- Runs AF-Cluster on one chain MSA (`.a3m`).
- Writes per-cluster `.a3m` files + a cluster summary TSV.

2. `scripts/make_boltz_jobs_from_clusters.py`
- Builds Boltz YAML jobs from top chain-A and chain-B clusters.
- Optionally adds template hints from an mmCIF file.

3. `scripts/run_afcluster_pipeline.sh`
- End-to-end launcher: cluster both chain MSAs, generate Boltz jobs, run Boltz.

4. `scripts/rank_extended.py`
- Ranks predicted structures by an extension proxy (`max CA-CA distance`).

## Requirements

- Python: `afcluster`, `pyyaml`, `gemmi`
- CLI: `boltz`

Install minimal deps:

```bash
pip install afcluster pyyaml gemmi
```

## Quick start

```bash
cd codex-AFCluster

python scripts/extract_chain_sequences_from_pdb.py \
  --pdb ../data/template_example/seed_090_frame_000.pdb \
  --outdir ./runs/integrin_ab_afcluster/seq

./scripts/run_afcluster_pipeline.sh \
  --chain-a-seq-file ./runs/integrin_ab_afcluster/seq/chain_A.seq \
  --chain-b-seq-file ./runs/integrin_ab_afcluster/seq/chain_B.seq \
  --chain-a-msa ../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/0/non_pairing.a3m \
  --chain-b-msa ../data/outputs/outputs_integrin_alpha5_beta1/integrin_alpha5_beta1/msa/1/non_pairing.a3m \
  --template-cif ../data/template_example/seed_090_frame_000.cif \
  --outdir ./runs/integrin_ab_afcluster
```

Then rank for extended conformers:

```bash
python scripts/rank_extended.py \
  --pred-root ./runs/integrin_ab_afcluster/boltz_outputs \
  --out ./runs/integrin_ab_afcluster/extended_rank.tsv
```

## Notes

- This is an exploration pipeline; it does not guarantee extended conformations.
- For larger sweeps, increase `--top-a`, `--top-b`, and Boltz diffusion samples.
- If your mmcif DB is only on cluster, run these scripts there with cluster paths.
