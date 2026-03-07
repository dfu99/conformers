# protenix-a5b1

A5B1 heterodimer + staged tag attachment pipeline.

## Data stream
- Sequence source: `data/a5b1/sequences/sequences_updated`
- Accepted heterodimer source (default): `data/runs/a5b1/protenix/outputs_integrin_alpha5_beta1/...`
- Workflow outputs: `data/runs/a5b1/staged_attachment/...`

## Why staged
Direct 4-chain prediction can place tags unrealistically in the core. This pipeline enforces a staged process:
1. Keep accepted A5B1 heterodimer as receptor reference.
2. Predict `A5B1 + SpyTag` (stage 1).
3. Predict `A5B1 + Streptavidin` (stage 2).
4. Align stage outputs back to the accepted receptor and merge ligand chains into one final structure.

## Scripts
- `scripts/setup_staged_attachment_workflow.py`:
  creates constraints/manifests and picks the heterodimer CIF.
- `scripts/build_staged_protenix_inputs.py`:
  builds stage-1 and stage-2 Protenix JSON inputs.
- `scripts/merge_staged_tagged_complex.py`:
  selects best stage outputs and writes one merged final CIF/PDB.
- `scripts/run_complete_tagged_pipeline.sh`:
  end-to-end runner.
- `scripts/submit_complete_tagged_pipeline_slurm.sh`:
  sbatch entrypoint for cluster.

## Cluster run
```bash
sbatch ~/scratch/conformers/pipelines/protenix-a5b1/scripts/submit_complete_tagged_pipeline_slurm.sh
```

## Key outputs
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.cif`
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.pdb`
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.merge_summary.json`

## Note on Protenix-Dock
`protenix-dock` is built for protein-ligand docking with ligand SDF inputs, not protein-protein tag attachment. For this A5B1+SpyTag+Streptavidin task, staged Protenix multichain inference + receptor-aligned merge is the compatible path in this repo.
