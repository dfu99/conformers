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
  selects stage outputs (`ranking`, `tail_distance`, or `hybrid`) and writes one merged final CIF/PDB.
- `scripts/run_complete_tagged_pipeline.sh`:
  end-to-end runner.
- `scripts/submit_complete_tagged_pipeline_slurm.sh`:
  sbatch entrypoint for cluster.
- `scripts/build_conjugates_first_protenix_inputs.py`:
  experimental branch input builder (`alpha+streptavidin`, `beta+spytag`, then `alpha+beta` docking).
- `scripts/run_conjugates_first_pipeline.sh`:
  experimental branch runner for conjugates-first strategy.

## Cluster run
```bash
sbatch ~/scratch/conformers/pipelines/protenix-a5b1/scripts/submit_complete_tagged_pipeline_slurm.sh
```

## Tail-aware sample selection
By default, merged output uses per-stage best `ranking_score` sample.

To prioritize tail-proximal placement during merge:
```bash
SELECTION_MODE=tail_distance \
STAGE1_LIGAND_ANCHOR_RESIDUE=10 \
STAGE2_LIGAND_ANCHOR_RESIDUE=1 \
bash pipelines/protenix-a5b1/scripts/run_complete_tagged_pipeline.sh
```

Hybrid mode (prefer high ranking among near-tail candidates):
```bash
SELECTION_MODE=hybrid \
STAGE1_MAX_TAIL_DISTANCE=20 \
STAGE2_MAX_TAIL_DISTANCE=30 \
bash pipelines/protenix-a5b1/scripts/run_complete_tagged_pipeline.sh
```

Notes:
- SpyTag reactive residue defaults to Asp10 from sequence `RGVPHIVMVDAYKRYK`.
- Tail residues are auto-read from `inputs/components.json` (`A:971`, `B:818` in current run).

## Experimental: Conjugates-first branch
Run:
```bash
bash pipelines/protenix-a5b1/scripts/run_conjugates_first_pipeline.sh
```

Design references used:
- `references/Summary.pdf`
- `references/Group Design PPT.pdf`

This branch currently performs:
1. `alpha5-Avi + streptavidin` conjugate inference
2. `beta1-SpyCatcher + spytag` conjugate inference
3. `alpha5-Avi + beta1-SpyCatcher` docking stage

Outputs summary JSON:
- `data/runs/a5b1/conjugates_first/outputs/final/conjugates_first_summary.json`

## Minimal remote control (submit/watch/fetch)
Use `scripts/pace_minimal.sh` to avoid repeated interactive SSH login loops.

### SSH auth (secure, non-interactive)
Do not store your GT account password for automation. Use an SSH key and keychain:

```bash
ssh-add --apple-use-keychain ~/.ssh/id_ed25519
ssh-copy-id -i ~/.ssh/id_ed25519.pub dfu71@login-phoenix.pace.gatech.edu
```

Optional SSH config alias:
```sshconfig
Host pace-phoenix
  HostName login-phoenix.pace.gatech.edu
  User dfu71
  IdentityFile ~/.ssh/id_ed25519
  IdentitiesOnly yes
  AddKeysToAgent yes
  UseKeychain yes
```

```bash
chmod +x pipelines/protenix-a5b1/scripts/pace_minimal.sh
```

Quick smoke check (tiny SLURM job + watch + fetch):
```bash
PACE_HOST='dfu71@login-phoenix.pace.gatech.edu' pipelines/protenix-a5b1/scripts/pace_minimal.sh check
PACE_HOST='dfu71@login-phoenix.pace.gatech.edu' pipelines/protenix-a5b1/scripts/pace_minimal.sh smoke 10 60
```
The smoke job writes `pace_smoke_<jobid>.log` and `pace_smoke_<jobid>.err`, sleeps for 60 seconds by default, and fetches both logs plus a small result file.

Submit real A5B1 tagged run, then watch and fetch:
```bash
PACE_HOST='dfu71@login-phoenix.pace.gatech.edu' pipelines/protenix-a5b1/scripts/pace_minimal.sh submit
PACE_HOST='dfu71@login-phoenix.pace.gatech.edu' pipelines/protenix-a5b1/scripts/pace_minimal.sh watch <job_id>
PACE_HOST='dfu71@login-phoenix.pace.gatech.edu' pipelines/protenix-a5b1/scripts/pace_minimal.sh fetch <job_id>
```

Defaults:
- Remote host alias: `pace` (override `PACE_HOST`).
- Remote repo root: `$HOME/scratch/conformers` (override `PACE_REMOTE_ROOT`).
- Submit script requests A100 80GB: `--partition=gpu-a100 --gres=gpu:A100:1 --constraint=A100-80GB`.

## Key outputs
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.cif`
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.pdb`
- `data/runs/a5b1/staged_attachment/outputs/final/a5b1_tagged_complete.merge_summary.json`

## Note on Protenix-Dock
`protenix-dock` is built for protein-ligand docking with ligand SDF inputs, not protein-protein tag attachment. For this A5B1+SpyTag+Streptavidin task, staged Protenix multichain inference + receptor-aligned merge is the compatible path in this repo.
