---
name: pace-slurm-ops
description: Use this skill when asked to submit, monitor, and fetch PACE cluster jobs (especially with ssh dfu71@login-phoenix.pace.gatech.edu, squeue/sacct status checks, and logs/results sync).
---

# PACE SLURM Ops

## Use This Skill When
- User asks to run jobs on PACE without repeated manual login steps.
- User asks for `squeue -u dfu71`, `sacct`, log checks, or artifact fetch from PACE.
- User asks for smoke validation of submit/watch/fetch workflow.

## Canonical Entry Point
- Use [`pipelines/protenix-a5b1/scripts/pace_minimal.sh`](../../../pipelines/protenix-a5b1/scripts/pace_minimal.sh) for `check`, `smoke`, `submit`, `watch`, and `fetch`.
- Default host alias in this repo is `pace-phoenix`.

## Preflight (Always)
1. Ensure VPN is connected.
2. Ensure key-based SSH auth works in non-interactive mode:
```bash
ssh -o BatchMode=yes pace-phoenix 'echo ok'
```
3. Verify remote SLURM availability:
```bash
PACE_HOST='pace-phoenix' pipelines/protenix-a5b1/scripts/pace_minimal.sh check
```

If auth fails with `Permission denied`, do one-time key enrollment:
```bash
ssh-copy-id -i ~/.ssh/id_ed25519.pub pace-phoenix
```
Do not automate password storage.

## Smoke Test Workflow (Hello World, 60s)
Run:
```bash
PACE_HOST='pace-phoenix' pipelines/protenix-a5b1/scripts/pace_minimal.sh smoke 5 60
```

During run, also verify user queue visibility:
```bash
ssh pace-phoenix "squeue -u dfu71"
```

Expected:
- Smoke submit prints `Smoke job submitted: <jobid>`.
- Watch loop prints `SQUEUE` lines and final `SACCT ... COMPLETED`.
- Fetch stage copies:
  - `logs/pace-smoke/pace_smoke_<jobid>.log`
  - `logs/pace-smoke/pace_smoke_<jobid>.err`
  - `data/runs/pace-smoke/smoke_<jobid>.txt`

## Real Pipeline Workflow
Submit:
```bash
PACE_HOST='pace-phoenix' pipelines/protenix-a5b1/scripts/pace_minimal.sh submit
```
Watch:
```bash
PACE_HOST='pace-phoenix' pipelines/protenix-a5b1/scripts/pace_minimal.sh watch <jobid>
```
Fetch:
```bash
PACE_HOST='pace-phoenix' pipelines/protenix-a5b1/scripts/pace_minimal.sh fetch <jobid>
```

## Operational Defaults
- PACE host alias: `pace-phoenix` (`dfu71@login-phoenix.pace.gatech.edu`).
- PACE account default for smoke jobs: `gts-yke8` (`PACE_ACCOUNT` override allowed).
- Smoke outputs are `.log` and `.err`; `.out` may be absent and is non-fatal.
- A5B1 submit script currently targets `--gres=gpu:RTX_6000:1`.

## Failure Handling
- `Could not resolve hostname` or timeout: VPN/routing issue.
- `Permission denied (publickey...)`: key not enrolled on remote.
- `--account option required`: missing account; set `PACE_ACCOUNT` or use default script behavior.
- If submit succeeds but watch fails, run:
```bash
ssh pace-phoenix "squeue -u dfu71"
ssh pace-phoenix "sacct -n -P -j <jobid> --format=JobIDRaw,State,ExitCode,Elapsed,End"
```
then run `fetch`.
