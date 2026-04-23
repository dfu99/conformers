# Code Style — Stub

**Canonical source:** `/home/dan/Documents/code/development/global/code-style.md`

All code style rules live there. This file is a pointer. **Do not edit it.**

When this project discovers a rule worth sharing system-wide, edit the canonical file (not this stub). That update propagates to every project on next read.

### For conformers, apply:

- Core Rules apply. **Project Overlay:** research project with GPU compute. Default CPU-only: no accidental `torch.cuda.is_available()` gates that silently grab the workstation GPU. PACE SLURM scripts follow `PACE.md` (account `gts-yke8`, pick the cheapest fitting GPU, always `module load cuda`). Training or evaluation scripts must emit a greppable final metrics summary; silent exit is not acceptable.

See the canonical file for full rule text, language-specific subsections, and the current changelog.
