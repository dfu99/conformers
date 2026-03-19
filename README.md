# conformers

Integrin conformer generation and validation pipelines for α5β1 and αVβ3.

## Pipelines

| Pipeline | Target | Method |
|----------|--------|--------|
| `pipelines/protenix-a5b1` | α5β1 tagged assembly | Staged Protenix inference + receptor-aligned merge |
| `pipelines/protenix-avb3-template` | αVβ3 template-guided | Protenix with template registration |
| `pipelines/boltz` | Heterodimer exploration | BoltzGen structure generation |
| `pipelines/afcluster` | MSA-clustered variants | AFCluster → BoltzGen |
| `pipelines/avb3-conformers` | αVβ3 conformer library | Steered MD → pseudo-AFM images |

## Data layout

```
data/
├── a5b1/           # α5β1 sequences and inputs
├── avb3/           # αVβ3 templates and inputs
└── runs/           # All pipeline outputs
    ├── a5b1/       #   Protenix, staged attachment
    ├── avb3/       #   Protenix, conformers, validation
    ├── boltz/      #   BoltzGen runs
    └── afcluster/  #   AFCluster runs
```

## Requirements

- **Protenix**: `venv_protenix` (do not share with Boltz)
- **Boltz/AFCluster**: `venv_boltz` (do not share with Protenix)
- **GPU**: A100 80GB for Protenix and BoltzGen on full integrins; RTX 6000 for MD
