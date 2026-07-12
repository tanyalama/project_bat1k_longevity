# project_bat1k_longevity

Compiled by Tanya M. Lama (tlama@smith.edu)

Workflows supporting the manuscript
*"Bat adaptations in mitochondrial damage recognition revealed through comparative multi-omics."*

## Overview

Genome assemblies and annotations for 23 bat and 7 non-bat species were used to
test which genes and gene families show signatures of adaptation associated with the
hallmarks of aging. This repository includes several workflows including the structural
modeling of bat-specific FBXO7 variants.

## Analyses

Each analysis is contained within a subfolder under [`analyses/`](analyses/), including a README
documenting the objective, inputs, pipeline, outputs, and scripts required for reproducibility.

| Analysis | Method | Documentation |
|---|---|---|
| Gene family expansion / contraction | OrthoFinder + CAFE5 | [`analyses/cafe/`](analyses/cafe/) |
| Annotation completeness (BUSCO) | compleasm | [`analyses/compleasm/`](analyses/compleasm/) |
| Signatures of selection | HyPhy (aBSREL, MEME, RELAX, BGM) | [`analyses/hyphy/`](analyses/hyphy/) |
| Comparative transcriptomics | EVE (Expression Variance and Evolution) | [`analyses/eve/`](analyses/eve/) |
| FBXO7–PINK1–PSMF1 variant modeling | stage 2 Boltz-2 protocol + PyRosetta ΔΔG | [`analyses/fbxo7_ternary_modeling/`](analyses/fbxo7_ternary_modeling/) |
| Western blot quantification | Welch *t*-tests + ggplot2 boxplots | [`analyses/western_blots/`](analyses/western_blots/) |
| Branch length vs. genes under selection | Negative binomial regression | [count2branches](https://github.com/lmdavalos/count2branches) |

## Repository layout

```
analyses/
  cafe/                     Gene family expansion / contraction (OrthoFinder + CAFE5)
  compleasm/                Annotation completeness (compleasm / BUSCO)
  hyphy/                    Signatures of selection (HyPhy)
  eve/                      Comparative transcriptomics (EVE)
  fbxo7_ternary_modeling/   Structural modeling of bat FBXO7 variants
  western_blots/            Western blot densitometry (mouse vs. Myotis bat)
```
