# project_bat1k_longevity

Workflows supporting the manuscript
*"Bat adaptations in mitochondrial damage recognition revealed through comparative multi-omics."*

## Overview

Genome assemblies and annotations for 23 bat and 7 non-bat species are used to
test which genes and gene families show signatures of adaptation associated with the
hallmarks of aging. The repository collects the analysis pipelines and the structural
modeling of bat-specific FBXO7 variants.

## Analyses

Each analysis lives in its own subfolder under [`analyses/`](analyses/), with a README
documenting objective, inputs, pipeline steps, outputs, and the scripts alongside it.

| Analysis | Method | Documentation |
|---|---|---|
| Gene family expansion / contraction | OrthoFinder + CAFE5 | [`analyses/cafe/`](analyses/cafe/) |
| Annotation completeness (BUSCO) | compleasm | [`analyses/compleasm/`](analyses/compleasm/) |
| Signatures of selection | HyPhy (aBSREL, MEME, RELAX, BGM) | [`analyses/hyphy/`](analyses/hyphy/) |
| Comparative transcriptomics | EVE (Expression Variance and Evolution) | [`analyses/eve/`](analyses/eve/) |
| FBXO7–PINK1–PSMF1 variant modeling | Boltz-2 + PyRosetta ΔΔG | [`analyses/fbxo7_ternary_modeling/`](analyses/fbxo7_ternary_modeling/) |
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

Each analysis folder contains a `README.md` and the scripts referenced in it.
