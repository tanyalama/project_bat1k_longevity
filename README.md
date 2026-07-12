# project_bat1k_longevity

Comparative multi-omics of bat lifespan evolution, supporting the manuscript
*"Bat adaptations in mitochondrial damage recognition revealed through comparative multi-omics."*

## Overview

Genome annotations for 23 bat and 7 non-bat species (generated with TOGA1) are used to
test which genes and gene families show signatures of adaptation associated with the
hallmarks of aging. The repository collects the analysis pipelines and the structural
modeling of bat-specific FBXO7 variants.

## Analyses

| Analysis | Method | Documentation |
|---|---|---|
| Gene family expansion / contraction | OrthoFinder + CAFE5 | [`scripts/cafe.md`](scripts/cafe.md) |
| Annotation completeness (BUSCO) | compleasm | [`scripts/compleasm.md`](scripts/compleasm.md) |
| Signatures of selection | HyPhy (aBSREL, MEME, RELAX, BGM) | [`scripts/hyphy_tools.md`](scripts/hyphy_tools.md) |
| Gene enrichment | gprofiler2, NIH DAVID | — |
| Branch length vs. genes under selection | Negative binomial regression | [count2branches](https://github.com/lmdavalos/count2branches) |
| FBXO7–PINK1–PSMF1 variant modeling | Boltz-2 + PyRosetta ΔΔG | [`fbxo7_ternary_modeling/`](fbxo7_ternary_modeling/) |

## Repository layout

```
scripts/            Analysis walkthroughs (CAFE, compleasm, HyPhy)
fbxo7_ternary_modeling/
                    Structural modeling of bat FBXO7 variants (see its own README)
```

Each `scripts/*.md` documents the objective, inputs, pipeline steps, and outputs for one
analysis. Full scripts referenced in the walkthroughs are added to the repository alongside
their documentation.
