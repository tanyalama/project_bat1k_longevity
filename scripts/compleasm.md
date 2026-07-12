# Annotation completeness (BUSCO) — compleasm

By Blair Bentley

## Objective

Generate BUSCO-style completeness scores with
[compleasm](https://github.com/huangnengCSU/compleasm) for the bat1k longevity dataset.
Inputs are **protein** FASTAs (longest transcript per gene), so scores reflect *annotation*
completeness rather than genome completeness.

## Pipeline

1. **Prepare inputs.** For each of 30 species, extract query protein sequences, strip `-`
   and `X` characters and blank lines, and reformat headers to include `gene=` for the
   OrthoFinder `primary_transcript.py` tool.
2. **Keep the longest transcript per protein** with OrthoFinder's `primary_transcript.py`.
3. **Run compleasm** (`compleasm protein -l mammalia -t 16`) as a SLURM array (one task per
   species) against the mammalia lineage database (`compleasm download mammalia`).
4. **Summarize in R.** Parse the per-species logs into a single table of Single (S),
   Duplicated (D), Fragmented (F), Missing (M), and Complete (C = S+D) scores.
5. **Plot** stacked completeness bars per species (`ggpubr`, viridis palette).

## Software

compleasm (conda), OrthoFinder (`primary_transcript.py`), R (`ggpubr`, `viridis`).

## Key outputs

- `all_species_compleasm_stats.txt` — S/D/F/M/C scores per species
- `all_species_compleasm.pdf` — stacked-bar completeness figure

## Acknowledgements

Thanks to Bill Thomas for suggesting compleasm.

## Scripts

Scripts live in [`compleasm/`](compleasm/):

- `compleasm_array.sh` — SLURM array running compleasm (one task per species)
- `compleasm_summarize_plot.R` — parse logs into a stats table and plot completeness bars
