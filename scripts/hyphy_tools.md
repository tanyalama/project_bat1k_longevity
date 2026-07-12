# Signatures of selection — HyPhy (aBSREL, MEME, RELAX)

By Nikki Paulat, Tanya Lama, Bill Thomas, Graham Hughes, and others, with guidance from
Sergei Kosakovsky Pond.

## Objective

Test for signatures of selection across the 23-bat + 7-non-bat species alignment for the
bat1k longevity project:

- **aBSREL** — branch-site test for episodic selection anywhere in the alignment / tree
  (exploratory scan of all genes).
- **MEME** — site-level episodic selection, run on genes flagged by aBSREL.
- **RELAX** — test for intensification vs. relaxation of selection.

All analyses use `--srv Yes` (site-to-site synonymous rate variation).

## Pipeline

1. **Prepare alignments.** Standardize extensions to `.fas`, drop empty alignments and
   those with <5 species, clean headers to species abbreviations (`REFERENCE` → `hg38`),
   and replace unsupported characters (`X`→`N`, MACSE frameshift `!`→`n`).
2. **Build the species tree.** Concatenate gene alignments with FASconCAT-G, then run
   IQ-TREE (partitioned, `-B 1000`) to select models and estimate the phylogeny. Reroot and
   drop bootstraps → `supermatrix_partition_iqtree_reroot_no_bootstraps.tree`.
   *([`fasconcat.sh`](hyphy/fasconcat.sh), [`slurm_iqtree.sh`](hyphy/slurm_iqtree.sh))*
3. **Clip trees per gene.** Because not every gene is present in every species, prune the
   species tree to each gene's taxa → one `.tfl` tree per gene.
   *([`modclip.py`](hyphy/modclip.py), [`cliptrees.R`](hyphy/cliptrees.R))*
4. **Tag foreground branches** with `{FG}` for branch-specific tests.
5. **Run aBSREL** as SLURM array jobs, batched at ≤1500 genes each (`--branches FG`).
   *([`template_slurm_absrel_array.sh`](hyphy/template_slurm_absrel_array.sh))*
6. **Run MEME** on genes with significant aBSREL results.
   *([`slurm_meme_array.sh`](hyphy/slurm_meme_array.sh))*
7. **Run RELAX** with test (`FG1`) vs. reference (`FG2`) branch sets.
   *([`slurm_relax_array1.sh`](hyphy/slurm_relax_array1.sh))*
8. **Parse results** — extract gene ID and p-value per test; collect genes with positive
   selection signal.

## Software

- HyPhy ([GitHub](https://github.com/veg/hyphy)), IQ-TREE, FASconCAT-G (v1.06.1)
- R (`ape`, `phangorn`) for tree clipping

## Notes

- aBSREL accepts rooted or unrooted trees (it unroots internally) and performs within-gene
  multiple-test correction.
- Foreground tags are ignored unless the job passes `--branches <tag>`.
- Batch array jobs are capped to stay under the cluster's queue limit (max 2000 jobs);
  the final batch's array length is the remainder.

## Key outputs

- `*_absrel_output.txt`, `*_MEME_output.txt`, `*_RELAX_output.txt` — per-gene results
- `absrel_positive_results.txt` — genes with significant episodic selection

## Scripts

Scripts live in [`hyphy/`](hyphy/):

- **Tree building** — `fasconcat.sh`, `slurm_iqtree.sh`
- **Per-gene tree clipping** — `modclip.py`, `cliptrees.R`
- **Selection tests (SLURM arrays)** — `template_slurm_absrel_array.sh`,
  `slurm_batch_hyphy.sh`, `slurm_meme_array.sh`, `slurm_relax_array1.sh`
- **Results parsing / plotting** — `padj.py`, `add_selection_antPal2.sh`,
  `relax_volcano_plot.R`
