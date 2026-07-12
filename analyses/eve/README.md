# Comparative transcriptomics — Expression Variance and Evolution (EVE)

Lineage-specific shifts in gene expression tested with the **EVE** model
(Expression Variance and Evolution; [Rohlfs & Nielsen 2015, *Syst. Biol.* 64:695–708](https://doi.org/10.1093/sysbio/syv042)).

## Objective

EVE models gene expression as a quantitative trait in a phylogenetic comparative
framework. It estimates β, the ratio of within-species (intraspecific) to between-species
(interspecific) expression variance, and tests whether a gene's expression along a focal
branch is consistent with stabilizing selection (a single optimum across the tree) or with
a branch-specific shift in optimum (diversifying / directional selection). Here the focal
branch is *Myotis myotis*, and the test asks which of 6,443 single-copy orthologs show an
evolutionarily adaptive expression shift unique to that lineage.

## Inputs (`data/`)

| File | Description |
|---|---|
| `6443_longevity_tpm.txt` | Expression matrix: header line = gene count (6443); each row = gene symbol + TPM values across all samples of all species (28 columns) |
| `renamed_longevity_tree.nh` | Species phylogeny, 8 taxa, Newick with branch lengths (numeric leaf labels 1–8; leaf 8 = *Myotis myotis*) |
| `longevity_num_reorder.txt` | Replicate counts per species in tree-leaf order (`4 3 3 4 3 4 4 3`, summing to the 28 expression columns) |
| `10_genelist/longevity_genelist.txt` | 6,443 gene symbols, in the same row order as the TPM matrix — used to re-label the numeric EVE output |

Single-copy orthologs were identified with OrthoFinder; per-species TPM matrices were merged
on the OrthoGroup identifier (see the manuscript Methods and `../cafe/` for the ortholog
pipeline). The molecular-clock mammalian phylogeny
([Álvarez-Carretero et al. 2022](https://doi.org/10.1038/s41586-021-04341-1)) was pruned to
the eight species used here.

## Pipeline

1. **Run EVE** — [`run_EVE_6443_longevity.sh`](scripts/run_EVE_6443_longevity.sh) calls the `EVEmodel` binary in branch-test mode:
   ```
   EVEmodel -O -o 7 -n 6443 \
     -t data/renamed_longevity_tree.nh \
     -i data/longevity_num_reorder.txt \
     -d data/6443_longevity_tpm.txt \
     -f _6443_longevity -v 10
   ```
   `-O` runs the one-θ vs. two-θ (branch-shift) test; `-o 7` places the second optimum so
   that a shift is tested on species 8 (*Myotis myotis*); `-n` is the gene count. For each
   gene EVE fits a single-optimum OU model (H₀, stabilizing selection) and a two-optimum
   model with a separate optimum on the focal branch (H₁, diversifying selection), and
   reports the likelihood-ratio test statistic.

2. **Convert branch-shift output** — [`convert_output_branch.sh`](scripts/convert_output_branch.sh) reshapes the raw LRT file
   (`BSThetaTestLRTs_6443_longevity.res`) into a per-gene table: it splits the
   space-delimited LRTs onto one line per gene, converts out of scientific notation, and
   pastes the values alongside the gene-symbol list, producing `result` (gene ⇥ LRT).
   ```
   bash convert_output_branch.sh \
     -f BSThetaTestLRTs_6443_longevity.res -s 8 -g longevity_genelist.txt -o result
   ```

3. **Significance & candidate genes** — [`eve_branchshift_significance.R`](scripts/eve_branchshift_significance.R)
   converts each gene's branch-shift LRT statistic to a p-value under the chi-square null
   (df = 1), applies a Bonferroni correction (p_adj < 0.05), and writes the significant-gene
   list (`result_6443_longevity_sig`) — the set reported in the manuscript. Downstream: a
   dropout test (removing bat expression data) validates lineage specificity, and candidate
   gene sets are annotated for KEGG / GO pathways with DAVID.

## Outputs (`results/`)

| File | Description |
|---|---|
| `oneThetaMLs_6443_longevity.res` | Per-gene log-likelihood of the single-optimum (H₀) model |
| `twoBSThetaMLs_6443_longevity.res` | Per-gene log-likelihood of the two-optimum branch-shift (H₁) model |
| `oneThetaMLparams_6443_longevity.res` | Maximum-likelihood parameter estimates, one-θ model |
| `twoBSThetaMLparams_6443_longevity.res` | Maximum-likelihood parameter estimates, two-θ branch-shift model |
| `BSThetaTestLRTs_6443_longevity.res` | Per-gene likelihood-ratio test statistics (H₀ vs. H₁) |
| `result` | Converted per-gene table: gene symbol ⇥ branch-shift LRT |

## Visualization

- [`eve_branchshift_significance.R`](scripts/eve_branchshift_significance.R) - converts the
  branch-shift LRT statistics to Bonferroni-corrected p-values and extracts the
  significant-gene list, with an empirical-CDF diagnostic of the LRT distribution.
- [`eve_tpm_dotplot.R`](scripts/eve_tpm_dotplot.R) - per-gene TPM dot plots across the
  comparative species panel (*Myotis* fibroblasts highlighted), for genes of interest
  emerging from the EVE screen.

## Software

- **EVEmodel** — the EVE release binary (`EVE_release/EVEmodel`);
  method: [Rohlfs & Nielsen 2015](https://doi.org/10.1093/sysbio/syv042),
  reference implementation [github.com/bnielsen1995/EVEmodel](https://github.com/bnielsen1995/EVEmodel).
- Standard Unix text tools (`tr`, `awk`, `paste`) for output conversion.
- Run on a SLURM cluster (`-p cpu`, 40 tasks/node; the 6,443-gene one-θ-vs-two-θ analysis
  completed in ~10,700 s).

## Contributors

Comparative transcriptomics (EVE): Blair Bentley, Tanya Lama, William Thomas.
