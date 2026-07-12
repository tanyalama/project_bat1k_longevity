# Gene family expansion and contraction — OrthoFinder + CAFE

By Tanya Lama & Nicole Paulat

## Objective

Analyze gene family (protein family) expansion and contraction across 24 species using
CAFE5, starting from TOGA1 annotations.

## Pipeline

1. **Extract query proteins.** From each species' TOGA `prot.fasta`, keep only query
   sequences (drop the human reference) → `query_{species}.pep`.
   *([`grep_query.py`](cafe/grep_query.py))*
2. **Drop one-to-one orthologs.** Using each species' `orthology_classification.tsv`,
   pad gene IDs to 5 digits, list `one2one` transcripts, and remove them from the protein
   FASTAs → `{species}_final.pep`.
   *([`pad_names.py`](cafe/pad_names.py), [`get_one2ones.sh`](cafe/get_one2ones.sh), [`parse_paralogs.py`](cafe/parse_paralogs.py))*
3. **Keep the longest transcript per gene**, then re-add gene symbols to headers from an
   Ensembl v104 Biomart transcript→symbol table → `{species}_longest_transcript_annotated.faa`.
   *([`get_longest_transcripts.py`](cafe/get_longest_transcripts.py), [`get_faas_gene_names.py`](cafe/get_faas_gene_names.py))*
4. **Run OrthoFinder** on the annotated longest-transcript FASTAs to define orthogroups.
   *([`orthofinder.sh`](cafe/orthofinder.sh) — `orthofinder -a 4 -f input_faas -o results`)*
5. **Curate orthogroups for CAFE.** Deduplicate all-vs-all similarity pairs, merge
   orthogroups >80% similar (via connected components), then drop orthogroups that are
   very large (>450 genes), highly variable (>100 gene range across species), or
   predominantly olfactory receptors.
   *([`3a_clean_pairs.py`](cafe/3a_clean_pairs.py), [`calculate_minimum_mergers.py`](cafe/calculate_minimum_mergers.py), [`perform_min_mergers.py`](cafe/perform_min_mergers.py),
   [`6a_filter_orthogroups.py`](cafe/6a_filter_orthogroups.py), [`remove_olfactory_orthogroups.py`](cafe/remove_olfactory_orthogroups.py))*
6. **Convert to CAFE input** (orthogroup × per-species gene counts, tab-separated).
   *([`7_convert_into_cafe_format.py`](cafe/7_convert_into_cafe_format.py))*
7. **Run CAFE5.** Estimate the error model (`-e`), then fit models with 1–3 lambdas on an
   ultrametric tree (topology matching `supermatrix_iqtree.contree`, minus hg38). Select
   the best model by comparing final likelihoods across runs.
   *([`sbatch_cafe1.sh`](cafe/sbatch_cafe1.sh))*
8. **Visualize** expansions/contractions with CafePlotter.

## Filtering rules for orthogroups

- Drop orthogroups with >450 genes total.
- Drop orthogroups with >100 gene-count variation between the highest and lowest species.
- Merge orthogroups sharing >80% similarity.
- Drop olfactory-receptor-dominated orthogroups.

## Software

OrthoFinder, CAFE5 (`bioconda::cafe`), CafePlotter, Biopython, pandas. Run on a SLURM
cluster (`module load conda/latest`).

## Key outputs

- `cafe_input_w_or.tsv` / `cafe_input_no_or.tsv` — CAFE count matrices (with / without
  olfactory receptors)
- `Base_results.txt`, `Base_error_model.txt` — per-run CAFE model fit
- CafePlotter trees annotated with significant expansions/contractions

## References

- [CAFE5](https://github.com/hahnlab/CAFE5) and its
  [tutorial](https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md)

## Scripts

All scripts referenced above live in [`cafe/`](cafe/):

- **Protein extraction / filtering** — `grep_query.py`, `pad_names.py`, `get_one2ones.sh`,
  `parse_paralogs.py`, `get_longest_transcripts.py`, `get_faas_gene_names.py`
- **OrthoFinder** — `orthofinder.sh`
- **Orthogroup curation** — `3a_clean_pairs.py`, `grep_cleaner_pairs.sh`,
  `calculate_minimum_mergers.py`, `perform_min_mergers.py`, `6a_filter_orthogroups.py`,
  `remove_olfactory_orthogroups.py`
- **CAFE input + run** — `7_convert_into_cafe_format.py`, `sbatch_cafe1.sh`
