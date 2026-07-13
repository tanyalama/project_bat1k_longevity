#!/usr/bin/env Rscript
# ==============================================================================
# EVE TPM dot plots — per-gene expression across species
# ------------------------------------------------------------------------------
# Companion visualization to the EVE branch-shift analysis (analyses/eve/).
# Plots per-sample transcripts-per-million (TPM) for any gene of interest,
# read directly from the comparative TPM matrix, with Myotis (Bat) highlighted.
#
# The TPM matrix (all_species_TPM_byGeneID_clean.txt) has one row per gene:
#   geneID  <sample_1>  <sample_2>  ...  <sample_N>
# Species of each sample is parsed from the sample-name prefix (bosTau, canFam,
# ...); the three fibroblast replicates (fibro_MMY63*) are the Myotis samples.
#
# Usage:
#   Rscript eve_tpm_dotplot.R <TPM_matrix> <GENE> [GENE2 GENE3 ...]
#
# Example:
#   Rscript eve_tpm_dotplot.R all_species_TPM_byGeneID_clean.txt SLC9A1
#   Rscript eve_tpm_dotplot.R all_species_TPM_byGeneID_clean.txt CDK5 SLC9A1 GABARAP
#
# One PNG is written per gene:  eve_tpm_dotplot_<GENE>.png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(viridis)
})

# ---- arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript eve_tpm_dotplot.R <TPM_matrix> <GENE> [GENE2 ...]")
}
matrix_path <- args[1]
genes       <- args[-1]

# ---- species order + palette (Bat first) -------------------------------------
species_order  <- c("Bat", "Cat", "Dog", "Horse", "Pig", "Cow", "Mouse", "Human")
species_colors <- c("#2FAC7D", viridis(7, option = "D", begin = 0, end = 0.55))
names(species_colors) <- species_order

# Map a sample-name prefix to its species label.
label_species <- function(sample) {
  dplyr::case_when(
    grepl("bosTau", sample) ~ "Cow",
    grepl("canFam", sample) ~ "Dog",
    grepl("equCab", sample) ~ "Horse",
    grepl("felCat", sample) ~ "Cat",
    grepl("homSap", sample) ~ "Human",
    grepl("musMus", sample) ~ "Mouse",
    grepl("MMY63",  sample) ~ "Bat",
    grepl("susScr", sample) ~ "Pig",
    TRUE ~ "Other"
  )
}

# ---- read matrix -------------------------------------------------------------
mat <- read.delim(matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- names(mat)[1]                       # first column holds gene IDs
sample_cols <- names(mat)[-1]

# ---- build the per-sample dot plot for one gene ------------------------------
tpm_dotplot_for_gene <- function(gene) {
  hit <- mat[mat[[gene_col]] == gene, , drop = FALSE]
  if (nrow(hit) == 0) {
    warning(sprintf("Gene '%s' not found in %s — skipping.", gene, matrix_path))
    return(invisible(NULL))
  }
  if (nrow(hit) > 1) {
    warning(sprintf("Gene '%s' has %d rows; using the first.", gene, nrow(hit)))
    hit <- hit[1, , drop = FALSE]
  }

  # named TPM vector (sample -> TPM) straight from the matrix row
  tpm <- as.numeric(hit[1, sample_cols])
  names(tpm) <- sample_cols

  df <- enframe(tpm, name = "sample", value = "TPM") |>
    mutate(species = factor(label_species(sample), levels = species_order))

  p <- ggplot(df, aes(x = species, y = TPM, color = species)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
    scale_color_manual(values = species_colors) +
    theme_classic(base_size = 12) +
    labs(x = "Species", y = "Transcripts per Million (TPM)", title = gene) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  out <- sprintf("eve_tpm_dotplot_%s.png", gene)
  ggsave(out, plot = p, width = 6, height = 5, dpi = 300)
  message(sprintf("Wrote %s  (Bat mean TPM = %.1f, n = %d)",
                  out,
                  mean(df$TPM[df$species == "Bat"]),
                  sum(df$species == "Bat")))
  invisible(p)
}

# ---- render every requested gene ---------------------------------------------
for (g in genes) tpm_dotplot_for_gene(g)
