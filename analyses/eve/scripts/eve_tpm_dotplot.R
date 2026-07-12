# ==============================================================================
# EVE TPM dot plots — per-gene expression across species
# ------------------------------------------------------------------------------
# Companion visualization to the EVE branch-shift analysis (analyses/eve/).
# Plots per-sample transcripts-per-million (TPM) for a gene of interest across
# the comparative species panel, with Myotis (Bat) highlighted first.
#
# Each gene is a named numeric vector (sample -> TPM). The species of each
# sample is parsed from the sample-name prefix (bosTau, canFam, ...); the three
# fibroblast replicates (fibro_MMY63*) are the Myotis samples.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tibble)
library(viridis)

# Fixed species display order (Bat first) and matching palette.
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

# Build the per-sample dot plot for one gene's named TPM vector.
tpm_dotplot <- function(tpm, gene_label) {
  df <- enframe(tpm, name = "sample", value = "TPM") |>
    mutate(species = factor(label_species(sample), levels = species_order))

  ggplot(df, aes(x = species, y = TPM, color = species)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
    scale_color_manual(values = species_colors) +
    theme_classic(base_size = 12) +
    labs(x = "Species", y = "Transcripts per Million (TPM)", title = gene_label) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ------------------------------------------------------------------------------
# Gene 1 (unlabeled in the source notebook — rename to the intended gene symbol)
# ------------------------------------------------------------------------------
gene1_tpm <- c(
  bosTau_SRR17288213 = 25.258, bosTau_SRR17288214 = 20.933, bosTau_SRR17288215 = 27.774, bosTau_SRR17288216 = 22.284,
  canFam_SRR17870687 = 28.606, canFam_SRR17870688 = 25.647, canFam_SRR17870691 = 26.118, canFam_SRR17870692 = 24.095,
  equCab_SRR18906497 = 20.459, equCab_SRR18906498 = 20.324, equCab_SRR18906499 = 20.362, equCab_SRR18906500 = 12.048,
  felCat_SRR24448362 = 25.312, felCat_SRR24448363 = 27.29,  felCat_SRR24448364 = 25.026,
  homSap_SRR23630190 = 32.185, homSap_SRR23630193 = 43.391, homSap_SRR23630200 = 52.779, homSap_SRR23630208 = 38.282,
  musMus_SRR24530105 = 29.1,   musMus_SRR24530106 = 36.764, musMus_SRR24530107 = 32.464,
  fibro_MMY6321_1.1_S16 = 151.96, fibro_MMY6321_1.2_S17 = 150.072, fibro_MMY6321_1.3_S18 = 148.618,
  susScr_SRR25461976 = 47.474, susScr_SRR25461977 = 40.391, susScr_SRR25461978 = 38.291
)

# ------------------------------------------------------------------------------
# CDK5
# ------------------------------------------------------------------------------
cdk5_tpm <- c(
  bosTau_SRR17288213 = 19.499, bosTau_SRR17288214 = 20.894, bosTau_SRR17288215 = 20.548, bosTau_SRR17288216 = 19.006,
  canFam_SRR17870687 = 25.319, canFam_SRR17870688 = 24.961, canFam_SRR17870691 = 23.778, canFam_SRR17870692 = 20.782,
  equCab_SRR18906497 = 34.409, equCab_SRR18906498 = 27.833, equCab_SRR18906499 = 30.726, equCab_SRR18906500 = 15.845,
  felCat_SRR24448362 = 10.924, felCat_SRR24448363 = 11.698, felCat_SRR24448364 = 9.425,
  homSap_SRR23630190 = 33.152, homSap_SRR23630193 = 35.053, homSap_SRR23630200 = 47.45,  homSap_SRR23630208 = 29.14,
  musMus_SRR24530105 = 23.195, musMus_SRR24530106 = 32.974, musMus_SRR24530107 = 26.913,
  fibro_MMY6321_1.1_S16 = 145.962, fibro_MMY6321_1.2_S17 = 136.643, fibro_MMY6321_1.3_S18 = 145.209,
  susScr_SRR25461976 = 25.418, susScr_SRR25461977 = 21.218, susScr_SRR25461978 = 26.89
)

# Render and save both panels.
ggsave("eve_tpm_dotplot_gene1.png", plot = tpm_dotplot(gene1_tpm, "Gene 1"),
       width = 6, height = 5, dpi = 300)
ggsave("eve_tpm_dotplot_CDK5.png",  plot = tpm_dotplot(cdk5_tpm, "CDK5"),
       width = 6, height = 5, dpi = 300)
