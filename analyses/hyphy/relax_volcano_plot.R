# ==============================================================================
# RELAX selection volcano plot
# ------------------------------------------------------------------------------
# Visualizes the HyPhy RELAX screen: per-gene selection-intensity parameter K
# against its FDR-corrected p-value. K > 1 = intensified selection, K < 1 =
# relaxed selection. A curated set of genes is highlighted and labelled.
#
# Input: a tab-delimited table with columns  GeneID  K  pFDR
# Usage: edit INPUT_TSV below (or pass it as the first argument to Rscript).
# ==============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

INPUT_TSV <- "gene_k_values_v2.tsv"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) INPUT_TSV <- args[1]

# Genes to highlight and label.
intensified_genes <- c("FOXP3", "SLX4", "XPA", "BCCIP", "WDR45", "ARID1B", "TRIM65", "SMCHD1")
relaxed_genes     <- c("PIK3CD", "CLSPN", "UMOD", "ERCC1", "EME1", "COA8",
                       "HUWE1", "IL1RL1", "TNIP3", "ZMYM6", "KAT6A")

data <- read.table(INPUT_TSV, header = TRUE, sep = "\t") |>
  filter(K > 0) |>                       # log2 undefined for K <= 0
  mutate(
    log2K = log2(K),
    Selection = case_when(
      pFDR < 0.05 & K > 1 ~ "Intensified",
      pFDR < 0.05 & K < 1 ~ "Relaxed",
      TRUE                ~ "Not significant"
    ),
    LabelGroup = case_when(
      GeneID %in% intensified_genes ~ "Intensified",
      GeneID %in% relaxed_genes     ~ "Relaxed",
      TRUE                          ~ NA_character_
    ),
    Label = ifelse(!is.na(LabelGroup), GeneID, NA_character_)
  )

volcano <- ggplot(data, aes(x = log2K, y = pFDR)) +
  # Background: non-highlighted genes as faint hollow points.
  geom_point(
    data = subset(data, is.na(LabelGroup)),
    shape = 21, fill = "gray70", color = "gray70",
    size = 2, stroke = 0.5, alpha = 0.3
  ) +
  # Foreground: highlighted genes, colored by selection direction.
  geom_point(
    data = subset(data, !is.na(LabelGroup)),
    aes(color = LabelGroup), size = 3
  ) +
  geom_text_repel(
    data = subset(data, !is.na(LabelGroup)),
    aes(label = Label),
    size = 3, max.overlaps = Inf,
    segment.color = "black", color = "black"
  ) +
  scale_color_manual(values = c(
    "Relaxed"     = "#32648EFF",   # viridis blue
    "Intensified" = "#B8DE29FF"    # viridis green
  )) +
  geom_vline(xintercept = 0,    linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  scale_y_reverse(limits = c(1, -0.10)) +   # most-significant genes at top
  labs(x = expression(log[2]~K), y = "p (FDR-corrected)", color = "Selection")

ggsave("relax_volcano.png", plot = volcano, width = 7, height = 6, dpi = 300)
