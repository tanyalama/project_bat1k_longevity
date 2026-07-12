# ==============================================================================
# EVE branch-shift significance — identify genes with shifted expression
# ------------------------------------------------------------------------------
# Downstream of the EVE branch-shift run (see run_EVE_6443_longevity.sh and
# convert_output_branch.sh in this folder). Takes the per-gene likelihood-ratio
# test (LRT) statistics from the two-theta-vs-one-theta branch-shift test on the
# Myotis myotis longevity lineage, converts them to p-values under the chi-square
# null (df = 1), applies a Bonferroni correction across all tested genes, and
# writes out the significant-gene list.
#
# This is the step that produced the significant longevity-associated gene set
# (6443 Myotis genes) reported in the manuscript.
#
# Input : a two-column, whitespace-delimited table with NO header:
#           column 1 = gene ID
#           column 2 = branch-shift LRT statistic
#         (the `result` file emitted by convert_output_branch.sh)
# Output: <prefix>_pval_matrix  — gene, raw p, Bonferroni-corrected p
#         <prefix>_sig          — subset with corrected p < 0.05
#
# Note: dataframes here were originally named *SOR (from an earlier Sorex
# araneus project). They have been renamed generically; the linked data is the
# Myotis myotis EVE result.
# ==============================================================================

library(readr)

# ---- Parameters --------------------------------------------------------------
INPUT   <- "result"                 # EVE branch-shift LRT output (gene, LRT)
PREFIX  <- "result_6443_longevity"  # output filename prefix
ALPHA   <- 0.05                     # significance threshold (corrected p)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) INPUT  <- args[1]
if (length(args) >= 2) PREFIX <- args[2]

# ---- Read LRT statistics -----------------------------------------------------
eve_result <- read_table(INPUT, col_names = FALSE)
colnames(eve_result)[1:2] <- c("gene", "LRT")

# ---- Convert LRT -> p-value under chi-square (df = 1), then Bonferroni --------
# Each branch-shift LRT statistic is compared to a chi-square distribution with
# 1 degree of freedom (one extra theta parameter under the two-theta model).
eve_result$p_raw <- pchisq(eve_result$LRT, df = 1, lower.tail = FALSE)
eve_result$p_bonferroni <- p.adjust(eve_result$p_raw, method = "bonferroni")

pval_matrix <- eve_result[, c("gene", "p_raw", "p_bonferroni")]
sig_genes   <- subset(pval_matrix, p_bonferroni < ALPHA)

# ---- Diagnostic: empirical CDF of the LRT statistics -------------------------
png(paste0(PREFIX, "_LRT_ecdf.png"), width = 700, height = 600)
plot(ecdf(eve_result$LRT),
     main = "EVE branch-shift LRT — empirical CDF",
     xlab = "LRT statistic", ylab = "Fn(x)")
dev.off()

# ---- Write results -----------------------------------------------------------
write.table(pval_matrix, paste0(PREFIX, "_pval_matrix"),
            na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(sig_genes, paste0(PREFIX, "_sig"),
            na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

cat(sprintf("Tested %d genes; %d significant at Bonferroni p < %.2f\n",
            nrow(pval_matrix), nrow(sig_genes), ALPHA))
