#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Build the gene-level significance table for the hsa04137 log2FC figure.
#
# One row per GENE (not per KEGG node) with |log2FC| >= sig_thresh, giving:
#   entrez, geneSymbol, log2fc, hex
# where `hex` is the viridis color that gene's own value maps to on the scale.
#
# IMPORTANT — node vs gene: KEGG collapses gene families into a single box
# (e.g. the "RAB7A" node maps {RAB7A=7879, RAB7B=338382}). pathview labels the
# box with a representative gene but COLORS it by the aggregate of all mapped
# members, so a box's displayed color can differ from an individual member's
# per-gene hex below. This table is gene-level and authoritative for gene
# identity + value; see the pathway PNG for the aggregated node coloring.
# ---------------------------------------------------------------------------
# Usage:  Rscript sig_genes_table.R
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))
suppressMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })

# ----------------------------- CONFIG --------------------------------------
data_file  <- "fc_data.tsv"
sig_thresh <- 2
vmin       <- NULL   # NULL -> data min; else numeric
vmax       <- NULL   # NULL -> data max
n_bins     <- 20
# 20-bin viridis ramp (purple -> teal -> yellow), matching the pathway figure.
ramp <- strsplit(paste0(
  "#440154,#40115A,#3C2160,#383167,#34416D,#315073,#2D6079,#297080,#258086,",
  "#21908C,#21908C,#399A81,#52A375,#6AAD6A,#83B75E,#9BC053,#B4CA47,#CCD43C,",
  "#E5DD30,#FDE725"), ",")[[1]]
# ---------------------------------------------------------------------------

df <- read.delim(data_file, stringsAsFactors = FALSE, na.strings = c("NA","","NaN"))
entrez <- as.character(df[[1]])
lfc    <- as.numeric(df[[2]])

if (is.null(vmin)) vmin <- min(lfc, na.rm = TRUE)
if (is.null(vmax)) vmax <- max(lfc, na.rm = TRUE)

# Per-gene color: clamp to [vmin,vmax], cut into n_bins, index the ramp.
# Matches pathview's node.color binning (seq of n_bins+1 cuts, include.lowest).
hex_for <- function(v) {
  if (is.na(v)) return(NA_character_)
  vv <- max(min(v, vmax), vmin)
  cuts <- seq(vmin, vmax, length.out = n_bins + 1)
  idx <- as.integer(cut(vv, cuts, include.lowest = TRUE))
  toupper(ramp[idx])
}

keep <- !is.na(lfc) & abs(lfc) >= sig_thresh
sym  <- mapIds(org.Hs.eg.db, keys = entrez[keep], column = "SYMBOL",
               keytype = "ENTREZID", multiVals = "first")

tab <- data.frame(
  entrez     = entrez[keep],
  geneSymbol = unname(sym),
  log2fc     = lfc[keep],
  hex        = vapply(lfc[keep], hex_for, character(1)),
  stringsAsFactors = FALSE
)
tab <- tab[order(-tab$log2fc), ]

# --- verification: every row must match fc_data.tsv exactly ---------------
truth <- setNames(lfc, entrez)
stopifnot(all(abs(tab$log2fc - truth[tab$entrez]) < 1e-9))
stopifnot(all(abs(tab$log2fc) >= sig_thresh))

write.csv(tab, "sig_genes_table.csv", row.names = FALSE)
cat(sprintf("Wrote sig_genes_table.csv: %d significant genes (|log2FC| >= %g); verified against %s\n",
            nrow(tab), sig_thresh, data_file))
print(tab, row.names = FALSE)
