#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Visualize log2 fold change on KEGG pathway hsa04137 (Mitophagy - animal)
# using pathview, with a viridis color scale
#   purple (#440154) = downregulated  ->  yellow (#FDE725) = upregulated
# ---------------------------------------------------------------------------
# Usage:  Rscript mitophagy_pathview.R
# Edit the CONFIG block below to change data file, pathway, palette, or limits.

# Pick up the workspace-local library (holds GenomeInfoDbData, which the
# conda pathview build does not bundle).
if (dir.exists("rlib")) .libPaths(c(normalizePath("rlib"), .libPaths()))

suppressMessages(library(pathview))

# ----------------------------- CONFIG --------------------------------------
data_file   <- "fc_data.tsv"   # tab-separated: column 1 = Entrez gene ID, column 2 = log2fc
pathway_id  <- "04137"          # KEGG pathway number (hsa04137 = Mitophagy - animal)
species     <- "hsa"            # KEGG species code (human)
out_suffix  <- "log2fc"         # tag added to output filenames

# Viridis anchor colors: low -> mid -> high
col_low     <- "#440154"        # viridis purple  (downregulated)
col_mid     <- "#21908C"        # viridis teal    (~zero)
col_high    <- "#FDE725"        # viridis yellow  (upregulated)

# Color-scale limits. Controls how log2fc values map onto the viridis gradient.
#   "asymmetric" -> use the data's true min/max as the two ends (nothing clamps;
#                   the teal midpoint lands at the CENTER of the range, not at 0).
#   "symmetric"  -> +/- limit_abs, so log2fc = 0 maps to the mid (teal) color.
scale_mode  <- "asymmetric"     # "asymmetric" or "symmetric"
limit_abs   <- 3                # only used when scale_mode == "symmetric"
n_bins      <- 20               # number of discrete color steps in the gradient

# Significance threshold: genes with |log2fc| below this are considered
# not significant and are grayed out (set to 0 to color everything).
sig_thresh  <- 2
na_color    <- "gray80"         # color for non-significant / missing genes
# ---------------------------------------------------------------------------

# Load data ------------------------------------------------------------------
df <- read.delim(data_file, header = TRUE, stringsAsFactors = FALSE,
                 na.strings = c("NA", "", "NaN"))
fc <- as.numeric(df[[2]])
names(fc) <- as.character(df[[1]])
message(sprintf("Loaded %d genes (%d with values, %d NA)",
                length(fc), sum(!is.na(fc)), sum(is.na(fc))))

# Gray out non-significant genes: |log2fc| below threshold -> NA (na_color).
n_gray_thresh <- sum(abs(fc) < sig_thresh, na.rm = TRUE)
fc[abs(fc) < sig_thresh] <- NA
message(sprintf("Grayed out %d genes with |log2fc| < %.2f (plus %d missing); %d remain colored",
                n_gray_thresh, sig_thresh, sum(is.na(df[[2]])), sum(!is.na(fc))))

# Resolve color-scale limit --------------------------------------------------
# raw_vals = all numeric log2fc (before graying) so the scale spans the true range.
raw_vals <- as.numeric(df[[2]])
if (scale_mode == "asymmetric") {
  gene_limit <- c(min(raw_vals, na.rm = TRUE), max(raw_vals, na.rm = TRUE))
  message(sprintf("Color scale (asymmetric): [%.2f, %.2f]; teal midpoint at %.2f",
                  gene_limit[1], gene_limit[2], mean(gene_limit)))
} else {
  if (is.null(limit_abs)) limit_abs <- max(abs(raw_vals), na.rm = TRUE)
  gene_limit <- limit_abs
  message(sprintf("Color scale (symmetric): +/- %.2f (0 = teal; outside clamped)", limit_abs))
}

# Render ---------------------------------------------------------------------
pv <- pathview(
  gene.data    = fc,
  pathway.id   = pathway_id,
  species      = species,
  gene.idtype  = "entrez",
  out.suffix   = out_suffix,
  kegg.native  = TRUE,          # render on the native KEGG PNG diagram
  same.layer   = FALSE,         # draw gene labels on a separate layer (cleaner)
  low          = list(gene = col_low,  cpd = col_low),
  mid          = list(gene = col_mid,  cpd = col_mid),
  high         = list(gene = col_high, cpd = col_high),
  limit        = list(gene = gene_limit, cpd = gene_limit),
  bins         = list(gene = n_bins,    cpd = n_bins),
  na.col       = na_color
)

# Add a "not significant" gray swatch to the legend ---------------------------
# pathview's color key only draws the gradient, so we stamp a gray box + label
# beneath it to explain the grayed-out (non-significant / missing) genes.
out_png <- sprintf("%s%s.%s.png", species, pathway_id, out_suffix)
suppressMessages(library(png))
img <- readPNG(out_png)
h <- dim(img)[1]; w <- dim(img)[2]

png(out_png, width = w, height = h)
op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot(NA, xlim = c(0, w), ylim = c(0, h), type = "n", axes = FALSE, xlab = "", ylab = "")
rasterImage(img, 0, 0, w, h)

# Swatch geometry (top-right, just below the gradient key). Tweak if the key moves.
sw_x0 <- w * 0.695; sw_x1 <- sw_x0 + (w * 0.028)   # gray box
sw_y1 <- h * 0.905; sw_y0 <- sw_y1 - (h * 0.020)   # y measured from bottom
rect(sw_x0, sw_y0, sw_x1, sw_y1, col = na_color, border = "black", lwd = 1)
text(sw_x1 + w * 0.006, (sw_y0 + sw_y1) / 2,
     sprintf("Not significant (|log2FC| < %g)", sig_thresh),
     adj = c(0, 0.5), cex = 1.1)
par(op)
invisible(dev.off())

message("Done. Output written to: ", out_png)
