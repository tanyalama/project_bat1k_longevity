## Load libraries
library(ggplot2)
library(dplyr)
library(viridis)      # for viridis color palette
library(scales)        # for log10 axis formatting if needed

## Load data
data <- read.table("relax_selection_summary.txt", header = TRUE, sep = "\t")
#data <- read.table("relax_selection_summary.txt", header = TRUE, sep = "\t")

## Log2 transform K, but only if K > 0
data <- data %>%
  filter(K > 0) %>%  # Avoid taking log of zero or negative
  mutate(
    log2K = log2(K),
    Selection = case_when(
      pFDR < 0.05 & K > 1 ~ "Intensified",
      pFDR < 0.05 & K < 1 ~ "Relaxed",
      TRUE ~ "Not significant"
    ),
    Selection = factor(Selection, levels = c("Relaxed", "Not significant", "Intensified"))
  )

## Pick your genes of interest
intensified_genes <- c("FOXP3", "SLX4", "XPA", "BCCIP", "WDR45", "ARID1B", "TRIM65", "SMCHD1")
relaxed_genes <- c("PIK3CD", "CLSPN", "UMOD", "ERCC1", "EME1", "COA8", "HUWE1", "IL1RL1", "TNIP3", "ZMYM6", "KAT6A")

## Create your LabelGroup
data$LabelGroup <- case_when(
  data$GeneID %in% intensified_genes ~ "Intensified",
  data$GeneID %in% relaxed_genes ~ "Relaxed",
  TRUE ~ NA_character_
)

data$Label <- ifelse(!is.na(data$LabelGroup), data$GeneID, NA)
  
## Pick your color palette 
show_col(viridis_pal()(20))
show(viridis_pal()(20))


## Final figure with labels
ggplot(data, aes(x = log2K, y = pFDR)) +
    # Gray hollow background points (non-highlighted)
    geom_point(
        data = subset(data, is.na(LabelGroup)), 
        shape = 21, fill = "gray70", color = "gray70", 
        size = 2, stroke = 0.5, alpha = 0.3
    ) +
    
    # Colored solid points (highlighted genes only)
    geom_point(
        data = subset(data, !is.na(LabelGroup)), 
        aes(color = LabelGroup),
        size = 3
    ) +
    
    # Black labels for highlighted genes
    geom_text_repel(
        data = subset(data, !is.na(LabelGroup)),
        aes(label = Label),
        size = 3,
        max.overlaps = Inf,
        segment.color = "black",
        color = "black"
    ) +
    
    # Custom viridis-inspired colors
    scale_color_manual(values = c(
        "Relaxed" = "#32648EFF",        # viridis blue
        "Intensified" = "#B8DE29FF"     # viridis bright green
    )) +
    
    # Reference lines
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    
    # Layout and theme
    theme_minimal(base_size = 14) +
    coord_cartesian(ylim = c(1, -0.10)) +
    scale_y_reverse(limits = c(1, -0.10)) 
  ```
  
We also included the entrez summary in Supplementary Table 12 for each gene, and a note regarding whether each gene is associated with aging hallmarks or a priori lists (e.g., antagonistic pleiotropies of aging).

You can now move foward with investigations of invidividual genes

## Part 9: alphaFold & alphaMissense
### 1. Visualization with pymissense
Already installed in /bin
Usage is:
