# Mitophagy pathway log2 fold-change visualization (hsa04137)

Visualizes RNA-seq log2 fold changes on the KEGG **Mitophagy – animal**
pathway (`hsa04137`) using [`pathview`](https://bioconductor.org/packages/pathview/),
with a viridis color scale (purple = downregulated, yellow = upregulated).

## Contents

```
mitophagy_pathway_lfc/
├── scripts/
│   ├── mitophagy_pathview.R   # renders the pathway diagram colored by log2FC
│   ├── sig_genes_table.R      # builds the gene-level significant-gene table
│   └── make_legend.py         # renders the standalone color-scale legend
├── data/
│   ├── fc_data.tsv            # input: Entrez gene ID + log2fc (tab-separated)
│   └── sig_genes_table.csv    # output: significant genes + symbol + HEX color
└── figures/
    ├── hsa04137.log2fc.png    # the colored pathway diagram
    └── legend.png             # standalone Arial legend
```

## Method

1. **Input** (`data/fc_data.tsv`): one row per gene, `<Entrez ID>\t<log2fc>`.
   Genes with `NA` values are treated as non-significant.
2. **Coloring** (`scripts/mitophagy_pathview.R`):
   - Viridis anchors — low `#440154` (purple) → mid `#21908C` (teal) →
     high `#FDE725` (yellow).
   - **Significance threshold**: genes with `|log2FC| < sig_thresh` (default 2)
     are grayed out (`gray80`), as are genes with no value.
   - **Color scale**: `scale_mode = "asymmetric"` uses the data's true min/max
     (here −7.08 to +13.08); `"symmetric"` uses `± limit_abs` with 0 at the
     teal midpoint.
   - A gray "not significant" swatch is stamped onto the legend as a
     post-processing step (pathview's native key has no NA slot).
3. **Legend** (`scripts/make_legend.py`): standalone key matching the diagram,
   all text in Arial ≥ 8 pt, titled "log fold change".

All tunable parameters live in the `CONFIG` block at the top of each script.

## Reproducing

**R (pathway diagram):**
```bash
Rscript scripts/mitophagy_pathview.R
```
Requires R with `pathview` and the Bioconductor annotation data packages
`org.Hs.eg.db` and `GenomeInfoDbData`. `pathview` downloads the KEGG KGML/PNG
for `hsa04137` at runtime (needs network access to rest.kegg.jp / kegg.jp).

**Python (legend):**
```bash
python scripts/make_legend.py
```
Requires `matplotlib` and numpy; uses the system Arial font if available.

## Notes

- The 5 genes with `NA` values and the genes with `|log2FC| < 2` render gray.
- 11 genes have `|log2FC| >= 2`; see `data/sig_genes_table.csv` (built by
  `scripts/sig_genes_table.R`) for their Entrez IDs, symbols, values, and the
  per-gene HEX codes. That script verifies every row's value against
  `fc_data.tsv` before writing.
- On the asymmetric scale, the single largest value (CALCOCO2, +13.08) and the
  next (NLRX1, +10.56) share the top yellow bins, while the +2 to +4 genes
  cluster near the teal midpoint (+3.0).

### Node coloring vs. gene identity (important)

KEGG collapses gene families into a single box. For example the box labeled
**RAB7A** maps two Entrez IDs, `{RAB7A=7879, RAB7B=338382}`, and the box labeled
**ATG9A** maps `{ATG9A=79065, ATG9B=285973}`. pathview draws each box with its
representative gene's *label* but colors it by the *aggregate* of all mapped
members. So a strongly-colored box can be driven by a non-representative family
member: the deep-purple "RAB7A" box is colored by **RAB7B (338382, −7.08)** —
RAB7A itself (7879) is +0.44 and not significant.

Consequently:
- `data/sig_genes_table.csv` is **gene-level** — it is authoritative for gene
  identity and value (verified against `fc_data.tsv`).
- The pathway PNG shows **node-level aggregated** coloring.
- When a colored box interests you, check the gene table for which family member
  actually carries the signal rather than reading the box label directly.
