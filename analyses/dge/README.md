# Mitophagy pathway log2 fold-change visualization (hsa04137)

Overlays *Myotis myotis* vs *Mus musculus* fibroblast RNA-seq log2 fold changes
on the KEGG **Mitophagy – animal** pathway (`hsa04137`) with a **viridis**
colour scale: purple = downregulated, teal = 0, green → yellow = upregulated.

Version 2 (2026-10). Replaces the pathview/R rendering (v1) with a
dependency-light Python script that reproduces pathview's native-KEGG
rendering step without Bioconductor annotation packages.

## Contents

```
dge/
├── scripts/
│   ├── mitophagy_kegg_overlay.py   # v2: fetches KGML + PNG from KEGG REST, recolours gene boxes,
│   │                               #     writes figure, legend and tables (run this)
│   ├── mitophagy_pathview.R        # v1 (superseded): pathview rendering, asymmetric scale
│   ├── sig_genes_table.R           # v1 (superseded): gene table for the pathview figure
│   └── make_legend.py              # v1 (superseded): standalone legend for the pathview figure
├── data/
│   ├── fc_data.tsv                 # original input: Entrez gene ID + log2fc (87 genes)
│   ├── fc_data_padj.tsv            # input: entrez, symbol, log2fc, padj (padj from Supp Table 12)
│   ├── sig_genes_table.csv         # output: 11 significant genes + symbol + HEX colour
│   └── node_table.csv              # output: one row per KEGG gene box (77) — members, driver, colour
├── kegg/                           # cached KEGG downloads (hsa04137.xml, hsa04137.png, hsa_genes.tsv)
└── figures/
    ├── hsa04137.log2fc.viridis.pdf # the coloured pathway (vector boxes/text over the KEGG raster)
    ├── hsa04137.log2fc.viridis.png # same, 300 dpi
    ├── legend.pdf / legend.png     # standalone colour key (Arial)
    └── hsa04137.log2fc.png         # v1 figure (superseded; asymmetric scale, teal midpoint at +3)
```

## Method

1. **Input** (`data/fc_data_padj.tsv`): `entrez`, `symbol`, `log2fc`, `padj` for the
   87 genes in `fc_data.tsv`. `log2fc` is the DESeq2 log2FoldChange (identical to
   `fc_data.tsv`); `padj` comes from manuscript **Supplementary Table 12** (M. myotis vs
   house mouse fibroblasts). Table 12 lists only genes with |log2FC| ≥ 2, so padj is NA
   for genes below the LFC threshold; those fail the significance test regardless. The
   11 genes present in both files have identical log2FC values (checked to 1e-4).
2. **Significance**: a gene is significant only if **|log2FC| ≥ 2 and padj < 0.05**.
   The script stops if a gene passes the LFC filter but has no padj.
3. **KEGG map**: the `hsa04137` KGML and base PNG are fetched from `rest.kegg.jp`
   (cached in `kegg/`). There are 77 gene boxes covering 105 Entrez IDs.
4. **Box states** (KEGG collapses gene families into one box):
   - *significant*: at least one member is significant. The box is coloured by the
     significant member with the largest |log2FC| (the *driver*) and labelled with
     that gene's symbol. No box here has more than one significant member.
   - *measured, not significant*: grey `#CCCCCC`.
   - *not measured*: white. No member is in the dataset, or every member's log2FC is NA.
5. **Colour scale**: viridis, split at the threshold. Nothing between −2 and +2 gets a
   colour; that interval is a grey band in the key.
   - Down, −7 → −2: viridis 0.00–0.20 (deep purple → violet).
   - Up, +2 → +7: viridis 0.62–1.00 (green → yellow).
   - Values beyond ±7 are clamped: CALCOCO2 (+13.08) and NLRX1 (+10.56) show yellow,
     RAB7B (−7.08) shows purple. The key marks the ends `≤ −7` / `≥ 7`.
6. **Typography**: all overlaid text is Arial, embedded as TrueType so the PDF is
   editable. KEGG's own pathway text is part of the raster.

All tunable parameters (`LFC_THRESH`, `PADJ_THRESH`, `LIMIT`, `DOWN_RANGE`/`UP_RANGE`, dataset caption) live
in the `CONFIG` block at the top of `scripts/mitophagy_kegg_overlay.py`.

## Reproducing

```bash
python scripts/mitophagy_kegg_overlay.py
```

Requires Python 3 with `matplotlib`, `numpy`, `pandas`, `pillow`, and network
access to `rest.kegg.jp` on the first run (downloads are cached in `kegg/`).

## Results

Significant genes (|log2FC| ≥ 2 and padj < 0.05; all 11 are on the map):

| Entrez | Symbol | log2FC | padj | Colour |
|---|---|---|---|---|
| 10241 | CALCOCO2 | +13.08 | 8.5e-28 | `#fde725` (clamped) |
| 79671 | NLRX1 | +10.56 | 2.5e-18 | `#fde725` (clamped) |
| 10133 | OPTN | +5.21 | 2.3e-174 | `#a5db36` |
| 4580 | MTX1 | +4.58 | 8.4e-123 | `#84d44b` |
| 5071 | PRKN | +3.82 | 0.0149 | `#63cb5f` |
| 1460 | CSNK2B | +3.23 | 8.0e-173 | `#4ac16d` |
| 84557 | MAP1LC3A | +2.34 | 3.3e-24 | `#2eb37c` |
| 7316 | UBC | +2.25 | 1.3e-147 | `#2cb17e` |
| 65018 | PINK1 | +2.22 | 6.6e-18 | `#2ab07f` |
| 285973 | ATG9B | −4.96 | 1.4e-70 | `#481d6f` |
| 338382 | RAB7B | −7.08 | 2.8e-57 | `#440154` (clamped) |

Box counts: 15 significant (PRKN appears three times on the map; MAP1LC3A/LC3 and
RAB7B twice each), 54 measured but not significant (grey), 8 not measured (white).
PRKN is the only gene whose call depends on the padj cutoff: it passes at 0.05 and
would fail at 0.01.

### Changes from v1

- **Significance (v2.1)**: genes are coloured only when |log2FC| ≥ 2 **and** padj < 0.05. The key has a grey band over −2…+2, and unmeasured genes are white rather than grey.
- **Scale**: v1 used an asymmetric scale (−7.08 … +13.08) whose teal midpoint
  fell at +3.0, so PINK1 (+2.2), UBC, MAP1LC3A, CSNK2B and PRKN (+3.8) all
  rendered teal — visually indistinguishable from "no change". v2 centres the
  scale on 0, and v2.1 removes colour from the non-significant interval entirely.
- **Labels**: v1 labelled boxes with KEGG's representative gene (e.g. "RAB7A",
  "ATG9A", "MTX2", "CSNK2A1") even when a different family member carried
  the colour. v2 labels coloured boxes with the driver gene (RAB7B, ATG9B,
  MTX1, CSNK2B); `data/node_table.csv` lists every member of every box.
- **Toolchain**: no R / pathview / `org.Hs.eg.db` dependency.
