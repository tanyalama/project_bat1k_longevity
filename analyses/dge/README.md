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
│   ├── fc_data.tsv                 # input: Entrez gene ID + log2fc (tab-separated, 87 genes)
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

1. **Input** (`data/fc_data.tsv`): one row per gene, `<Entrez ID>\t<log2fc>`.
   `NA` values are treated as not measured.
2. **KEGG map**: `hsa04137` KGML and base PNG are fetched from
   `rest.kegg.jp` (cached in `kegg/`). Each KGML `gene` entry gives a box
   position and the Entrez IDs it collapses (77 boxes, 105 Entrez IDs; 84 of
   the 87 input genes are on the map).
3. **Significance threshold**: genes with `|log2FC| < 2` are not significant
   and rendered grey (`#CCCCCC`), as are unmeasured genes.
4. **Box colouring**: a box is coloured by the significant member with the
   largest |log2FC| (its *driver*) and **labelled with the driver's symbol**,
   so the label on a coloured box always names the gene carrying the signal.
   In this dataset no box has more than one significant member, so the
   aggregation rule never has to arbitrate. Uncoloured boxes are labelled
   with the first member that has an HGNC symbol.
5. **Colour scale**: matplotlib `viridis`, **symmetric and centred on 0**
   (`±7`, clamped). 0 is the teal midpoint, so every significant up-gene lands
   in the green→yellow half and every down-gene in the purple half. Three
   values exceed the limit and are clamped to the end colours: CALCOCO2
   (+13.08), NLRX1 (+10.56) → yellow; RAB7B (−7.08) → purple. The in-figure
   key and `figures/legend.*` mark the ends as `≤ −7` / `≥ 7`.
6. **Typography**: Arial for all overlaid text; PDF text is embedded as
   TrueType (editable). KEGG's own pathway text is part of the raster.

All tunable parameters (`SIG_THRESH`, `LIMIT`, `CMAP`, dataset caption) live
in the `CONFIG` block at the top of `scripts/mitophagy_kegg_overlay.py`.

## Reproducing

```bash
python scripts/mitophagy_kegg_overlay.py
```

Requires Python 3 with `matplotlib`, `numpy`, `pandas`, `pillow`, and network
access to `rest.kegg.jp` on the first run (downloads are cached in `kegg/`).

## Results

Significant genes (`|log2FC| ≥ 2`; all 11 are on the map):

| Entrez | Symbol | log2FC | Colour |
|---|---|---|---|
| 10241 | CALCOCO2 | +13.08 | `#fde725` (clamped) |
| 79671 | NLRX1 | +10.56 | `#fde725` (clamped) |
| 10133 | OPTN | +5.21 | `#aadc32` |
| 4580 | MTX1 | +4.58 | `#8bd646` |
| 5071 | PRKN | +3.82 | `#69cd5b` |
| 1460 | CSNK2B | +3.23 | `#54c568` |
| 84557 | MAP1LC3A | +2.34 | `#35b779` |
| 7316 | UBC | +2.25 | `#34b679` |
| 65018 | PINK1 | +2.22 | `#32b67a` |
| 285973 | ATG9B | −4.96 | `#46337f` |
| 338382 | RAB7B | −7.08 | `#440154` (clamped) |

15 boxes are coloured (PRKN appears three times on the map, MAP1LC3A/LC3 and
RAB7B twice); 54 boxes are measured but not significant; 8 boxes contain no
measured gene.

### Changes from v1

- **Scale**: v1 used an asymmetric scale (−7.08 … +13.08) whose teal midpoint
  fell at +3.0, so PINK1 (+2.2), UBC, MAP1LC3A, CSNK2B and PRKN (+3.8) all
  rendered teal — visually indistinguishable from "no change". v2 centres the
  scale on 0.
- **Labels**: v1 labelled boxes with KEGG's representative gene (e.g. "RAB7A",
  "ATG9A", "MTX2", "CSNK2A1") even when a different family member carried
  the colour. v2 labels coloured boxes with the driver gene (RAB7B, ATG9B,
  MTX1, CSNK2B); `data/node_table.csv` lists every member of every box.
- **Toolchain**: no R / pathview / `org.Hs.eg.db` dependency.
