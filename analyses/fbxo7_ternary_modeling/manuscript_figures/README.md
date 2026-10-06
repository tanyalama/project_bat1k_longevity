# Manuscript figures: FBXO7 Rosetta interface ddG and AlphaMissense (v2)

Decoy-level Rosetta analysis and AlphaMissense lookup for the final four-species FBXO7 site-sets. This folder is self-contained:
every figure and supplementary table below regenerates from `data/` by running the scripts in numeric order
(`python scripts/00_...py` through `05_...py`; any working directory). Requires numpy, pandas, scipy, matplotlib.

It supersedes the single-decoy results in `../results/fbxo7_ddg_combined.csv` and `../scripts/plot_ddg_heatmap.py` (kept, not deleted).

## Final site-sets (human FBXO7, UniProt Q9Y3I1)

| Lineage | Substitutions |
|---|---|
| *Myotis myotis* | T19E; T47A; A64T; S109P; Q127E; F146V; G409I |
| *Myotis nigricans* | T19E; T47A; S109P; Q127E; F146V; G409I |
| *Desmodus rotundus* | T19E; T47A; S109P; S110H; Q127E; F146V; D191G; L290P; G409R |
| *Diphylla ecaudata* | T19E; T47E; S109P; S110C; Q127E; F146V; L290P; E292R; G409R |

P17Q, Q119R and I151V were modeled in Rosetta but are in no final site-set; they appear only in Tables S16/S18.

## Figures

| File (`figures/`) | Content | Script |
|---|---|---|
| `Fig3B_alphamissense_Myotis` | Main Fig 3B: AlphaMissense mean pathogenicity per residue, *M. myotis* + *M. nigricans* sites | `03_make_AlphaMissense_figs.py` |
| `Fig3D_rosetta_ddG_Myotis_viridis` | Main Fig 3D: Rosetta ddG heatmap, *Myotis* sites | `02_make_Fig3D_Myotis_heatmap.py` |
| `FigS_alphamissense_four_species` | Supplement: AlphaMissense, all 14 final-set variants | `03_make_AlphaMissense_figs.py` |
| `FigS_ddG_heatmap_four_species` | Supplement: ddG heatmap, 14 variants + 4 species site-sets | `04_make_FigS_heatmap_four_species.py` |
| `FigS9_rosetta_ddG_decoys_viridis` | Supplement: per-decoy ddG rainplot (180 x 100 mm) | `05_make_FigS_rainplot.py` |

Supplement figure numbers are placeholders. Colour: viridis green = stabilizing, purple = destabilizing, grey = not distinguishable from WT.

## Data and tables (`data/`)

| File | Description |
|---|---|
| `all_decoys_interface_metrics.tsv` | Input. Relax decoys (25 per sample per interface; WT + 21 variants/site-sets) with `dG_separated` and interface metrics |
| `TableS18_rosetta_ddG_decoy_stats.csv` | Supp Table 18 (42 rows): median ddG vs WT, bootstrap 95% CI, Mann-Whitney p, BH-adjusted p, call. Written by `01_...py` |
| `AlphaMissense_Q9Y3I1.csv` | AlphaMissense scores for all FBXO7 missense variants (AlphaFold DB `AF-Q9Y3I1-F1`; Cheng et al. 2023, CC BY 4.0) |
| `SuppTable16_AlphaMissense_updated.csv` | Supp Table 16 (18 variants). Written by `00_...py` |

## Methods summary

- **ddG** = median `dG_separated` of the 25 decoys of a sample minus the median of the 25 WT decoys, per interface (REU; negative = stabilizing).
- **Interval and test**: bootstrap (10,000 resamples of n = 25) of the difference of medians gives the 95% CI; two-sided Mann-Whitney U against WT; Benjamini-Hochberg across all 42 tests.
- **Call**: significant only if BH-adjusted p < 0.05, the CI excludes 0, and |ddG| > 2.0 REU. The 2.0 REU floor is fixed for consistency across table and figures.
  A WT-vs-WT resampling estimate (97.5th percentile of the median difference) is 2.36 REU at the PINK1 and 2.14 REU at the PSMF1 interface;
  thresholds of 2.0, 2.3 and 2.4 give the same 34 significant calls.
- **Baseline band (rainplot, PINK1 panel)**: range of medians of sites outside the PINK1- and PSMF1-binding regions (Ubl and both CDK6 segments; species site-sets excluded), -10.2 to -5.9 REU.
- **Regions** (UniProt Q9Y3I1): Ubl/Parkin-binding 1-88, PINK1-binding 92-129, CDK6-binding 129-169, PSMF1/dimerization 180-324, CDK6-binding 381-522.
- **AlphaMissense**: mean pathogenicity = mean of the 19 substitution scores at a residue; "likely benign" < 0.34. Parkin was not part of the ternary model, so effects of the Ubl-region variants on Parkin binding were not assessed.
