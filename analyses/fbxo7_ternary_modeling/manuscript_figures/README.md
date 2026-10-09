# Manuscript figures: FBXO7 Rosetta interface ddG and AlphaMissense (v2)

Decoy-level Rosetta analysis and AlphaMissense lookup for the final four-species FBXO7 site-sets. This folder is self-contained:
every figure and supplementary table below regenerates from `data/` by running the scripts in numeric order
(`python scripts/00_...py` through `07_...py`; any working directory). Requires numpy, pandas, scipy, matplotlib; the Fig 3C scripts (`06_`, `07_`) also need PyMOL (conda `pymol-open-source`) and Pillow.

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
| `Fig3C_ternary_surface` | Main Fig 3C: FBXO7 (chain A, grey), PINK1 kinase domain (chain B, blue) and PSMF1 (chain C, teal) as opaque surfaces. Yellow patches mark the wild-type FBXO7 residue at each *M. myotis* position that is visible from this view (T19, T47, S109, Q119, Q127); A64 and G409 face away and are named in the footnote; F146V is left out by choice. Canvas 135 x 100 mm is a placeholder until the layout slot is fixed | `06_make_Fig3C_render.py`, `07_compose_Fig3C_labels.py` |
| `FigS_alphamissense_four_species` | Supplement: AlphaMissense, all 14 final-set variants | `03_make_AlphaMissense_figs.py` |
| `FigS_ddG_heatmap_four_species` | Supplement: ddG heatmap, 14 variants + 4 species site-sets | `04_make_FigS_heatmap_four_species.py` |
| `FigS9_rosetta_ddG_decoys_viridis` | Supplement: per-decoy ddG rainplot (180 x 100 mm) | `05_make_FigS_rainplot.py` |

Supplement figure numbers are placeholders. Colour: viridis green = stabilizing, purple = destabilizing, grey = not distinguishable from WT. In `FigS_alphamissense_four_species` the labels and lines are blue and the variant-score dots yellow so they cannot be read as those calls; Fig 3B (*Myotis*) keeps purple labels and green dots.

## Data and tables (`data/`)

| File | Description |
|---|---|
| `all_decoys_interface_metrics.tsv` | Input. Relax decoys (25 per sample per interface; WT + 21 variants/site-sets) with `dG_separated` and interface metrics |
| `TableS18_rosetta_ddG_decoy_stats.csv` | Supp Table 18 (42 rows): median ddG vs WT, bootstrap 95% CI, Mann-Whitney p, BH-adjusted p, call. Written by `01_...py` |
| `AlphaMissense_Q9Y3I1.csv` | AlphaMissense scores for all FBXO7 missense variants (AlphaFold DB `AF-Q9Y3I1-F1`; Cheng et al. 2023, CC BY 4.0) |
| `SuppTable16_AlphaMissense_updated.csv` | Supp Table 16 (17 single substitutions: the 14 in the final site-sets plus P17Q, Q119R, I151V). Written by `00_...py` |
| `FBXO7_PINK1kinase_PSMF1_forcedtocif_model_0.cif` | Input to Fig 3C. Ternary model: chain A FBXO7 (522 residues, identical to UniProt Q9Y3I1), chain B PINK1 kinase domain, chain C PSMF1. The B-factor column appears to hold per-residue pLDDT. How the model was generated is not recorded here |
| `Fig3C_anchors.json` | Written by `06_...py`: visible pixels and image position of each *M. myotis* site patch, from a flat-colour pass with the same camera; read by `07_...py` |
| `Fig3C_chainmask.png` | Written by `06_...py`: per-pixel chain map of the render (0 background, 1 FBXO7, 2 PINK1, 3 PSMF1); read by `07_...py` |

## Methods summary

- **ddG** = median `dG_separated` of the 25 decoys of a sample minus the median of the 25 WT decoys, per interface (REU; negative = stabilizing).
- **Interval and test**: bootstrap (10,000 resamples of n = 25) of the difference of medians gives the 95% CI; two-sided Mann-Whitney U against WT; Benjamini-Hochberg across all 42 tests.
- **Call**: significant only if BH-adjusted p < 0.05, the CI excludes 0, and |ddG| > 2.0 REU. The 2.0 REU floor is fixed for consistency across table and figures.
  A WT-vs-WT resampling estimate (97.5th percentile of the median difference) is 2.36 REU at the PINK1 and 2.14 REU at the PSMF1 interface;
  thresholds of 2.0, 2.3 and 2.4 give the same 34 significant calls.
- **PINK1-interface background**: all ten substitutions outside the PINK1-binding segment read -3 to -10 REU at the PINK1 interface, which suggests an offset between the wild-type and variant models rather than site-specific effects (cause to be confirmed). The rainplot no longer draws a band for it; interpret PINK1-interface values relative to this background (Q127E, -16.3, and S110C, +6.8, are the PINK1-segment substitutions clearly beyond it).
- **Regions** (UniProt Q9Y3I1): Ubl/Parkin-binding 1-88, PINK1-binding 92-129, CDK6-binding 129-169, PSMF1/dimerization 180-324, CDK6-binding 381-522.
- **AlphaMissense**: mean pathogenicity = mean of the scores of all possible substitutions at a residue (19 per position; not a variant count); "likely benign" < 0.34. Parkin was not part of the ternary model, so effects of the Ubl-region variants on Parkin binding were not assessed.
- **Fig 3C**: PyMOL all-surface render of the ternary model with PINK1 above PSMF1 (rotation 90 degrees about the vertical and -30 degrees about the horizontal axis, chosen from a scan of 168 orientations for the visible area of the three PINK1-interface site patches). A flat-colour pass with the same camera measures the visible pixels of each site patch and chain (`Fig3C_anchors.json`, `Fig3C_chainmask.png`); `07_...py` uses them to place labels and leader lines on the white background without overlaps. F146V is deliberately left uncoloured and unlabelled (`--omit F146V`, the default). Keep `surface_quality` 1 and `ray_trace_mode` 0 as set in `06_...py`; outlined ray mode with higher surface quality did not finish in testing.
