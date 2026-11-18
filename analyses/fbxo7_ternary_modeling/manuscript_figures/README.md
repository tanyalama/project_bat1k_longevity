# FBXO7 manuscript figures

Scripts and inputs for the Figure 3 panels and their supplementary variants. Every script resolves its
own paths from `__file__`, so each can be run from any working directory:

    python scripts/<name>.py

`TableS18_rosetta_ddG_decoy_stats.csv` is the single source of truth for every ddG value and call shown
in a figure. No figure recomputes a median, a confidence interval or a call; they are read from that
table, so a figure cannot disagree with it. The table is derived by `01_make_TableS18_rosetta_stats.py`
from `all_decoys_interface_metrics.tsv` (25 Rosetta decoys per variant per interface). The noise floor
is the fixed +/- 2.0 REU stated in the manuscript Methods; cells within it, or whose 95% CI spans zero,
are drawn grey.

## Two versions of Figure 3

Each of the four panels exists in two versions, distinguished by filename suffix.

| suffix | content | panel sizes (W x H mm) |
|---|---|---|
| `_myotis_psmf1` | the ten positions variable between *Myotis myotis* and *M. nigricans*, plus three positions in the UniProt PSMF1-binding region (191, 290, 292) where both *Myotis* retain the human residue | A 100x110, B 110x70, C 90x100, D 90x88 |
| `_four_species` | all 17 modelled substitutions at 14 positions, across all four focal species | A 104x110, B 110x78, C 90x100, D 90x112 |

Positions 191, 290 and 292 lie in the PSMF1-binding region as UniProt annotates it (Q9Y3I1, 180-324),
but **none of the three contacts PSMF1 in the co-folded model**: minimum heavy-atom distances to chain C
are 16.4, 14.4 and 20.8 A. D191 sits 2.4 A from PINK1 instead. The FBXO7-PSMF1 interface in the model is
formed by other residues (104 FBXO7 residues within 5 A of chain C; 31 of them inside 180-324). Treat the
three as region members, not interface contacts.

## Scripts, in run order

| script | output |
|---|---|
| `00_make_SuppTable16_AlphaMissense.py` | Supplementary Table 16 (AlphaMissense scores) |
| `01_make_TableS18_rosetta_stats.py` | Table S18: per-variant ddG medians, CIs and calls |
| `02_make_Fig3D_Myotis_heatmap.py` | Fig 3D, original *Myotis*-only panel (10 rows, 90x70 mm) |
| `03_make_AlphaMissense_figs.py` | every AlphaMissense panel: both Fig 3B versions and the 180 mm supplementary strips |
| `04_make_FigS_heatmap_four_species.py` | supplementary four-species ddG heatmap |
| `05_make_FigS_rainplot.py` | Fig S9 per-decoy ddG rainplot |
| `06_make_Fig3C_render.py` | PyMOL surface render + chain/site mask. `--sites "19:T19E 191:D191 ..."` overrides the default site set |
| `07_compose_Fig3C_labels.py` | places the leader-line labels over the render |
| `08_make_Fig3A_alignment.py` | both Fig 3A alignment panels (`myotis_psmf1`, `four_species` or `both`) |
| `09_make_Fig3D_myotis_psmf1.py` | Fig 3D, *Myotis* + PSMF1-region version |
| `10_make_Fig3D_four_species.py` | Fig 3D, four-species version |
| `11_make_Fig3C_versions.sh` | runs 06 then 07 twice, once per version, renaming the intermediates |

`06` and `11` need PyMOL; everything else needs only numpy, pandas, matplotlib and gemmi.

## Inputs

- `data/alignment/` holds the 30-species FBXO7 codon alignment, the species tree, the display-name table
  and the manual residue overrides used by `08`. `species_names.tsv` is the same file as the copy in
  `analyses/phylo_LQ_timetree/data/`; keep the two in step.
- `data/FBXO7_PINK1kinase_PSMF1_forcedtocif_model_0.cif` is the Boltz ternary model used by `06`.
- `data/Fig3C_anchors_*.json` and `data/Fig3C_chainmask_*.png` are written by `06` and consumed by `07`.
- `data/FBXO7_site_residues_*.csv` and `data/Fig3A_alignment_QC_*.txt` are written by `08`; the QC file
  records any disagreement between the alignment and the modelled variant lists.

## Known discrepancy

The alignment gives *Desmodus rotundus* and *Diphylla ecaudata* a T at position 151 and an E at 119,
neither of which appears in their modelled site-sets. This is recorded in both
`Fig3A_alignment_QC_*.txt` files and is unresolved.
