# Phylogeny and longevity-quotient figure (103 bats + 8 outgroups)

This folder builds the figure that sets the four focal species in their phylogenetic and
life-history context. The canvas is 180 x 190 mm.

| Panel | Content |
|---|---|
| a | T2 timetree from Morales et al. 2026 (*Nature* 658:141). Clades without focal species are collapsed to families. Bands show the 95% HPD of each node age, and numbers mark nodes with bootstrap support below 100%. |
| b | Longevity quotient (LQ) for each row (open circles = outgroups) |
| c | Body mass vs maximum lifespan: 582 non-flying mammals (open triangles), the 40 study bats with lifespan data (filled circles) and the line expected for non-flying mammals. Ribbons link the four focal species from panel b. |
| d | Focal species cards: *Myotis myotis* vs *M. nigricans*, and *Desmodus rotundus* vs *Diphylla ecaudata* |
| e | Hallmarks of aging banner, in López-Otín order |

![figure](figures/phylo_LQ_hallmarks_180x190.png)

## Node uncertainty (panel a)

Nothing was re-estimated: the support and dating values all come from the authors' own outputs.

- **Dating.** `T2_timetree.tre` is the mean of the 77 MCMCTree runs on coding sequences in
  `Divergence_Time_estimates/Caster_tree/coding_sequences` (all nodes agree within 0.03 Myr).
  Each band runs from the mean of the per-run lower 95% HPD bounds to the mean of the upper bounds. The median band width is about 12 Myr.
  - The CSV also gives the between-run spread of mean ages, with a median width of about 1.7 Myr. This measures how well the runs agree with each other, not how uncertain the dates are, so it is not what is plotted.
  - Bands that extend past the left edge of the plot are clipped and end in an arrowhead. This affects the root, whose HPD runs to about 124 Myr.
- **Topology.** Bootstrap values (`bs1`/1000) and quartet support (`q1`) are taken from the
  neutral-window CASTER tree `multiZhg38.concatNeutralWin_500_COV80.scored_U2.tre`, which has the same topology as T2.
  - 104 of the 108 nodes that can be scored have 100% bootstrap support. The two root nodes cannot be scored because that tree is unrooted.
  - The panel shows the two visible nodes below 100%: anc_Emballonuroidea_Vespertilionoidea (89.9%, q1 0.34) and anc_Scrotifera (68.8%).
  - The other two nodes below 100% fall inside collapsed clades: anc_N087 within *Myotis* (96.8%) and anc_N024 within *Hipposideros* (99.9%).

| Node | Age (Myr) | 95% HPD | Bootstrap | q1 |
|---|---|---|---|---|
| anc_Chiroptera | 66.8 | 64.7–69.6 | 100 | 0.73 |
| anc_Phyllostomidae | 40.4 | 33.5–48.0 | 100 | 0.91 |
| anc_Desmodontinae | 27.0 | 16.8–35.8 | 100 | 0.84 |
| anc_Myotis | 16.2 | 11.6–22.6 | 100 | 0.96 |

## Scripts

| Script | What it does |
|---|---|
| `scripts/01_extract_node_uncertainty.py` | Parses the MCMCTree `FigTree.tre` files and the CASTER support tree and writes `data/T2_node_support_age_uncertainty.csv`. It must be run where the Morales et al. outputs are unpacked; on Unity that is `/scratch4/workspace/tlama_smith_edu-project_longevity/dating/`. The command is in the script header. |
| `scripts/02_make_phylo_LQ_figure.py` | Builds the figure from `data/` and writes `figures/phylo_LQ_hallmarks_180x190.{pdf,png}`. It takes no arguments and is run from this folder. |

`02_make_phylo_LQ_figure.py` was assembled from the interactive session that produced the
approved figure, and it reproduces that figure pixel for pixel.

Requirements: Python 3 with numpy, pandas, matplotlib, dendropy, pillow, scipy and openpyxl. The figure uses the Arial font.

## Data

| File | Content |
|---|---|
| `data/T2_timetree.labeled.nwk` | T2 timetree with a stable label on every internal node (anc_Chiroptera, anc_Myotis, …); ages in units of 100 Myr |
| `data/T2_internal_nodes_retention.tsv` | Internal nodes, their tips and their T2 ages |
| `data/T2_node_support_age_uncertainty.csv` | Output of script 01: per-node age, HPD, between-run spread, neutral-site age, bootstrap, q1 |
| `data/bat103_max_lifespan.xlsx` | Maximum lifespan, body mass and LQ for the 111 taxa, with sources (sheet `MLS_LQ_table`). LQ = MLS / 10^(0.537576 + 0.189234·log10 mass in g). |
| `data/data_lifespan_body_mammals.csv` | Lifespan and body mass for the 583-species mammal comparison set (ln units). The one bat in it, *Lasionycteris noctivagans*, is not in the study and is excluded. |
| `data/images/phylopic_whale_93d5c0be.png`, `phylopic_mouse_fe09db98.png` | PhyloPic silhouettes, both CC0: whale by Jonathan Pople (filed under *Cetotherium riabinini*); mouse by Ben Murrell (*Mus musculus*). The mouse marks the smallest mammal in the set, *Mus minutoides*. |

**Images not committed.** The bat paintings and the hallmarks graphic are artwork, so they
are not committed to this public repository. Place them in `data/images/` before running script 02:
`Myotis_myotis.jpg`, `Myotis_nigricans.jpg`, `Desmodus_rotundus.jpg`, `Diphylla_ecaudata.jpg`, `hallmarks_of_aging_v1.png`.
Credits for the artwork go in the figure caption.
