# Western blot quantification: mouse vs. *Myotis* bat

Densitometric quantification of mitophagy- and ubiquitin-pathway proteins in
primary dermal fibroblasts, comparing 3 mouse and 9 *Myotis* individuals. Each
protein's band intensity was normalized to a loading control (α-tubulin or
GAPDH, as noted per target), then compared between species with an unpaired
(Welch) two-sample *t*-test.

## Samples

| Species | Individuals |
|---|---|
| Mouse (*n* = 3) | M1A, M4D, M5E |
| *Myotis* bat (*n* = 9) | MM6309A/B/C, MM6313A/B/C, MMY6678A/B/E |

## Proteins quantified

Normalized per α-tubulin (unless noted): NDP52, LC3B, P62, OPTN, PRKN, PINK1,
PINK1_v2, BNIP3. Normalized per GAPDH: K48-linked Ub, K63-linked Ub, Total Ub.
Phospho-Ub normalized per α-tubulin.

## Running

```bash
Rscript western_blots.R
```

Requires R with `tidyverse`, `viridis`, and `patchwork`. Densitometry values are
embedded in the script (one vector per protein), so no external input file is
needed. Produces:

- `western_blots_panel.png` — per-protein boxplots (mouse vs. bat), 3×4 panel
- `western_blots_ttests.csv` — tidy table: per-protein group means, *t*, df,
  *p*-value, significance stars, and which species is higher

## Result summary

Five proteins differ significantly between species (Welch *t*-test):

| Protein | *t* | df | *p* | Higher in |
|---|---|---|---|---|
| NDP52 | 10.1 | 8.00 | 7.8e-06 (***) | Bat |
| LC3B | 8.15 | 8.07 | 3.6e-05 (***) | Bat |
| K63-linked Ub | 7.61 | 10.0 | 1.8e-05 (***) | Bat |
| Total Ub | 4.78 | 9.05 | 9.8e-04 (***) | Bat |
| PRKN | 3.34 | 9.28 | 8.3e-03 (**) | Bat |

The remaining proteins (P62, OPTN, PINK1, PINK1_v2, BNIP3, Phospho-Ub,
K48-linked Ub) show no significant species difference.

Significance codes: `***` *p* < 0.001, `**` *p* < 0.01, `*` *p* < 0.05.
