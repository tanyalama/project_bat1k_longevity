#!/usr/bin/env bash
# The two Figure 3C versions. Run from the manuscript_figures folder, in an environment with PyMOL.
# Each pass writes figures/Fig3C_render.png + data/Fig3C_{anchors.json,chainmask.png}; those are renamed
# after the labels are composed, so the second pass starts from a clean slate.
set -e
S=scripts; F=figures; D=data

# version B -- M. myotis sites plus the three PSMF1-region positions (191, 290, 292)
python $S/06_make_Fig3C_render.py --sites 19:T19E 47:T47A 64:A64T 109:S109P 119:Q119R 127:Q127E 146:F146V \
                                          191:D191 290:L290 292:E292 409:G409I --omit F146V
python $S/07_compose_Fig3C_labels.py --out Fig3C_ternary_myotis_psmf1 --size_mm 90 100
mv $F/Fig3C_render.png $F/Fig3C_render_myotis_psmf1.png
mv $D/Fig3C_anchors.json $D/Fig3C_anchors_myotis_psmf1.json
mv $D/Fig3C_chainmask.png $D/Fig3C_chainmask_myotis_psmf1.png

# version A -- every position carrying a modelled substitution in any of the four focal species
python $S/06_make_Fig3C_render.py --sites 17:P17Q 19:T19E 47:T47A/E 64:A64T 109:S109P 110:S110C/H 119:Q119R \
                                          127:Q127E 146:F146V 151:I151V 191:D191G 290:L290P 292:E292R 409:G409I/R --omit
python $S/07_compose_Fig3C_labels.py --out Fig3C_ternary_four_species --size_mm 90 100
mv $F/Fig3C_render.png $F/Fig3C_render_four_species.png
mv $D/Fig3C_anchors.json $D/Fig3C_anchors_four_species.json
mv $D/Fig3C_chainmask.png $D/Fig3C_chainmask_four_species.png
