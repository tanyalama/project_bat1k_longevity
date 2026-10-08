"""Supp Table 16: AlphaMissense scores for FBXO7 (human, UniProt Q9Y3I1) variants.
AlphaMissense scores are a lookup, not a rerun: data/AlphaMissense_Q9Y3I1.csv is the AlphaFold DB per-protein table
(https://alphafold.ebi.ac.uk/files/AF-Q9Y3I1-F1-aa-substitutions.csv; Cheng et al. 2023, CC BY 4.0).
Position mean = mean of the 19 substitution scores at that residue. Classes: LBen = benign (<0.34), Amb = ambiguous, LPath = pathogenic (>0.564).
Variants = the 17 single substitutions modelled with Rosetta (Table S18): the 14 in the final site-sets plus P17Q, Q119R and I151V (in no site-set)."""
import pathlib, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; DATA = HERE.parent / 'data'
am = pd.read_csv(DATA / 'AlphaMissense_Q9Y3I1.csv')
am['pos'] = am.protein_variant.str[1:-1].astype(int); am['wt'] = am.protein_variant.str[0]; am['mut'] = am.protein_variant.str[-1]
meanpos = am.groupby('pos').am_pathogenicity.mean(); cls = {'LBen': 'benign', 'Amb': 'ambiguous', 'LPath': 'pathogenic'}
SETS = {'Myotis myotis': ['T19E','T47A','A64T','S109P','Q127E','F146V','G409I'], 'Myotis nigricans': ['T19E','T47A','S109P','Q127E','F146V','G409I'],
        'Desmodus rotundus': ['T19E','T47A','S109P','S110H','Q127E','F146V','D191G','L290P','G409R'],
        'Diphylla ecaudata': ['T19E','T47E','S109P','S110C','Q127E','F146V','L290P','E292R','G409R']}
ROSETTA = ['P17Q','T19E','T47A','T47E','A64T','S109P','S110H','S110C','Q119R','Q127E','F146V','I151V','D191G','L290P','E292R','G409R','G409I']
allv = sorted(set(ROSETTA) | {v for s in SETS.values() for v in s}, key=lambda v: (int(v[1:-1]), v))
rows = []
for v in allv:
    r = am[am.protein_variant == v].iloc[0]
    rows.append(dict(Variant=v, Position=r.pos, WT_residue=r.wt, Variant_residue=r.mut, Mean_pathogenicity_at_position=round(meanpos[r.pos], 4),
                     Variant_pathogenicity=r.am_pathogenicity, AlphaMissense_class=cls[r.am_class],
                     Species_site_sets=', '.join(s for s, l in SETS.items() if v in l) or 'none (Rosetta-modelled only)'))
pd.DataFrame(rows).to_csv(DATA / 'SuppTable16_AlphaMissense_updated.csv', index=False)
