"""AlphaMissense plots for FBXO7 (human Q9Y3I1). Main Fig 3B: M. myotis + M. nigricans sites. Supplement: all final four-species sites.
Inputs: AlphaMissense_Q9Y3I1.csv (AlphaFold DB per-protein table), SuppTable16_AlphaMissense_updated.csv"""
import pathlib, os
HERE = pathlib.Path(__file__).resolve().parent if '__file__' in globals() else pathlib.Path('scripts').resolve()
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True); TSV_PATH = str(DATA / 'all_decoys_interface_metrics.tsv')

import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'pdf.fonttype':42,'axes.linewidth':0.6})
am = pd.read_csv(str(DATA / 'AlphaMissense_Q9Y3I1.csv')); am['pos'] = am.protein_variant.str[1:-1].astype(int)
mp = am.groupby('pos').am_pathogenicity.mean()
T16 = pd.read_csv(str(DATA / 'SuppTable16_AlphaMissense_updated.csv')).set_index('Variant')
SETS = {'Myotis myotis':['T19E','T47A','A64T','S109P','Q127E','F146V','G409I'], 'Myotis nigricans':['T19E','T47A','S109P','Q127E','F146V','G409I'],
        'Desmodus rotundus':['T19E','T47A','S109P','S110H','Q127E','F146V','D191G','L290P','G409R'],
        'Diphylla ecaudata':['T19E','T47E','S109P','S110C','Q127E','F146V','L290P','E292R','G409R']}
vir = plt.get_cmap('viridis'); PUR = vir(0.0); GRN = plt.get_cmap('viridis')(0.72)

def am_plot(variants, place, stem, short_labels, lab_col=PUR, dot_col=GRN):
    by_pos = {}
    for v in variants: by_pos.setdefault(int(v[1:-1]), []).append(v)
    fig, ax = plt.subplots(figsize=(6.2, 3.1))
    ax.plot(mp.index, mp.values, color='#555555', lw=1.0, zorder=2)
    ax.axhspan(0, 0.34, color='#EEEEEE', zorder=0); ax.text(521, 0.02, 'AlphaMissense "likely benign" < 0.34', ha='right', va='bottom', fontsize=5.5, color='#777777')
    for pos, vs in sorted(by_pos.items()):
        lab = (vs[0][0] + str(pos)) if short_labels else (vs[0][:-1] + '/'.join(v[-1] for v in vs))
        lx, y = place[pos]
        ax.plot([pos, pos], [0, 1.0], color=lab_col, lw=0.9, ls=(0, (3, 2)), zorder=1, clip_on=False)
        ax.plot([pos, lx], [1.0, y - 0.015], color=lab_col, lw=0.5, clip_on=False)
        ax.text(lx, y, lab, ha='center', va='bottom', fontsize=6.5, color=lab_col, clip_on=False)
        for v in vs: ax.scatter([pos], [T16.loc[v, 'Variant_pathogenicity']], s=16, color=dot_col, edgecolor='black', lw=0.5, zorder=4)
    ax.scatter([], [], s=16, color=dot_col, edgecolor='black', lw=0.5, label='variant score')
    ax.plot([], [], color='#555555', lw=1.0, label='mean over all possible substitutions at the position')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), fontsize=6, frameon=False, ncol=2)
    ax.set_xlim(1, 522); ax.set_ylim(0, 1.0); ax.set_xticks([1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 522])
    ax.set_xlabel('Residue sequence number', fontsize=7.5); ax.set_ylabel('Mean pathogenicity', fontsize=7.5); ax.tick_params(labelsize=6.5)
    for k in ['top', 'right']: ax.spines[k].set_visible(False)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.72, bottom=0.27)
    fig.savefig(str(FIGS / stem) + '.pdf', bbox_inches='tight'); fig.savefig(str(FIGS / stem) + '.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    return fig

myo = sorted(set(SETS['Myotis myotis']) | set(SETS['Myotis nigricans']), key=lambda v: int(v[1:-1]))
allv = sorted({v for s in SETS.values() for v in s}, key=lambda v: (int(v[1:-1]), v))
am_plot(myo, {19:(19,1.34), 47:(47,1.20), 64:(76,1.34), 109:(100,1.20), 127:(127,1.34), 146:(155,1.20), 409:(409,1.34)}, 'Fig3B_alphamissense_Myotis', True)
am_plot(allv, {19:(19,1.34), 47:(50,1.20), 64:(80,1.34), 109:(92,1.20), 110:(118,1.34), 127:(140,1.20), 146:(165,1.34), 191:(195,1.20),
               290:(275,1.34), 292:(312,1.20), 409:(409,1.34)}, 'FigS_alphamissense_four_species', False, lab_col=vir(0.30), dot_col=vir(1.0))
print(myo); print(allv)
