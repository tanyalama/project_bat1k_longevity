"""Fig S9 (180 mm x 100 mm): per-decoy Rosetta interface ddG 'rainplot' (viridis two-colour). 
Inputs: TSV_PATH (all_decoys_interface_metrics.tsv) and TableS18_rosetta_ddG_decoy_stats.csv (calls, medians, noise floor).
Colour = call from Table S18: viridis green = stabilizing, viridis purple = destabilizing, grey = within WT noise."""
import pathlib, os
HERE = pathlib.Path(__file__).resolve().parent if '__file__' in globals() else pathlib.Path('scripts').resolve()
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True); TSV_PATH = str(DATA / 'all_decoys_interface_metrics.tsv')

import re, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

META_GREY = globals().get('META_GREY', '#7F7F7F')
try: apply_house_style()
except NameError: plt.rcParams['font.family'] = 'Arial'

d = pd.read_csv(TSV_PATH, sep='\t')
t = pd.read_csv(str(DATA / 'TableS18_rosetta_ddG_decoy_stats.csv'))
callcol = [c for c in t.columns if c.startswith('Call')][0]
NOISE = 2.0  # Fixed for consistency across figures
nm = {'Myotis myotis site-set':'MyoMyosites_rerun','Myotis nigricans site-set':'MyoNigsites_rerun',
      'Desmodus rotundus site-set':'DesRotsites_rerun','Diphylla ecaudata site-set':'DipEcasites_rerun'}
t['sample'] = t['Variant'].map(lambda v: nm.get(v, v))
cols = ['FBXO7-PINK1', 'FBXO7-PSMF1']
call = t.set_index(['sample', 'Interface'])[callcol]; med = t.set_index(['sample', 'Interface'])['ddG median (REU)']

groups = [("Ubl/Parkin-binding\n(1–88)", ['T19E','T47A','T47E','A64T']), ("PINK1-binding\n(92–129)", ['S109P','S110H','S110C','Q127E']),
          ("CDK6-binding\n(129–169)", ['F146V']), ("PSMF1-binding\n(180–324)", ['D191G','L290P','E292R']),
          ("CDK6-binding\n(381–522)", ['G409I','G409R']), ("Species site-sets", list(nm.values()))]
pretty = {'MyoMyosites_rerun':'$\\it{Myotis\\ myotis}$','MyoNigsites_rerun':'$\\it{Myotis\\ nigricans}$',
          'DesRotsites_rerun':'$\\it{Desmodus\\ rotundus}$','DipEcasites_rerun':'$\\it{Diphylla\\ ecaudata}$'}
order = [s for _, ss in groups for s in ss]
assert set(order) <= set(t['sample']), set(order) - set(t['sample'])

vir = plt.get_cmap('viridis')
COL = {'stabilizing': mpl.colors.to_hex(vir(0.72)), 'destabilizing': mpl.colors.to_hex(vir(0.0)), 'ns': '#9A9A9A'}
rng = np.random.default_rng(1)

# ---- fixed 18 cm wide x 10 cm high canvas; axes placed in mm so the saved file is exactly this size ----
MM = 1 / 25.4; W, H = 180, 100
L, R, B, T, GAP = 23, 22, 13.5, 5.5, 7                      # margins (mm): y labels | group labels | x label + legend | titles
pw = (W - L - R - GAP) / 2; ph = H - B - T
fig = plt.figure(figsize=(W * MM, H * MM))
axes = [fig.add_axes([L / W, B / H, pw / W, ph / H])]
axes.append(fig.add_axes([(L + pw + GAP) / W, B / H, pw / W, ph / H], sharey=axes[0]))
FS_T, FS_TICK, FS_LAB, FS_SMALL = 7, 6, 6.5, 5.5
for ax, lab in zip(axes, cols):
    wmed = np.median(d[(d['sample'] == 'WT') & (d.label == lab)]['dG_separated'].values)
    ax.axvspan(-NOISE, NOISE, color='#E4E4E4', zorder=0)
    ax.axvline(0, color='black', lw=0.6, ls='--', zorder=1)
    for i, s in enumerate(order):
        v = d[(d['sample'] == s) & (d.label == lab)]['dG_separated'].values - wmed
        c = COL[call[(s, lab)]]; faded = call[(s, lab)] == 'ns'
        p = ax.violinplot([v], positions=[i], orientation='horizontal', widths=0.85, showextrema=False)
        for b in p['bodies']: b.set_facecolor(c); b.set_edgecolor('none'); b.set_alpha(0.18 if faded else 0.32)
        ax.scatter(v, np.full_like(v, i) + rng.uniform(-.2, .2, len(v)), s=2.2, color=c, alpha=0.45 if faded else 0.8, lw=0, zorder=2)
        ax.plot([np.median(v)] * 2, [i - .38, i + .38], color='black', lw=1.0, zorder=3)
    y0 = -.5
    for name, ss in groups:
        y1 = y0 + len(ss)
        if y0 > -.5: ax.axhline(y0, color=META_GREY, lw=0.5, zorder=1)
        if lab == 'FBXO7-PSMF1':
            ax.plot([1.02, 1.02], [y0 + .15, y1 - .15], color=META_GREY, lw=0.7, clip_on=False, transform=ax.get_yaxis_transform())
            ax.text(1.04, (y0 + y1) / 2, name, ha='left', va='center', fontsize=FS_SMALL,
                    color=META_GREY, clip_on=False, linespacing=0.95, transform=ax.get_yaxis_transform())
        y0 = y1
    ax.set_title(lab.replace('-', '–') + ' interface', loc='left', fontsize=FS_T, pad=2)
    ax.set_xlabel('ΔΔG (REU): decoy − WT median', fontsize=FS_LAB, labelpad=1.5)
    ax.tick_params(labelsize=FS_TICK, length=2, width=0.5, pad=1.5)
    for k in ['top', 'right']: ax.spines[k].set_visible(False)
    for k in ['left', 'bottom']: ax.spines[k].set_linewidth(0.5)
    ax.margins(x=0.03); ax.set_ylim(len(order) - .5, -0.7)
axes[1].tick_params(labelleft=False)
axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels([pretty.get(s, s) for s in order], fontsize=FS_TICK)
handles = [Patch(facecolor=COL['stabilizing'], alpha=0.6, label='stabilizing'), Patch(facecolor=COL['destabilizing'], alpha=0.6, label='destabilizing'),
           Patch(facecolor=COL['ns'], alpha=0.4, label='not distinguishable from WT'), Patch(facecolor='#E4E4E4', label=f'±{NOISE} REU WT-vs-WT noise'),
           Line2D([0], [0], color='black', lw=1.0, label='median of 25 decoys')]
fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=FS_SMALL, bbox_to_anchor=(0.5, 0.0),
           handlelength=1.2, handleheight=0.8, columnspacing=1.0, borderaxespad=0.3)
mpl.rcParams['savefig.bbox'] = 'standard'                # keep the exact 180 x 100 mm canvas (no tight cropping)
fig.savefig(str(FIGS / 'FigS9_rosetta_ddG_decoys_viridis.pdf'), bbox_inches=None)
fig.savefig(str(FIGS / 'FigS9_rosetta_ddG_decoys_viridis.png'), dpi=300, bbox_inches=None)
