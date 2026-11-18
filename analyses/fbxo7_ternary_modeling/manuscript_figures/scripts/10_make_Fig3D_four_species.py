"""Fig 3D (all four species, viridis; 90 x 112 mm): medians and calls read from TableS18_rosetta_ddG_decoy_stats.csv (single source of truth)."""
import pathlib, os
HERE = pathlib.Path(__file__).resolve().parent if '__file__' in globals() else pathlib.Path('scripts').resolve()
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True); TSV_PATH = str(DATA / 'all_decoys_interface_metrics.tsv')

import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'pdf.fonttype':42,'axes.linewidth':0.6})
G = '#888888'
t = pd.read_csv(str(DATA / 'TableS18_rosetta_ddG_decoy_stats.csv')); cc = [c for c in t.columns if c.startswith('Call')][0]
NOISE = float(cc.split('floor ')[1].split(' REU')[0]); assert NOISE == 2.0, NOISE
nm = {'MyoMyosites_rerun':'Myotis myotis site-set','MyoNigsites_rerun':'Myotis nigricans site-set','DesRotsites_rerun':'Desmodus rotundus site-set','DipEcasites_rerun':'Diphylla ecaudata site-set'}
groups = [("Ubl/Parkin-binding (1–88)", ['P17Q','T19E','T47A','T47E','A64T']), ("PINK1-binding (92–129)", ['S109P','S110H','S110C','Q119R','Q127E']),
          ("CDK6-binding (129–169)", ['F146V','I151V']), ("PSMF1-binding (180–324)", ['D191G','L290P','E292R']),
          ("CDK6-binding (381–522)", ['G409I','G409R']), ("Species site-sets", list(nm))]
order = [s for _, ss in groups for s in ss]; cols = ['FBXO7-PINK1','FBXO7-PSMF1']
key = lambda s: nm.get(s, s)
M = np.array([[t[(t.Variant==key(s))&(t.Interface==c)]['ddG median (REU)'].item() for c in cols] for s in order])
S = np.array([[t[(t.Variant==key(s))&(t.Interface==c)][cc].item()!='ns' for c in cols] for s in order])
vir = plt.get_cmap('viridis'); cmap = LinearSegmentedColormap.from_list('v2', [vir(0.72), (1,1,1,1), vir(0.0)])
vmax = float(np.ceil(np.abs(M).max())); norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
lab = {'MyoMyosites_rerun':'$\\it{M.\\ myotis}$ (8 sites)','MyoNigsites_rerun':'$\\it{M.\\ nigricans}$ (8 sites)','DesRotsites_rerun':'$\\it{D.\\ rotundus}$ (9 sites)','DipEcasites_rerun':'$\\it{D.\\ ecaudata}$ (9 sites)'}
fmt = lambda v: f"{v:.0f}" if abs(v) >= 10 else f"{v:.1f}"
# ---- fixed 90 mm wide x 70 mm high canvas (right column of the 180 x 190 mm Fig 3; axes placed in mm, saved without tight cropping) ----
MM = 1 / 25.4; W, H = 90, 112
L, B, T, HW = 25, 8, 11, 29                     # mm: row-label margin | footnote | headers | heatmap width
CAX_X = 78.9                                     # mm: left edge of the colourbar (the drawing is laid out to fill the 90 mm width)
HH = H - B - T
FS_CELL, FS_TICK, FS_SMALL = 6.5, 6, 5.5
fig = plt.figure(figsize=(W * MM, H * MM)); ax = fig.add_axes([L / W, B / H, HW / W, HH / H])
ax.imshow(np.where(S, M, np.nan), cmap=cmap, norm=norm, aspect='auto')
for i in range(M.shape[0]):
    for j in range(2):
        if not S[i,j]: ax.add_patch(mpl.patches.Rectangle((j-.5,i-.5),1,1,facecolor='#E3E3E3',edgecolor='none'))
        rgba = cmap(norm(M[i,j])) if S[i,j] else (.89,.89,.89,1); lum = .2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
        ax.text(j,i,fmt(M[i,j]),ha='center',va='center',fontsize=FS_CELL,color='white' if lum<.45 else ('#555555' if not S[i,j] else 'black'))
ax.set_xticks(np.arange(-.5,2,1),minor=True); ax.set_yticks(np.arange(-.5,len(order),1),minor=True); ax.grid(which='minor',color='white',lw=1.0); ax.tick_params(which='minor',length=0)
ax.set_xticks([0,1]); ax.set_xticklabels(['PINK1\ninterface','PSMF1\ninterface'],fontsize=FS_TICK); ax.xaxis.tick_top()
fig.text((L + HW / 2) / W, (H - 1.2) / H, 'FBXO7 bound to:', ha='center', va='top', fontsize=FS_TICK)
ax.set_yticks(range(len(order))); ax.set_yticklabels([lab.get(s,s) for s in order],fontsize=FS_TICK); ax.tick_params(length=0,pad=2)
for sp in ax.spines.values(): sp.set_visible(False)
y0 = -.5
for name, ss in groups:
    y1 = y0+len(ss); label = name.replace(' (','\n(') if '(' in name else name.replace(' site','\nsite')
    if y0 > -.5: ax.axhline(y0,color='black',lw=.6)
    ax.plot([1+1.0/HW]*2,[y0+.12,y1-.12],color=G,lw=.8,clip_on=False,transform=ax.get_yaxis_transform())
    ax.text(1+2.2/HW,(y0 + y1) / 2,label,ha='left',va='center',fontsize=FS_SMALL,color=G,clip_on=False,linespacing=0.95,transform=ax.get_yaxis_transform()); y0 = y1
ax.set_xlim(-.5,1.5); ax.set_ylim(len(order)-.5,-.5)
cax = fig.add_axes([CAX_X / W, 40 / H, 3 / W, 34 / H]); cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_ticks([-15,-10,-5,0,5,10,15]); cb.ax.tick_params(labelsize=FS_SMALL,length=1.5,width=0.5,pad=1.5); cb.outline.set_linewidth(0.5); cb.set_label('ΔΔG (REU), median of 25 decoys',fontsize=FS_SMALL,labelpad=1.5)
cb.ax.text(.5,1.04,'destabilizing',transform=cb.ax.transAxes,ha='center',va='bottom',fontsize=FS_SMALL,color=G); cb.ax.text(.5,-.04,'stabilizing',transform=cb.ax.transAxes,ha='center',va='top',fontsize=FS_SMALL,color=G)
fig.text(2 / W, 1.5 / H, f'grey: within WT noise (|ΔΔG| ≤ {NOISE} REU or 95% CI spans 0). Site-sets are the\ncombined substitutions of each species, modelled together.',fontsize=FS_SMALL,color='#555555',va='bottom',linespacing=1.15)
mpl.rcParams['savefig.bbox'] = 'standard'
fig.savefig(str(FIGS / 'Fig3D_heatmap_four_species.pdf'),bbox_inches=None); fig.savefig(str(FIGS / 'Fig3D_heatmap_four_species.png'),dpi=300,bbox_inches=None)
print(pd.DataFrame(M,index=order,columns=cols).round(2).assign(PINK1_sig=S[:,0],PSMF1_sig=S[:,1]).to_string())
