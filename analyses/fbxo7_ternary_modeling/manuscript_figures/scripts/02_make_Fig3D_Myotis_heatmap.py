"""Fig 3D (Myotis only, viridis): medians and calls read from TableS18_rosetta_ddG_decoy_stats.csv (single source of truth)."""
import pathlib, os
HERE = pathlib.Path(__file__).resolve().parent if '__file__' in globals() else pathlib.Path('scripts').resolve()
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True); TSV_PATH = str(DATA / 'all_decoys_interface_metrics.tsv')

import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'pdf.fonttype':42,'axes.linewidth':0.6})
G = '#888888'
t = pd.read_csv(str(DATA / 'TableS18_rosetta_ddG_decoy_stats.csv')); cc = [c for c in t.columns if c.startswith('Call')][0]
NOISE = float(cc.split('floor ')[1].split(' REU')[0]); assert NOISE == 2.0, NOISE
nm = {'MyoMyosites_rerun':'Myotis myotis site-set','MyoNigsites_rerun':'Myotis nigricans site-set'}
groups = [("Ubl/Parkin-binding (1–88)", ['T19E','T47A','A64T']), ("PINK1-binding (92–129)", ['S109P','Q127E']),
          ("CDK6-binding (129–169)", ['F146V']), ("CDK6-binding (381–522)", ['G409I']), ("Species site-sets", list(nm))]
order = [s for _, ss in groups for s in ss]; cols = ['FBXO7-PINK1','FBXO7-PSMF1']
key = lambda s: nm.get(s, s)
M = np.array([[t[(t.Variant==key(s))&(t.Interface==c)]['ddG median (REU)'].item() for c in cols] for s in order])
S = np.array([[t[(t.Variant==key(s))&(t.Interface==c)][cc].item()!='ns' for c in cols] for s in order])
vir = plt.get_cmap('viridis'); cmap = LinearSegmentedColormap.from_list('v2', [vir(0.72), (1,1,1,1), vir(0.0)])
vmax = float(np.ceil(np.abs(M).max())); norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
lab = {'MyoMyosites_rerun':'$\\it{M.\\ myotis}$ (all 7 sites)','MyoNigsites_rerun':'$\\it{M.\\ nigricans}$ (6 sites)','A64T':'A64T †'}
fmt = lambda v: f"{v:.0f}" if abs(v) >= 10 else f"{v:.1f}"
fig, ax = plt.subplots(figsize=(5.0, 3.8)); ax.imshow(np.where(S, M, np.nan), cmap=cmap, norm=norm, aspect='auto')
for i in range(M.shape[0]):
    for j in range(2):
        if not S[i,j]: ax.add_patch(mpl.patches.Rectangle((j-.5,i-.5),1,1,facecolor='#E3E3E3',edgecolor='none'))
        rgba = cmap(norm(M[i,j])) if S[i,j] else (.89,.89,.89,1); lum = .2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
        ax.text(j,i,fmt(M[i,j]),ha='center',va='center',fontsize=7,color='white' if lum<.45 else ('#555555' if not S[i,j] else 'black'))
ax.set_xticks(np.arange(-.5,2,1),minor=True); ax.set_yticks(np.arange(-.5,len(order),1),minor=True); ax.grid(which='minor',color='white',lw=1.2); ax.tick_params(which='minor',length=0)
ax.set_xticks([0,1]); ax.set_xticklabels(['PINK1\ninterface','PSMF1\ninterface'],fontsize=6.5); ax.xaxis.tick_top()
ax.text(0.5,-1.55,'FBXO7 bound to:',ha='center',va='bottom',fontsize=6.5)
ax.set_yticks(range(len(order))); ax.set_yticklabels([lab.get(s,s) for s in order],fontsize=7); ax.tick_params(length=0,pad=4)
for sp in ax.spines.values(): sp.set_visible(False)
y0 = -.5
for name, ss in groups:
    y1 = y0+len(ss)
    if y0 > -.5: ax.axhline(y0,color='black',lw=.8)
    ax.plot([1.56,1.56],[y0+.12,y1-.12],color=G,lw=.9,clip_on=False)
    ax.text(1.64,(y0 + y1) / 2,name if len(ss)==1 else name.replace(' (','\n('),ha='left',va='center',fontsize=6,color=G,clip_on=False,linespacing=1.0); y0 = y1
ax.set_xlim(-.5,1.5); ax.set_ylim(len(order)-.5,-.5)
cax = fig.add_axes([0.80,0.20,0.028,0.40]); cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_ticks([-15,-10,-5,0,5,10,15]); cb.ax.tick_params(labelsize=6,length=2); cb.set_label('ΔΔG (REU)\nmedian of 25 decoys',fontsize=6.3)
cb.ax.text(.5,1.05,'destabilizing',transform=cb.ax.transAxes,ha='center',va='bottom',fontsize=6,color=G); cb.ax.text(.5,-.05,'stabilizing',transform=cb.ax.transAxes,ha='center',va='top',fontsize=6,color=G)
fig.text(0.30,0.025,f'grey: within WT noise (|ΔΔG| ≤ {NOISE} REU or 95% CI spans 0)\n† $\\it{{M.\\ myotis}}$ only; all other sites are shared by both species',fontsize=5.8,color='#555555',va='bottom')
fig.subplots_adjust(left=0.29,right=0.54,top=0.83,bottom=0.12)
fig.savefig(str(FIGS / 'Fig3D_rosetta_ddG_Myotis_viridis.pdf'),bbox_inches='tight'); fig.savefig(str(FIGS / 'Fig3D_rosetta_ddG_Myotis_viridis.png'),dpi=300,bbox_inches='tight',pad_inches=0.05)
print(pd.DataFrame(M,index=order,columns=cols).round(2).assign(PINK1_sig=S[:,0],PSMF1_sig=S[:,1]).to_string())
