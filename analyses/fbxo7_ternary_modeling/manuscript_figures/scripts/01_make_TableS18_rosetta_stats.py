"""Decoy-level Rosetta ΔΔG statistics -> Table S18 (figures are made by make_Fig3D_Myotis_viridis.py, make_FigS_heatmap_four_species.py, make_FigS9_decoys_viridis.py).
Input: all_decoys_interface_metrics.tsv (22 samples x 25 decoys x 2 interfaces)."""
import pathlib, os
HERE = pathlib.Path(__file__).resolve().parent if '__file__' in globals() else pathlib.Path('scripts').resolve()
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True); TSV_PATH = str(DATA / 'all_decoys_interface_metrics.tsv')

import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats

d = pd.read_csv(TSV_PATH, sep='\t')
rng = np.random.default_rng(1)
cols = ['FBXO7-PINK1', 'FBXO7-PSMF1']
wt = {lab: d[(d['sample']=='WT') & (d.label==lab)]['dG_separated'].values for lab in cols}

# noise floor: fixed at 2.0 REU for consistency across all figures and tables
NOISE = 2.0

rows = []
for s in sorted(d['sample'].unique()):
    if s == 'WT': continue
    for lab in cols:
        v = d[(d['sample']==s) & (d.label==lab)]['dG_separated'].values; w = wt[lab]
        boot = np.median(rng.choice(v,(10000,25)),1) - np.median(rng.choice(w,(10000,25)),1)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        U, p = stats.mannwhitneyu(v, w, alternative='two-sided')
        rows.append(dict(sample=s, interface=lab, n=len(v), ddG_median=np.median(v)-np.median(w), ddG_mean=v.mean()-w.mean(),
                         ci_lo=lo, ci_hi=hi, MWU_p=p, var_sd=v.std(ddof=1)))
res = pd.DataFrame(rows)
res['padj'] = stats.false_discovery_control(res.MWU_p, method='bh')
res['sig'] = (res.padj < 0.05) & ((res.ci_hi < 0) | (res.ci_lo > 0)) & (res.ddG_median.abs() > NOISE)
res['call'] = np.where(~res.sig, 'ns', np.where(res.ddG_median < 0, 'stabilizing', 'destabilizing'))

groups = [("Ubl/Parkin-binding (1–88)", ['P17Q','T19E','T47A','T47E','A64T']),
          ("PINK1-binding (92–129)", ['S109P','S110H','S110C','Q119R','Q127E']),
          ("CDK6-binding (129–169)", ['F146V','I151V']),
          ("PSMF1-binding (180–324)", ['D191G','L290P','E292R']),
          ("CDK6-binding (381–522)", ['G409I','G409R']),
          ("Species site-sets", ['MyoMyosites_rerun','MyoNigsites_rerun','DesRotsites_rerun','DipEcasites_rerun'])]
pretty = {'MyoMyosites_rerun':'$\\it{Myotis\\ myotis}$','MyoNigsites_rerun':'$\\it{Myotis\\ nigricans}$',
          'DesRotsites_rerun':'$\\it{Desmodus\\ rotundus}$','DipEcasites_rerun':'$\\it{Diphylla\\ ecaudata}$'}
nm = {'MyoMyosites_rerun':'Myotis myotis site-set','MyoNigsites_rerun':'Myotis nigricans site-set',
      'DesRotsites_rerun':'Desmodus rotundus site-set','DipEcasites_rerun':'Diphylla ecaudata site-set'}
order = [s for _, ss in groups for s in ss]
assert set(order) == set(res['sample'])
M = res.pivot(index='sample', columns='interface', values='ddG_median').loc[order, cols]
S = res.pivot(index='sample', columns='interface', values='sig').loc[order, cols].astype(bool)
distal_rows = groups[0][1] + groups[2][1] + groups[4][1]  # Ubl + both CDK6 segments (outside PINK1/PSMF1-binding)

# ---------- Table S18 ----------
tab = res.copy(); tab['variant'] = tab['sample'].map(lambda s: nm.get(s, s))
REGION = {**{s:'Ubl/Parkin-binding (1–88)' for s in ['P17Q','T19E','T47A','T47E','A64T']},
          **{s:'PINK1-binding (92–129)' for s in ['S109P','S110H','S110C','Q119R','Q127E']},
          **{s:'CDK6-binding (129–169)' for s in ['F146V','I151V']},
          **{s:'PSMF1-binding (180–324)' for s in ['D191G','L290P','E292R']},
          **{s:'CDK6-binding (381–522)' for s in ['G409I','G409R']},
          **{s:'Species site-set' for s in ['MyoMyosites_rerun','MyoNigsites_rerun','DesRotsites_rerun','DipEcasites_rerun']}}
grp = REGION; tab['region'] = tab['sample'].map(grp)
tab['order'] = tab['sample'].map({s: i for i, s in enumerate(order)}); tab = tab.sort_values(['order', 'interface'])
tab = tab[['variant','region','interface','n','ddG_median','ci_lo','ci_hi','ddG_mean','var_sd','MWU_p','padj','call']]
tab.columns = ['Variant','FBXO7 region','Interface','n decoys','ddG median (REU)','95% CI low','95% CI high','ddG mean (REU)',
               'SD across decoys','Mann-Whitney p','BH-adjusted p','Call (noise floor 2.0 REU)']
tab.round(4).to_csv(str(DATA / 'TableS18_rosetta_ddG_decoy_stats.csv'), index=False)
