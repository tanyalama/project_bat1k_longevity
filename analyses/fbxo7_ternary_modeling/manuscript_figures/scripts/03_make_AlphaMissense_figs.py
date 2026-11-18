"""AlphaMissense plots for FBXO7 (human Q9Y3I1). Main Fig 3B: M. myotis substitutions. Supplement: all sites in the four species composites.
Inputs: AlphaMissense_Q9Y3I1.csv (AlphaFold DB per-protein table), SuppTable16_AlphaMissense_updated.csv"""
import pathlib
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
HERE = pathlib.Path(__file__).resolve().parent
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True)
from matplotlib.patches import Patch
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'pdf.fonttype':42,'axes.linewidth':0.6})
am = pd.read_csv(DATA / 'AlphaMissense_Q9Y3I1.csv'); am['pos'] = am.protein_variant.str[1:-1].astype(int)
mp = am.groupby('pos').am_pathogenicity.mean()
T16 = pd.read_csv(DATA / 'SuppTable16_AlphaMissense_updated.csv').set_index('Variant')
SETS = {'Myotis myotis':['T19E','T47A','A64T','S109P','Q119R','Q127E','F146V','G409I'], 'Myotis nigricans':['P17Q','T19E','T47A','S109P','Q127E','F146V','I151V','G409I'],
        'Desmodus rotundus':['T19E','T47A','S109P','S110H','Q127E','F146V','D191G','L290P','G409R'],
        'Diphylla ecaudata':['T19E','T47E','S109P','S110C','Q127E','F146V','L290P','E292R','G409R']}
vir = plt.get_cmap('viridis'); PUR = vir(0.0); GRN = plt.get_cmap('viridis')(0.72)

def _hit(seg, box, n=80):
    (x1, y1), (x2, y2) = seg
    return any(box[0] <= x1 + (x2 - x1) * t <= box[1] and box[2] <= y1 + (y2 - y1) * t <= box[3] for t in np.linspace(0, 1, n))
def _cross(a, b):
    def o(p, q, r): return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0])
    return o(a[0], a[1], b[0]) * o(a[0], a[1], b[1]) < 0 and o(b[0], b[1], a[0]) * o(b[0], b[1], a[1]) < 0
def auto_place(ax, entries, rows=(1.17, 1.30, 1.43, 1.56), pad=5, fs=6.5, row_cost=25):
    """Label layout by depth-first search: labels go left to right, each tried at the free (row, x) slots nearest its position first;
    a slot is free when its box and leader line touch no earlier label or leader. Backtracks when a later label has no free slot."""
    r = ax.figure.canvas.get_renderer(); inv = ax.transData.inverted(); sizes = []
    for pos, lab in entries:
        t = ax.text(0, 0, lab, fontsize=fs); bb = t.get_window_extent(r); t.remove()
        (x0, y0), (x1, y1) = inv.transform((bb.x0, bb.y0)), inv.transform((bb.x1, bb.y1)); sizes.append((x1 - x0, y1 - y0))
    def free(pos, lx, y, w, h, placed):
        box = (lx - w / 2 - pad / 2, lx + w / 2 + pad / 2, y, y + h); seg = ((pos, 1.0), (lx, y - 0.015))
        return None if any(not (box[1] < p[0][0] or box[0] > p[0][1] or box[3] < p[0][2] or box[2] > p[0][3]) or _hit(seg, p[0], 40) or _hit(p[1], box, 40) or _cross(seg, p[1]) for p in placed) else (box, seg)
    def solve(k, placed, out):
        if k == len(entries): return out
        pos, lab = entries[k]; w, h = sizes[k]
        cands = sorted(((abs(dx) + row_cost * ri, pos + dx, y) for ri, y in enumerate(rows) for dx in range(-120, 121, 3) if -8 <= pos + dx - w / 2 and pos + dx + w / 2 <= 530), key=lambda c: c[0])
        for _, lx, y in cands:
            f = free(pos, lx, y, w, h, placed)
            if f:
                res = solve(k + 1, placed + [f], {**out, pos: (lx, y)})
                if res: return res
        return None
    sol = solve(0, [], {}); assert sol, 'no collision-free label layout found'
    return sol

def am_plot(variants, place, stem, short_labels, lab_col=PUR, dot_col=GRN, size_mm=None, lay=None, rows=None, marks=None, legend_note=None, footer=False, yticks=None, row_cost=25):
    by_pos = {}
    for v in variants: by_pos.setdefault(int(v[1:-1]), []).append(v)
    FS = dict(lab=6.5, ax=7.5, tick=6.5, leg=6, note=5.5)
    if size_mm:   # fixed canvas in mm: the axes box is placed in mm and the file is saved without tight cropping
        Wm, Hm = size_mm; FS = dict(lab=6, ax=6.5, tick=6, leg=5.5, note=5.5); Lm, Rm, Bm, AHm = lay or (12, 2.5, 12.5, 44)
        fig = plt.figure(figsize=(Wm / 25.4, Hm / 25.4)); ax = fig.add_axes([Lm / Wm, Bm / Hm, (Wm - Lm - Rm) / Wm, AHm / Hm])
    else:
        fig, ax = plt.subplots(figsize=(6.2, 3.1)); fig.subplots_adjust(left=0.10, right=0.98, top=0.72, bottom=0.27)
    ax.set_xlim(1, 522); ax.set_ylim(0, 1.0)
    ax.plot(mp.index, mp.values, color='#555555', lw=1.0, zorder=2)
    ax.axhspan(0, 0.34, color='#EEEEEE', zorder=0)
    if not footer: ax.text(521, 0.02, 'AlphaMissense "likely benign" < 0.34', ha='right', va='bottom', fontsize=FS['note'], color='#777777')
    labels = {pos: ((vs[0][0] + str(pos)) if short_labels else (vs[0][:-1] + '/'.join(v[-1] for v in vs))) for pos, vs in by_pos.items()}
    labels = {p: l + (marks or {}).get(p, '') for p, l in labels.items()}   # dagger marks: M. myotis-only / M. nigricans-only substitutions (as in Fig 3D)
    if place is None: place = auto_place(ax, sorted(labels.items()), rows=(rows or ((1.17, 1.30, 1.43) if size_mm else (1.17, 1.30, 1.43, 1.56))), fs=FS['lab'], row_cost=row_cost)
    for pos, vs in sorted(by_pos.items()):
        lab = labels[pos]; lx, y = place[pos]
        ax.plot([pos, pos], [0, 1.0], color=lab_col, lw=0.9, ls=(0, (3, 2)), zorder=1, clip_on=False)
        ax.plot([pos, lx], [1.0, y - 0.015], color=lab_col, lw=0.5, clip_on=False)
        ax.text(lx, y, lab, ha='center', va='bottom', fontsize=FS['lab'], color=lab_col, clip_on=False)
        for v in vs: ax.scatter([pos], [T16.loc[v, 'Variant_pathogenicity']], s=16, color=dot_col, edgecolor='black', lw=0.5, zorder=4)
    ax.scatter([], [], s=16, color=dot_col, edgecolor='black', lw=0.5, label='variant score')
    ax.plot([], [], color='#555555', lw=1.0, label='mean over all possible substitutions at the position')
    if legend_note: ax.plot([], [], ' ', label=legend_note)
    if size_mm and not footer: fig.legend(*ax.get_legend_handles_labels(), loc='lower center', bbox_to_anchor=(0.5, 0.4 / Hm), fontsize=FS['leg'], frameon=False, ncol=(3 if legend_note else 2), borderaxespad=0, handletextpad=0.4, columnspacing=1.5)
    elif not size_mm: ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), fontsize=FS['leg'], frameon=False, ncol=2)
    ax.set_xlim(1, 522); ax.set_ylim(0, 1.0); ax.set_xticks([1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 522])
    ax.set_xlabel('Residue sequence number', fontsize=FS['ax'], labelpad=1.5 if size_mm else 4, loc='left' if footer else 'center'); ax.set_ylabel('Mean pathogenicity', fontsize=FS['ax'], labelpad=2 if size_mm else 4); ax.tick_params(labelsize=FS['tick'], pad=2 if size_mm else 3.5)
    for k in ['top', 'right']: ax.spines[k].set_visible(False)
    if yticks is not None: ax.set_yticks(yticks)
    if footer:   # one footer row: axis label at the left, legend at the right, on the same line (for very short canvases)
        fig.canvas.draw(); bb = ax.xaxis.label.get_window_extent(); yc = (bb.y0 + bb.y1) / 2 / fig.bbox.height
        hh, ll = ax.get_legend_handles_labels(); hh.append(Patch(facecolor='#EEEEEE', edgecolor='#BBBBBB', lw=0.4)); ll.append('AlphaMissense "likely benign" < 0.34')   # the grey band, named in the legend instead of printed on the data
        fig.legend(hh, ll, loc='center right', bbox_to_anchor=((Wm - Rm) / Wm, yc), fontsize=FS['leg'], frameon=False, ncol=3, borderaxespad=0, handletextpad=0.4, columnspacing=1.5, handlelength=1.6)
    if size_mm: fig.savefig(str(FIGS / stem) + '.pdf', bbox_inches=None); fig.savefig(str(FIGS / stem) + '.png', dpi=300, bbox_inches=None)
    else:
        fig.subplots_adjust(left=0.10, right=0.98, top=0.72, bottom=0.27)
        fig.savefig(str(FIGS / stem) + '.pdf', bbox_inches='tight'); fig.savefig(str(FIGS / stem) + '.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    return fig

myo = sorted(SETS['Myotis myotis'], key=lambda v: int(v[1:-1]))   # main Fig 3 is focused on the M. myotis substitutions
allv = sorted({v for s in SETS.values() for v in s}, key=lambda v: (int(v[1:-1]), v))
# --- Fig 3B, wide version (180 x 60 mm): the 10 positions that vary in the two Myotis (as in Fig 3A / 3D) ---
myo10 = sorted(SETS['Myotis myotis'] + ['P17Q', 'I151V'], key=lambda v: int(v[1:-1]))
AH = 38.5; ROWS = tuple(1 + (2.6 + 2.8 * k) / AH for k in range(2))        # label rows 2.6 / 5.4 mm above the axes box; 9 mm top margin
FIG3B_WIDE = am_plot(myo10, None, 'Fig3B_alphamissense_Myotis_180x60', True, lab_col=vir(0.30), dot_col=vir(1.0), size_mm=(180, 60), lay=(12, 2.5, 12.5, AH), rows=ROWS,
                     marks={17: '\u2021', 64: '\u2020', 119: '\u2020', 151: '\u2021'}, legend_note='\u2020 M. myotis only; \u2021 M. nigricans only')
print(myo10)
print(myo); print(allv)
assert len(allv) == 17 and {'P17Q', 'Q119R', 'I151V'} <= set(allv), allv
# --- supplementary strips, 180 mm wide: all four-species variants (17 variants at 14 positions: the four composites incl. P17Q, Q119R, I151V) ---
# Axis label (left) and legend (right) share one footer row; label rows: 2 if a collision-free layout exists, else 3.
def strip(H, stem, yticks, first=1.8, step=2.6, n_opts=(2, 3), row_cost=25, top_extra=0.0):
    """Label rows start `first` mm above the axes and are `step` mm apart; the axes height is what is left of H. row_cost: how many residues of
    sideways shift one extra label row is worth (a low value keeps labels close to their own dashed line and the leaders steeper)."""
    for n_rows in n_opts:
        top = 0.6 + (first + step * (n_rows - 1)) + 2.1 + top_extra; ah = H - 8.0 - top      # mm: top margin incl. label rows; 8 mm footer below the axes
        try:
            plt.close('all')
            return am_plot(allv, None, stem, False, lab_col=vir(0.30), dot_col=vir(1.0), size_mm=(180, H), lay=(12, 2.5, 8.0, ah),
                           rows=tuple(1 + (first + step * k) / ah for k in range(n_rows)), footer=True, yticks=yticks, row_cost=row_cost), ah, n_rows
        except AssertionError as e: print('rows', n_rows, 'failed:', e)
    raise RuntimeError('no label layout found')
assert len(allv) == 17 and {'P17Q', 'Q119R', 'I151V'} <= set(allv), allv
FIG_S30, AH30, NR30 = strip(30, 'FigS_alphamissense_four_species_180x30', [0, 0.5, 1.0])
FIG_S50, AH50, NR50 = strip(50, 'FigS_alphamissense_four_species_180x50', [0, 0.2, 0.4, 0.6, 0.8, 1.0], first=3.0, step=2.6, n_opts=(2, 3), row_cost=8, top_extra=1.5)
print('180x30: label rows', NR30, 'axes mm', round(AH30, 1), '| 180x50: label rows', NR50, 'axes mm', round(AH50, 1))

PSMF1_SITES = ['D191G', 'L290P', 'E292R']          # PSMF1-interface substitutions; both Myotis retain the human residue at these positions

# --- version B: the Myotis panel plus the three PSMF1-interface positions ---
myo_psmf1 = sorted(set(SETS['Myotis myotis']) | set(SETS['Myotis nigricans']) | set(PSMF1_SITES), key=lambda v: int(v[1:-1]))
am_plot(myo_psmf1, None, 'Fig3B_alphamissense_myotis_psmf1', True, lab_col=vir(0.30), dot_col=vir(1.0), size_mm=(110, 70),
        marks={17: '\u2021', 64: '\u2020', 119: '\u2020', 151: '\u2021', 191: '\u00a7', 290: '\u00a7', 292: '\u00a7'},
        legend_note=None)   # 110 mm is too narrow for a third legend column; the dagger marks are defined in the caption

# --- version A: every substitution in the four focal species ---
AH4 = 44.0; ROWS4 = tuple(1 + (2.6 + 3.0 * k) / AH4 for k in range(2))     # two label rows above the axes box
am_plot(allv, None, 'Fig3B_alphamissense_four_species', True, lab_col=vir(0.30), dot_col=vir(1.0), size_mm=(110, 78),
        lay=(12, 2.5, 16.0, AH4), rows=ROWS4,
        legend_note=None)
