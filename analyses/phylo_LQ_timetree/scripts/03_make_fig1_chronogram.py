"""Figure 1 of the bat-longevity manuscript (180 x 190 mm).

Panel A  time-calibrated phylogeny of the 30 species with genome assemblies, drawn as a chronogram
         (branch lengths in Myr, root at the left, the present at the right) over a geologic time
         scale in the official ICS/CGMW epoch colours.
Panel B  longevity quotient (LQ) per species, log scale, aligned row-for-row with panel A.
Panel C  maximum lifespan against body mass for 89 bats and 582 other mammals, with the non-flying
         mammal expectation line used to compute LQ; ribbons link the four focal species to panel B.
Panel D  the two long/short-lived species pairs.
Panel E  the nine hallmarks of aging (third-party graphic; see README).

Run from anywhere:  python scripts/03_make_fig1_chronogram.py
Requires numpy, pandas, matplotlib, pillow and openpyxl. Two third-party image inputs are NOT in
this repository (see README): data/images/hallmarks_of_aging.png and the four bat paintings in
data/images/. The script stops with a clear message if they are missing.
"""
import json, re, collections, pathlib, sys
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image

HERE = pathlib.Path(__file__).resolve().parent
DATA, FIGS, IMG = HERE.parent / 'data', HERE.parent / 'figures', HERE.parent / 'data' / 'images'
FIGS.mkdir(exist_ok=True)
OUT_STEM = str(FIGS / 'fig1_chronogram_180x190')
PATHS = {
    'timetree': DATA / 'forcedtop_longevity_timetree.tre',   # dated tree (relaxed clock, fixed topology) -> panel A
    'tree30':   DATA / 'tree30.nwk',                         # IQ-TREE supermatrix tree, for the topology cross-check only
    'names':    DATA / 'species_names.tsv',                  # assembly id -> display name, bat/outgroup
    'ics':      DATA / 'ics_intervals.json',                 # ICS/CGMW interval colours (see 00_export_ics_intervals.R)
    'lq':       DATA / 'bat103_max_lifespan.xlsx',           # MLS_LQ_table sheet: lifespan, mass, LQ
    'bats_all': DATA / 'bat_lifespan_all_bats.csv',          # every bat with a lifespan record -> panel C
    'mm':       DATA / 'data_lifespan_body_mammals.csv',     # 582 non-bat mammals -> panel C
    'whale':    IMG / 'phylopic_whale_93d5c0be.png',         # CC0 PhyloPic
    'mouse':    IMG / 'phylopic_mouse_fe09db98.png',          # CC0 PhyloPic
    'hall':     IMG / 'hallmarks_of_aging_v1.png',           # third-party, not redistributed (same file script 02 uses)
    'Myotis myotis':     IMG / 'Myotis_myotis.jpg',           # third-party paintings, not redistributed
    'Myotis nigricans':  IMG / 'Myotis_nigricans.jpg',
    'Desmodus rotundus': IMG / 'Desmodus_rotundus.jpg',
    'Diphylla ecaudata': IMG / 'Diphylla_ecaudata.jpg',
}
_missing = [str(p.relative_to(HERE.parent)) for p in PATHS.values() if not pathlib.Path(p).exists()]
if _missing:
    sys.exit('Missing input files (see README for the ones not redistributed here):\n  ' + '\n  '.join(_missing))
PATHS = {k: str(v) for k, v in PATHS.items()}

def apply_house_style(sizes=(9, 8, 7)):
    """House figure defaults, inlined so this script has no dependency outside PyPI."""
    t, l, k = sizes
    mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                         'font.size': t, 'axes.titlesize': t, 'axes.labelsize': t, 'legend.fontsize': l,
                         'xtick.labelsize': k, 'ytick.labelsize': k, 'axes.linewidth': 0.6, 'patch.linewidth': 0.6,
                         'lines.linewidth': 1.2, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                         'xtick.major.size': 3.0, 'ytick.major.size': 3.0, 'axes.spines.right': False,
                         'axes.spines.top': False, 'axes.titlelocation': 'left', 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'figure.dpi': 200, 'savefig.dpi': 300})

def parse_nwk(s):
    s = re.sub(r'\[[^\]]*\]', '', s.strip().rstrip(';')); pos = 0
    def node():
        nonlocal pos
        if s[pos] == '(':
            pos += 1; kids = [node()]
            while s[pos] == ',': pos += 1; kids.append(node())
            assert s[pos] == ')'; pos += 1
            m = re.match(r"([^:,();]*)(?::([\d.eE+-]+))?", s[pos:]); pos += m.end()
            return {'name': None, 'label': m.group(1) or None, 'bl': float(m.group(2) or 0), 'kids': kids}
        m = re.match(r"([^:,();]+)(?::([\d.eE+-]+))?", s[pos:]); pos += m.end()
        return {'name': m.group(1), 'label': None, 'bl': float(m.group(2) or 0), 'kids': []}
    return node()
def tips(d): return [d['name']] if not d['kids'] else [t for k in d['kids'] for t in tips(k)]
def all_splits(d, allt=None):
    allt = frozenset(tips(d)) if allt is None else allt; out = {}
    def rec(n):
        ts = frozenset(tips(n))
        if n['kids']:
            for k in n['kids']: rec(k)
            if 1 < len(ts) < len(allt) - 1: out[min(ts, allt - ts, key=lambda s: (len(s), sorted(s)))] = n['label']
    rec(d); return out

ICS = json.load(open(PATHS['ics']))          # official ICS/CGMW interval colours, exported from the R deeptime package
apply_house_style(sizes=(9, 8, 7))
# ---------------- data ----------------
T30 = parse_nwk(open(PATHS['tree30']).read()); t30 = tips(T30); allt = frozenset(t30)
TT = parse_nwk(open(PATHS['timetree']).read())                 # dated tree: topology AND node ages for panel a
assert set(tips(TT)) == allt, 'time tree taxa differ from tree30'
_rt = {}
def _walk(n, acc=0.0):
    acc += n['bl']; _rt[frozenset(tips(n))] = acc
    for k in n['kids']: _walk(k, acc)
_walk(TT)
ROOTAGE = max(_rt[frozenset([t])] for t in t30)
assert ROOTAGE - min(_rt[frozenset([t])] for t in t30) < 1e-4, 'the time tree is not ultrametric'
AGE = {s: ROOTAGE - d for s, d in _rt.items()}
RF = len(set(all_splits(T30)) ^ set(all_splits(TT)))                     # topology difference vs the substitution tree
names = pd.read_csv(PATHS['names'], sep='\t'); disp = dict(zip(names.id, names.display_name)); grp = dict(zip(names.id, names.group)); inv = {v: k for k, v in disp.items()}
splits30 = all_splits(TT)    # panel a now follows the dated tree
xl = pd.ExcelFile(PATHS['lq']); lq = xl.parse('MLS_LQ_table'); tbl = lq.set_index('tip_code')
PROXY = {'HLdesRot8A': 'HLdesRot15A', 'HLdipEca1': 'HLdipEca2A', 'HLtadBra2': 'HLtadBra3', 'mm10': 'mm39', 'bosTau9': 'HLbosTau10'}   # same species, other assembly
def trow(i):
    c = PROXY.get(i, i); return tbl.loc[c] if c in tbl.index else None
EX = pd.read_csv(PATHS['bats_all']); EXLQ = dict(zip(EX.species, EX.LQ)); LQ = {}
for i in t30:
    if grp[i] != 'outgroup': LQ[i] = float(EXLQ.get(disp[i], np.nan))      # bats: extended lifespan table (103-bat table + user databases + AnAge adult weight)
    else:
        r = trow(i); LQ[i] = float(r.LQ) if r is not None and pd.notna(r.LQ) else np.nan
mm = pd.read_csv(PATHS['mm']); mm['is_bat'] = mm.species.eq('Lasionycteris_noctivagans'); nb = mm[~mm.is_bat]; assert len(nb) == 582
b = EX[['species', 'body_mass_g', 'MLS_years']].dropna().reset_index(drop=True)   # every bat with a maximum lifespan and a body mass, analysed here or not
FOC = {"Myotis myotis": ("#1F4E79",), "Myotis nigricans": ("#9E2B25",), "Desmodus rotundus": ("#4E8FC0",), "Diphylla ecaudata": ("#D4786E",)}   # blue = long-lived, red = short-lived (matches the synthesis panel)
LONG = {"Myotis myotis", "Desmodus rotundus"}; FOCID = {'Myotis myotis': 'HLmyoMyo6', 'Myotis nigricans': 'HLmyoNig1', 'Desmodus rotundus': 'HLdesRot8A', 'Diphylla ecaudata': 'HLdipEca1'}
F = lq[lq.species.isin(FOC)].set_index('species'); assert set(F.index) == set(FOC)
imgs = {}
for s in FOC:
    im = Image.open(PATHS[s]).convert('RGB'); im.thumbnail((500, 500)); imgs[s] = np.asarray(im)
def sil(path, col="#4D4D4D"):
    a = np.asarray(Image.open(path).convert("LA")).astype(float) / 255
    rgb = np.array(mpl.colors.to_rgb(col)); out = np.zeros(a.shape[:2] + (4,)); out[..., :3] = rgb; out[..., 3] = a[..., 1] * (1 - a[..., 0]) if a[..., 0].mean() < 0.5 else a[..., 1]
    return out
WH = sil(PATHS['whale']); MO = sil(PATHS['mouse'])
HP = PATHS['hall']; HI = None
# ---------------- hallmark icons (as in the 103-bat figure) ----------------
HP=PATHS['hall']
HI=Image.open(HP).convert("RGB"); pass
sx=HI.size[0]/1560; sy=HI.size[1]/1100
boxes={"Genomic instability":((575,215,762,332),(585,130,740,205)),
 "Telomere attrition":((800,210,1008,362),(860,130,1000,205)),
 "Epigenetic alterations":((925,378,1092,490),(1055,318,1220,392)),
 "Loss of proteostasis":((972,500,1152,578),(1098,560,1285,636)),
 "Deregulated nutrient sensing":((878,688,1100,822),(955,805,1140,920)),
 "Altered intercellular communication":((708,752,876,892),(680,900,912,1018)),
 "Mitochondrial dysfunction":((510,718,672,820),(408,818,618,895)),
 "Stem cell exhaustion":((462,508,618,666),(278,553,445,630)),
 "Cellular senescence":((432,388,652,478),(342,312,522,390))}

HA=np.asarray(Image.open(HP).convert("RGBA")).copy()
S=lambda b: (int(b[0]*sx),int(b[1]*sy),int(b[2]*sx),int(b[3]*sy))
labcol={}
for k,(ib,lb) in boxes.items():
    x0,y0,x1,y1=S(lb); reg=HA[y0:y1,x0:x1]; px=reg[reg[...,3]>200][:,:3]
    labcol[k]=mpl.colors.to_hex(np.median(px,0)/255)
    reg[...,3]=0  # erase label text from source
grey=lambda a: (np.abs(a[...,0].astype(int)-a[...,1])<14)&(np.abs(a[...,1].astype(int)-a[...,2])<14)&(a[...,0]>140)&(a[...,0]<225)
icons={}
for k,(ib,lb) in boxes.items():
    x0,y0,x1,y1=S(ib); c=HA[y0:y1,x0:x1].copy()
    if k not in ("Deregulated nutrient sensing",): c[grey(c),3]=0
    else:
        g=grey(c); h,w=g.shape; m=np.zeros_like(g); m[:int(h*0.25),int(w*0.8):]=True; c[g&m,3]=0
    ys,xs__=np.where(c[...,3]>30); c=c[ys.min():ys.max()+1,xs__.min():xs__.max()+1]
    icons[k]=c
pass

c=icons["Epigenetic alterations"]; h,w=c.shape[:2]; c[:int(h*0.2),int(w*0.85):,3]=0; icons["Epigenetic alterations"]=c
c=icons["Deregulated nutrient sensing"]; h,w=c.shape[:2]; g=grey(c); c[:int(h*0.3),int(w*0.55):,3]=np.where(g[:int(h*0.3),int(w*0.55):],0,c[:int(h*0.3),int(w*0.55):,3]); icons["Deregulated nutrient sensing"]=c
ORDER=["Genomic instability","Telomere attrition","Epigenetic alterations","Loss of proteostasis","Deregulated nutrient sensing","Mitochondrial dysfunction","Cellular senescence","Stem cell exhaustion","Altered intercellular communication"]
LBL={"Genomic instability":"Genomic\ninstability","Telomere attrition":"Telomere\nattrition","Epigenetic alterations":"Epigenetic\nalterations","Loss of proteostasis":"Loss of\nproteostasis",
     "Deregulated nutrient sensing":"Deregulated\nnutrient sensing","Mitochondrial dysfunction":"Mitochondrial\ndysfunction","Cellular senescence":"Cellular\nsenescence","Stem cell exhaustion":"Stem cell\nexhaustion","Altered intercellular communication":"Altered intercellular\ncommunication"}
B0,H0=0.20,0.80

from scipy import ndimage
c=icons["Epigenetic alterations"]; h,w=c.shape[:2]
ys,xs_=np.where(c[...,3]>30); 
# drop small connected blobs in the top-right quadrant

lab_,nl=ndimage.label(c[...,3]>30); sizes=ndimage.sum(np.ones_like(lab_),lab_,range(1,nl+1))
for i_,s_ in enumerate(sizes,1):
    yy,xx=np.where(lab_==i_)
    if s_<0.002*h*w and xx.mean()>0.8*w and yy.mean()<0.35*h: c[lab_==i_,3]=0
icons["Epigenetic alterations"]=c

c=icons["Epigenetic alterations"]; h,w=c.shape[:2]
lab_,nl=ndimage.label(c[...,3]>30); sizes=ndimage.sum(np.ones_like(lab_),lab_,range(1,nl+1))
for i_,s_ in enumerate(sizes,1):
    yy,xx=np.where(lab_==i_)
    if xx.min()>0.85*w and yy.max()<0.4*h: c[lab_==i_,3]=0; print("removed",s_)
ys,xs__=np.where(c[...,3]>30); icons["Epigenetic alterations"]=c[ys.min():ys.max()+1,xs__.min():xs__.max()+1]


# ---------------- rooted cladogram of tree30 ----------------
og = frozenset(i for i in t30 if grp[i] == 'outgroup'); bats = allt - og
# ---- panel a is built directly from the dated tree, using its own root (Euarchontoglires | the rest) ----
kids = {}; AGEC = {}
def _struct(n):
    c_ = frozenset(tips(n)); AGEC[c_] = AGE[c_]
    kids[c_] = [frozenset(tips(k)) for k in n['kids']]
    for k in n['kids']: _struct(k)
_struct(TT)
cl = sorted(kids, key=len)
assert all(len(v) in (0, 2) for v in kids.values()), 'the dated tree is not bifurcating'
assert bats in kids, 'Chiroptera is not a clade in the dated tree'
_sk = lambda c_: (0 if c_ <= bats else 1, -len(c_), min(disp[x] for x in c_))
ROWS = []
def emit(c_):
    if len(c_) == 1: ROWS.append(next(iter(c_))); return
    for k in sorted(kids[c_], key=_sk): emit(k)
emit(allt); Y = {i: r for r, i in enumerate(ROWS)}; rowof = dict(FOCID)
Yc = {}
def meas(c_):
    if len(c_) == 1: Yc[c_] = Y[next(iter(c_))]; return
    for k in kids[c_]: meas(k)
    Yc[c_] = float(np.mean([Yc[k] for k in kids[c_]]))
meas(allt)
# ---- x = time elapsed from the root (chronogram) ----
D_ = {c_: ROOTAGE - AGE[c_] for c_ in cl}
assert all(abs(D_[frozenset([t])] - ROOTAGE) < 1e-4 for t in t30)
XMAX = ROOTAGE
PINK = '#E58FA0'; PINKD = '#A63D5A'; DARKL = '#3A3A3A'; BLUE = '#3B528B'; GREY = '#7F7F7F'
def build_tree(ax, W_mm=70.2, label_mm=35.5, ticks=(100, 80, 60, 40, 20, 0)):
    """Chronogram: x = millions of years before present, root at the left, the present at the right.
    Every tip ends at the present, so tip labels align without leader lines."""
    unit = (W_mm - label_mm) / (XMAX * 1.02); mmx = lambda mm_: mm_ / unit
    xlab = XMAX + mmx(1.2)
    ax.set_xlim(-mmx(1.2), xlab + mmx(label_mm))
    for c_ in cl:
        if len(c_) == 1: continue
        x = D_[c_]; ys = [Yc[k] for k in kids[c_]]
        ax.plot([x, x], [min(ys), max(ys)], color=DARKL, lw=0.7, solid_capstyle='butt', zorder=2)
        for k in kids[c_]: ax.plot([x, D_[k]], [Yc[k]] * 2, color=DARKL, lw=0.7, solid_capstyle='butt', zorder=2)
    ax.plot([-mmx(0.8), 0], [Yc[allt]] * 2, color=DARKL, lw=0.7, zorder=2)
    for i in ROWS:
        isf = any(i == v for v in FOCID.values()); nm = disp[i]
        col = (FOC[nm][0] if isf else ('#6E6E6E' if grp[i] == 'outgroup' else '#1A1A1A'))
        ax.text(xlab, Y[i], nm, fontsize=7, style='italic', fontweight='bold' if isf else 'normal', color=col, ha='left', va='center')
    cb = bats; ax.plot([D_[cb]], [Yc[cb]], 'o', ms=3.6, color=BLUE, mec='white', mew=0.4, zorder=5)
    # ---- geologic time scale: official ICS/CGMW colours (deeptime), epochs then periods ----
    y0 = len(ROWS) - 0.2
    for rows, hgt, fs in ((ICS['epochs'], 1.9, 6.0),):           # epochs only; the period tier was redundant
        yb = y0
        for iv in rows:
            xa, xb = ROOTAGE - min(iv['max_age'], ROOTAGE), ROOTAGE - iv['min_age']
            if xb <= 0 or xa >= ROOTAGE: continue
            ax.add_patch(Rectangle((xa, yb), xb - xa, hgt, fc=iv['color'], ec='white', lw=0.35, clip_on=False, zorder=4))
            box_mm = (xb - xa) * unit                                   # pick the longest label that fits the box
            rnd = ax.figure.canvas.get_renderer()
            for lbl in (iv['name'], iv['abbr'], None):
                if lbl is None: break
                t_ = ax.text(0, 0, lbl, fontsize=fs)
                w_mm = t_.get_window_extent(rnd).width / ax.figure.dpi * 25.4
                t_.remove()
                if w_mm <= box_mm - 0.25: break
            if lbl: ax.text((xa + xb) / 2, yb + hgt / 2, lbl, fontsize=fs, ha='center', va='center',
                            color=iv['lab_color'], clip_on=False, zorder=5)
    ytop = y0 + 1.6
    for t_ in ticks:                                                     # age ticks under the scale
        xt = ROOTAGE - t_
        ax.plot([xt, xt], [ytop, ytop + 0.3], color='#7F7F7F', lw=0.5, clip_on=False, zorder=4)
        ax.text(xt, ytop + 1.35, f'{t_:g}', fontsize=7, ha='center', va='center', color='#4D4D4D', clip_on=False)
    ax.text(ROOTAGE / 2, ytop + 2.9, 'Million years ago', fontsize=7.5, ha='center', va='center', color='#4D4D4D', clip_on=False)
    ax.set_ylim(len(ROWS) + 4.6, -0.6); ax.axis('off')
def tree_axes(fig, rect):
    ax = fig.add_axes(rect); build_tree(ax); return ax
def lq_strip(ax):
    ax.axvline(1, color='#BDBDBD', lw=0.6, ls=(0, (2, 2)), zorder=0); rng = np.random.default_rng(1)
    for i in ROWS:
        v = LQ[i]; y = Y[i]
        if np.isnan(v): ax.text(0.34, y, 'n.d.', fontsize=6, color='#9A9A9A', va='center', ha='left'); continue
        nm = disp[i]
        if nm in FOC: ax.plot(v, y, 'o', ms=4, color=FOC[nm][0], mec='white', mew=0.4, zorder=4)
        elif grp[i] == 'outgroup': ax.plot(v, y, 'o', ms=2.4, mfc='white', mec=GREY, mew=0.5)
        else: ax.plot(v, y, 'o', ms=2.4, mfc=GREY, mec=GREY, mew=0.5)
    ax.set_xscale('log'); ax.set_xlim(0.3, 12); ax.set_xticks([0.5, 1, 2, 5, 10]); ax.set_xticklabels(['0.5', '1', '2', '5', '10']); ax.minorticks_off()
    ax.set_ylim(len(ROWS) + 4.6, -0.6); ax.set_yticks([])
    for s in ('left', 'right', 'top'): ax.spines[s].set_visible(False)
    ax.set_xlabel('Longevity quotient (LQ)', fontsize=7.5)
lq_strip2 = lq_strip
def ribbon(ov, x0, y0, x1, y1, col, w=0.0055, alpha=0.45):
    xm = (x0 + x1) / 2
    verts = [(x0, y0 + w), (xm, y0 + w), (xm, y1 + w * 0.5), (x1, y1 + w * 0.5), (x1, y1 - w * 0.5), (xm, y1 - w * 0.5), (xm, y0 - w), (x0, y0 - w), (x0, y0 + w)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ov.add_patch(PathPatch(Path(verts, codes), fc=col, ec='none', alpha=alpha, zorder=1))
toF = lambda ax, x, y: ax.figure.transFigure.inverted().transform(ax.transData.transform((x, y)))
PANELS = {'A': (0.013, 0.978, 'Phylogeny'), 'B': (0.424, 0.978, 'Longevity quotient'), 'C': (0.648, 0.978, 'Lifespan ~ body size'),
          'D': (0.628, 0.543, 'Species pairs'), 'E': (0.013, 0.148, 'Hallmarks of Aging')}
def panel_letter(fig, l):
    x, y, title = PANELS[l]
    fig.text(x, y, l, fontsize=8.5, fontweight='bold', ha='center', va='center', color='#1A1A1A',
             bbox=dict(boxstyle='circle,pad=0.32', fc='white', ec='#1A1A1A', lw=0.8))
    fig.text(x + 0.020, y, title, fontsize=9.5, fontweight='bold', ha='left', va='center', color='#1A1A1A')

def render2(fname):
    fig=plt.figure(figsize=(180/25.4,190/25.4)); R=lambda r:[r[0],B0+r[1]*H0,r[2],r[3]*H0]; _ft=fig.text; ftext=lambda x,y,*a,**k: _ft(x,B0+y*H0,*a,**k)
    axT=tree_axes(fig,R([0.035,0.08,0.39,0.87])); axL=fig.add_axes(R([0.43,0.08,0.155,0.87])); lq_strip2(axL)
    for ln in axT.lines:
        if ln.get_color()=="#E8692C": ln.set_color("#3B528B")
    for tx in axT.texts:
        if tx.get_text()=="Chiroptera": tx.set_color("#3B528B")
    ov=fig.add_axes([0,0,1,1],facecolor="none"); ov.set_xlim(0,1); ov.set_ylim(0,1); ov.axis("off"); ov.set_zorder(0.5)
    axS=fig.add_axes(R([0.655,0.50,0.28,0.40])); axS.set_zorder(1); axS.patch.set_alpha(0)
    mx=np.exp(nb.ln_body); my=np.exp(nb.ln_lifespan)
    xs_=np.logspace(0.5,7.9,60); yfit=10**(0.537576+0.189234*np.log10(xs_))
    axS.plot(mx,my,"^",ms=2.4,mfc="none",mec="#B5B5B5",mew=0.45,zorder=1,label=f"other mammals (n = {len(nb)})")
    axS.plot(xs_,yfit,color="#9A9A9A",lw=1.0,zorder=2)
    axS.plot(b.body_mass_g,b.MLS_years,"o",ms=2.8,color="#4D4D4D",mec="white",mew=0.3,zorder=3,label=f"bats with lifespan data (n = {len(b)})")
    for pair in (("Myotis myotis","Myotis nigricans"),("Desmodus rotundus","Diphylla ecaudata")):
        pp=F.loc[list(pair)]; axS.plot(pp.body_mass_g,pp.MLS_years,color=FOC[pair[0]][0],lw=0.9,zorder=4)
    for s,r in F.iterrows(): axS.plot(r.body_mass_g,r.MLS_years,"o",ms=5.5,color=FOC[s][0],mec="white",mew=0.5,zorder=5)
    hum=mm[mm.species=="Homo_sapiens"].iloc[0]; nmr=mm[mm.species=="Heterocephalus_glaber"].iloc[0]
    axS.annotate("Human",(np.exp(hum.ln_body),np.exp(hum.ln_lifespan)),xytext=(-4,0),textcoords="offset points",fontsize=7,ha="right",va="center",color="#4D4D4D")
    axS.plot(np.exp(nmr.ln_body),np.exp(nmr.ln_lifespan),"^",ms=4.2,mfc="white",mec="#4D4D4D",mew=0.8,zorder=7)
    axS.plot(np.exp(hum.ln_body),np.exp(hum.ln_lifespan),"^",ms=4.2,mfc="white",mec="#4D4D4D",mew=0.8,zorder=7)
    axS.annotate("Naked mole-rat",(np.exp(nmr.ln_body),np.exp(nmr.ln_lifespan)),xytext=(400,75),textcoords="data",fontsize=7,ha="left",va="center",color="#4D4D4D",
                 arrowprops=dict(arrowstyle="-|>",lw=0.5,color="#6E6E6E",mutation_scale=5,shrinkA=1,shrinkB=2.5))
    wh=mm.loc[mm.ln_body.idxmax()]; mo=mm.loc[mm.ln_body.idxmin()]
    axS.add_artist(AnnotationBbox(OffsetImage(WH,zoom=0.022),(np.exp(wh.ln_body),np.exp(wh.ln_lifespan)),xybox=(-14,10),boxcoords="offset points",frameon=False,box_alignment=(1.0,0.5),zorder=6))
    axS.add_artist(AnnotationBbox(OffsetImage(MO,zoom=0.011),(np.exp(mo.ln_body),np.exp(mo.ln_lifespan)),xybox=(0,-13),boxcoords="offset points",frameon=False,box_alignment=(0.5,0.5),zorder=6))
    axS.set_xscale("log"); axS.set_yscale("log"); axS.minorticks_off()
    axS.set_xticks([10,1e3,1e5,1e7]); axS.set_xticklabels(["10 g","1 kg","100 kg","10 t"])
    axS.set_yticks([1,3,10,30,100]); axS.set_yticklabels(["1","3","10","30","100"]); axS.set_xlim(2,1e8); axS.set_ylim(1.2,150)
    axS.yaxis.tick_right(); axS.yaxis.set_label_position("right")
    axS.spines["left"].set_visible(False); axS.spines["top"].set_visible(False); axS.spines["right"].set_visible(True)
    axS.set_xlabel("Body mass",fontsize=7.5); axS.set_ylabel("Maximum lifespan (years)",fontsize=7.5)
    axS.legend(loc="lower right",frameon=False,fontsize=7,handletextpad=0.2,borderaxespad=0.1,markerscale=1.4)
    for s in FOC:
        x0,y0=toF(axL,12,Y[rowof[s]]); x1,y1=toF(axS,F.loc[s,"body_mass_g"],F.loc[s,"MLS_years"])
        ribbon(ov,x0,y0,x1-0.004,y1,FOC[s][0],alpha=0.5)
    ov2=fig.add_axes(R([0,0,1,1]),facecolor="none"); ov2.set_xlim(0,1); ov2.set_ylim(0,1); ov2.axis("off"); ov2.set_zorder(2)
    PAIRS = [("Vespertilionidae", "Myotis myotis", "Myotis nigricans", 0.235), ("Phyllostomidae", "Desmodus rotundus", "Diphylla ecaudata", 0.020)]
    IMW, IMH = 0.072, 0.145; XC = 0.787                                  # image width/height and the head-to-head meeting point
    FACE = {"Myotis myotis": "L", "Myotis nigricans": "L", "Desmodus rotundus": "R", "Diphylla ecaudata": "R"}   # direction each portrait faces as supplied
    for fam, sl, ss, cy in PAIRS:
        ov2.text(0.628, cy + 0.152, fam, fontsize=8, fontweight="bold", color=FOC[sl][0], va="bottom")
        for s, side in ((sl, -1), (ss, +1)):
            r = F.loc[s]; c = FOC[s][0]; long_ = s in LONG
            ix = XC - IMW - 0.002 if side < 0 else XC + 0.002
            im_ = imgs[s][:, ::-1] if FACE[s] == ("L" if side < 0 else "R") else imgs[s]   # left member must face right, right member left
            ia = fig.add_axes(R([ix, cy, IMW, IMH])); ia.imshow(im_); ia.axis("off"); ia.set_zorder(2)
            tx = XC - IMW - 0.010 if side < 0 else XC + IMW + 0.010
            ha = "right" if side < 0 else "left"
            ov2.text(tx, cy + 0.130, s, fontsize=7.5, style="italic", va="top", ha=ha)
            ov2.text(tx, cy + 0.104, ("long-lived" if long_ else "short-lived"), fontsize=6.5, va="top", ha=ha, color="#6E6E6E")
            ov2.text(tx, cy + 0.078, f"LQ {r.LQ:.2f}", fontsize=9, fontweight="bold", va="top", ha=ha, color=c)
            ov2.text(tx, cy + 0.046, f"{r.MLS_years:g} years \u00b7 {r.body_mass_g:.3g} g", fontsize=6.5, va="top", ha=ha, color="#4D4D4D")
            bx = tx - 0.085 if side < 0 else tx
            ov2.add_patch(mpl.patches.Rectangle((bx, cy + 0.018), 0.085, 0.006, color=c, lw=0))
    for l in ("A", "B", "C", "D"): panel_letter(fig, l)
    bn=fig.add_axes([0.0,0.0,1.0,0.185]); bn.axis("off"); bn.set_xlim(0,1); bn.set_ylim(0,1)
    bn.plot([0.02,0.98],[0.97,0.97],color="#D0D0D0",lw=0.6)
    panel_letter(fig, "E")
    W=0.96/9
    for i,k in enumerate(ORDER):
        cx=0.02+W*(i+0.5); im_=icons[k]; h,w=im_.shape[:2]
        boxw=W*0.82*180; boxh=0.085*190  # mm
        sc=min(boxw/w,boxh/h); ww=w*sc/180; hh=h*sc/190
        ia=fig.add_axes([cx-ww/2,0.072+(0.085-hh)/2,ww,hh]); ia.imshow(im_,interpolation="antialiased"); ia.axis("off")
        fig.text(cx,0.062,LBL[k],fontsize=6.5,ha="center",va="top",color=labcol[k],fontweight="bold",linespacing=1.05)
    mpl.rcParams['savefig.bbox'] = 'standard'; mpl.rcParams['savefig.pad_inches'] = 0   # exact 180 x 190 mm canvas (house style sets tight cropping)
    fig.savefig(fname+".png",dpi=250,bbox_inches=None); fig.savefig(fname+".pdf",bbox_inches=None)
    return fig

fig = render2(OUT_STEM)
