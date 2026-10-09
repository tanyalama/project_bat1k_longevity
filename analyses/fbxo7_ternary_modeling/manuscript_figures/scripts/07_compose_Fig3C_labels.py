#!/usr/bin/env python3
"""Fig 3C labels: crops the all-surface render, fits it to a canvas of given size, and adds chain labels and one yellow label per visible
M. myotis site, each with a leader line to the measured location of the site patch (Fig3C_anchors.json). Label positions are searched
automatically so that labels sit on the white background, do not overlap, and leader lines do not cross.
Run (after 06_make_Fig3C_render.py; any working directory):  python scripts/07_compose_Fig3C_labels.py [--size_mm W H]
Reads figures/Fig3C_render.png, data/Fig3C_anchors.json, data/Fig3C_chainmask.png; writes figures/Fig3C_ternary_surface.png and .pdf.
The canvas size (default 135 x 100 mm) is a placeholder until the layout slot is fixed.
"""
import argparse, json, math
import numpy as np
from PIL import Image, ImageFilter
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import pathlib
HERE = pathlib.Path(__file__).resolve().parent                      # scripts/
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True)
ap = argparse.ArgumentParser()
ap.add_argument('--out', default='Fig3C_ternary_surface')
ap.add_argument('--size_mm', nargs=2, type=float, default=[135, 100]); ap.add_argument('--dpi', type=int, default=400)
ap.add_argument('--min_px', type=int, default=30, help='sites with fewer visible pixels (in the 2400-px render) are not labelled on the figure')
a = ap.parse_args()
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'pdf.fonttype': 42, 'ps.fonttype': 42})
A = json.load(open(DATA / 'Fig3C_anchors.json')); W0, H0 = A['size']
img = Image.open(FIGS / 'Fig3C_render.png').convert('RGBA'); lab = Image.open(DATA / 'Fig3C_chainmask.png'); assert img.size == (W0, H0) == lab.size
al = np.array(img)[..., 3] > 16; ys, xs = np.nonzero(al); bx0, bx1, by0, by1 = xs.min(), xs.max() + 1, ys.min(), ys.max() + 1
Wc, Hc = [int(round(v / 25.4 * a.dpi)) for v in a.size_mm]
mmpx = a.dpi / 25.4                                     # canvas pixels per mm
mt, mb, mx = 7 * mmpx, 13 * mmpx, 4 * mmpx             # top / bottom (footnote) / side margins
cw, ch = bx1 - bx0, by1 - by0
s = min((Wc - 2 * mx) / cw, (Hc - mt - mb) / ch); nw, nh = int(round(cw * s)), int(round(ch * s))
px0, py0 = int(round((Wc - nw) / 2)), int(round(mt + (Hc - mt - mb - nh) / 2))
mol = img.crop((bx0, by0, bx1, by1)).resize((nw, nh), Image.LANCZOS); labc = np.array(lab.crop((bx0, by0, bx1, by1)).resize((nw, nh), Image.NEAREST))
canvas = Image.new('RGBA', (Wc, Hc), (255, 255, 255, 0)); canvas.paste(mol, (px0, py0), mol)
tf = lambda x, y: ((x - bx0) * s + px0, (y - by0) * s + py0)
sil = np.zeros((Hc, Wc), bool); sil[py0:py0 + nh, px0:px0 + nw] = np.array(mol)[..., 3] > 16
occ = np.array(Image.fromarray((sil * 255).astype(np.uint8)).filter(ImageFilter.MaxFilter(2 * int(0.012 * Wc) + 1))) > 0
cx_mol, cy_mol = px0 + nw / 2, py0 + nh / 2

fig = plt.figure(figsize=(Wc / a.dpi, Hc / a.dpi), dpi=a.dpi); ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(np.array(canvas), extent=(0, Wc, Hc, 0), interpolation='nearest')
ax.set_xlim(0, Wc); ax.set_ylim(Hc, 0); ax.axis('off'); fig.canvas.draw(); rend = fig.canvas.get_renderer()

# ---- label specifications ----
TXT = dict(fbxo7='#6b6b6b', pink1='#3c3c96', psmf1='#3a7a7a'); NAME = dict(fbxo7='FBXO7', pink1='PINK1', psmf1='PSMF1')
def chain_anchor(k, direction):
    yy, xx = np.nonzero(labc == k); d = np.array(direction) / np.linalg.norm(direction); sc = (xx - cx_mol) * d[0] + (yy - cy_mol) * d[1]
    top = np.argsort(-sc)[:300]; mx_, my_ = xx[top].mean(), yy[top].mean(); i = top[int(np.argmin((xx[top] - mx_) ** 2 + (yy[top] - my_) ** 2))]
    return float(xx[i]), float(yy[i]), tuple(d)
specs = []
for o, k, dr in [('pink1', 2, (0.7, -0.7)), ('psmf1', 3, (0.2, 1.0)), ('fbxo7', 1, (1.0, 0.6))]:
    x, y, d = chain_anchor(k, dr); specs.append(dict(kind='chain', key=o, text=NAME[o], anchor=(x, y), pref=math.degrees(math.atan2(d[1], d[0]))))
skipped = []
for v, d in A['sites'].items():
    if d.get('omitted'): continue                    # left off the figure on purpose (render --omit); no label, no footnote
    if d['n'] < a.min_px: skipped.append((v, d['n'])); continue
    ax_, ay_ = tf(d['x'], d['y']); ang = math.degrees(math.atan2(ay_ - cy_mol, ax_ - cx_mol)); specs.append(dict(kind='site', key=v, text=v, anchor=(ax_, ay_), pref=ang))
for sp in specs:
    if sp['kind'] == 'chain': sp['art'] = ax.text(0, 0, sp['text'], ha='center', va='center', fontsize=8, fontweight='bold', color=TXT[sp['key']], zorder=6, bbox=dict(boxstyle='round,pad=0.30,rounding_size=0.45', fc='white', ec=A['colours'][sp['key']], lw=0.8))
    else: sp['art'] = ax.text(0, 0, sp['text'], ha='center', va='center', fontsize=7, fontweight='bold', color='black', zorder=6, bbox=dict(boxstyle='round,pad=0.28,rounding_size=0.45', fc=A['colours']['site'], ec='none'))
fig.canvas.draw(); rend = fig.canvas.get_renderer()
for sp in specs: bb = sp['art'].get_bbox_patch().get_window_extent(rend); sp['w'], sp['h'] = bb.width, bb.height

def seg_inter(p1, p2, p3, p4):
    c = lambda a_, b_, c_: (c_[1] - a_[1]) * (b_[0] - a_[0]) - (b_[1] - a_[1]) * (c_[0] - a_[0])
    return c(p1, p3, p4) * c(p2, p3, p4) < 0 and c(p1, p2, p3) * c(p1, p2, p4) < 0
def seg_box(p, q, box):
    for t in np.linspace(0, 1, 40):
        x, y = p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])
        if box[0] <= x <= box[2] and box[1] <= y <= box[3]: return True
    return False
placed = []; pad = 0.012 * Wc
def search(sp):
    best = None; ax_, ay_ = sp['anchor'] if sp['kind'] == 'site' else (px0 + sp['anchor'][0], py0 + sp['anchor'][1])
    sp['anchor'] = (ax_, ay_)
    for dev in range(0, 181, 10):
        for sgn in ((1,) if dev in (0, 180) else (1, -1)):
            ang = math.radians(sp['pref'] + sgn * dev)
            for r in np.arange(0.03 * Wc, 0.45 * Wc, 0.012 * Wc):
                cx, cy = ax_ + r * math.cos(ang), ay_ + r * math.sin(ang); box = (cx - sp['w'] / 2 - pad, cy - sp['h'] / 2 - pad, cx + sp['w'] / 2 + pad, cy + sp['h'] / 2 + pad)
                if box[0] < 0.01 * Wc or box[2] > 0.99 * Wc or box[1] < 0.01 * Hc or box[3] > Hc - 0.11 * Hc: continue
                if occ[int(box[1]):int(box[3]), int(box[0]):int(box[2])].any(): continue
                if any(not (box[2] < q['box'][0] or box[0] > q['box'][2] or box[3] < q['box'][1] or box[1] > q['box'][3]) for q in placed): continue
                P, Q = (ax_, ay_), (cx, cy)
                if any(seg_inter(P, Q, q['P'], q['Q']) for q in placed) or any(seg_box(P, Q, q['box']) for q in placed): continue
                if any(seg_box(q['P'], q['Q'], box) for q in placed): continue
                ts = np.linspace(0, 1, 60); inside = np.mean([occ[min(Hc - 1, int(ay_ + t * (cy - ay_))), min(Wc - 1, int(ax_ + t * (cx - ax_)))] for t in ts]) * r
                cost = r + 2.5 * inside + 0.8 * dev * (0.01 * Wc)
                if best is None or cost < best[0]: best = (cost, cx, cy, box)
    return best
for sp in sorted(specs, key=lambda q: (q['kind'] != 'chain')):       # chain labels first, then the sites
    b = search(sp); assert b is not None, 'no free position for ' + sp['text']
    _, cx, cy, box = b; sp['pos'] = (cx, cy); placed.append(dict(box=box, P=sp['anchor'], Q=(cx, cy), sp=sp))
for q in placed:
    sp = q['sp']; col = A['colours'][sp['key']] if sp['kind'] == 'chain' else '#222222'; ax_, ay_ = sp['anchor']
    ax.plot([ax_, sp['pos'][0]], [ay_, sp['pos'][1]], color=col, lw=0.6, zorder=4, solid_capstyle='round')
    ax.add_patch(Circle((ax_, ay_), 0.0032 * Wc, fc=col if sp['kind'] == 'chain' else 'black', ec='white', lw=0.4, zorder=5))
    sp['art'].set_position(sp['pos'])
note = 'Not visible from this side: ' + ', '.join(v for v, n in skipped if n == 0) if any(n == 0 for _, n in skipped) else ''
fn = [l for l in [note] if l]
small = [v for v, d in A['sites'].items() if 0 < d['n'] < 200 and not d.get('omitted')]
if small: fn.append(f"{', '.join(small)}: almost fully buried, only a sliver of surface shows")
ax.text(mx, Hc - 0.035 * Hc, '\n'.join(fn), ha='left', va='bottom', fontsize=5.5, color='#555555', linespacing=1.3, zorder=6)
fig.canvas.draw(); rend = fig.canvas.get_renderer()
# ---- QC ----
tb = [(sp['text'], sp['art'].get_bbox_patch().get_window_extent(rend)) for sp in specs]
inside = [t for t, b in tb if b.x0 < 0 or b.x1 > Wc or b.y0 < 0 or b.y1 > Hc]
over = [(t1, t2) for i, (t1, b1) in enumerate(tb) for t2, b2 in tb[i + 1:] if b1.overlaps(b2)]
print('canvas mm', round(Wc / mmpx, 2), 'x', round(Hc / mmpx, 2), '| molecule scale', round(s, 3), '| labels', [t for t, _ in tb], '| outside:', inside, '| overlapping:', over, '| skipped (n px):', skipped, '| footnote:', fn)
fig.savefig(str(FIGS / (a.out + '.png')), dpi=a.dpi, transparent=False, facecolor='white'); fig.savefig(str(FIGS / (a.out + '.pdf')), dpi=a.dpi, facecolor='white')
