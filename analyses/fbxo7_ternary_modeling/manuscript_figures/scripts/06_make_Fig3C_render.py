#!/usr/bin/env python3
"""Fig 3C render: FBXO7 (chain A) - PINK1 kinase domain (chain B) - PSMF1 (chain C), all-surface, transparent background.
Orientation: PINK1 above PSMF1 (canonical frame from the chain centroids), then rotated by THETA about the vertical axis and PHI about the
horizontal axis. The same camera is used for (1) a flat-colour mask render, from which the visible pixels of each M. myotis site and of each
chain are measured, and (2) the final ray-traced render. Outputs: <out>_render.png (RGBA) and <out>_anchors.json.
Run (any working directory; no arguments needed):  python scripts/06_make_Fig3C_render.py
Reads data/FBXO7_PINK1kinase_PSMF1_forcedtocif_model_0.cif; writes figures/Fig3C_render.png (RGBA), data/Fig3C_anchors.json and data/Fig3C_chainmask.png.
Defaults reproduce the manuscript figure: theta 90, phi -30, 2400 x 1800 px, F146V left uncoloured (--omit F146V).
Requires PyMOL (conda: pymol-open-source), numpy, Pillow. Keep surface_quality 1 and ray_trace_mode 0: outlined ray mode with higher surface quality did not finish in testing.
"""
import argparse, json, os, tempfile, time
import numpy as np
from PIL import Image
import pymol
from pymol import cmd

import pathlib
HERE = pathlib.Path(__file__).resolve().parent                      # scripts/
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; FIGS.mkdir(exist_ok=True)
ap = argparse.ArgumentParser()
ap.add_argument('--cif', default=str(DATA / 'FBXO7_PINK1kinase_PSMF1_forcedtocif_model_0.cif'))
ap.add_argument('--theta', type=float, default=90.0); ap.add_argument('--phi', type=float, default=-30.0)
ap.add_argument('--w', type=int, default=2400); ap.add_argument('--h', type=int, default=1800)
ap.add_argument('--quality', type=int, default=1); ap.add_argument('--antialias', type=int, default=2); ap.add_argument('--buffer', type=float, default=1.0)
ap.add_argument('--omit', nargs='*', default=['F146V'], help='variant names (e.g. F146V) to leave uncoloured and unlabelled in the final figure; they are still measured')
a = ap.parse_args()
T0 = time.time(); log = lambda m: print('[%5.1fs] %s' % (time.time() - T0, m), flush=True)

SITES = {19: 'T19E', 47: 'T47A', 64: 'A64T', 109: 'S109P', 119: 'Q119R', 127: 'Q127E', 146: 'F146V', 409: 'G409I'}   # M. myotis, human FBXO7 numbering
CHAIN = {'fbxo7': 'A', 'pink1': 'B', 'psmf1': 'C'}
COL = {'fbxo7': '#bcbcbc', 'pink1': '#4f4fa6', 'psmf1': '#5f9696', 'site': '#f5e51b'}
MASKRGB = {'fbxo7': (150, 150, 150), 'pink1': (20, 20, 120), 'psmf1': (20, 120, 120)}
SITERGB = {19: (255, 0, 0), 47: (0, 255, 0), 64: (0, 0, 255), 109: (255, 0, 255), 119: (0, 255, 255), 127: (255, 128, 0), 146: (128, 0, 255), 409: (255, 255, 0)}

pymol.finish_launching(['pymol', '-cq'])
cmd.reinitialize(); cmd.load(a.cif, 'cplx')
for o, ch in CHAIN.items(): cmd.create(o, f'cplx and chain {ch}')
cmd.delete('cplx')
X0 = {o: np.array(cmd.get_coords(o), dtype=float) for o in CHAIN}
up = X0['pink1'].mean(0) - X0['psmf1'].mean(0); up /= np.linalg.norm(up)
ref = np.array([1.0, 0, 0]) if abs(up[0]) < .9 else np.array([0, 1.0, 0]); xax = ref - up * (ref @ up); xax /= np.linalg.norm(xax); zax = np.cross(xax, up)
Rc = np.vstack([xax, up, zax]); origin = np.vstack(list(X0.values())).mean(0)
Ry = lambda t: np.array([[np.cos(np.radians(t)), 0, np.sin(np.radians(t))], [0, 1, 0], [-np.sin(np.radians(t)), 0, np.cos(np.radians(t))]])
Rx = lambda t: np.array([[1, 0, 0], [0, np.cos(np.radians(t)), -np.sin(np.radians(t))], [0, np.sin(np.radians(t)), np.cos(np.radians(t))]])
R = Rx(a.phi) @ Ry(a.theta) @ Rc
for o in X0: cmd.load_coords(((X0[o] - origin) @ R.T).astype('float32'), o, state=1)
cmd.reset(); cmd.zoom('all', complete=1, buffer=a.buffer)
site_sel = 'fbxo7 and resi ' + '+'.join(str(p) for p, v in SITES.items() if v not in a.omit)   # yellow patches in the final render
for o, h in COL.items(): cmd.set_color('c_' + o, [int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)])
for o, rgb in MASKRGB.items(): cmd.set_color('k_' + o, [v / 255 for v in rgb])
for p, rgb in SITERGB.items(): cmd.set_color(f'k_{p}', [v / 255 for v in rgb])

# ---- (1) mask render: flat colours, no lighting, no anti-aliasing ----
for k, v in [('ray_shadows', 0), ('antialias', 0), ('ambient', 1.0), ('direct', 0.0), ('reflect', 0.0), ('specular', 0), ('spec_reflect', 0.0), ('depth_cue', 0),
             ('ray_opaque_background', 1), ('transparency', 0), ('surface_quality', a.quality), ('orthoscopic', 1), ('ray_trace_mode', 0)]: cmd.set(k, v)
cmd.bg_color('black'); cmd.hide('everything'); cmd.show('surface', 'fbxo7 or pink1 or psmf1')
for o in CHAIN: cmd.color('k_' + o, o)
for p in SITES: cmd.color(f'k_{p}', f'fbxo7 and resi {p}')
f = tempfile.mktemp(suffix='.png'); cmd.png(f, a.w, a.h, ray=1); M = np.array(Image.open(f).convert('RGB')).astype(int); os.remove(f); log('mask render done')
lab = np.zeros(M.shape[:2], np.uint8)
for k, o in enumerate(['fbxo7', 'pink1', 'psmf1'], 1): lab[np.abs(M - np.array(MASKRGB[o])).sum(-1) < 30] = k
for p, rgb in SITERGB.items(): lab[np.abs(M - np.array(rgb)).sum(-1) < 40] = 1          # site patches belong to FBXO7
Image.fromarray(lab).save(str(DATA / 'Fig3C_chainmask.png'))   # 0 background, 1 FBXO7, 2 PINK1, 3 PSMF1

def pick(rgb, tol=40):
    m = np.abs(M - np.array(rgb)).sum(-1) < tol; n = int(m.sum())
    if not n: return dict(n=0)
    yy, xx = np.nonzero(m); cx, cy = xx.mean(), yy.mean(); i = int(np.argmin((xx - cx) ** 2 + (yy - cy) ** 2))
    return dict(n=n, cx=float(cx), cy=float(cy), x=int(xx[i]), y=int(yy[i]))
out = dict(size=[a.w, a.h], theta=a.theta, phi=a.phi, sites={SITES[p]: dict(pos=p, omitted=SITES[p] in a.omit, **pick(rgb)) for p, rgb in SITERGB.items()},
           chains={o: pick(rgb, 30) for o, rgb in MASKRGB.items()}, colours=COL)

# ---- (2) final render ----
for k, v in [('ray_shadows', 1), ('antialias', a.antialias), ('ambient', 0.5), ('direct', 0.5), ('reflect', 0.45), ('specular', 0.1), ('spec_reflect', 0.15),
             ('shininess', 15), ('ray_opaque_background', 0), ('ray_trace_mode', 0)]: cmd.set(k, v)
cmd.bg_color('white'); cmd.hide('everything'); cmd.show('surface', 'fbxo7 or pink1 or psmf1')
for o in CHAIN: cmd.color('c_' + o, o)
cmd.color('c_site', site_sel)
cmd.png(str(FIGS / 'Fig3C_render.png'), a.w, a.h, ray=1, dpi=300); log('final render done')
json.dump(out, open(DATA / 'Fig3C_anchors.json', 'w'), indent=1)
for k, v in out['sites'].items(): print(k, v['n'], end=' | ')
print(); cmd.quit()
