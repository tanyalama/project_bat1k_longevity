#!/usr/bin/env python3
"""
Overlay RNA-seq log2 fold changes on the KEGG "Mitophagy - animal" map (hsa04137)
with a viridis colour scale: purple = downregulated, yellow = upregulated.

This is a dependency-light re-implementation of pathview(kegg.native = TRUE,
same.layer = FALSE): the KGML and base PNG are fetched from KEGG REST, every
gene box is recoloured from the log2FC of its mapped Entrez members, and a label
is drawn on top. No Bioconductor annotation packages are needed.

Usage:   python scripts/mitophagy_kegg_overlay.py
Inputs:  data/fc_data.tsv            (<Entrez ID>\t<log2fc>, header row)
Outputs: figures/hsa04137.log2fc.viridis.{pdf,png}   pathway diagram (PDF = vector text/boxes)
         figures/legend.{pdf,png}                      standalone colour key
         data/node_table.csv                           one row per KEGG gene box
         data/sig_genes_table.csv                      gene-level significant genes + HEX
Requires: matplotlib, numpy, pandas, pillow; network access to rest.kegg.jp.
"""
import os, sys, urllib.request, xml.etree.ElementTree as ET
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.font_manager as fm
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize, to_hex, to_rgb
from matplotlib.cm import ScalarMappable
from PIL import Image

# ----------------------------- CONFIG --------------------------------------
PATHWAY     = "hsa04137"
DATA_FILE   = "data/fc_data.tsv"
OUT_STEM    = f"figures/{PATHWAY}.log2fc.viridis"
SIG_THRESH  = 2.0        # |log2FC| below this -> grey (not significant)
LIMIT       = 7.0        # colour scale is symmetric: -LIMIT (purple) .. 0 (teal) .. +LIMIT (yellow);
                         # values beyond +/-LIMIT are clamped to the end colours
CMAP        = "viridis"
GRAY        = "#CCCCCC"  # R gray80: not significant / not measured
DATASET     = "M. myotis vs M. musculus fibroblasts, log2FC"
FONT_TTFS   = ("/System/Library/Fonts/Supplemental/Arial.ttf", "/Library/Fonts/Arial.ttf",
               "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
# ---------------------------------------------------------------------------

root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(root_dir)
os.makedirs("figures", exist_ok=True); os.makedirs("kegg", exist_ok=True)

def fetch(url, dst):
    if not os.path.exists(dst):
        urllib.request.urlretrieve(url, dst)
    return dst

kgml = fetch(f"https://rest.kegg.jp/get/{PATHWAY}/kgml", f"kegg/{PATHWAY}.xml")
png  = fetch(f"https://rest.kegg.jp/get/{PATHWAY}/image", f"kegg/{PATHWAY}.png")
glist = fetch("https://rest.kegg.jp/list/hsa", "kegg/hsa_genes.tsv")

# Entrez -> symbol (KEGG list: hsa:ID \t type \t position \t "SYM1, SYM2; description")
sym = {}
for line in open(glist):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 4 and ";" in p[3]:
        s = p[3].split(";")[0].split(",")[0].strip()
        if s and " " not in s:
            sym[p[0].split(":")[1]] = s

# log2FC input
fc = pd.read_csv(DATA_FILE, sep="\t", dtype={0: str}, na_values=["NA", "", "NaN"])
fc.columns = ["entrez", "log2fc"]; fc["log2fc"] = pd.to_numeric(fc.log2fc, errors="coerce")
fcmap = dict(zip(fc.entrez, fc.log2fc))

# KGML gene boxes -> node value. KEGG collapses gene families into one box; a box
# is coloured by the significant member with the largest |log2FC| ("driver") and
# labelled with that driver's symbol, otherwise with KEGG's own box label.
rows = []
for e in ET.parse(kgml).getroot().findall("entry"):
    g = e.find("graphics")
    if e.get("type") != "gene" or g.get("type") != "rectangle":
        continue
    members = [m.split(":")[1] for m in e.get("name").split()]
    vals = {m: fcmap.get(m, np.nan) for m in members}
    measured = {m: v for m, v in vals.items() if pd.notna(v)}
    sig = {m: v for m, v in measured.items() if abs(v) >= SIG_THRESH}
    drv = max(sig, key=lambda m: abs(sig[m])) if sig else None
    kegg_label = g.get("name").split(",")[0].strip().rstrip(".")
    # uncoloured boxes: first member with a clean (non-LOC) symbol, else KEGG's label
    clean = [sym[m] for m in members if m in sym and not sym[m].startswith("LOC")]
    default_label = clean[0] if clean else kegg_label
    rows.append(dict(entry_id=e.get("id"), kegg_label=kegg_label,
                     label=sym.get(drv, kegg_label) if drv else default_label,
                     members=";".join(members), n_members=len(members),
                     n_measured=len(measured), n_sig=len(sig),
                     sig_members=";".join(f"{sym.get(m, m)}({m}):{v:.2f}" for m, v in sig.items()),
                     driver_entrez=drv, driver_symbol=sym.get(drv) if drv else None,
                     node_log2fc=sig[drv] if drv else np.nan,
                     x=float(g.get("x")), y=float(g.get("y")),
                     w=float(g.get("width")), h=float(g.get("height"))))
nd = pd.DataFrame(rows)
assert len(nd) > 0, "no gene boxes parsed from KGML"
multi = nd[nd.n_sig > 1]
if len(multi):
    print(f"WARNING: {len(multi)} boxes have >1 significant member (coloured by max |log2FC|):",
          multi[["kegg_label", "sig_members"]].to_string(index=False), file=sys.stderr)

# colours
cmap = plt.get_cmap(CMAP); norm = Normalize(-LIMIT, LIMIT, clip=True)
colour = lambda v: GRAY if pd.isna(v) else to_hex(cmap(norm(v)))
nd["hex"] = nd.node_log2fc.map(colour)
nd["clamped"] = nd.node_log2fc.abs() > LIMIT

# fonts (Arial if available)
for cand in FONT_TTFS:
    if os.path.exists(cand):
        fm.fontManager.addfont(cand); plt.rcParams["font.family"] = "Arial"; break
plt.rcParams["pdf.fonttype"] = 42

# --- pathway figure ----------------------------------------------------------
base = Image.open(png).convert("RGB"); W, H = base.size
fig = plt.figure(figsize=(W / 72, H / 72), dpi=72)          # 1 pt == 1 KEGG pixel
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.imshow(base, extent=(0, W, H, 0), interpolation="none"); ax.set_xlim(0, W); ax.set_ylim(H, 0)
fig.canvas.draw(); rend = fig.canvas.get_renderer()
for r in nd.itertuples():
    ax.add_patch(Rectangle((r.x - r.w / 2, r.y - r.h / 2), r.w, r.h,
                           facecolor=r.hex, edgecolor="black", lw=0.8, zorder=2))
    lum = np.dot(to_rgb(r.hex), [0.299, 0.587, 0.114])
    t = ax.text(r.x, r.y + 0.5, r.label, ha="center", va="center", fontsize=10,
                color="white" if lum < 0.45 else "black", zorder=3)
    for fs in np.arange(10, 5.5, -0.25):                      # shrink label to fit box
        t.set_fontsize(fs)
        if t.get_window_extent(rend).width <= r.w - 3:
            break
# replace KEGG footer with provenance
ax.add_patch(Rectangle((2, H - 36), 160, 34, facecolor="white", edgecolor="none", zorder=2))
ax.text(10, H - 24, f"Data on KEGG graph {PATHWAY}", fontsize=8, va="center", zorder=3)
ax.text(10, H - 12, DATASET, fontsize=8, va="center", zorder=3)

def draw_key(fig, cax_rect, sw_rect, txt_xy, fs=9):
    cax = fig.add_axes(cax_rect)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([-LIMIT, -LIMIT / 2, 0, LIMIT / 2, LIMIT])
    cb.set_ticklabels([f"≤ −{LIMIT:g}", f"−{LIMIT/2:g}", "0", f"{LIMIT/2:g}", f"≥ {LIMIT:g}"])
    cb.ax.tick_params(labelsize=fs, length=3, width=0.6); cb.outline.set_linewidth(0.6)
    cb.set_label("log$_2$ fold change", fontsize=fs, labelpad=3)
    sax = fig.add_axes(sw_rect); sax.axis("off")
    sax.add_patch(Rectangle((0, 0), 1, 1, facecolor=GRAY, edgecolor="black", lw=0.6))
    fig.text(*txt_xy, f"Not significant (|log$_2$FC| < {SIG_THRESH:g}) or not measured",
             fontsize=fs, va="center")

draw_key(fig, [0.72, 0.945, 0.255, 0.022], [0.72, 0.868, 0.022, 0.022], (0.748, 0.879))
fig.savefig(OUT_STEM + ".pdf", facecolor="white")
fig.savefig(OUT_STEM + ".png", dpi=300, facecolor="white")
plt.close(fig)

# --- standalone legend -------------------------------------------------------
lf = plt.figure(figsize=(3.6, 1.3), dpi=300)
draw_key(lf, [0.08, 0.58, 0.84, 0.16], [0.08, 0.12, 0.05, 0.16], (0.15, 0.20))
lf.savefig("figures/legend.pdf", bbox_inches="tight", facecolor="white")
lf.savefig("figures/legend.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(lf)

# --- tables --------------------------------------------------------------------
nd.to_csv("data/node_table.csv", index=False)
all_members = {m for s in nd.members for m in s.split(";")}
sig = fc[fc.log2fc.abs() >= SIG_THRESH].copy()
sig["geneSymbol"] = sig.entrez.map(sym); sig["hex"] = sig.log2fc.map(colour)
sig["clamped"] = sig.log2fc.abs() > LIMIT; sig["in_map"] = sig.entrez.isin(all_members)
sig = sig.sort_values("log2fc", ascending=False)[["entrez", "geneSymbol", "log2fc", "hex", "clamped", "in_map"]]
# verify every table value against the input file before writing
assert all(np.isclose(sig.log2fc, [fcmap[e] for e in sig.entrez]))
sig.to_csv("data/sig_genes_table.csv", index=False)

print(f"{len(nd)} gene boxes; {nd.node_log2fc.notna().sum()} coloured, "
      f"{((nd.n_sig == 0) & (nd.n_measured > 0)).sum()} measured-but-NS, {(nd.n_measured == 0).sum()} unmeasured")
print(f"{len(fc)} input genes ({fc.log2fc.notna().sum()} with values); {len(sig)} significant; "
      f"{sig.in_map.sum()} of them on the map; {sig.clamped.sum()} clamped at |log2FC| > {LIMIT:g}")
print("wrote", OUT_STEM + ".pdf/.png, figures/legend.pdf/.png, data/node_table.csv, data/sig_genes_table.csv")
