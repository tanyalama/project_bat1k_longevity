#!/usr/bin/env python3
"""
Overlay RNA-seq differential expression (Myotis myotis vs Mus musculus fibroblasts)
on the KEGG "Mitophagy - animal" map (hsa04137).

Only genes passing BOTH thresholds (|log2FC| >= LFC_THRESH and padj < PADJ_THRESH)
are coloured: downregulated -> viridis purples, upregulated -> viridis greens/yellow.
Measured genes that fail either threshold are grey; genes absent from the dataset
are white. The colour key has a grey band over (-LFC_THRESH, +LFC_THRESH), so no
colour is ever assigned to a value inside the non-significant range.

This is a dependency-light re-implementation of pathview(kegg.native = TRUE,
same.layer = FALSE): KGML and base PNG come from KEGG REST, every gene box is
repainted from its mapped Entrez members, and a label is drawn on top.

Usage:   python scripts/mitophagy_kegg_overlay.py
Inputs:  data/fc_data_padj.tsv   entrez, symbol, log2fc, padj (tab-separated).
         log2fc = DESeq2 log2FoldChange (same values as data/fc_data.tsv);
         padj from Supplementary Table 12, which lists only |log2FC| >= 2 genes,
         so padj is NA for genes below the LFC threshold (they fail regardless).
Outputs: figures/hsa04137.log2fc.viridis.{pdf,png}   pathway diagram
         figures/legend.{pdf,png}                      standalone colour key
         data/node_table.csv                           one row per KEGG gene box
         data/sig_genes_table.csv                      gene-level significant genes + HEX
Requires: matplotlib, numpy, pandas, pillow; network access to rest.kegg.jp (first run).
"""
import os, sys, urllib.request, xml.etree.ElementTree as ET
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.font_manager as fm
from matplotlib.patches import Rectangle
from matplotlib.colors import to_hex, to_rgb
from PIL import Image

# ----------------------------- CONFIG --------------------------------------
PATHWAY     = "hsa04137"
DATA_FILE   = "data/fc_data_padj.tsv"
OUT_STEM    = f"figures/{PATHWAY}.log2fc.viridis"
LFC_THRESH  = 2.0        # significant requires |log2FC| >= this ...
PADJ_THRESH = 0.05       # ... AND padj < this
LIMIT       = 7.0        # |log2FC| at which the colour ramps reach their ends; beyond is clamped
CMAP        = "viridis"
DOWN_RANGE  = (0.00, 0.20)   # viridis fraction used for downregulated: -LIMIT -> -LFC_THRESH (purples)
UP_RANGE    = (0.62, 1.00)   # viridis fraction used for upregulated:  +LFC_THRESH -> +LIMIT (green -> yellow)
GRAY        = "#CCCCCC"  # measured, not significant
WHITE       = "#FFFFFF"  # not in dataset
DATASET     = "M. myotis vs M. musculus fibroblasts (DESeq2)"
# hyphy RELAX tags, limited to the genes named in the Fig S3 caption (FBXO7 is not on hsa04137).
# symbol -> (badge text, selection regime, badge anchor relative to box: (dx, dy) in KEGG px from box centre)
SELECTION_TAGS = {
    "OPA1":  ("k > 1", "intensified", (44, 0)),     # RELAX k = 48.56, LRT p = 0.0153
    "HUWE1": ("k < 1", "relaxed",     (-46, -11)),  # RELAX, Supp Table 8 (BH p < 0.0001)
}
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

kgml  = fetch(f"https://rest.kegg.jp/get/{PATHWAY}/kgml", f"kegg/{PATHWAY}.xml")
png   = fetch(f"https://rest.kegg.jp/get/{PATHWAY}/image", f"kegg/{PATHWAY}.png")
glist = fetch("https://rest.kegg.jp/list/hsa", "kegg/hsa_genes.tsv")

# Entrez -> symbol (KEGG list: hsa:ID \t type \t position \t "SYM1, SYM2; description")
sym = {}
for line in open(glist):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 4 and ";" in p[3]:
        s = p[3].split(";")[0].split(",")[0].strip()
        if s and " " not in s:
            sym[p[0].split(":")[1]] = s

# ---- input + significance ------------------------------------------------------
fc = pd.read_csv(DATA_FILE, sep="\t", dtype={"entrez": str}, na_values=["NA", "", "NaN"])
fc["log2fc"] = pd.to_numeric(fc.log2fc, errors="coerce")
fc["padj"] = pd.to_numeric(fc.padj, errors="coerce")
fc["measured"] = fc.log2fc.notna()
fc["significant"] = (fc.log2fc.abs() >= LFC_THRESH) & (fc.padj < PADJ_THRESH)
# a gene passing the LFC filter must have a padj, otherwise significance is undetermined
undetermined = fc[(fc.log2fc.abs() >= LFC_THRESH) & fc.padj.isna()]
assert undetermined.empty, f"genes with |log2FC| >= {LFC_THRESH} but no padj:\n{undetermined}"
gene = fc.set_index("entrez")

# ---- colour mapping ---------------------------------------------------------------
cmap = plt.get_cmap(CMAP)
def colour(v, sig, measured):
    if not measured:
        return WHITE
    if not sig:
        return GRAY
    f = min((abs(v) - LFC_THRESH) / (LIMIT - LFC_THRESH), 1.0)   # 0 at threshold, 1 at LIMIT
    if v > 0:
        return to_hex(cmap(UP_RANGE[0] + f * (UP_RANGE[1] - UP_RANGE[0])))
    return to_hex(cmap(DOWN_RANGE[1] - f * (DOWN_RANGE[1] - DOWN_RANGE[0])))

# ---- KGML gene boxes -> node state --------------------------------------------------
# KEGG collapses gene families into one box. A box is coloured by its significant
# member with the largest |log2FC| (the "driver") and labelled with that gene;
# grey if any member is measured but none significant; white if none measured.
rows = []
for e in ET.parse(kgml).getroot().findall("entry"):
    g = e.find("graphics")
    if e.get("type") != "gene" or g.get("type") != "rectangle":
        continue
    members = [m.split(":")[1] for m in e.get("name").split()]
    present = [m for m in members if m in gene.index]
    measured = [m for m in present if gene.at[m, "measured"]]
    sig = [m for m in measured if gene.at[m, "significant"]]
    drv = max(sig, key=lambda m: abs(gene.at[m, "log2fc"])) if sig else None
    kegg_label = g.get("name").split(",")[0].strip().rstrip(".")
    clean = [sym[m] for m in members if m in sym and not sym[m].startswith("LOC")]
    state = "significant" if sig else ("measured_ns" if measured else "not_measured")
    val = gene.at[drv, "log2fc"] if drv else np.nan
    rows.append(dict(entry_id=e.get("id"), kegg_label=kegg_label,
                     label=sym.get(drv, kegg_label) if drv else (clean[0] if clean else kegg_label),
                     state=state, members=";".join(members), n_members=len(members),
                     n_measured=len(measured), n_sig=len(sig),
                     sig_members=";".join(f"{sym.get(m, m)}({m}):{gene.at[m, 'log2fc']:.2f}" for m in sig),
                     driver_entrez=drv, driver_symbol=sym.get(drv) if drv else None,
                     node_log2fc=val, node_padj=gene.at[drv, "padj"] if drv else np.nan,
                     hex=colour(val, bool(sig), bool(measured)),
                     x=float(g.get("x")), y=float(g.get("y")),
                     w=float(g.get("width")), h=float(g.get("height"))))
nd = pd.DataFrame(rows)
assert len(nd) > 0, "no gene boxes parsed from KGML"
nd["clamped"] = nd.node_log2fc.abs() > LIMIT
if (nd.n_sig > 1).any():
    print("WARNING: boxes with >1 significant member (coloured by max |log2FC|):\n",
          nd.loc[nd.n_sig > 1, ["kegg_label", "sig_members"]].to_string(index=False), file=sys.stderr)

# ---- fonts --------------------------------------------------------------------------
for cand in FONT_TTFS:
    if os.path.exists(cand):
        fm.fontManager.addfont(cand); plt.rcParams["font.family"] = "Arial"; break
plt.rcParams["pdf.fonttype"] = 42

# ---- colour key: down ramp | grey band | up ramp, plus white "not measured" ----------
def draw_key(fig, bar_rect, sw_rect, fs=9):
    ax = fig.add_axes(bar_rect)
    xs = np.linspace(-LIMIT, LIMIT, 561)
    for a, b in zip(xs[:-1], xs[1:]):
        mid = (a + b) / 2
        sig = abs(mid) >= LFC_THRESH
        ax.axvspan(a, b, color=colour(mid, sig, True), lw=0)
    for t in (-LFC_THRESH, LFC_THRESH):
        ax.axvline(t, color="black", lw=0.6)
    ax.set_xlim(-LIMIT, LIMIT); ax.set_yticks([])
    ticks = [-LIMIT, -LFC_THRESH, 0, LFC_THRESH, LIMIT]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"≤ −{LIMIT:g}", f"−{LFC_THRESH:g}", "0", f"{LFC_THRESH:g}", f"≥ {LIMIT:g}"])
    ax.tick_params(labelsize=fs, length=3, width=0.6)
    for s in ax.spines.values(): s.set_linewidth(0.6)
    ax.set_xlabel("log$_2$ fold change (padj < %g)" % PADJ_THRESH, fontsize=fs, labelpad=3)
    # swatches
    items = [(GRAY, f"Not significant (|log$_2$FC| < {LFC_THRESH:g} or padj ≥ {PADJ_THRESH:g})"),
             (WHITE, "Not measured (absent from dataset or log$_2$FC = NA)")]
    x0, y0, sw, sh, gap = sw_rect
    lx = x0 + sw * 1.35 + (0.035 if sw >= 0.03 else 0.0)   # label x; clears the wider badge
    for i, (c, txt) in enumerate(items):
        y = y0 - i * gap
        sax = fig.add_axes([x0, y, sw, sh]); sax.axis("off")
        sax.add_patch(Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="black", lw=0.6))
        sax.set_xlim(0, 1); sax.set_ylim(0, 1)
        fig.text(lx, y + sh / 2, txt, fontsize=fs, va="center")
    seen = []
    for txt, regime, _ in SELECTION_TAGS.values():
        if (txt, regime) in seen: continue
        seen.append((txt, regime))
        y = y0 - (len(items) + len(seen) - 1) * gap
        fig.text(x0 + sw / 2, y + sh / 2, txt, ha="center", va="center", fontsize=fs * 0.85,
                 fontweight="bold", bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6",
                                              facecolor="white", edgecolor="black", lw=0.9))
        fig.text(lx, y + sh / 2,
                 f"{regime.capitalize()} selection (hyphy RELAX)", fontsize=fs, va="center")

# ---- pathway figure -----------------------------------------------------------------
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
from matplotlib.patches import FancyBboxPatch
def badge(ax, cx, cy, txt, fs=7.5, zorder=4):
    ax.text(cx, cy, txt, ha="center", va="center", fontsize=fs, fontweight="bold", zorder=zorder + 1,
            bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6", facecolor="white",
                      edgecolor="black", lw=0.9))
nd["selection"] = None
for symb, (txt, regime, (dx, dy)) in SELECTION_TAGS.items():
    hit = nd[nd.label == symb]
    assert len(hit) >= 1, f"selection tag gene {symb} not found on map"
    for i, r in hit.iterrows():
        ax.add_patch(Rectangle((r.x - r.w / 2, r.y - r.h / 2), r.w, r.h, facecolor="none",
                               edgecolor="black", lw=2.2, zorder=3.5))     # emphasised outline
        badge(ax, r.x + dx, r.y + dy, txt)
        nd.at[i, "selection"] = f"{regime} ({txt})"
ax.add_patch(Rectangle((2, H - 36), 200, 34, facecolor="white", edgecolor="none", zorder=2))
ax.text(10, H - 24, f"Data on KEGG graph {PATHWAY}", fontsize=8, va="center", zorder=3)
ax.text(10, H - 12, DATASET, fontsize=8, va="center", zorder=3)
draw_key(fig, [0.70, 0.948, 0.285, 0.020], (0.70, 0.880, 0.020, 0.016, 0.0215))
fig.savefig(OUT_STEM + ".pdf", facecolor="white")
fig.savefig(OUT_STEM + ".png", dpi=300, facecolor="white")
plt.close(fig)

# ---- standalone legend --------------------------------------------------------------
lf = plt.figure(figsize=(4.2, 2.1), dpi=300)
draw_key(lf, [0.06, 0.76, 0.88, 0.10], (0.06, 0.45, 0.045, 0.08, 0.13))
lf.savefig("figures/legend.pdf", bbox_inches="tight", facecolor="white")
lf.savefig("figures/legend.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(lf)

# ---- tables -------------------------------------------------------------------------
nd.to_csv("data/node_table.csv", index=False)
on_map = {m for s in nd.members for m in s.split(";")}
sig = fc[fc.significant].copy()
sig["geneSymbol"] = sig.entrez.map(sym).fillna(sig.symbol)
sig["hex"] = [colour(v, True, True) for v in sig.log2fc]
sig["clamped"] = sig.log2fc.abs() > LIMIT
sig["in_map"] = sig.entrez.isin(on_map)
sig = sig.sort_values("log2fc", ascending=False)[["entrez", "geneSymbol", "log2fc", "padj", "hex", "clamped", "in_map"]]
sig.to_csv("data/sig_genes_table.csv", index=False)

print(f"{len(nd)} gene boxes: " + ", ".join(f"{k} {v}" for k, v in nd.state.value_counts().items()))
print(f"{len(fc)} input genes ({fc.measured.sum()} measured); {fc.significant.sum()} significant "
      f"(|log2FC| >= {LFC_THRESH:g} & padj < {PADJ_THRESH:g}); {sig.in_map.sum()} on map; {sig.clamped.sum()} clamped")
