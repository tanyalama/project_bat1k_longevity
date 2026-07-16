#!/usr/bin/env python3
"""
Standalone viridis color-scale legend for the hsa04137 mitophagy log2FC figure.

Renders the gradient key (matching pathview's 20-bin viridis ramp over the data
range) plus a gray "not significant" swatch, titled "log fold change".
All text is Arial, >= 8 pt.

Usage:  python make_legend.py
Output: legend.png (300 dpi)

Edit the CONFIG block to change the range, threshold label, or font size.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
import numpy as np
import os

# ----------------------------- CONFIG --------------------------------------
VMIN, VMAX   = -7.08, 13.08     # color-scale range (match the pathway figure)
MID          = 3.0              # midpoint tick (center of asymmetric range)
SIG_THRESH   = 2                # significance threshold shown in the gray label
GRAY_HEX     = "#CCCCCC"        # matches R gray80 used on non-significant nodes
FONT_SIZE    = 9                # base font size in pt (must be >= 8)
TITLE        = "log fold change"
OUT_PNG      = "legend.png"
# 20-bin viridis ramp (purple -> teal -> yellow), as produced by
# pathview:::colorpanel2(20, low="#440154", mid="#21908C", high="#FDE725")
COLORS = ("#440154,#40115A,#3C2160,#383167,#34416D,#315073,#2D6079,#297080,"
          "#258086,#21908C,#21908C,#399A81,#52A375,#6AAD6A,#83B75E,#9BC053,"
          "#B4CA47,#CCD43C,#E5DD30,#FDE725").split(",")
# ---------------------------------------------------------------------------

# Register Arial if present on the system (falls back to default otherwise).
for cand in ("/System/Library/Fonts/Supplemental/Arial.ttf",
             "/Library/Fonts/Arial.ttf",
             "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"):
    if os.path.exists(cand):
        fm.fontManager.addfont(cand)
        plt.rcParams["font.family"] = "Arial"
        break
plt.rcParams["pdf.fonttype"] = 42   # keep text editable in vector exports

n = len(COLORS)
bounds = np.linspace(VMIN, VMAX, n + 1)
cmap = ListedColormap(COLORS)
norm = BoundaryNorm(bounds, cmap.N)

fig = plt.figure(figsize=(3.6, 1.35), dpi=300)

# Gradient colorbar (top)
cax = fig.add_axes([0.08, 0.55, 0.74, 0.16])
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                  orientation="horizontal", spacing="proportional")
cb.set_ticks([VMIN, MID, VMAX])
cb.set_ticklabels([f"{VMIN:.1f}", f"{MID:.1f}", f"{VMAX:.1f}"])
cb.ax.tick_params(labelsize=FONT_SIZE, length=3, width=0.6)
cb.outline.set_linewidth(0.6)
cb.set_label(TITLE, fontsize=FONT_SIZE, labelpad=4)

# Gray "not significant" swatch (bottom)
sax = fig.add_axes([0.08, 0.10, 0.045, 0.16])
sax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor=GRAY_HEX,
                            edgecolor="black", linewidth=0.6))
sax.set_xlim(0, 1); sax.set_ylim(0, 1); sax.axis("off")
fig.text(0.14, 0.18, f"Not significant (|log2FC| < {SIG_THRESH})",
         fontsize=FONT_SIZE, ha="left", va="center")

fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Wrote {OUT_PNG} (base font {FONT_SIZE} pt, Arial)")
