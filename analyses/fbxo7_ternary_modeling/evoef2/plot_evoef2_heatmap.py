#!/usr/bin/env python3
"""
EvoEF2 ddG heatmap for the FBXO7 ternary complex (FBXO7-PINK1kinase-PSMF1).

Per-variant interface binding ddG (EvoEF2 energy units) for the two FBXO7
interfaces (FBXO7-PINK1 = chains A_B; FBXO7-PSMF1 = chains A_C), plus the four
bat species multi-mutant site-sets, as a diverging heatmap.

EvoEF2 protocol: RepairStructure -> ComputeBinding(WT) -> BuildMutant ->
ComputeBinding(mutant); ddG = dG_bind(mut) - dG_bind(WT). EvoEF2 binding energy
is more negative for stronger binding, so:
  * stabilizing ddG (negative) -> green   (viridis @ 0.60, #22a784)
  * destabilizing ddG (positive) -> purple (viridis @ 0.18, #443a83)
  * |ddG| <= FLOOR (0.5) -> white; labels in the band are printed grey.
EvoEF2 has no InterfaceAnalyzer-style published noise floor; the 0.5-unit band
is a conservative label-greying threshold below EvoEF2's discrimination.
Color limit is +/- 6 EvoEF2 units (covers the table max |ddG| ~5.2).
The diverging scale is centered on the semantic zero (ddG = 0).

Input : data/evoef2_ddg_combined.csv (long; columns sample, interface, ddG)
Output: fbxo7_evoef2_ddg_heatmap.png (300 dpi)
Usage : python plot_evoef2_heatmap.py [input_csv] [output_png]

Author: FBXO7 ternary modeling analysis, project_bat1k_longevity.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans",
                                   "Arimo", "DejaVu Sans"]
mpl.rcParams["mathtext.fontset"] = "custom"
mpl.rcParams["mathtext.rm"] = "Arial"
mpl.rcParams["mathtext.it"] = "Arial:italic"
mpl.rcParams["axes.unicode_minus"] = False

FLOOR = 0.5     # |ddG| <= this -> white cell, grey label (EvoEF2 discrimination)
LIM   = 6.0     # symmetric color limit (EvoEF2 units)

INPUT  = sys.argv[1] if len(sys.argv) > 1 else "data/evoef2_ddg_combined.csv"
OUTPUT = sys.argv[2] if len(sys.argv) > 2 else "fbxo7_evoef2_ddg_heatmap.png"

GROUPS = [
    ("PINK1 binding\ndomain", [
        ("S109P", "S109P", False), ("S110H", "S110H", False),
        ("S110C", "S110C", False), ("Q127E", "Q127E", False)]),
    ("PSMF1 binding\ndomain", [
        ("D191G", "D191G", False), ("L290P", "L290P", False),
        ("E292R", "E292R", False)]),
    ("Other FBXO7\nsites", [
        ("T19E", "T19E", False), ("T47A", "T47A", False),
        ("T47E", "T47E", False), ("A64T", "A64T", False),
        ("F146V", "F146V", False), ("G409R", "G409R", False),
        ("G409I", "G409I", False)]),
    ("Species\nsite-sets", [
        ("Myotis_myotis",    "Myotis myotis",     True),
        ("Myotis_nigricans", "Myotis nigricans",  True),
        ("Desmodus_rotundus","Desmodus rotundus", True),
        ("Diphylla_ecaudata","Diphylla ecaudata", True)]),
]

def load_matrix(path):
    df = pd.read_csv(path)
    piv = df.pivot(index="sample", columns="interface", values="ddG")
    labels, italic, M, spans = [], [], [], []
    i = 0
    for gname, members in GROUPS:
        start = i
        for key, disp, ital in members:
            if key not in piv.index:
                raise KeyError(f"sample {key!r} missing from {path}")
            labels.append(disp); italic.append(ital)
            M.append([piv.loc[key, "A_B"], piv.loc[key, "A_C"]])
            i += 1
        spans.append((gname, start, i - 1))
    return labels, italic, np.array(M, float), spans

def build_cmap():
    green  = tuple(mpl.colormaps["viridis"](0.60))
    purple = tuple(mpl.colormaps["viridis"](0.18))
    b = FLOOR / (2 * LIM)
    lo, hi = 0.5 - b, 0.5 + b
    cmap = LinearSegmentedColormap.from_list(
        "vir_div", [(0.0, green), (lo, (1, 1, 1, 1)),
                    (hi, (1, 1, 1, 1)), (1.0, purple)])
    return cmap, Normalize(vmin=-LIM, vmax=LIM), green, purple

def main():
    labels, italic, M, spans = load_matrix(INPUT)
    nrow = len(labels)
    cmap, norm, green, purple = build_cmap()

    fig, ax = plt.subplots(figsize=(4.6, 6.4))
    img = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")

    for i in range(nrow):
        for j in range(2):
            v = M[i, j]
            if abs(v) <= FLOOR:
                tc = "#8a8a8a"
            else:
                rgba = cmap(norm(v))
                lum = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
                tc = "white" if lum < 0.5 else "#222"
            ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                    color=tc, fontsize=8.5)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["FBXO7\u2013PINK1", "FBXO7\u2013PSMF1"])
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.set_yticks(range(nrow)); ax.set_yticklabels(labels)
    for t, it in zip(ax.get_yticklabels(), italic):
        if it:
            t.set_fontstyle("italic")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-.5, nrow, 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.5)
    ax.tick_params(which="minor", length=0)

    xb = 1.62
    for name, a, b2 in spans:
        ax.annotate("", xy=(xb, a - 0.4), xytext=(xb, b2 + 0.4),
                    xycoords=("data", "data"), annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color="#555", lw=1.2))
        ax.text(xb + 0.10, (a + b2) / 2, name, rotation=270,
                va="center", ha="left", fontsize=7.5, color="#444")

    cb = fig.colorbar(img, ax=ax, fraction=0.045, pad=0.30,
                      ticks=[-6, -4, -2, 0, 2, 4, 6])
    cb.set_label("\u0394\u0394G (EvoEF2 units)"); cb.outline.set_visible(False)
    cb.ax.text(0.5, 1.02, "destabilizing", transform=cb.ax.transAxes,
               ha="center", va="bottom", fontsize=6.5, color=purple)
    cb.ax.text(0.5, -0.02, "stabilizing", transform=cb.ax.transAxes,
               ha="center", va="top", fontsize=6.5, color=green)

    ax.set_title("Predicted binding free energy (EvoEF2)",
                 fontsize=10, pad=22, loc="left")
    fig.savefig(OUTPUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUTPUT}  ({nrow} variants x 2 interfaces)")

if __name__ == "__main__":
    main()
