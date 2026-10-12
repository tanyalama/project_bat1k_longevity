#!/usr/bin/env python3
"""Build the 103-bat phylogeny / longevity-quotient figure (180 x 190 mm canvas).

Panels: a, T2 timetree (Morales et al. 2026, Nature 658:141) collapsed to families, with
95% HPD node-age bands and bootstrap support < 100%; b, longevity quotient (LQ) per row;
c, body mass vs maximum lifespan (582 non-flying mammals + 40 study bats) with ribbons to
the four focal species; d, focal-species cards; e, hallmarks-of-aging banner.

Run from the folder that contains data/ and figures/:
    python scripts/02_make_phylo_LQ_figure.py
Writes figures/phylo_LQ_hallmarks_180x190.pdf and .png.
Requires numpy, pandas, matplotlib, dendropy, pillow, scipy, openpyxl.
Inputs are listed in README.md. Node uncertainty comes from
scripts/01_extract_node_uncertainty.py (data/T2_node_support_age_uncertainty.csv).
Note: assembled from the interactive session that produced the approved figure; it
reproduces that figure pixel-for-pixel.
"""
import os, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib as mpl, matplotlib.pyplot as plt
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
PATHS={'130b90bd-6878-4397-b04a-00d2a0ec6a5d': 'data/T2_internal_nodes_retention.tsv', '4d77b5b9-9b73-4fc2-8587-cc44100061fe': 'data/bat103_max_lifespan.xlsx', '5187b13d-bd02-438c-adbb-f8e3917fb1b5': 'data/T2_timetree.labeled.nwk', '58c75f9b-ce2a-4712-b80a-61e8af592f2a': 'data/images/hallmarks_of_aging_v1.png', '688f2ef8-a0a9-49c3-85bb-815fa8762225': 'data/data_lifespan_body_mammals.csv', '6ac0e1d4-27a3-419a-811a-804104713f62': 'data/images/Desmodus_rotundus.jpg', '8472d540-5e52-41ec-b43d-9d8f658c8a29': 'data/images/Myotis_nigricans.jpg', 'b178d461-9954-4918-98bc-10f70cd6e259': 'data/images/Diphylla_ecaudata.jpg', 'df706c40-0150-400c-bb0e-6ff055c6db61': 'data/images/Myotis_myotis.jpg'}
def artifact_path(v): return PATHS[v]

# ---- house figure style (Arial; publication rules)
META_GREY='#888888'
MAGMA_END=0.85
DIVERGING_ANCHORS=['#2C3E8C', '#7B9FD4', '#F2F2F2', '#F4665C', '#FEB77E']
def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    """Set matplotlib rcParams for publication-grade output. Call once before plotting.

    This sets mechanics (role-mapped font-size ladder, outward ticks, frameless
    legends, 300-dpi save, Type-42 embedded fonts) — not a house aesthetic.
    Frame, font and the size ladder are parameters.

    frame : 'open' (bottom+left spines, default) | 'boxed' (all four) | 'none'
    font  : sans-serif family name; None = system default sans-serif
    sizes : (base, secondary, tick) — titles/axis-labels, legend/annotation, ticks
    grid  : whether to draw axes.grid (default False)
    """
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")

    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    """§5.7: bold panel letter outside top-left of axes. case ∈ {'lower','upper'}."""
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


def set_frame(ax, style="open"):
    """§3: set spine visibility on an existing axes. style ∈ {'open','boxed','none'}."""
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def house_cmap(name="magma_house", end=None):
    """This user's continuous colormap: magma truncated before its pale-yellow
    high end (default end=0.85, top colour ~#FEB77E). Registers the result
    under `name` so it can also be referenced as a string, and returns the
    Colormap object."""
    import numpy as np
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    if end is None:
        end = MAGMA_END
    base = mpl.colormaps["magma"]
    cmap = LinearSegmentedColormap.from_list(name, base(np.linspace(0.0, end, 256)))
    try:
        mpl.colormaps.register(cmap, name=name, force=True)
    except Exception:
        pass
    return cmap


def house_diverging_cmap(name="magma_diverging"):
    """This user's DIVERGING colormap, for scaled or zero-centred data
    (z-scores, log ratios, fold changes, differences, residuals): deep blue
    -> near-white -> the same peach the truncated magma tops out at. Use this
    instead of magma whenever zero is meaningful, and always centre it (e.g.
    matplotlib TwoSlopeNorm(vcenter=0), or symmetric vmin/vmax) -- an
    off-centre diverging map puts the neutral colour on a non-zero value.
    Registers under `name` and returns the Colormap object."""
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(name, DIVERGING_ANCHORS)
    try:
        mpl.colormaps.register(cmap, name=name, force=True)
    except Exception:
        pass
    return cmap


def apply_house_style(*args, **kwargs):
    """Call figure-style's apply_figure_style() then apply this user's house
    defaults on top: Arial font family, magma truncated at end=0.85
    (registered as "magma_house") as the default continuous colormap, and the
    house diverging map (registered as "magma_diverging") for scaled or
    zero-centred data. Returns whatever apply_figure_style returns."""
    import matplotlib as mpl
    result = None
    try:
        result = apply_figure_style(*args, **kwargs)
    except NameError:
        # figure-style not loaded in this kernel; still apply house defaults
        pass
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    house_cmap()
    house_diverging_cmap()          # registers "magma_diverging" for scaled data
    mpl.rcParams["image.cmap"] = "magma_house"
    return result


# ---- recovered block 30
import pandas as pd
X=artifact_path("4d77b5b9-9b73-4fc2-8587-cc44100061fe")
xl=pd.ExcelFile(X); print(xl.sheet_names)
for s in xl.sheet_names[:4]:
    d=xl.parse(s); print("==",s,d.shape); print(d.head(3).to_string()[:1500]); print(list(d.columns))
nodes=pd.read_csv(artifact_path("130b90bd-6878-4397-b04a-00d2a0ec6a5d"),sep="\t")

# ---- recovered block 31
import dendropy, numpy as np, re
lq=xl.parse("MLS_LQ_table")
nwk=open(artifact_path("5187b13d-bd02-438c-adbb-f8e3917fb1b5")).read()
t=dendropy.Tree.get(data=nwk,schema="newick",preserve_underscores=True,rooting="force-rooted")
age=dict(zip(nodes.node_label,nodes.age_T2_timetree*100))
for n in t.preorder_node_iter():
    n.age = 0.0 if n.is_leaf() else age[n.label]
info=lq.set_index("tip_code")
focal={"HLmyoMyo6":"long","HLmyoNig1":"short","HLdesRot15A":"long","HLdipEca2A":"short"}
fam=lambda tip: info.loc[tip,"family"]
grp=lambda tip: info.loc[tip,"family"] if info.loc[tip,"group"]!="outgroup" else "OG:"+info.loc[tip,"species"]
# collapse: maximal clades with no focal tip and single family (outgroups: never merged across species)
def tips(n): return [l.taxon.label for l in n.leaf_iter()]
collapsed=[]  # (node, label)
def walk(n):
    ts=tips(n)
    fams={grp(x) for x in ts}
    if not (set(ts)&set(focal)) and len(fams)==1:
        collapsed.append(n); return
    for k in n.child_nodes(): walk(k)
walk(t.seed_node)
rows=[]
for n in collapsed:
    ts=tips(n); f=grp(ts[0]); famsize=(lq.family==fam(ts[0])).sum() if not f.startswith("OG:") else 1
    rows.append((f,len(ts),famsize,[info.loc[x,"species"] for x in ts][:3]))

# ---- recovered block 32
sp=lambda x: info.loc[x,"species"]
code_of={v:k for k,v in info.species.items()}
def tipset(n): return frozenset(l.taxon.label for l in n.leaf_iter())
clades={tipset(n):n for n in t.postorder_node_iter()}
def by_family(f): return frozenset(info.index[info.family==f])
groups=[]
# outgroups
groups+= [("Elephas maximus",frozenset(["HLeleMax2"]),"og"),("Euarchontoglires (human, mouse)",frozenset(["hg38","mm39"]),"og"),
          ("Sorex fumeus",frozenset(["HLsorFum1"]),"og"),("Fereuungulata (4)",frozenset(["HLequAsi3","HLbosTau10","HLmanPen4","HLneoVis2"]),"og")]
vesp=by_family("Vespertilionidae"); myo=frozenset(x for x in vesp if sp(x).startswith("Myotis"))
phy=by_family("Phyllostomidae"); des=frozenset(code_of[s] for s in ["Desmodus rotundus","Diaemus youngii","Diphylla ecaudata"])
special={"Vespertilionidae":[("Myotis",myo,"focal"),("other Vespertilionidae",vesp-myo,"bat")],
         "Phyllostomidae":[("Macrotus waterhousii",frozenset([code_of["Macrotus waterhousii"]]),"bat"),("Micronycteris megalotis",frozenset([code_of["Micronycteris megalotis"]]),"bat"),
                           ("Desmodontinae",des,"focal"),("other Phyllostomidae",phy-des-{code_of["Macrotus waterhousii"],code_of["Micronycteris megalotis"]},"bat")]}
for f in lq.loc[lq.group!="outgroup","family"].unique():
    if f in special: groups+=special[f]
    else: groups.append((f,by_family(f),"bat"))
bad=[g[0] for g in groups if g[1] not in clades]

# ---- recovered block 33
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Polygon
apply_house_style(sizes=(7,6,5.5))
gmap={g[1]:g for g in groups}
t.ladderize(ascending=False)
order=[]
def emit(n):
    ts=tipset(n)
    if ts in gmap: order.append((n,gmap[ts])); return
    for k in n.child_nodes(): emit(k)
emit(t.seed_node)
Y={}
for i,(n,g) in enumerate(order): Y[n]=i
def ypos(n):
    if n in Y: return Y[n]
    ys=[ypos(k) for k in n.child_nodes()]; Y[n]=np.mean(ys); return Y[n]
ypos(t.seed_node)
FOC={"Myotis myotis":("#3B0F70","Myotis pair"),"Myotis nigricans":("#9C7FC4","Myotis pair"),
     "Desmodus rotundus":("#B5367A","Desmodontinae pair"),"Diphylla ecaudata":("#F08FB0","Desmodontinae pair")}
LONG={"Myotis myotis","Desmodus rotundus"}
GREY="#9A9A9A"; DARK="#4D4D4D"
def build_tree(ax):
    # branches
    def draw(n):
        if n in dict((o[0],1) for o in order): return
        for k in n.child_nodes():
            ax.plot([n.age,n.age],[Y[n],Y[k]],color=DARK,lw=0.6,solid_capstyle="butt")
            xend=k.age
            ax.plot([n.age,xend],[Y[k],Y[k]],color=DARK,lw=0.6)
            draw(k)
    ax.plot([t.seed_node.age+4,t.seed_node.age],[Y[t.seed_node]]*2,color=DARK,lw=0.6)
    draw(t.seed_node)
    for n,(lab,ts,kind) in order:
        y=Y[n]; ntip=len(ts)
        col={"og":"#BDBDBD","bat":"#7F7F7F","focal":"#2B2B2B"}[kind]
        if ntip>1:
            h=min(0.42,0.12+0.03*ntip)
            ax.add_patch(Polygon([[n.age,y],[0,y-h],[0,y+h]],closed=True,fc=col,ec="none",alpha=0.9 if kind!="og" else 0.8))
        txt=f"{lab} ({ntip})" if ntip>1 and "(" not in lab else lab
        ital = lab.split()[0][0].isupper() and len(lab.split())>=2 and lab.split()[1][0].islower() and "(" not in lab
        ax.text(-1.5,y,txt,va="center",ha="left",fontsize=6,style="italic" if (ital or lab=="Myotis") else "normal",
                fontweight="bold" if kind=="focal" else "normal",color="#1A1A1A" if kind!="og" else "#6E6E6E")
    ax.set_xlim(t.seed_node.age+5,-38); ax.set_ylim(len(order)-0.4,-0.6)
    for s in ("left","right","top"): ax.spines[s].set_visible(False)
    ax.set_yticks([]); ax.set_xticks([90,60,30,0]); ax.set_xlabel("Million years ago",fontsize=6)
    ax.spines["bottom"].set_bounds(t.seed_node.age,0)
    chi=[n for n in t.preorder_node_iter() if n.label=="anc_Chiroptera"][0]
    ax.plot(chi.age,Y[chi],marker="o",ms=2.5,color="#E8692C",zorder=5)
    ax.text(chi.age+1,Y[chi]-0.55,"Chiroptera",fontsize=5.5,ha="right",color="#E8692C")
def lq_strip(ax):
    ax.axvline(1,color="#BDBDBD",lw=0.6,ls=(0,(2,2)),zorder=0)
    rng=np.random.default_rng(1)
    for n,(lab,ts,kind) in order:
        d=lq[lq.tip_code.isin(ts)&lq.LQ.notna()]
        for _,r in d.iterrows():
            if r.species in FOC: continue
            ax.plot(r.LQ,Y[n]+rng.uniform(-0.18,0.18),"o",ms=2.2,mfc="white" if kind=="og" else GREY,mec=GREY,mew=0.5)
        for _,r in d[d.species.isin(FOC)].iterrows():
            ax.plot(r.LQ,Y[n],"o",ms=4,color=FOC[r.species][0],mec="white",mew=0.4,zorder=4)
        if len(d)==0: ax.text(0.33,Y[n],"–",fontsize=5.5,color="#BDBDBD",va="center",ha="center")
    ax.set_xscale("log"); ax.set_xlim(0.3,12); ax.set_xticks([0.5,1,2,5,10]); ax.set_xticklabels(["0.5","1","2","5","10"])
    ax.set_ylim(len(order)-0.4,-0.6); ax.set_yticks([])
    for s in ("left","right","top"): ax.spines[s].set_visible(False)
    ax.set_xlabel("Longevity quotient (LQ)",fontsize=6)
toF=lambda ax,x,y: fig.transFigure.inverted().transform(ax.transData.transform((x,y)))
rowof={s:[n for n,g in order if code_of[s] in g[1]][0] for s in FOC}

# ---- recovered block 34
def lq_strip(ax):
    ax.axvline(1,color="#BDBDBD",lw=0.6,ls=(0,(2,2)),zorder=0)
    rng=np.random.default_rng(1)
    for n,(lab,ts,kind) in order:
        d=lq[lq.tip_code.isin(ts)&lq.LQ.notna()]
        for _,r in d.iterrows():
            if r.species in FOC: continue
            ax.plot(r.LQ,Y[n]+rng.uniform(-0.18,0.18),"o",ms=2.2,mfc="white" if kind=="og" else GREY,mec=GREY,mew=0.5)
        for _,r in d[d.species.isin(FOC)].iterrows():
            ax.plot(r.LQ,Y[n],"o",ms=4,color=FOC[r.species][0],mec="white",mew=0.4,zorder=4)
    ax.set_xscale("log"); ax.set_xlim(0.3,12); ax.set_xticks([0.5,1,2,5,10]); ax.set_xticklabels(["0.5","1","2","5","10"])
    ax.set_ylim(len(order)-0.4,-0.6); ax.set_yticks([])
    for s in ("left","right","top"): ax.spines[s].set_visible(False)
    ax.set_xlabel("Longevity quotient (LQ)",fontsize=6)
def tree_axes(fig,rect,xl=-62):
    ax=fig.add_axes(rect); build_tree(ax); ax.set_xlim(t.seed_node.age+5,xl); return ax
def ribbon(ov,x0,y0,x1,y1,col,w=0.0055,alpha=0.45):
    xm=(x0+x1)/2
    verts=[(x0,y0+w),(xm,y0+w),(xm,y1+w*0.5),(x1,y1+w*0.5),(x1,y1-w*0.5),(xm,y1-w*0.5),(xm,y0-w),(x0,y0-w),(x0,y0+w)]
    codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
    ov.add_patch(PathPatch(Path(verts,codes),fc=col,ec="none",alpha=alpha,zorder=1))
F=lq[lq.species.isin(FOC)].set_index("species")

# ---- recovered block 35
groups=[(("Euarchontoglires (2)" if g[0].startswith("Euarch") else g[0]),g[1],g[2]) for g in groups]
gmap={g[1]:g for g in groups}; order=[(n,gmap[tipset(n)]) for n,_ in order]
def lq_strip2(ax): lq_strip(ax); ax.minorticks_off()
# ---------- A (polished) ----------

# ---- recovered block 37
from PIL import Image
imgs={}
for s_,vid in [("Myotis myotis","df706c40-0150-400c-bb0e-6ff055c6db61"),("Myotis nigricans","8472d540-5e52-41ec-b43d-9d8f658c8a29"),
              ("Desmodus rotundus","6ac0e1d4-27a3-419a-811a-804104713f62"),("Diphylla ecaudata","b178d461-9954-4918-98bc-10f70cd6e259")]:
    im=Image.open(artifact_path(vid)).convert("RGB"); im.thumbnail((500,500)); imgs[s_]=np.asarray(im)
from PIL import Image
FOC={"Myotis myotis":("#440154",),"Myotis nigricans":("#9A82C2",),"Desmodus rotundus":("#21918C",),"Diphylla ecaudata":("#8ED3CC",)}
def render(fname):
    fig=plt.figure(figsize=(180/25.4,130/25.4))
    axT=tree_axes(fig,[0.02,0.08,0.40,0.87]); axL=fig.add_axes([0.43,0.08,0.11,0.87]); lq_strip2(axL)
    # focal clade wedges coloured by pair
    for p in axT.patches:
        pass
    axS=fig.add_axes([0.655,0.50,0.30,0.40])
    axS.plot(xs,10**(0.537576+0.189234*np.log10(xs)),color="#BDBDBD",lw=0.8)
    axS.text(4e4,10**(0.537576+0.189234*np.log10(4e4))*0.68,"expected for\nnon-flying mammals",fontsize=5.5,color="#8A8A8A",ha="left",va="top")
    axS.plot(b.body_mass_g,b.MLS_years,"o",ms=2.6,mfc=GREY,mec=GREY,mew=0.4,alpha=0.8,label=f"bats (n = {len(b)})")
    axS.plot(o.body_mass_g,o.MLS_years,"o",ms=2.6,mfc="white",mec=GREY,mew=0.6,label=f"non-bat mammals (n = {len(o)})")
    for pair in (("Myotis myotis","Myotis nigricans"),("Desmodus rotundus","Diphylla ecaudata")):
        pp=F.loc[list(pair)]; axS.plot(pp.body_mass_g,pp.MLS_years,color=FOC[pair[0]][0],lw=0.9,zorder=3)
    for s,r in F.iterrows(): axS.plot(r.body_mass_g,r.MLS_years,"o",ms=5.5,color=FOC[s][0],mec="white",mew=0.5,zorder=5)
    axS.set_xscale("log"); axS.set_yscale("log"); axS.minorticks_off()
    axS.set_xticks([10,1e3,1e5,1e7]); axS.set_xticklabels(["10 g","1 kg","100 kg","10 t"])
    axS.set_yticks([3,10,30,100]); axS.set_yticklabels(["3","10","30","100"]); axS.set_xlim(2,2e7); axS.set_ylim(2.5,130)
    axS.yaxis.tick_right(); axS.yaxis.set_label_position("right")
    axS.spines["left"].set_visible(False); axS.spines["top"].set_visible(False); axS.spines["right"].set_visible(True)
    axS.set_xlabel("Body mass",fontsize=6); axS.set_ylabel("Maximum lifespan (years)",fontsize=6)
    axS.legend(loc="lower right",frameon=False,fontsize=5.5,handletextpad=0.2,borderaxespad=0.1)
    fig.text(0.655,0.955,"Long-lived bats outlive their body size",fontsize=7,ha="left",va="top")
    ov=fig.add_axes([0,0,1,1],facecolor="none"); ov.set_xlim(0,1); ov.set_ylim(0,1); ov.axis("off"); ov.set_zorder(-1)
    for s in FOC:
        x0,y0=toF(axL,12,Y[rowof[s]]); x1,y1=toF(axS,F.loc[s,"body_mass_g"],F.loc[s,"MLS_years"])
        ribbon(ov,x0,y0,x1-0.004,y1,FOC[s][0],alpha=0.5)
    # pair cards with portraits
    ov2=fig.add_axes([0,0,1,1],facecolor="none"); ov2.set_xlim(0,1); ov2.set_ylim(0,1); ov2.axis("off")
    cols={"Myotis myotis":(0.625,0.235),"Myotis nigricans":(0.625,0.06),"Desmodus rotundus":(0.815,0.235),"Diphylla ecaudata":(0.815,0.06)}
    ov2.text(0.625,0.415,"Vespertilionidae",fontsize=6,fontweight="bold",color=FOC["Myotis myotis"][0])
    ov2.text(0.815,0.415,"Phyllostomidae: vampire bats",fontsize=6,fontweight="bold",color=FOC["Desmodus rotundus"][0])
    for s,(cx,cy) in cols.items():
        r=F.loc[s]; c=FOC[s][0]; long_=s in LONG
        ia=fig.add_axes([cx,cy,0.075,0.15]); ia.imshow(imgs[s]); ia.axis("off")
        ov2.text(cx+0.08,cy+0.125,s,fontsize=6,style="italic",va="top")
        ov2.text(cx+0.08,cy+0.095,("long-lived" if long_ else "short-lived"),fontsize=5.5,va="top",color="#6E6E6E")
        ov2.text(cx+0.08,cy+0.068,f"LQ {r.LQ:.2f}",fontsize=7,fontweight="bold",va="top",color=c if long_ else "#5A5A5A")
        ov2.text(cx+0.08,cy+0.035,f"{r.MLS_years:g} years · {r.body_mass_g:.3g} g",fontsize=5.5,va="top",color="#4D4D4D")
        ov2.add_patch(mpl.patches.Rectangle((cx+0.08,cy+0.005),0.085,0.006,color=c,lw=0))
    for a,l in ((axT,"a"),(axL,"b"),(axS,"c")): panel_letter(a,l)
    fig.text(0.61,0.43,"d",fontsize=8,fontweight="bold")
    # colour the two focal wedges by pair
    for n,(lab,ts,kind) in order:
        if lab in ("Myotis","Desmodontinae"):
            colr=FOC["Myotis myotis" if lab=="Myotis" else "Desmodus rotundus"][0]
            for p in axT.patches:
                yy=p.get_xy()[0][1]
                if abs(yy-Y[n])<1e-6: p.set_facecolor(colr)
            for tx in axT.texts:
                if tx.get_text().startswith(lab): tx.set_color(colr)
    fig.savefig(fname+".png",dpi=250); fig.savefig(fname+".pdf")
    return fig

# ---- recovered block 42
toF=lambda ax,x,y: ax.figure.transFigure.inverted().transform(ax.transData.transform((x,y)))

# ---- recovered block 44
mm=pd.read_csv(artifact_path("688f2ef8-a0a9-49c3-85bb-815fa8762225"))

# ---- recovered block 45
mm["genus"]=mm.species.str.split("_").str[0]
batgen=set(lq.loc[lq.group!="outgroup","species"].str.split().str[0])
extra={"Pteropus","Eidolon","Epomophorus","Epomops","Macroglossus","Nyctimene","Dobsonia","Megaloglossus","Syconycteris","Chaerephon","Nyctinomops","Promops","Cynomops","Lasionycteris","Lasiurus","Nycticeius","Scotophilus","Plecotus","Barbastella","Nyctalus","Pipistrellus","Vespertilio","Murina","Kerivoula","Tylonycteris","Hypsugo","Scotoecus","Chalinolobus","Nyctophilus","Vespadelus","Mormopterus","Cheiromeles","Otomops","Molossops","Glossophaga","Anoura","Lonchophylla","Sturnira","Stenoderma","Ametrida","Chiroderma","Enchisthenes","Dermanura","Vampyrum","Chrotopterus","Lophostoma","Tonatia","Mimon","Lampronycteris","Pteronotus","Noctilio","Saccopteryx","Peropteryx","Balantiopteryx","Cormura","Emballonura","Coleura","Nycteris","Megaderma","Cardioderma","Lavia","Macroderma","Rhinopoma","Hipposideros","Asellia","Triaenops","Cloeotis","Rhinonicteris","Mystacina","Natalus","Thyroptera","Furipterus","Myzopoda","Craseonycteris","Antrozous","Corynorhinus","Idionycteris","Euderma","Miniopterus","Tadarida","Molossus","Eumops","Mops","Desmodus","Diphylla","Diaemus","Artibeus","Carollia","Phyllostomus","Rousettus","Myotis","Eptesicus","Rhinolophus","Taphozous","Macrotus","Leptonycteris","Uroderma","Platyrrhinus","Vampyressa","Rhinophylla","Choeronycteris","Erophylla","Brachyphylla","Monophyllus","Trachops","Micronycteris","Glyphonycteris","Trinycteris","Lionycteris","Lonchorhina","Macrophyllum","Centurio","Ectophylla","Cynopterus","Eonycteris","Hypsignathus","Aselliscus","Doryrhina","Rhynchonycteris","Cistugo","Mormoops","Ia","Scotomanes","Harpiocephalus","Otonycteris","Arielulus","Glischropus","Laephotis","Neoromicia","Nycticeinops","Rhogeessa","Baeodon","Cynopterus","Penthetor","Ptenochirus","Thoopterus","Haplonycteris","Aproteles","Harpyionycteris","Pteralopex","Melonycteris","Notopteris","Acerodon","Styloctenium","Lissonycteris","Myonycteris","Plerotes","Scotonycteris","Casinycteris","Megaerops","Balionycteris","Aethalops","Alionycteris","Otopteropus","Latidens","Sphaerias","Paranyctimene","Dyacopterus","Hypsugo","Falsistrellus","Mirostrellus","Perimyotis","Parastrellus","Submyotodon","Kerivoula","Phoniscus","Murina","Harpiola","Lasiurus","Dasypterus","Aeorestes","Nyctalus","Scoteanax","Scotorepens","Thainycteris"}
mm["is_bat"]=mm.genus.isin(batgen|extra)
ours=set(lq.species.str.replace(" ","_"))

# ---- recovered block 46
nb=mm[~mm.is_bat]
X=np.c_[np.ones(len(nb)),nb.ln_body]; beta,res,_,_=np.linalg.lstsq(X,nb.ln_lifespan,rcond=None)
resid=nb.ln_lifespan-X@beta; s=np.sqrt((resid**2).sum()/(len(nb)-2))

# ---- recovered block 51
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
def sil(path,col="#4D4D4D"):
    a=np.asarray(Image.open(path).convert("LA")).astype(float)/255
    rgb=np.array(mpl.colors.to_rgb(col)); out=np.zeros(a.shape[:2]+(4,)); out[...,:3]=rgb; out[...,3]=a[...,1]*(1-a[...,0]) if a[...,0].mean()<0.5 else a[...,1]
    return out
WH=sil("data/images/phylopic_whale_93d5c0be.png"); MO=sil("data/images/phylopic_mouse_fe09db98.png")

# ---- recovered block 62
t=dendropy.Tree.get(data=nwk,schema="newick",preserve_underscores=True,rooting="force-rooted")
for n in t.preorder_node_iter(): n.age = 0.0 if n.is_leaf() else age[n.label]
t.ladderize(ascending=False)
clades={tipset(n):n for n in t.postorder_node_iter()}
order=[]; 
def emit(n):
    ts=tipset(n)
    if ts in gmap: order.append((n,gmap[ts])); return
    for k in n.child_nodes(): emit(k)
emit(t.seed_node); Y={}
for i,(n,g) in enumerate(order): Y[n]=i
ypos(t.seed_node)
rowof={s:[n for n,g in order if code_of[s] in g[1]][0] for s in FOC}

# ---- recovered block 63
HP=artifact_path("58c75f9b-ce2a-4712-b80a-61e8af592f2a")
HI=Image.open(HP).convert("RGB"); print(HI.size)
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
icons={};cols_={}

# ---- recovered block 64
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

# ---- recovered block 65
c=icons["Epigenetic alterations"]; h,w=c.shape[:2]; c[:int(h*0.2),int(w*0.85):,3]=0; icons["Epigenetic alterations"]=c
c=icons["Deregulated nutrient sensing"]; h,w=c.shape[:2]; g=grey(c); c[:int(h*0.3),int(w*0.55):,3]=np.where(g[:int(h*0.3),int(w*0.55):],0,c[:int(h*0.3),int(w*0.55):,3]); icons["Deregulated nutrient sensing"]=c
ORDER=["Genomic instability","Telomere attrition","Epigenetic alterations","Loss of proteostasis","Deregulated nutrient sensing","Mitochondrial dysfunction","Cellular senescence","Stem cell exhaustion","Altered intercellular communication"]
LBL={"Genomic instability":"Genomic\ninstability","Telomere attrition":"Telomere\nattrition","Epigenetic alterations":"Epigenetic\nalterations","Loss of proteostasis":"Loss of\nproteostasis",
     "Deregulated nutrient sensing":"Deregulated\nnutrient sensing","Mitochondrial dysfunction":"Mitochondrial\ndysfunction","Cellular senescence":"Cellular\nsenescence","Stem cell exhaustion":"Stem cell\nexhaustion","Altered intercellular communication":"Altered intercellular\ncommunication"}
B0,H0=0.20,0.80
import re as _re

# ---- recovered block 72
c=icons["Epigenetic alterations"]; h,w=c.shape[:2]
ys,xs_=np.where(c[...,3]>30); 
# drop small connected blobs in the top-right quadrant
from scipy import ndimage
lab_,nl=ndimage.label(c[...,3]>30); sizes=ndimage.sum(np.ones_like(lab_),lab_,range(1,nl+1))
for i_,s_ in enumerate(sizes,1):
    yy,xx=np.where(lab_==i_)
    if s_<0.002*h*w and xx.mean()>0.8*w and yy.mean()<0.35*h: c[lab_==i_,3]=0
icons["Epigenetic alterations"]=c

# ---- recovered block 73
c=icons["Epigenetic alterations"]; h,w=c.shape[:2]
lab_,nl=ndimage.label(c[...,3]>30); sizes=ndimage.sum(np.ones_like(lab_),lab_,range(1,nl+1))
for i_,s_ in enumerate(sizes,1):
    yy,xx=np.where(lab_==i_)
    if xx.min()>0.85*w and yy.max()<0.4*h: c[lab_==i_,3]=0; print("removed",s_)
ys,xs__=np.where(c[...,3]>30); icons["Epigenetic alterations"]=c[ys.min():ys.max()+1,xs__.min():xs__.max()+1]

# ---- node-age uncertainty (95% HPD bands) and bootstrap support in panel a
sup=pd.read_csv("data/T2_node_support_age_uncertainty.csv").set_index("node")
inside=set()
for n,(lab,ts,kind) in order:
    for d in n.preorder_iter():
        if d is not n: inside.add(d)
vis=[n for n in t.preorder_node_iter() if not n.is_leaf() and n not in inside]
HPD_COL="#B5654A"
def _band(ax,x0,x1,y,h=0.22):
    ax.add_patch(Polygon([[x0,y-h],[x1,y-h],[x1,y+h],[x0,y+h]],closed=True,fc=HPD_COL,ec="none",alpha=0.38,zorder=0.5))
def add_node_uncertainty(ax):
    xmax=t.seed_node.age+4.5
    for n in vis:
        lo,hi=sup.loc[n.label,["hpd_lo_mean","hpd_hi_mean"]]
        _band(ax,lo,min(hi,xmax),Y[n])
        if hi>xmax:   # band clipped at the plot edge: arrowhead marks the continuation
            ax.add_patch(Polygon([[xmax,Y[n]-0.22],[xmax+1.6,Y[n]],[xmax,Y[n]+0.22]],closed=True,fc=HPD_COL,ec="none",alpha=0.38,zorder=0.5))
        bs=sup.loc[n.label,"bootstrap_pct"]
        if bs<100:
            ax.text(n.age+0.8,Y[n]-0.32,f"{bs:.0f}",fontsize=5,ha="right",va="bottom",color="#4D4D4D",zorder=6)
    kx=t.seed_node.age+2; ky=11.0
    _band(ax,kx-5,kx,ky)
    ax.text(kx-6.5,ky,"node age,\n95% HPD",fontsize=5,va="center",ha="left",color="#4D4D4D",linespacing=1.0)
    ax.text(kx-2.5,ky+1.35,"90",fontsize=5,va="center",ha="center",color="#4D4D4D")
    ax.text(kx-6.5,ky+1.35,"bootstrap %\n(blank = 100)",fontsize=5,va="center",ha="left",color="#4D4D4D",linespacing=1.0)
_build_tree_orig=build_tree
def build_tree(ax):
    _build_tree_orig(ax); add_node_uncertainty(ax)
    ax.set_xlim(t.seed_node.age+6.5,-38)

# ---- final render
from scipy import ndimage
m_=lq[lq.MLS_years.notna()&lq.body_mass_g.notna()]; b=m_[m_.group!="outgroup"]; assert len(b)==40
def render2(fname):
    fig=plt.figure(figsize=(180/25.4,190/25.4)); R=lambda r:[r[0],B0+r[1]*H0,r[2],r[3]*H0]; _ft=fig.text; ftext=lambda x,y,*a,**k: _ft(x,B0+y*H0,*a,**k)
    axT=tree_axes(fig,R([0.035,0.08,0.39,0.87])); axL=fig.add_axes(R([0.43,0.08,0.11,0.87])); lq_strip2(axL)
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
    axS.plot(b.body_mass_g,b.MLS_years,"o",ms=2.8,color="#4D4D4D",mec="white",mew=0.3,zorder=3,label=f"bats in this study (n = {len(b)})")
    for pair in (("Myotis myotis","Myotis nigricans"),("Desmodus rotundus","Diphylla ecaudata")):
        pp=F.loc[list(pair)]; axS.plot(pp.body_mass_g,pp.MLS_years,color=FOC[pair[0]][0],lw=0.9,zorder=4)
    for s,r in F.iterrows(): axS.plot(r.body_mass_g,r.MLS_years,"o",ms=5.5,color=FOC[s][0],mec="white",mew=0.5,zorder=5)
    hum=mm[mm.species=="Homo_sapiens"].iloc[0]; nmr=mm[mm.species=="Heterocephalus_glaber"].iloc[0]
    axS.annotate("Human",(np.exp(hum.ln_body),np.exp(hum.ln_lifespan)),xytext=(-4,0),textcoords="offset points",fontsize=5.5,ha="right",va="center",color="#4D4D4D")
    axS.plot(np.exp(nmr.ln_body),np.exp(nmr.ln_lifespan),"^",ms=4.2,mfc="white",mec="#4D4D4D",mew=0.8,zorder=7)
    axS.plot(np.exp(hum.ln_body),np.exp(hum.ln_lifespan),"^",ms=4.2,mfc="white",mec="#4D4D4D",mew=0.8,zorder=7)
    axS.annotate("Naked mole-rat",(np.exp(nmr.ln_body),np.exp(nmr.ln_lifespan)),xytext=(400,75),textcoords="data",fontsize=5.5,ha="left",va="center",color="#4D4D4D",
                 arrowprops=dict(arrowstyle="-|>",lw=0.5,color="#6E6E6E",mutation_scale=5,shrinkA=1,shrinkB=2.5))
    wh=mm.loc[mm.ln_body.idxmax()]; mo=mm.loc[mm.ln_body.idxmin()]
    axS.add_artist(AnnotationBbox(OffsetImage(WH,zoom=0.022),(np.exp(wh.ln_body),np.exp(wh.ln_lifespan)),xybox=(-14,10),boxcoords="offset points",frameon=False,box_alignment=(1.0,0.5),zorder=6))
    axS.add_artist(AnnotationBbox(OffsetImage(MO,zoom=0.011),(np.exp(mo.ln_body),np.exp(mo.ln_lifespan)),xybox=(0,-13),boxcoords="offset points",frameon=False,box_alignment=(0.5,0.5),zorder=6))
    axS.set_xscale("log"); axS.set_yscale("log"); axS.minorticks_off()
    axS.set_xticks([10,1e3,1e5,1e7]); axS.set_xticklabels(["10 g","1 kg","100 kg","10 t"])
    axS.set_yticks([1,3,10,30,100]); axS.set_yticklabels(["1","3","10","30","100"]); axS.set_xlim(2,1e8); axS.set_ylim(1.2,150)
    axS.yaxis.tick_right(); axS.yaxis.set_label_position("right")
    axS.spines["left"].set_visible(False); axS.spines["top"].set_visible(False); axS.spines["right"].set_visible(True)
    axS.set_xlabel("Body mass",fontsize=6); axS.set_ylabel("Maximum lifespan (years)",fontsize=6)
    axS.legend(loc="lower right",frameon=False,fontsize=5.5,handletextpad=0.2,borderaxespad=0.1,markerscale=1.4)
    ftext(0.70,0.955,"Long-lived bats outlive their body size",fontsize=7,ha="left",va="top")
    ftext(0.675,0.955,"c",fontsize=8,fontweight="bold",va="top")
    for s in FOC:
        x0,y0=toF(axL,12,Y[rowof[s]]); x1,y1=toF(axS,F.loc[s,"body_mass_g"],F.loc[s,"MLS_years"])
        ribbon(ov,x0,y0,x1-0.004,y1,FOC[s][0],alpha=0.5)
    ov2=fig.add_axes(R([0,0,1,1]),facecolor="none"); ov2.set_xlim(0,1); ov2.set_ylim(0,1); ov2.axis("off"); ov2.set_zorder(2)
    cols={"Myotis myotis":(0.598,0.215),"Myotis nigricans":(0.598,0.035),"Desmodus rotundus":(0.80,0.215),"Diphylla ecaudata":(0.80,0.035)}
    ov2.text(0.598,0.395,"Vespertilionidae",fontsize=6,fontweight="bold",color=FOC["Myotis myotis"][0])
    ov2.text(0.80,0.395,"Phyllostomidae: vampire bats",fontsize=6,fontweight="bold",color=FOC["Desmodus rotundus"][0])
    ftext(0.578,0.395,"d",fontsize=8,fontweight="bold")
    for s,(cx,cy) in cols.items():
        r=F.loc[s]; c=FOC[s][0]; long_=s in LONG
        ia=fig.add_axes(R([cx-0.01,cy,0.085,0.15])); ia.imshow(imgs[s]); ia.axis("off"); ia.set_zorder(2)
        tx=cx+0.088
        ov2.text(tx,cy+0.135,s,fontsize=6,style="italic",va="top")
        ov2.text(tx,cy+0.108,("long-lived" if long_ else "short-lived"),fontsize=5.5,va="top",color="#6E6E6E")
        ov2.text(tx,cy+0.082,f"LQ {r.LQ:.2f}",fontsize=7,fontweight="bold",va="top",color=c if long_ else "#5A5A5A")
        ov2.text(tx,cy+0.050,f"{r.MLS_years:g} years · {r.body_mass_g:.3g} g",fontsize=5.5,va="top",color="#4D4D4D")
        ov2.add_patch(mpl.patches.Rectangle((tx,cy+0.022),0.085,0.006,color=c,lw=0))
    for a,l in ((axT,"a"),(axL,"b")): panel_letter(a,l)
    for n,(lab,ts,kind) in order:
        if lab in ("Myotis","Desmodontinae"):
            colr=FOC["Myotis myotis" if lab=="Myotis" else "Desmodus rotundus"][0]
            for p in axT.patches:
                if abs(p.get_xy()[0][1]-Y[n])<1e-6: p.set_facecolor(colr)
            for tx_ in axT.texts:
                if tx_.get_text().startswith(lab): tx_.set_color(colr)
    bn=fig.add_axes([0.0,0.0,1.0,0.185]); bn.axis("off"); bn.set_xlim(0,1); bn.set_ylim(0,1)
    bn.plot([0.02,0.98],[0.97,0.97],color="#D0D0D0",lw=0.6)
    fig.text(0.02,0.176,"e",fontsize=8,fontweight="bold",va="top"); fig.text(0.045,0.176,"Hallmarks of aging",fontsize=7,va="top")
    W=0.96/9
    for i,k in enumerate(ORDER):
        cx=0.02+W*(i+0.5); im_=icons[k]; h,w=im_.shape[:2]
        boxw=W*0.82*180; boxh=0.085*190  # mm
        sc=min(boxw/w,boxh/h); ww=w*sc/180; hh=h*sc/190
        ia=fig.add_axes([cx-ww/2,0.072+(0.085-hh)/2,ww,hh]); ia.imshow(im_,interpolation="antialiased"); ia.axis("off")
        fig.text(cx,0.062,LBL[k],fontsize=5.5,ha="center",va="top",color=labcol[k],fontweight="bold",linespacing=1.05)
    fig.savefig(fname+".png",dpi=250,bbox_inches=None); fig.savefig(fname+".pdf",bbox_inches=None)
    return fig

plt.close("all"); fig=render2("figures/phylo_LQ_hallmarks_180x190")
