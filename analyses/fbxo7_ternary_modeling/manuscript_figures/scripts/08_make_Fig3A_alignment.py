"""Fig 3A: FBXO7 protein alignment at the variable sites (human UniProt Q9Y3I1 numbering) on the species tree.
Row order, names and clade colours come from data/alignment/species_names.tsv; residues from the 30-species TOGA/MACSE/HmmCleaner codon alignment
(translated here), with explicit per-cell replacements for masked codons in data/alignment/site_overrides.tsv. Changing SITES redraws the figure.
Colouring (read off the original hand-made figure): fixed residue palette; the commonest residue of a column is left unshaded when it is at least
UNSHADE_FREQ of the sequences, except columns in ALWAYS_SHADE. Only the rows in SHADE_ONLY are coloured (set to None to colour all rows). Also writes the residue table and a QC report (masked cells, overrides, and
substitutions seen in the four focal species that are not in the Rosetta composites listed in Table S16)."""
import pathlib, itertools, collections, io
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from Bio import Phylo
HERE = pathlib.Path(__file__).resolve().parent if '__file__' in globals() else pathlib.Path('scripts').resolve()
DATA, FIGS = HERE.parent / 'data', HERE.parent / 'figures'; AD = DATA / 'alignment'; FIGS.mkdir(exist_ok=True)
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'pdf.fonttype': 42})

# ------------------------------ configuration ------------------------------
import sys
VERSION = sys.argv[1] if len(sys.argv) > 1 else 'both'        # 'myotis_psmf1' | 'four_species' | 'both'
H = 110                                                        # canvas height in mm (shared by both versions)
CONF = {
 # the ten positions variable in the two Myotis, plus three PSMF1-region positions (191, 290, 292) where
 # BOTH Myotis retain the human residue -- those substitutions occur only in the two vampire bats
 'myotis_psmf1': dict(SITES=[17, 19, 47, 64, 109, 119, 127, 146, 151, 191, 290, 292, 409],
                      BOXES=[('HLmyoMyo6', 'HLmyoNig1')],
                      SHADE_ONLY=['HLmyoMyo6', 'HLmyoNig1'], W=100),
 # every position carrying a modelled substitution in any of the four focal species
 'four_species': dict(SITES=[17, 19, 47, 64, 109, 110, 119, 127, 146, 151, 191, 290, 292, 409],
                      BOXES=[('HLmyoMyo6', 'HLmyoNig1'), ('HLdesRot8A', 'HLdipEca1')],
                      SHADE_ONLY=['HLmyoMyo6', 'HLmyoNig1', 'HLdesRot8A', 'HLdipEca1'], W=104)}

def render(VERSION):
    SITES, BOXES, SHADE_ONLY, W = (CONF[VERSION][k] for k in ('SITES', 'BOXES', 'SHADE_ONLY', 'W'))

    REF = 'hg38'
    OUTGROUP = ['hg38', 'mm10', 'felCat9', 'canFam4', 'equCab3', 'susScr11', 'bosTau9']
    OUTGROUP_TOPOLOGY = '((hg38,mm10),((felCat9,canFam4),(equCab3,(susScr11,bosTau9))))'   # accepted mammal phylogeny; the file tree places horse with human/mouse (71% support)
    FOREGROUND = ['HLsacLep1', 'HLmyoMyo6']          # branches below the MRCA of these tips (Yangochiroptera) are drawn blue
    FOCAL = {'HLmyoMyo6': 'Myotis myotis', 'HLmyoNig1': 'Myotis nigricans', 'HLdesRot8A': 'Desmodus rotundus', 'HLdipEca1': 'Diphylla ecaudata'}
    UNSHADE_FREQ, ALWAYS_SHADE = 0.5, {127}
    UNSHADE_REF = {191, 290, 292}   # PSMF1-interface positions added for context: leave the human residue unshaded, so shading there means 'differs from human'
    PALETTE = {  # observed in the original figure: A C E F G H I L M P Q R S T V; assumed by analogy (same colour family): D K N W Y
        'A': '#306eb8', 'V': '#377df7', 'I': '#4fa0f4', 'L': '#79a7d7', 'F': '#3880e6', 'M': '#22549a', 'W': '#3770c8', 'C': '#e68080', 'G': '#e69a4d', 'P': '#cdcc42',
        'S': '#43982a', 'T': '#5dca45', 'Q': '#6eed5c', 'N': '#6eed5c', 'H': '#52b4b4', 'Y': '#52b4b4', 'E': '#984098', 'D': '#984098', 'R': '#ea4129', 'K': '#ea4129'}

    # ------------------------------ residues ------------------------------
    names = pd.read_csv(AD / 'species_names.tsv', sep='\t'); ids = names.id.tolist(); disp = dict(zip(names.id, names.display_name))
    CODON = {a + b + c: aa for (a, b, c), aa in zip(itertools.product('TCAG', repeat=3), 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
    IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'N': 'ACGT'}
    def aa_of(c):
        if c == '---': return '-'
        if any(x not in IUPAC for x in c): return 'X'
        s = {CODON[''.join(p)] for p in itertools.product(*(IUPAC[x] for x in c))}; return s.pop() if len(s) == 1 else 'X'
    seqs, k = {}, None
    for l in open(AD / 'FBXO7_30sp_codon_aln.fasta'):
        l = l.strip()
        if l.startswith('>'): k = l[1:]; seqs[k] = []
        elif k: seqs[k].append(l)
    prot = {k: ''.join(aa_of(''.join(v).upper()[i:i + 3]) for i in range(0, len(''.join(v)) - 2, 3)) for k, v in seqs.items()}
    assert set(ids) == set(prot), set(ids) ^ set(prot)
    col, n = {}, 0
    for i, a in enumerate(prot[REF]):
        if a != '-': n += 1; col[n] = i
    R = pd.DataFrame({p: [prot[i][col[p]] for i in ids] for p in SITES}, index=ids)
    masked = [(i, p) for i in ids for p in SITES if R.loc[i, p] == 'X']
    ov = pd.read_csv(AD / 'site_overrides.tsv', sep='\t'); applied = []
    for r in ov.itertuples():
        if r.human_pos in SITES and R.loc[r.species_id, r.human_pos] in ('X', '-'): R.loc[r.species_id, r.human_pos] = r.residue; applied.append((r.species_id, r.human_pos, r.residue))
    R.rename(index=disp).to_csv(FIGS.parent / 'data' / f'FBXO7_site_residues_{VERSION}.csv', index_label='species')

    # ------------------------------ tree ------------------------------
    tr = Phylo.read(io.StringIO(open(AD / 'tree30.nwk').read()), 'newick'); tr.root_with_outgroup(OUTGROUP)
    def to_d(c): return {'name': c.name, 'kids': [to_d(x) for x in c.clades]}
    def parse(s):
        pos = 0
        def node():
            nonlocal pos
            if s[pos] == '(':
                pos += 1; kids = [node()]
                while s[pos] == ',': pos += 1; kids.append(node())
                pos += 1; return {'name': None, 'kids': kids}
            j = pos
            while s[j] not in ',()': j += 1
            nm = s[pos:j]; pos = j; return {'name': nm, 'kids': []}
        return node()
    tips = lambda d: [d['name']] if not d['kids'] else [t for k in d['kids'] for t in tips(k)]
    T0 = to_d(tr.root); og_kids = [k for k in T0['kids'] if set(tips(k)) <= set(OUTGROUP)]; rest = [k for k in T0['kids'] if all(k is not o for o in og_kids)]
    assert og_kids and rest, 'outgroup clade not found at the root after rerooting'
    T = {'name': None, 'kids': [parse(OUTGROUP_TOPOLOGY)] + rest}; assert set(tips(T)) == set(ids)   # root = (outgroup, bats)
    pref = {i: r for r, i in enumerate(ids)}
    def order(d):
        d['kids'] = sorted((order(k) for k in d['kids']), key=lambda k: min(pref[t] for t in tips(k))); return d
    T = order(T); assert tips(T) == ids, 'row order in species_names.tsv is not compatible with the tree topology'
    def height(d): d['h'] = 0 if not d['kids'] else 1 + max(height(k) for k in d['kids']); return d['h']
    height(T)
    def find_mrca(d, tset):
        for k in d['kids']:
            if tset <= set(tips(k)): return find_mrca(k, tset)
        return d
    fg = find_mrca(T, set(FOREGROUND))
    def mark(d, on):
        d['fg'] = on
        for k in d['kids']: mark(k, on)
    mark(T, False); mark(fg, True)

    # ------------------------------ shading ------------------------------
    white = {}
    for p in SITES:
        vals = [v for v in R[p] if v not in ('X', '-')]; top, cnt = collections.Counter(vals).most_common(1)[0]
        white[p] = top if (cnt / len(vals) >= UNSHADE_FREQ and p not in ALWAYS_SHADE) else None
        if p in UNSHADE_REF: white[p] = R.loc[REF, p]

    # ------------------------------ draw (canvas in mm) ------------------------------
    FS_RES, FS_NAME, FS_HEAD = 5.8, 5.8, 6.5
    TOP, BOT, PAD, TREE_W, GAP = 9.0, 1.0, 1.5, 12.0, 1.5
    ROW = (H - TOP - BOT) / len(ids)                                     # row pitch (mm) follows from the canvas height
    _f = plt.figure(figsize=(W / 25.4, H / 25.4)); _r = _f.canvas.get_renderer()
    NAME_W = max(_f.text(0, 0, disp[i], fontsize=FS_NAME, style='italic').get_window_extent(_r).width for i in ids) / _f.dpi * 25.4; plt.close(_f)
    XT0, XT1 = PAD, PAD + TREE_W                                         # tree
    XN = XT1 + GAP; XC = XN + NAME_W + GAP                               # names, then the residue columns
    CW = (W - PAD - XC) / len(SITES)                                     # column width (mm)
    assert CW >= 3.6, f'columns too narrow ({CW:.2f} mm); shorten TREE_W or the names'
    fig = plt.figure(figsize=(W / 25.4, H / 25.4)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis('off')
    yrow = {i: TOP + ROW * (r + .5) for r, i in enumerate(ids)}
    BLUE, BLACK = '#1f3fe0', 'black'
    def draw_y(d):
        if d['kids']: d['y'] = np.mean([draw_y(k) for k in d['kids']])
        else: d['y'] = yrow[d['name']]
        return d['y']
    def pos_x(d): return XT1 - (XT1 - XT0) * d['h'] / T['h'] if d['kids'] else XT1
    def branches(d):
        d['x'] = pos_x(d)
        for k in d['kids']:
            branches(k); c = BLUE if k['fg'] else BLACK
            ax.plot([d['x'], k['x']], [k['y'], k['y']], color=c, lw=0.7, solid_capstyle='butt')
        if d['kids']:
            for k in d['kids']: ax.plot([d['x']] * 2, [k['y'], d['y']], color=BLUE if k['fg'] else BLACK, lw=0.7, solid_capstyle='butt')
    draw_y(T); branches(T)
    lum = lambda h: (lambda c: .2126 * c[0] + .7152 * c[1] + .0722 * c[2])(mpl.colors.to_rgb(h))
    for j, p in enumerate(SITES):
        x0 = XC + CW * j
        ax.text(x0 + CW / 2, TOP - 5.2, str(p), ha='center', va='center', fontsize=FS_HEAD, fontweight='bold')
        ax.add_patch(Polygon([[x0 + CW / 2 - 1.1, TOP - 3.1], [x0 + CW / 2 + 1.1, TOP - 3.1], [x0 + CW / 2, TOP - 1.5]], closed=True, color='black', lw=0))
        for i in ids:
            a = R.loc[i, p]; y = yrow[i]; shaded = a in PALETTE and a != white[p] and (SHADE_ONLY is None or i in SHADE_ONLY)
            if shaded: ax.add_patch(Rectangle((x0 + .2, y - ROW / 2), CW - .4, ROW, facecolor=PALETTE[a], edgecolor='none'))
            fc = 'white' if shaded and lum(PALETTE[a]) < .45 else ('#9a9a9a' if a in ('X', '-') else 'black')
            ax.text(x0 + CW / 2, y + .05, a, ha='center', va='center', fontsize=FS_RES, fontweight='bold', color=fc)
    for i in ids: ax.text(XN, yrow[i], disp[i], ha='left', va='center', fontsize=FS_NAME, style='italic')
    for a_, b_ in BOXES: ax.add_patch(Rectangle((XN - 0.8, yrow[a_] - ROW / 2), XC + CW * len(SITES) - XN + 1.2, yrow[b_] - yrow[a_] + ROW, fill=False, edgecolor='#d7191c', lw=0.8))
    fig.savefig(FIGS / f'Fig3A_alignment_sites_{VERSION}.pdf'); fig.savefig(FIGS / f'Fig3A_alignment_sites_{VERSION}.png', dpi=300)

    # ------------------------------ QC report ------------------------------
    L = ['Fig 3A alignment QC', f'sites: {SITES}', f'alignment: {len(ids)} sequences, {len(prot[REF])} codons; human numbering via {REF}', '',
         f'masked/ambiguous codons at the sites: {len(masked)}']
    L += [f'  {disp[i]} ({i}) position {p}: ' + (f'replaced by {dict(((a, b), c) for a, b, c in applied)[(i, p)]} from site_overrides.tsv' if (i, p) in {(a, b) for a, b, _ in applied} else 'still X (drawn as X)') for i, p in masked]
    t16 = pd.read_csv(DATA / 'SuppTable16_AlphaMissense_updated.csv'); L += ['', 'substitutions vs human at SITES in the four focal species, compared with the Rosetta composites (Table S16 Species_site_sets):']
    for i, nm in FOCAL.items():
        obs = {f"{R.loc[REF, p]}{p}{R.loc[i, p]}" for p in SITES if R.loc[i, p] not in ('X', '-') and R.loc[i, p] != R.loc[REF, p]}
        mod = {v for v in t16[t16.Species_site_sets.str.contains(nm, regex=False)].Variant if int(v[1:-1]) in SITES}   # only positions shown in the figure
        L.append(f'  {nm}: in alignment, not modelled: {sorted(obs - mod) or "none"} | modelled, not in alignment: {sorted(mod - obs) or "none"}')
    L += ['', 'unshaded residue per column (commonest, >= %.0f%% of sequences): ' % (100 * UNSHADE_FREQ) + ', '.join(f'{p}:{white[p] or "-"}' for p in SITES)]
    open(FIGS.parent / 'data' / f'Fig3A_alignment_QC_{VERSION}.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))


for _v in (['myotis_psmf1', 'four_species'] if VERSION == 'both' else [VERSION]):
    render(_v)
