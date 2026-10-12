#!/usr/bin/env python3
"""Per-node age uncertainty and topological support for the T2 timetree (Morales et al. 2026).

Nothing is re-estimated: all values come from the authors' own outputs.
  * Node ages and 95% HPD: the 77 MCMCTree runs in
    Divergence_Time_estimates/Caster_tree/coding_sequences/bin_*/divtime/mcmc/FigTree.tre
    (Caster_tree = T2 topology; ages in units of 100 Myr). T2_timetree.tre is the mean of
    these 77 runs (checked: max difference 0.03 Myr).
  * Bootstrap (bs1, out of 1000) and quartet support (q1): the neutral-window CASTER tree
    tree_T2/trees.neutralWindows/multiZhg38.concatNeutralWin_500_COV80.scored_U2.tre
    (same topology as T2; unrooted, so the two root nodes have no value).
  * Mean ages from the 104 neutral-site (random-window) runs are given for comparison.

Output columns (ages in Myr): node, n_tips, age_mean, hpd_lo_mean, hpd_hi_mean (mean of the
per-run 95% HPD bounds; plotted as bands), hpd_lo_min, hpd_hi_max, between_bin_q025/q975
(spread of run means; measures agreement between runs, not dating uncertainty),
neutral_random_mean, bootstrap_pct, q1, tips.

Usage:
  python scripts/01_extract_node_uncertainty.py DIVTIME_DIR SUPPORT_TREE LABELED_TREE OUT_CSV
On smith (Divergence_Time_estimates.tar.gz extracted there):
  python scripts/01_extract_node_uncertainty.py \\
     /scratch4/workspace/tlama_smith_edu-project_longevity/dating/Divergence_Time_estimates \\
     /scratch4/workspace/tlama_smith_edu-project_longevity/tree_T2/trees.neutralWindows/multiZhg38.concatNeutralWin_500_COV80.scored_U2.tre \\
     data/T2_timetree.labeled.nwk data/T2_node_support_age_uncertainty.csv
Requires numpy and pandas.
"""
import glob, re, sys
import numpy as np, pandas as pd


def parse_newick(s):
    """Return (internal_nodes, root); each internal node is (tipset, height, annotation text)."""
    s = re.sub(r"\s+(?=[(\[])|(?<=[,(])\s+", "", s.strip().rstrip(";"))
    i = 0
    nodes = []

    def node():
        nonlocal i
        kids = []
        if s[i] == "(":
            i += 1
            while True:
                kids.append(node())
                if s[i] == ",":
                    i += 1
                    continue
                if s[i] == ")":
                    i += 1
                    break
        m = re.match(r"\s*('[^']*'|[^:,()\[\s]*)\s*(\[[^\]]*\])?\s*(?::\s*([\d.eE+-]+))?", s[i:])
        i += m.end()
        name, ann, bl = m.group(1), m.group(2), float(m.group(3) or 0)
        if s[i:i + 1] == "[":
            k = s.index("]", i)
            ann = (ann or "") + s[i:k + 1]
            i = k + 1
        if kids:
            tips = frozenset().union(*[k[0] for k in kids])
            h = max(k[1] + k[2] for k in kids)
            nodes.append((tips, h, (name or "") + (ann or "")))
            return (tips, h, bl)
        return (frozenset([name.strip("'")]), 0.0, bl)

    root = node()
    return nodes, root


def tree_text(path):
    t = open(path).read()
    m = re.search(r"UTREE\s+\d+\s*=\s*(.*;)", t)
    return m.group(1) if m else t


def main(divdir, support_tree, labeled_tree, out_csv):
    lab_nodes, lab_root = parse_newick(open(labeled_tree).read())
    ALL = lab_root[0]
    canon = lambda ts: ts if len(ts) <= len(ALL) / 2 else frozenset(ALL - ts)
    label = {canon(ts): a.strip("'") for ts, h, a in lab_nodes}
    real = {a.strip("'"): ts for ts, h, a in lab_nodes}

    recs = []
    for grp in ["coding_sequences", "neutral_sites/random"]:
        files = sorted(glob.glob(f"{divdir}/Caster_tree/{grp}/bin_*/divtime/mcmc/FigTree.tre"))
        assert files, f"no FigTree.tre under {grp}"
        for f in files:
            ns, rt = parse_newick(tree_text(f))
            assert rt[0] == ALL, f"tip set differs: {f}"
            for ts, h, a in ns:
                m = re.search(r"95%=\{([\d.]+),\s*([\d.]+)\}", a)
                lo, hi = (map(float, m.groups()) if m else (np.nan, np.nan))
                recs.append((grp, canon(ts), h, lo, hi))
    R = pd.DataFrame(recs, columns=["grp", "clade", "age", "lo", "hi"])
    assert R.clade.isin(label.keys()).all(), "a dated clade is not in the labelled T2 tree"

    C = R[R.grp == "coding_sequences"].groupby("clade")
    A = pd.DataFrame({"age_mean": C.age.mean(), "hpd_lo_mean": C.lo.mean(), "hpd_hi_mean": C.hi.mean(),
                      "hpd_lo_min": C.lo.min(), "hpd_hi_max": C.hi.max(),
                      "between_bin_q025": C.age.quantile(.025), "between_bin_q975": C.age.quantile(.975)}) * 100
    A["neutral_random_mean"] = R[R.grp == "neutral_sites/random"].groupby("clade").age.mean() * 100

    sn, sr = parse_newick(open(support_tree).read())
    assert sr[0] == ALL
    sup = {}
    for ts, h, a in sn:
        d = dict(re.findall(r"(\w+)=([\d.eE+-]+)", a))
        if d:
            sup[canon(ts)] = d
    A["bootstrap_pct"] = [float(sup[k]["bs1"]) / 10 if k in sup else np.nan for k in A.index]
    A["q1"] = [float(sup[k]["q1"]) if k in sup else np.nan for k in A.index]
    A["node"] = [label[k] for k in A.index]
    A["n_tips"] = [len(real[n]) for n in A.node]
    A["tips"] = [";".join(sorted(real[n])) for n in A.node]
    cols = ["node", "n_tips", "age_mean", "hpd_lo_mean", "hpd_hi_mean", "hpd_lo_min", "hpd_hi_max",
            "between_bin_q025", "between_bin_q975", "neutral_random_mean", "bootstrap_pct", "q1", "tips"]
    A[cols].sort_values("age_mean", ascending=False).round(3).to_csv(out_csv, index=False)
    print(f"{len(A)} nodes written to {out_csv}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
