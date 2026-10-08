#!/usr/bin/env python3
"""Percent missing data in the 30-species codon alignments (TOGA orthologs -> MACSE v2 exon alignment -> concatenated -> HmmCleaner).

Definition (identical to the calculation behind the manuscript numbers)
  * A position is "missing" if it is anything other than A, C, G or T after upper-casing. In these alignments that is N: HmmCleaner-masked
    positions and alignment gaps both appear as N (the only characters present are A, C, G, N, T).
  * Per-sequence missing fraction  = missing positions / sequence length (empty sequences are skipped).
  * Per-alignment missing fraction = total missing positions / total positions, pooled over all sequences in the alignment
    (a length-weighted figure; it is NOT the mean of the per-sequence fractions).
  * Summary statistics are taken over alignments (median, mean, interquartile range) and, separately, over sequences.
  * The interquartile range uses the index convention  sorted[n//4], sorted[3n//4]  (as in the original calculation). The interpolated
    ("inclusive") quartiles are also printed for comparison; they differ only in the last digit or two.
  * Only files ending in .fasta are read. The 7082_fastas directory also holds FASconCAT-G_v1.05.pl and fasconcat.sh; reading those as
    FASTA inflates the statistics, so they are excluded.

Standard library only (Python >= 3.6). Usage:
    python3 alignment_missingness.py --root /path/to/data --out qc_out
Defaults read <root>/final_alignment/fastas and <root>/7082_fastas; use --dirs to give other directories (label=path pairs allowed).
Outputs: <out>/alignment_missingness_per_alignment.csv, <out>/alignment_missingness_summary.txt
"""
import argparse, collections, csv, glob, os, statistics, sys

BASES = 'ACGT'


def read_fasta(path):
    """Return {header: uppercase sequence}. A repeated header replaces the earlier record (as in the original calculation)."""
    seqs, name = collections.OrderedDict(), None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: ''.join(v).upper() for k, v in seqs.items()}


def n_missing(seq):
    return len(seq) - sum(seq.count(b) for b in BASES)


def analyse(path, label):
    seqs = read_fasta(path)
    tot = miss = 0
    seq_fracs = []
    for name, s in seqs.items():
        if not s:
            continue
        m = n_missing(s)
        tot += len(s)
        miss += m
        seq_fracs.append((m / len(s), name))
    row = dict(set=label, alignment=os.path.basename(path), n_sequences=len(seqs), n_nonempty=len(seq_fracs),
               total_positions=tot, missing_positions=miss, frac_missing=(miss / tot) if tot else float('nan'),
               max_seq_frac=max((f for f, _ in seq_fracs), default=float('nan')),
               n_seq_gt50=sum(f > 0.5 for f, _ in seq_fracs), n_seq_gt80=sum(f > 0.8 for f, _ in seq_fracs))
    return row, seq_fracs


def pct(x):
    return '%.2f%%' % (100 * x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', default='/project/tlama_umass_edu/projects/project_bat1k_longevity/data')
    ap.add_argument('--dirs', nargs='*', default=None, help='directories of .fasta alignments, optionally label=path')
    ap.add_argument('--out', default='alignment_missingness_out')
    a = ap.parse_args()
    dirs = a.dirs or ['final_alignment/fastas=' + os.path.join(a.root, 'final_alignment', 'fastas'), '7082_fastas=' + os.path.join(a.root, '7082_fastas')]
    sets = []
    for d in dirs:
        label, path = d.split('=', 1) if '=' in d else (os.path.basename(os.path.normpath(d)), d)
        sets.append((label, path))
    os.makedirs(a.out, exist_ok=True)

    rows, seq_all, seq_hi = [], [], []
    for label, path in sets:
        files = sorted(glob.glob(os.path.join(path, '*.fasta')))
        print('%s: %d .fasta files' % (label, len(files)), file=sys.stderr)
        for f in files:
            row, sf = analyse(f, label)
            if not row['total_positions']:
                continue
            rows.append(row)
            seq_all.extend(x for x, _ in sf)
            seq_hi.extend((x, row['alignment'], n) for x, n in sf if x > 0.8)

    cols = ['set', 'alignment', 'n_sequences', 'n_nonempty', 'total_positions', 'missing_positions', 'frac_missing', 'max_seq_frac', 'n_seq_gt50', 'n_seq_gt80']
    with open(os.path.join(a.out, 'alignment_missingness_per_alignment.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (('%.6f' % r[k]) if isinstance(r[k], float) else r[k]) for k in cols})

    af = sorted(r['frac_missing'] for r in rows)
    n = len(af)
    L = ['Missing data (non-ACGT) in codon alignments', 'directories: ' + '; '.join('%s = %s' % s for s in sets), '']
    for label, _ in sets:
        v = sorted(r['frac_missing'] for r in rows if r['set'] == label)
        L.append('  %-22s %6d alignments, median %s, mean %s' % (label, len(v), pct(statistics.median(v)), pct(statistics.mean(v))))
    L += ['', 'PER ALIGNMENT (pooled over sequences; %d alignments)' % n,
          '  median %s | mean %s | IQR %s - %s (index quartiles sorted[n//4], sorted[3n//4])' % (
              pct(statistics.median(af)), pct(statistics.mean(af)), pct(af[n // 4]), pct(af[3 * n // 4])),
          '  interpolated quartiles (statistics.quantiles, inclusive): %s - %s' % tuple(pct(x) for x in (statistics.quantiles(af, n=4, method='inclusive')[0], statistics.quantiles(af, n=4, method='inclusive')[2])),
          '  percentiles 5/25/50/75/95/99: ' + ' / '.join(pct(af[min(n - 1, int(p * n))]) for p in (.05, .25, .5, .75, .95, .99)),
          '  maximum %s (%s)' % (pct(af[-1]), max(rows, key=lambda r: r['frac_missing'])['alignment']),
          '  alignments with >10%% missing: %d | >20%%: %d' % (sum(x > .1 for x in af), sum(x > .2 for x in af))]
    sq = sorted(seq_all)
    ns = len(sq)
    L += ['', 'PER SEQUENCE (%d sequences)' % ns,
          '  median %s | mean %s | maximum %s' % (pct(statistics.median(sq)), pct(statistics.mean(sq)), pct(sq[-1])),
          '  sequences with >50%% missing: %d | >70%%: %d | >80%%: %d' % (sum(x > .5 for x in sq), sum(x > .7 for x in sq), sum(x > .8 for x in sq))]
    for frac, aln, name in sorted(seq_hi, reverse=True):
        L.append('    >80%%: %s in %s = %s' % (name, aln, pct(frac)))
    few = [(r['alignment'], r['n_sequences']) for r in rows if r['n_sequences'] < 5]
    L += ['', 'ALIGNMENTS WITH FEWER THAN 5 SEQUENCES: %d' % len(few)] + ['    %s (%d sequences)' % x for x in few]
    open(os.path.join(a.out, 'alignment_missingness_summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
