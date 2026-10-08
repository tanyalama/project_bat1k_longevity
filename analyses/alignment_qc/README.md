# alignment_qc: missing data in the 30-species codon alignments

## Objective

Quantify how much missing data (gaps and masked positions) the 30-species codon alignments contain, for the Methods
statement on alignment completeness.

## Inputs

Codon alignments (`*.fasta`, 23 bat and 7 non-bat species) from the TOGA ortholog, MACSE v2 exon alignment and HmmCleaner
pipeline, in two directories under the project data root:

| Directory | Alignments |
|---|---|
| `final_alignment/fastas` | 10,893 |
| `7082_fastas` | 7,077 |

`7082_fastas` also contains `FASconCAT-G_v1.05.pl` and `fasconcat.sh`; they are not alignments and are ignored (only `*.fasta` is read).

## Pipeline

[`scripts/alignment_missingness.py`](scripts/alignment_missingness.py) (Python 3 standard library only):

```
python3 scripts/alignment_missingness.py --root <project data root> --out results
```

Definitions:

* A position is missing if it is anything other than A, C, G or T after upper-casing. In these files that is N: alignment gaps and
  HmmCleaner-masked positions both appear as N and cannot be told apart.
* Per-sequence missing fraction = missing positions / sequence length.
* Per-alignment missing fraction = total missing positions / total positions, pooled over all sequences in the alignment. It is
  length-weighted and is not the mean of the per-sequence fractions.
* The interquartile range uses the index convention `sorted[n//4]`, `sorted[3n//4]`; the interpolated quartiles are also reported and
  agree to two decimals.

## Outputs

* [`results/alignment_missingness_summary.txt`](results/alignment_missingness_summary.txt): summary statistics.
* [`results/alignment_missingness_per_alignment.csv`](results/alignment_missingness_per_alignment.csv): one row per alignment
  (set, file, sequences, positions, missing positions, missing fraction, maximum per-sequence fraction, sequences above 50% and 80%).

## Results

| | Value |
|---|---|
| Alignments | 17,970 (10,893 + 7,077) |
| Missing data per alignment, median | 0.27% |
| Missing data per alignment, mean | 1.06% |
| Missing data per alignment, interquartile range | 0.06-0.92% |
| Missing data per alignment, 95th / 99th percentile | 4.66% / 12.87% |
| Most incomplete alignment | 41.43% (`ENST00000303766.RBMY1F`) |
| Alignments above 10% / 20% missing | 285 / 53 |
| Sequences | 501,221 |
| Missing data per sequence, median / mean | 0.00% / 0.83% |
| Sequences above 50% / 70% / 80% missing | 379 / 60 / 3 |

By directory, the per-alignment median and mean are 0.29% and 0.88% for `final_alignment/fastas` and 0.24% and 1.33% for `7082_fastas`.

## Notes

* Three retained sequences are just above 80% missing: `HLmolMol2` in `LCE2C` (80.23%) and `felCat9` in `RP11-402P6.15` and `CXorf49`
  (80.08% each). The per-sequence cutoff therefore appears to be about 80%, not a strict "more than 80%".
* Twelve alignments in the analysed set have only 3 or 4 sequences (`ZNF559`, `RFPL3`, `IGFL3`, `BTNL3`, `LGALS14`, `CT47A7`,
  `CT47A5`, `IGHD6-25`, `NPIPB8`, `CXorf49B`, `ZNF790`, `OR8G5`).
* 17,970 alignment files were found on disk against 17,975 stated in the manuscript Methods; the difference of five is unexplained.
