# EvoEF2 interface ddG — FBXO7-PINK1-PSMF1 ternary complex

Physics-based binding ddG for FBXO7 bat-adaptive variants, computed with
[EvoEF2](https://github.com/tommyhuangthu/EvoEF2) (Huang lab) as an independent
cross-check on the Rosetta InterfaceAnalyzer ddG (`../scripts/plot_ddg_heatmap.py`).

## Interfaces
- **AB** — FBXO7 (chain A) x PINK1 kinase (chain B)
- **AC** — FBXO7 (chain A) x PSMF1 (chain C)

Both two-chain complexes were extracted from the FastRelax-relaxed WT ternary
model (`WT_relaxed_best.pdb`).

## Protocol
`RepairStructure` -> `ComputeBinding` (WT) -> `BuildMutant` (all variants in one
mutant file) -> `ComputeBinding` (each mutant). ddG = dG_bind(mut) - dG_bind(WT).

EvoEF2 binding energy is more negative for stronger binding, so **negative ddG =
stabilizing, positive ddG = destabilizing** — the same sign convention as the
Rosetta ddG heatmap in this project.

## Scope
- 14 single-site variants: T19E, T47A, T47E, A64T, S109P, S110H, S110C, Q127E,
  F146V, D191G, L290P, E292R, G409R, G409I
- 4 bat-species multi-mutant site-sets: *Myotis myotis*, *Myotis nigricans*,
  *Desmodus rotundus*, *Diphylla ecaudata*

## Files
- `run_evoef2_ddg.sh` — full pipeline (build complexes must pre-exist in
  `evoef2_run/{AB,AC}/complex.pdb`; EvoEF2 binary in `EvoEF2/`)
- `evoef2_mutants.txt` — 14 single-site mutations, EvoEF2 syntax `{WT}A{pos}{MUT};`
- `evoef2_species_mut.txt` — 4 species site-sets (comma-joined substitutions/line)
- `plot_evoef2_heatmap.py` — diverging heatmap (green=stabilizing, purple=
  destabilizing, white=|ddG|<=0.5), Arial, matches the Rosetta heatmap styling
- `evoef2_ddg_combined.csv` — long-format results (sample, interface, wt_bind,
  mut_bind, ddG, method)
- `fbxo7_evoef2_ddg_supp_table.csv` — wide supplementary table
- `fbxo7_evoef2_ddg_heatmap.png` — rendered figure

## Key result
Signal at the **FBXO7-PINK1** interface concentrates on the PINK1-pocket
residues: S109P (-4.0), S110C (-1.25), S110H (-1.22), Q127E (-0.21) stabilizing;
D191G (+4.58) destabilizing. Species combos: *Diphylla ecaudata* -5.21,
*Myotis myotis* / *M. nigricans* -4.21 (identical — they differ only at the
non-contact residue A64T), *Desmodus rotundus* -0.59.

The **FBXO7-PSMF1** interface reads 0.0 for every substitution: none of the 14
positions makes direct side-chain contact with PSMF1 in the relaxed complex.

## Method comparison caveat
EvoEF2 `ComputeBinding` + local `BuildMutant` repacking only registers residues
in direct interfacial contact, so distal substitutions read as exactly 0.0. This
is sparser than Rosetta FastRelax, which redistributes strain to distal sites.
Compare **rank and sign at contact residues**, not absolute magnitudes — the two
methods use different energy units and scales.
