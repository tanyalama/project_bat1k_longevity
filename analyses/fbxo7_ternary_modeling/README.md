# FBXO7 bat-variant modeling — FBXO7–PINK1–PSMF1 ternary complex

Structural modeling and interface energetics for bat-lineage FBXO7 variants in the
FBXO7–PINK1(kinase domain)–PSMF1 ternary complex, supporting the manuscript
*"Bat adaptations in mitochondrial damage recognition revealed through comparative multi-omics."*

## Overview

FBXO7 variants observed across four bat lineages were modeled onto a Boltz-2 ternary
complex template and scored for their effect on protein–protein binding at two interfaces:

- **A_B** — FBXO7–PINK1 (kinase domain)
- **A_C** — FBXO7–PSMF1

The workflow has two stages:

1. **Structure prediction (Boltz-2).** Each variant sequence (mutations encoded directly in
   chain A) is folded with Boltz-2 using the reference complex `model_0.cif` as a forced
   template (backbone-deviation threshold 5 Å) so all models share the same global
   architecture. Mutations are introduced *in the Boltz-2 sequence*, not in Rosetta.
   Ten diffusion samples are drawn per target; the best-scoring model (`model_0`, ranked by
   Boltz `confidence_score`) is carried forward.

2. **Interface energetics (Rosetta / PyRosetta).** Each Boltz model undergoes constrained
   FastRelax (Cα harmonic constraints, sigma = 0.5 Å, ref2015, 5 independent seeds; the
   lowest-energy decoy is kept), then `InterfaceAnalyzerMover` computes per-interface
   `dG_separated`. ddG = dG(variant) − dG(WT), with WT generated through the identical
   Boltz to Rosetta pipeline for a self-consistent reference.

## Variant set

Four bat lineages, each a combination of FBXO7 substitutions (positions on FBXO7, UniProt Q9Y3I1):

| Lineage | Substitutions |
|---|---|
| Myotis myotis | T19E; T47A; A64T; S109P; Q127E; F146V; G409I |
| Myotis nigricans | T19E; T47A; S109P; Q127E; F146V; G409I |
| Desmodus rotundus | T19E; T47A; S109P; S110H; Q127E; F146V; D191G; L290P; G409R |
| Diphylla ecaudata | T19E; T47E; S109P; S110C; Q127E; F146V; L290P; E292R; G409R |

Individual single substitutions (13 unique, excluding S109P which is analyzed separately):
T19E, T47A, T47E, A64T, S110H, S110C, Q127E, F146V, D191G, L290P, E292R, G409I, G409R.

## Repository layout

```
fbxo7_ternary_modeling/
  boltz_configs/       Boltz-2 YAML input for each target (WT + 17 variants);
                       chain A = (mutated) FBXO7, B = PINK1 kinase, C = PSMF1;
                       forced template model_0.cif, threshold 5.0
  scripts/
    relax_score_boltz_decoy.py   Constrained FastRelax (5 seeds) + InterfaceAnalyzer ddG
    score_interface_only.py      Interface analysis only (re-score a relaxed structure)
  structures/
    template/          FBXO7_PINK1kinase_PSMF1_model_0.cif  (stage-1 reference complex)
    boltz_models/      Best-scoring model_0 PDB for WT and each variant;
                       S109P provided under both selection strategies
                       (best-scoring model_0 and best-templated model_3)
  results/             Rosetta interface metrics (dG_separated, ddG, dSASA, unsat H-bonds,
                       packstat, nres) per interface
```

## Chains

| Chain | Protein | UniProt | Length |
|---|---|---|---|
| A | FBXO7 | Q9Y3I1 | 522 |
| B | PINK1 (kinase domain) | Q9BXM7 | 356 |
| C | PSMF1 | Q92530 | 270 |

## Software

- Boltz-2 (v2.2.1), run with `--use_msa_server --no_kernels --recycling_steps 3
  --diffusion_samples 10`. `--no_kernels` disables the fused cuequivariance CUDA kernels
  (numerically equivalent; used because that package was unavailable in the environment).
- PyRosetta (ref2015 score function), constrained FastRelax + InterfaceAnalyzerMover.

## Usage

```bash
# 1. Predict a structure (on a GPU node)
boltz predict boltz_configs/S109P.yaml --use_msa_server --no_kernels \
    --recycling_steps 3 --diffusion_samples 10 --output_format pdb

# 2. Relax + score interfaces (ddG vs WT reference)
python scripts/relax_score_boltz_decoy.py <model_0.pdb> <sample_name>
```

## Notes

- ddG is reported in Rosetta Energy Units (REU). Negative ddG = stabilizing.
- WT reference is produced by the same Boltz to Rosetta pipeline as the variants so that the
  dG(variant) − dG(WT) subtraction is protocol-consistent.
