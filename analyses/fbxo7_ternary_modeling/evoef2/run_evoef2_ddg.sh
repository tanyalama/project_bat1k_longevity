#!/bin/bash
# EvoEF2 interface ddG pipeline for the FBXO7-PINK1-PSMF1 ternary complex.
#
# Computes binding ddG for 14 single-site FBXO7 variants and 4 bat-species
# multi-mutant site-sets, at both interfaces:
#   AB = FBXO7 (chain A) - PINK1 kinase (chain B)
#   AC = FBXO7 (chain A) - PSMF1        (chain C)
#
# Protocol (Huang lab EvoEF2, github.com/tommyhuangthu/EvoEF2):
#   RepairStructure -> ComputeBinding(WT) -> BuildMutant(all in one file)
#   -> ComputeBinding(each mutant); ddG = dG_bind(mut) - dG_bind(WT).
#   EvoEF2 binding energy is more negative for stronger binding, so
#   negative ddG = stabilizing, positive ddG = destabilizing (same sign
#   convention as the Rosetta InterfaceAnalyzer ddG in this project).
#
# NOTE: EvoEF2 ComputeBinding + local BuildMutant repacking only registers
# residues in direct interfacial contact; distal substitutions read as 0.0.
# This is sparser than Rosetta FastRelax, which redistributes strain to
# distal sites. Compare rank/sign at contact residues, not absolute values.
#
# Inputs (per interface subdir evoef2_run/{AB,AC}/):
#   complex.pdb        two-chain complex extracted from WT_relaxed_best.pdb
# Mutant files (chain A, EvoEF2 syntax {WT}A{pos}{MUT}; one line per variant):
#   evoef2_mutants.txt         14 single-site variants
#   evoef2_species_mut.txt     4 species (comma-joined substitutions per line)
# Output: {AB,AC}_results.txt and evoef2_species_results.txt
#         (columns: variant iface wt_bind mut_bind ddG)
set -e
BASE=$(cd "$(dirname "$0")" && pwd)
EVO=$BASE/EvoEF2/EvoEF2
WORK=$BASE/evoef2_run
VARIANTS="T19E T47A T47E A64T S109P S110H S110C Q127E F146V D191G L290P E292R G409R G409I"
SPECIES="Myotis_myotis Myotis_nigricans Desmodus_rotundus Diphylla_ecaudata"

# ---- single-site variants -------------------------------------------------
run_single () {
  DIR=$1
  cd $WORK/$DIR
  cp $BASE/evoef2_mutants.txt mutants.txt
  [ -f complex_Repair.pdb ] || $EVO --command=RepairStructure --pdb=complex.pdb > repair.log 2>&1
  $EVO --command=ComputeBinding --pdb=complex_Repair.pdb > wt_bind.log 2>&1
  WT=$(grep '^Total' wt_bind.log | tail -1 | awk '{print $NF}')
  $EVO --command=BuildMutant --pdb=complex_Repair.pdb --mutant_file=mutants.txt > buildmut.log 2>&1
  echo "variant iface wt_bind mut_bind ddG" > ${DIR}_results.txt
  i=0
  for v in $VARIANTS; do
    i=$((i+1)); M=$(printf "complex_Repair_Model_%04d.pdb" $i)
    $EVO --command=ComputeBinding --pdb=$M > mb_${v}.log 2>&1
    MB=$(grep '^Total' mb_${v}.log | tail -1 | awk '{print $NF}')
    DD=$(python3 -c "print(round($MB-($WT),2))")
    echo "$v $DIR $WT $MB $DD" >> ${DIR}_results.txt
  done
  echo "=== $DIR single-site done ==="; cat ${DIR}_results.txt
}

# ---- species multi-mutant site-sets ---------------------------------------
run_species () {
  echo "species iface wt_bind mut_bind ddG" > $BASE/evoef2_species_results.txt
  for DIR in AB AC; do
    SUB=$WORK/$DIR/species; rm -rf "$SUB"; mkdir -p "$SUB"; cd "$SUB"
    cp ../complex_Repair.pdb .
    cp $BASE/evoef2_species_mut.txt species_mut.txt
    $EVO --command=ComputeBinding --pdb=complex_Repair.pdb > wt_bind.log 2>&1
    WT=$(grep '^Total' wt_bind.log | tail -1 | awk '{print $NF}')
    $EVO --command=BuildMutant --pdb=complex_Repair.pdb --mutant_file=species_mut.txt > buildmut.log 2>&1
    i=1
    for SP in $SPECIES; do
      M=$(printf "complex_Repair_Model_%04d.pdb" $i)
      $EVO --command=ComputeBinding --pdb="$M" > "mb_${SP}.log" 2>&1
      MB=$(grep '^Total' "mb_${SP}.log" | tail -1 | awk '{print $NF}')
      DD=$(python3 -c "print(round($MB-($WT),2))")
      echo "$SP $DIR $WT $MB $DD" >> $BASE/evoef2_species_results.txt
      i=$((i+1))
    done
  done
  echo "=== species done ==="; cat $BASE/evoef2_species_results.txt
}

run_single AB
run_single AC
run_species
echo "=== ALL DONE ==="
