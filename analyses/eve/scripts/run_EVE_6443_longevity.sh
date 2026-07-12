#!/bin/bash
#SBATCH -J 6443_longevity
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 4:00:00
#SBATCH --mem=10000  # Requested Memory
#SBATCH --ntasks-per-node=40

#running 6443_longevity_tpm.txt 
# -o 7 tests for shifts in species 8 (myoMyo)
/work/pi_tlama_smith_edu/bin/EVE_release/EVEmodel -O -o 7 -n 6443 -t data/renamed_longevity_tree.nh -i data/longevity_num_reorder.txt -d data/6443_longevity_tpm.txt -f _6443_longevity -v 10
