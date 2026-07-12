#!/bin/bash
# NOTE: run with the compleasm conda environment activated; download the mammalia DB first (compleasm download mammalia)
#SBATCH -J compleasm_protein_longest_allSpp
#SBATCH -o ./logs/all_spp/%x_%a_%A.log
#SBATCH -e ./logs/all_spp/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 02:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=64G  # Requested Memory
#SBATCH --array=1-30

SPP=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./all_spp.txt)
echo $SPP

compleasm protein -p /scratch3/workspace/bbentley_smith_edu-busco/refs/all_spp/cleaned/primary_transcripts/${SPP}.reform.query.prot.fasta -l mammalia -o compleasm/all_spp/${SPP} -t 16
