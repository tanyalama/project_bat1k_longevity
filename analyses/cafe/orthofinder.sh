#!/bin/bash
#SBATCH --output=ortho.out
#SBATCH --error=mv-ortho.err
#SBATCH --job-name=ortho
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --ntasks=20
#SBATCH --partition=cpu

module load conda/latest
conda activate compleasm

WORK_DIR=/scratch/workspace/nmcdonald1_access-ci_org-shared/nikki_work/cafe/6_orthofinder

orthofinder -a 4 -f $WORK_DIR/input_faas -o $WORK_DIR/results
