#!/bin/bash

#SBATCH --job-name=cafe1
#SBATCH --output=cafe1.out
#SBATCH --error=cafe1.out
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH -p cpu

module load conda/latest
conda activate cafe

# Run CAFE with cafe_input_w_or.tsv
echo "Running CAFE with OR and lambda1"
cafe5 -i cafe_input_w_or.tsv -t ultrametric_tree_clipped.nwk -eBase_error_model.txt -o cafe1_run1 -y tree_lambda1.txt
echo "Running CAFE with OR and lambda2"
cafe5 -i cafe_input_w_or.tsv -t ultrametric_tree_clipped.nwk -eBase_error_model.txt -o cafe1_run2 -y tree_lambda2.txt
echo "Running CAFE with OR and lambda3"
cafe5 -i cafe_input_w_or.tsv -t ultrametric_tree_clipped.nwk -eBase_error_model.txt -o cafe1_run3 -y tree_lambda3.txt
