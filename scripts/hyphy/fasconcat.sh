#! /bin/bash

#SBATCH --job-name=fasconcat
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --mem=200G #probably don't need this much, idk maybe 80G?
#SBATCH --time=04:00:00 #took 1.5 hours with 15.5k genes
#SBATCH -p cpu

#perl ./FASconCAT-G_v1.05.pl -n -p -l -s 
perl ./FASconCAT-G_v1.06.1.pl -n -p -l -s
