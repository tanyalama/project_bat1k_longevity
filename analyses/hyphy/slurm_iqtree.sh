#! /bin/bash
  
#SBATCH --job-name=iqtree-partitioned
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --mem=100G #15.5k genes needed min 47G
#SBATCH --time=24:00:00 #ran for 7 hours with 15.5k genes
#SBATCH -p cpu

iqtree -s *_supermatrix_partition.nex -p *_supermatrix_partition.txt --prefix supermatrix_partition_iqtree -B 1000 -T AUTO
