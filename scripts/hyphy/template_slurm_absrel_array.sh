#!/bin/bash
#SBATCH --job-name=absrel_array
#SBATCH --output=absrel_logs/absrel_%A_%a.out   # Log file per task (create 'logs' dir first)
#SBATCH --error=absrel_logs/absrel_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --array=0-1500         # Adjust based on number of lines in batch${i}.txt

# Load any necessary modules
module load conda/latest gcc/9.4.0
conda activate hyphy

# Read the input file corresponding to the SLURM array task ID
BATCH="batch<NUM>"
WORK_DIR=/myworkingdirectory/hyphy
LIST_FILE=$WORK_DIR/${BATCH}.txt
INPUT_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $LIST_FILE)

cd $WORK_DIR/$BATCH
#mkdir -p results/
echo "Processing file: $INPUT_FILE"

NAME="${INPUT_FILE%.exon*}"

hyphy absrel CPU=1 --pvalue 0.05 --alignment $INPUT_FILE --tree $INPUT_FILE.labeled.tfl --branches FG --srv Yes >> $WORK_DIR/absrel_hyp/$BATCH/${NAME}_absrel_output.txt
