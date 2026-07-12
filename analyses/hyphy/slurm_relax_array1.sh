#!/bin/bash
#SBATCH --job-name=relax_array
#SBATCH --output=logs/relax1_%A_%a.out   # Log file per task (create 'logs' dir first)
#SBATCH --error=logs/relax1_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --array=0-1500         # Adjust based on number of lines in alignment_list.txt

# Load any necessary modules
module load conda/latest gcc/9.4.0
conda activate hyphy

# Read the input file corresponding to the SLURM array task ID
BATCH="batch1"
WORK_DIR=/myworkingdirectory/hyphy
LIST_FILE=$WORK_DIR/${BATCH}.txt
INPUT_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $LIST_FILE)

mkdir -p $WORK_DIR/relax/$BATCH

cd $WORK_DIR/$BATCH
echo "Processing file: $INPUT_FILE"

NAME="${INPUT_FILE%.exon*}"

hyphy relax CPU=1 --pvalue 0.05 --alignment $INPUT_FILE --tree $INPUT_FILE.labeled1.tfl --srv Yes --test FG1 --reference FG2 >> $WORK_DIR/relax/$BATCH/${NAME}_RELAX_output.txt

## --output /work/###_RELAX_output.json
