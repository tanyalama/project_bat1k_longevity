#!/bin/bash
#SBATCH --job-name=meme_array
#SBATCH --output=logs/meme_%A_%a.out   # Log file per task (create 'logs' dir first)
#SBATCH --error=logs/meme_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --array=0-1166         # Adjust based on number of lines in alignment_list.txt

# Load any necessary modules
module load conda/latest gcc/9.4.0
conda activate hyphy

# Read the input file corresponding to the SLURM array task ID
WORK_DIR=/myworkingdirectory/hyphy/meme
LIST_FILE=$WORK_DIR/alignment_list.txt
INPUT_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $LIST_FILE)

cd $WORK_DIR/fastas
#mkdir -p results/
echo "Processing file: $INPUT_FILE"

NAME="${INPUT_FILE%.exon*}"

hyphy meme CPU=1 --pvalue 0.05 --alignment $INPUT_FILE --tree $INPUT_FILE.tfl --srv Yes >> $WORK_DIR/results/${NAME}_MEME_output.txt

## --output /work/###_MEME_output.json
