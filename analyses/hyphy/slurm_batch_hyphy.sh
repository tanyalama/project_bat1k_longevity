#!/bin/bash
#SBATCH --job-name=batch_hyphy
#SBATCH --output=batch_hyphy_%j.out
#SBATCH --error=batch_hyphy_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -p cpu

# === Arguments ===
# $1 - INPUT_LIST: Path to the input file listing all input files, one per line.

# === Configuration ===
MAX_REMAINING=500
JOB_LIST="$1"
PREV_JOB_ID=""

if [[ -z "$JOB_LIST" || ! -f "$JOB_LIST" ]]; then
    echo "Usage: $0 <job_list_file>"
    exit 1
fi

while IFS= read -r JOB_SCRIPT; do
    if [[ ! -f "$JOB_SCRIPT" ]]; then
        echo "Skipping missing job script: $JOB_SCRIPT"
        continue
    fi

    # Wait for previous job's array to drop below threshold
    if [[ -n "$PREV_JOB_ID" ]]; then
        echo "Waiting for fewer than $MAX_REMAINING total remaining tasks from job $PREV_JOB_ID..."

        # Get total number of array tasks from original job submission
        TOTAL_TASKS=$(sacct -j "$PREV_JOB_ID" --format=JobID,JobName%20,State --noheader -X | grep -E "^${PREV_JOB_ID}_[0-9]+" | wc -l)

        if [[ "$TOTAL_TASKS" -eq 0 ]]; then
            # fallback: try to parse from scontrol
            TOTAL_TASKS=$(scontrol show job "$PREV_JOB_ID" | grep -oP 'NumTasks=\K[0-9]+')
        fi

        [[ -z "$TOTAL_TASKS" || "$TOTAL_TASKS" -eq 0 ]] && TOTAL_TASKS=99999  # fail-safe

        while :; do
            # Count completed tasks
            COMPLETED_TASKS=$(sacct -j "$PREV_JOB_ID" --format=JobID,State --noheader -X \
                | grep -E "^${PREV_JOB_ID}_[0-9]+" | grep -E 'COMPLETED|CANCELLED|FAILED|TIMEOUT|BOOT_FAIL|PREEMPTED' | wc -l)

            REMAINING=$((TOTAL_TASKS - COMPLETED_TASKS))
            [[ $REMAINING -lt 0 ]] && REMAINING=0

            echo "  Remaining array tasks: $REMAINING"

            if (( REMAINING < MAX_REMAINING )); then
                echo "Below threshold; submitting next job."
                break
            else
                echo "Waiting 3 minutes..."
                sleep 180
            fi
        done
    fi

    # Submit the job script (which already includes #SBATCH --array)
    JOB_OUTPUT=$(sbatch "$JOB_SCRIPT")
    echo "Submitted: $JOB_SCRIPT → $JOB_OUTPUT"

    # Extract new job ID
    PREV_JOB_ID=$(echo "$JOB_OUTPUT" | awk '{print $4}')
    sleep 10
done < "$JOB_LIST"
