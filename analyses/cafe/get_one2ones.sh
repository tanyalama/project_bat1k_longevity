#!/bin/bash

# Loop over each TSV file in list_orthology_tsvs
while read -r file_path; do
    # Extract the characters 3-9 from the file name
    file_name=$(basename "$file_path")
    file_prefix=${file_name%_orthology_classification.tsv}

    # Extract the q_transcript values for rows where orthology_class == 'one2one'
    awk -F $'\t' '$5 == "one2one" { print $4 }' "$file_path" | sort -u > "one2one_${file_prefix}.txt"
done < list_orthology_tsvs
