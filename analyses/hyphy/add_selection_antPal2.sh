#!/bin/bash

# Input: BH-adjusted results file
INPUT_FILE="relax_results_bh_adjusted.txt"
OUTPUT_FILE="relax_results_bh_adjusted_with_selection.txt"

# Clear or create the output file
> "$OUTPUT_FILE"

# Read input file line-by-line
while IFS=$'\t' read -r gene_id gene_name p_value bh_adj; do
  # Handle header
  if [[ "$gene_id" == "gene_id" ]]; then
    echo -e "${gene_id}\t${gene_name}\t${p_value}\t${bh_adj}\tselection_type\tk_value" >> "$OUTPUT_FILE"
    continue
  fi

  selection=""
  k_value=""

  for file in output/*"${gene_id}"*output.txt; do
    if [[ -f "$file" ]]; then
      # Get the line that starts with '>Evidence for ...'
      line=$(grep '^>Evidence for' "$file")

      if [[ -n "$line" ]]; then
        selection=$(echo "$line" | grep -Eo 'intensification of selection|relaxation of selection' | head -n1)
      fi

      # Get the last K value line and extract the number after '='
      k_line=$(grep "Relaxation/intensification parameter (K) =" "$file" | tail -n1)
      if [[ -n "$k_line" ]]; then
        #k_value=$(echo "$k_line" | awk -F'= ' '{print $2}')
        k_value=$(echo "$k_line" | awk -F'= ' '{sub(/^[ \t]+/, "", $2); print $2}')
      fi

      break  # Stop after first matching file
    fi
  done

  echo -e "${gene_id}\t${gene_name}\t${p_value}\t${bh_adj}\t${selection}\t${k_value}" >> "$OUTPUT_FILE"

done < "$INPUT_FILE"

echo "Parsing completed. Results written to $OUTPUT_FILE"
