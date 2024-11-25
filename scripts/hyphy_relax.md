---
title: 'hyphy relax'
disqus: hackmd
---

Intensification or relaxation of selection with hyphy relax
===
By Tanya Lama

### Objective: 

### Acknowledgements
Thanks to Graham Hughes for running RELAX on each of our test/reference branch pairs

## Table of Contents

[TOC]

### Download the relax results from drive
Three sets: myoMyo myoNig; desRot dipEca; tadBra artJam

How many alignments were tested for each pair? (wc -l...)

### Pull the ENST.geneId and pvalue for each gene. Trying to determine how many of the total set of alignments tested have significantly relaxed or intensified selection.

Use parse.sh:
```
#!/bin/bash
#SBATCH --job-name=parse
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --partition=cpu
python3 parse.py
```

To execute parse.py:
```
import re
import glob
# Define the directory path and pattern for your files
file_pattern = "relax/*HLmyoMyo6_vs_HLmyoNig1_output.txt"
output_file = "myoMyo_relax_results_parsed.txt"
# Regular expressions to capture gene name and p-value
gene_regex = r"ENST[0-9]+\.([A-Za-z0-9]+)"
p_value_regex = r"\*\*p\s*=\s*([\d\.]+)\*\*"
# Open the output file for writing
with open(output_file, "w") as out_file:
    # Loop over each file matching the pattern
    for file_path in glob.glob(file_pattern):
        with open(file_path, "r") as file:
            content = file.read()
            
            # Search for gene name and p-value
            gene_match = re.search(gene_regex, file_path)
            p_value_match = re.search(p_value_regex, content)
            
            # If both are found, write to output file
            if gene_match and p_value_match:
                gene_name = gene_match.group(1)
                p_value = p_value_match.group(1)
                out_file.write(f"{gene_name}\t{p_value}\n")
```
### BH adjust pvalues 
Execute padj.py on a list of pvalues
```
import numpy as np
from statsmodels.stats.multitest import multipletests
# Input and output file names
input_file = "myo_pvalues.txt"  # Replace with your file containing p-values (one per line)
output_file = "bh_adjusted_pvalues.txt"
# Load p-values from the input file
with open(input_file, "r") as f:
    myo_pvalues = [float(line.strip()) for line in f if line.strip()]  # Read and convert to float
# Adjust p-values using the Benjamini-Hochberg procedure
bh_adjusted_pvalues = multipletests(myo_pvalues, method='fdr_bh')[1]
# Save the BH-adjusted p-values to the output file
with open(output_file, "w") as f:
    f.write("Original p-values,BH-adjusted p-values\n")
    for original, adjusted in zip(myo_pvalues, bh_adjusted_pvalues):
        f.write(f"{original},{adjusted}\n")
print(f"BH-adjusted p-values have been saved to {output_file}")
```
See bh_adjusted_pvalues.txt for output 

How many genes do you have left after BH correction?
n=262 total

### Parse selection intensity
For a list of geneIDs, parse whether the gene is under intensifying or relaxation of selection, print to file 
```
#!/bin/bash
# Define input file containing gene IDs
GENE_ID_FILE="myo_geneIDs"
# Define output file
OUTPUT_FILE="intensity_myo.txt"
# Clear or create the output file
> "$OUTPUT_FILE"
# Loop through each gene ID in the file
while read -r gene_id; do
  if [ -n "$gene_id" ]; then
    # Use a wildcard to find the file matching the gene_id
    for file in relax/*"${gene_id}".cleanLb_hmm_manual_RELAX_HLmyoMyo6_vs_HLmyoNig1_output.txt; do
      if [ -f "$file" ]; then
        # Run the grep command and format the output with the gene_id
        grep ">Evidence for *" "$file" | while read -r line; do
          echo -e "${gene_id}\t${line}" >> "$OUTPUT_FILE"
        done
      fi
    done
  fi
done < "$GENE_ID_FILE"
echo "Parsing completed. Results written to $OUTPUT_FILE"
```
See intensity_myo.txt for output
You should now have an excel sheet with geneID, original pvalue, BH adjusted pvalue, and selection intensity for padj<0.05 genes. 

We returned to the /relax folder to pull k= values for genes discussed in the text. We also included the entrez summary in Supplementary Table 12 for each gene, and a note regarding whether each gene is associated with aging hallmarks or a priori lists (e.g., antagonistic pleiotropies of aging).

You can now move foward with investigations of invidividual genes




