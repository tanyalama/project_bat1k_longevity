#!/usr/bin/env python
## this is a script which pads q_gene names with 00s such that ever name has 5 digits. e.g., reg_10 is now reg_00010. We can then search for the q_gene name itself in 1_longest_transcript.sh (next step) without worrying about confusion over reg_10 reg_100 and reg_1000 being treated as matching q_gene names 

import sys
import re
import pandas as pd

# Get the value for species from the command line argument
species = sys.argv[1]

def pad_gene_names(gene_names):
    padded_names = []
    for name in gene_names:
        if pd.isna(name):
            padded_names.append(name)  # Keep NaNs as-is
            continue
        name = str(name)
        match = re.match(r'(\D+)(\d+)', name)
        if match:
            prefix = match.group(1)
            number = match.group(2)
            padded_number = number.rjust(5, '0')
            padded_name = prefix + padded_number
            padded_names.append(padded_name)
        else:
            padded_names.append(name)
    return padded_names

# Read in your file as a DataFrame
df = pd.read_csv(f'{species}_orthology_classification.tsv', sep='\t')

# Apply the function to all values in column 3
df['q_gene'] = pad_gene_names(df['q_gene'])

# Write the modified DataFrame back to a new file
df.to_csv(f'./modified/{species}_orthology_classification.tsv', sep='\t', index=False)
