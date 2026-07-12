import numpy as np
from statsmodels.stats.multitest import multipletests

# Input and output file names
input_file = "relax_results_parsed.txt"   # From previous script
output_file = "relax_results_bh_adjusted.txt"   # Output file with BH-adjusted p-values

gene_ids = []
gene_names = []
p_values = []

# Load data, skipping the header
with open(input_file, "r") as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 3:
            gene_id, gene_name, p_val_str = parts
            try:
                p_val = float(p_val_str)
                gene_ids.append(gene_id)
                gene_names.append(gene_name)
                p_values.append(p_val)
            except ValueError:
                continue  # Skip if p-value is not a valid float

# Perform BH correction
bh_adjusted = multipletests(p_values, method='fdr_bh')[1]

# Write results
with open(output_file, "w") as f:
    f.write("gene_id\tgene_name\tp_value\tBH_adjusted_p_value\n")
    for gene_id, gene_name, orig_p, adj_p in zip(gene_ids, gene_names, p_values, bh_adjusted):
        f.write(f"{gene_id}\t{gene_name}\t{orig_p}\t{adj_p}\n")

print(f"BH-adjusted p-values saved to: {output_file}")
