python!
# Step 7: Convert merged_output_curated.csv to CAFE input format
# Now we need to convert the merged_output_curated.csv into CAFE format.
# CAFE requires an input file in a tab-separated format with gene IDs
# (in this case, Orthogroup IDs) in the first column and gene counts
# (number of copies) for each species in subsequent columns.

import pandas as pd

# Function to count genes in a cell or assign 0 to empty cells
def count_genes(cell_value):
    if pd.isna(cell_value) or cell_value.strip() == "":
        return "0"
    return str(len(cell_value.split(", ")))


# File paths
input_csv_file = 'filtered_merged_orthogroups.tsv'
output_cafe_file = 'cafe_input_file_w_olfactory.tsv'
#input_csv_file = 'no_olfactory_filtered_merged_orthogroups.tsv'
#output_cafe_file = 'cafe_input_file_no_olfactory.tsv'

# Read the CSV file
df = pd.read_csv(input_csv_file, sep='\t')

# Replace empty cells with "0" and count genes in non-empty cells
for col in df.columns[1:]:
    df[col] = df[col].fillna("").apply(count_genes)

# Prepare the CAFE input file in tab-separated format
with open(output_cafe_file, 'w') as f:
    f.write("Gene_ID\t" + "\t".join(df.columns[1:]) + "\n")
    for _, row in df.iterrows():
        f.write(row["Orthogroup"] + "\t" + "\t".join(row[1:]) + "\n")

print("CAFE input file has been created:", output_cafe_file)
