import re           # For regular expressions
import glob         # For file pattern matching
import subprocess   # To run the external R script

### Step 1: Extract species names from the tree file
# Initialize an empty list to hold species names
species = []

# Open the tree file (Newick format) for reading
with open("supermatrix_partition_iqtree_reroot_no_bootstraps.tree", "r") as tree_file:
    for line in tree_file:
        # Use regex to find all species names that match the pattern before a colon (e.g., Homo123:)
        #matches = re.findall(r'([A-Za-z]+[0-9]+):', line)
        #matches = re.findall(r'([A-Za-z]+[0-9]+[A-Za-z]*):', line)
        matches = re.findall(r'([A-Za-z][A-Za-z0-9_]*):', line)
        species.extend(matches)

# Remove duplicates to get a unique list of species
species = list(set(species))

### Check which species are present in each alignment file
# Open the output file where missing species info will be written
with open("details.txt", "w") as out_file:

    # Get a list of all alignment files with .fa, .fas, or .fasta extensions
    fasta_files = glob.glob("*.fa") + glob.glob("*.fas") + glob.glob("*.fasta")

    # Loop over each FASTA alignment file
    for fasta in fasta_files:
        gene = fasta.strip()  # Use the filename as the gene name
        out_file.write(gene)  # Write the gene name to the output file

        keep = []  # List to hold species found in the current alignment

        # Open and parse the FASTA file
        with open(fasta, "r") as f:
            for line in f:
                if line.startswith(">"):  # Lines starting with '>' contain species names
                    # Extract species name (first word after '>')
                    sp = line.strip().lstrip(">").split()[0]
                    keep.append(sp)

        # Compare the global species list with those present in this FASTA file
        for sp in species:
            if sp not in keep:
                # If the species is missing from the alignment, write it to the output
                out_file.write(f"|{sp}")
        
        # Add newline after processing each file
        out_file.write("\n")

### Step 3: Executes the R script named 'cliptrees.R'
#subprocess.run(["Rscript", "cliptrees.R"])
