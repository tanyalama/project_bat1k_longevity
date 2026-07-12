import sys
from Bio import SeqIO

# Check for command-line argument
if len(sys.argv) < 2:
    print("Please provide the species name as a command line argument")
    exit()

species = sys.argv[1]

input_file = f"{species}.prot.fasta"
output_file = f"query_{species}.pep"

# Filter records that contain both "PROT" and "QUERY" in the description
filtered_records = (
    record for record in SeqIO.parse(input_file, "fasta")
    if "PROT" in record.description and "QUERY" in record.description
)

# Write filtered records to output FASTA
with open(output_file, "w") as out_f:
    SeqIO.write(filtered_records, out_f, "fasta")

print(f"Filtered FASTA written to: {output_file}")
# The output of this script is a protein fasta called query_${species}.pep that retains ONLY protein query sequences
