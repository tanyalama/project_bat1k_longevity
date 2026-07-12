import sys
from Bio import SeqIO

species = sys.argv[1]

# ---------------------------------------------
# Step 1: Load exact transcript IDs from the one-to-one file
# Example: ENST00000355498.74
# ---------------------------------------------
with open(f'one2one_{species}.txt') as f:
    keep_ids = set(line.strip() for line in f if line.strip())

# ---------------------------------------------
# Step 2: Parse the FASTA file and classify entries
# ---------------------------------------------
input_fasta = f'../2_query_protein_fastas/query_{species}.pep'
matched_fasta = f'ensts_{species}.txt'
final_fasta = f'{species}_final.pep'

matched_records = []
unmatched_records = []

for record in SeqIO.parse(input_fasta, 'fasta'):
    # record.id will be "ENST00000355498.74"
    if record.id in keep_ids:
        matched_records.append(record)
    else:
        unmatched_records.append(record)

# ---------------------------------------------
# Step 3: Write the output FASTA files
# ---------------------------------------------
# Write matched entries to ensts_{species}.txt
SeqIO.write(matched_records, matched_fasta, 'fasta')
# Write remaining entries to {species}_final.pep
SeqIO.write(unmatched_records, final_fasta, 'fasta')
