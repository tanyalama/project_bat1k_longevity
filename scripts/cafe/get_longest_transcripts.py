import sys
import os
from collections import defaultdict
from Bio import SeqIO

species = sys.argv[1]

# Input files
orthology_file = f"../3_orthology_tsvs/modified/{species}_orthology_classification.tsv"
fasta_file = f"../4_one2ones/{species}_final.pep"

# Output directories
js_dir = f"{species}_js"
finaljs_dir = f"{species}_finalJs"
os.makedirs(js_dir, exist_ok=True)
os.makedirs(finaljs_dir, exist_ok=True)

# ---------------------------------------------
# Step 1: Parse orthology file to map genes → transcripts
# and record full rows for later filtering
# ---------------------------------------------
gene_to_transcripts = defaultdict(list)
gene_orthology_lines = defaultdict(list)  # For paralog filtering and orthology checks

with open(orthology_file) as f:
    for line in f:
        fields = line.strip().split()
        if len(fields) >= 5:
            gene_id = fields[2]  # column 3 (0-based)
            transcript_id = fields[3]  # column 4 (0-based)
            orthology_type = fields[4]  # column 5
            if gene_id != "None":
                gene_to_transcripts[gene_id].append(transcript_id)
                gene_orthology_lines[gene_id].append((line.strip(), orthology_type))

# Save gene IDs for record-keeping
with open(f"{species}_geneIDs", "w") as f:
    for gene in sorted(gene_to_transcripts.keys()):
        f.write(f"{gene}\n")

# ---------------------------------------------
# Step 2: Parse the FASTA file into a dict
# ---------------------------------------------
fasta_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# ---------------------------------------------
# Step 3: For each gene, find the longest associated transcript
# ---------------------------------------------
for gene_id, transcripts in gene_to_transcripts.items():
    longest_seq = None
    max_length = 0

    for transcript_id in transcripts:
        record = fasta_records.get(transcript_id)
        if record:
            seq_len = len(record.seq)
            if seq_len > max_length:
                max_length = seq_len
                longest_seq = record

    # Save the longest transcript record
    if longest_seq:
        output_path = os.path.join(finaljs_dir, gene_id)
        with open(output_path, "w") as out_f:
            SeqIO.write(longest_seq, out_f, "fasta")

    # Also write all transcript IDs to js_dir for reference
    with open(os.path.join(js_dir, f"x_{gene_id}"), "w") as out_f:
        for tid in transcripts:
            out_f.write(f"{tid}\n")

# ---------------------------------------------
# Step 4: Concatenate all longest transcripts (sorted by gene name)
# ---------------------------------------------
output_faa_dir = "longest_transcript_faas"
os.makedirs(output_faa_dir, exist_ok=True)

combined_faa_path = os.path.join(output_faa_dir, f"{species}_longest_transcript.faa")
with open(combined_faa_path, "w") as combined_out:
    for gene_file in sorted(os.listdir(finaljs_dir)):
        gene_path = os.path.join(finaljs_dir, gene_file)
        for record in SeqIO.parse(gene_path, "fasta"):
            SeqIO.write(record, combined_out, "fasta")

# ---------------------------------------------
# Step 5: Create list_of_paralogs_${species}
# Only include genes where orthology type ≠ "1:1"
# ---------------------------------------------
paralog_list_path = f"list_of_paralogs_{species}"
with open(paralog_list_path, "w") as paralog_out:
    for gene_file in sorted(os.listdir(finaljs_dir)):
        gene = gene_file  # file is named by gene ID
        # Look up original orthology entries
        for line, ortho_type in gene_orthology_lines.get(gene, []):
            if ortho_type != "1:1":
                paralog_out.write(f"{line}\n")

# ---------------------------------------------
# Step 6: Verify no one2one orthologs in finalJs or combined fasta
# ---------------------------------------------
one2one_genes_in_finaljs = []

for gene_file in os.listdir(finaljs_dir):
    gene = gene_file
    orthology_entries = gene_orthology_lines.get(gene, [])
    # If any orthology types == "1:1", collect for warning
    if any(ortho_type == "1:1" for _, ortho_type in orthology_entries):
        one2one_genes_in_finaljs.append(gene)

if one2one_genes_in_finaljs:
    print("WARNING: The following genes classified as 1:1 orthologs were found in the finalJs folder:")
    for g in one2one_genes_in_finaljs:
        print(f" - {g}")
    # Optional: stop execution if this is critical
    # raise ValueError("1:1 orthologs found in finalJs folder — please check filtering.")
else:
    print("SUCCESS: No 1:1 orthologs found in finalJs folder. Check passed.")
