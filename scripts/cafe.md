---
title: Gene Family Contraction and Expansion with OrthoFinder and CAFE - step by step instructions

---

# Gene Family Contraction and Expansion with OrthoFinder and CAFE - step by step instructions 
By Tanya Lama & Nicole Paulat

## Objective 
We will use CAFE to analyze protein family evolution (gene family expansion/contraction). For our comparative analyses, annotations have been generated using TOGA1.

## Table of Contents
[TOC]

## Background/Overview

## Set up environment
### Installing Orthofinder and dependencies on Unity
Create a conda environment with orthofinder and all dependencies. 

## Set up directories and files
Nucleotide and protein FASTA files are located in the /completed_togas folder

## Set up `1_raw_protein_fastas` files
We want to retain just the query protein sequences from each protein.fasta file, and drop all of the reference (human) sequences. 

Create a folder with all of your protein_fasta files (output from toga) called `1_raw_protein_fastas` and create a list called `list_protein_fastas`
```bash!
mkdir -p cafe/1_raw_protein_fastas
cd completed togas
for dir in */; do if [[ -f "$dir/prot.fasta" ]]; then cp "$dir/prot.fasta" ../cafe/1_raw_protein_fastas/${dir%/}.prot.fasta; fi; done
#check that all the files were copied properly (not empty), then make list
cd ../cafe/1_raw_protein_fastas
ls ./*.prot.fa > list_protein_fastas
```

Run `grep_query.py` over all species
```bash!
for i in *.prot.fasta; do NAME="${i%.prot.fasta}"; python ../grep_query.py $NAME; done
```
#### `grep_query.py`
```python!
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
```

### Set up `2_query_protein_fastas`
Move all of the outputs called `query_${species}.pep` to a folder called `2_query_protein_fastas`

```bash!
for i in *.pep; do mv $i ../2_query_protein_fastas/; done
```

#### What's in orthology.tsv? 
orthology.tsv classifies orthology relationships: 
* `one2many` -- a gene in the reference species is orthologous to multiple genes in the query species
* `one2one` -- (pairwise orthology) a gene in the reference species has only one ortholog in the query species
* `many2many` -- orthologs in the reference and the query species have undergone lineage-specific duplications
* `many2one` -- a gene in the query species is orthologous to multiple genes in the reference species
* `one2zero` -- a gene is found in human but there is no ortholog in the query species

## Set up 3_orthology_tsvs 
We will now create a folder called `3_orthology_tsvs/` and prepare a list of one2one orthologs for each species.

```bash!
cd ../completed_togas/
for i in */; do cp ${i}/orthology_classification.tsv ../cafe/3_orthology_tsvs/${i%/}_orthology_classification.tsv; done
for dir in */; do echo "${dir%/}" >> ../cafe/3_orthology_tsvs/species_list.txt; done
cd ../cafe/3_orthology_tsvs
mkdir modified
```

### Pad `q_gene` names with zeroes
In order to select the longest transcript (step 5) we need to compare the length of q_transcripts (ENSTs) grouped by geneID (q_gene). However, grepping reg_10 might compare the length of ENSTs also to reg_100, reg_1000, and reg_10000. Therefore we padded the q_gene names with 0s to ensure that all q_gene names have 5 digits. 

```bash!
cd orthology_tsvs
# for each species, run python pad_names.py <genSpe> (e.g. antPal2)
for i in $(cat species_list.txt); do python pad_names.py ${i}; done
```

#### `pad_names.py`
```python!
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
```

#### Create a list of files in the /modified folder
```bash
ls modified/*.tsv > list_orthology_tsvs
```

### Create `one2ones_${species}.txt`
We need a list of one2one **`q_transcripts`** for each species so that we can drop them from the protein fastas for CAFE. Loop over the files in `list_orthology_tsvs`, selecting only the `q_transcripts` which have an orthology_class listed as 'one2one'. Drop any duplicate q_transcripts, and save the result as `one2ones_${species}.txt` 

```bash!
./get_one2ones.sh
```
#### `get_one2ones.sh`
```bash!
#!/bin/bash

# Loop over each TSV file in list_orthology_tsvs
while read -r file_path; do
    # Extract the characters 3-9 from the file name
    file_name=$(basename "$file_path")
    file_prefix=${file_name%_orthology_classification.tsv}

    # Extract the q_transcript values for rows where orthology_class == 'one2one'
    awk -F $'\t' '$5 == "one2one" { print $4 }' "$file_path" | sort -u > "one2one_${file_prefix}.txt"
done < list_orthology_tsvs
```

#### one2ones
Move these `one2one_{species}.txt` files a folder called `4_one2ones/`
```bash!
for i in $(cat species_list.txt); do mv one2one_${i}.txt ../4_one2ones/; done 
```

### Drop one2one orthologs from `${species}.pep` files
We need to do the following: 
For each one2one ortholog (e.g. `ENST00000000233.22`), we need to find the matching header line in the protein fasta and drop the matching line AND the line after it (the protein data for that one2one ortholog). Example:

>ENST00000000233.22 | PROT | QUERY
MKKPRAAVGSRPRKQTSSREEREKHGFNSSQAKPTIPDAGQARRQEEVVLQASVSPYHQFRDTAEVTVFRDSLLSWYDLKKRDLPWRRQAEGEVDLDRRAYAVWVSEVMLQQTQVATVINYYTRWMQKWPTLQDLAGASLEEVNQLWAGLGYYSRGRRLQEGARKVVEELGGHVPRTAETLQRLLPGVGRYTAGAIASIAFGQVTGVVDGNVVRVLCRVRGIGADPRNTFVSQQLWSLAQQLVDPARPGDFNQAAMELGAIVCTPQHPHCSQCPVQSLCRAHQRVERTQLPASQSLPGSPDVEECAPNTGQCQLCAPPTQPWDQTLGVANFPRKASRKPPREECSATCVLEQPRAPGGARVLLVQRPNSGLLAGLWEFPSVSVHPSGRHQRQALLQELRSWAGPLPAARLRHLGQVVHVFSHIKLTYQVYGLALEGQAPVTGTAPGARWLTREEFHTAAVSTAMKKVFRVYEGQQSGTCKSSKRCQVSTASRRKKPRPGQQVLDSFFRPRIPTDAPRPNSTA-*

==**NOTE: Formatting changed depending on the TOGA version, so make sure that the `prot.fasta` files are all the same format (run `sed 's/\.[^\.]*//' proteinAlignments.fasta > prot.fasta` if format is `>ENST00000000233.ARF5.22 | PROT | QUERY`)**==

To do this, loop parse_paralogs.py over species list
```bash!
cd 4_one2ones
for i in $(cat species_list.txt); do python ../parse_paralogs.py ${i}; done
```
#### `parse_paralogs.py`
```python!
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
```

The final product is a protein fasta without any one2one orthologs, `{species}_final.pep`, as well as a list file of one2one orthologs, `ensts_{species].txt`

### Filter for longest transcript
As a last step, we will retain only the longest transcript (one ENST for each `q_gene`). 

Loop over species list with `get_longest_transcripts.py`
For 24 species this will take about 20-30 minutes

```bash!
for i in $(cat species_list.txt); do echo ${i}; python ../get_longest_transcripts.py ${i}; done
```

#### `get_longest_transcripts.py`
```python!
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
```

This script produces multiple output files:
* `${species}_geneIDs` which is any gene in `${species}_final.pep` that is not one2none in `${species}_orthology.tsv`
* `${species}_js/x_{q_gene}` directory of output files of longest transcript id (ENST) per gene (ENSG)
* `${species}_finalJs/{q_gene}` directory of output protein .faa files of longest transcipt (ENST) for each ENSG
* `longest_transcripts_faas` directory of concatenated protein .faa files for each species

The files needed are in `longest_transcripts_faas/`; the rest are intermediate files for checking results if needed.

#### remove_species.sh
These js and finalJs folders are really big. You may want to remove them, but don't do it on the login node. Execute 3_remove_species.sh or 4_remove_js_finalJs.sh to remove a single species folders or ALL the js and finalJs folders. 

### Add Gene Symbols to `{species}_longest_transcripts.faa`
Since there are formatting differences between TOGA versions, some `prot.fasta` files didn't have ENSTXXX.GeneSym.XX, and those that did were converted to ENSTXXX.XX. But to make our lives easier, we want gene names for Orthofinder.

To re-add gene names to the headers, use `get_faas_gene_names.py`
This requires that you make a list file of all `{species}_longest_transcripts.faa` files and download a .tsv with transcript ID and Gene Name from Biomart via the archived [Ensembl version 104](https://may2021.archive.ensembl.org/biomart/martview/) from the Ensembl version 104 Human gene dataset.

```bash!
conda activate compleasm
mkdir -p labeled_input_faas
cd input_faas/
ls -1 > longest_transcripts_file_list.txt
python ../get_faas_gene_names.py ../ensembl_104_transcript_to_symbol.tsv longest_transcript_file_list.txt
mv *_annotated.faa ../labeled_input_faas
```
#### `get_faas_gene_names.py`
```python!
import sys
import os
from Bio import SeqIO

def load_mapping(tsv_path):
    mapping = {}
    with open(tsv_path) as f:
        header = next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            enst_id, gene_name = parts[0], parts[1]
            base_id = enst_id.split('.')[0]
            mapping[base_id] = gene_name
    return mapping

def strip_version(enst_id):
    return enst_id.split('.')[0]

def get_version(enst_id):
    return enst_id.split('.')[1] if '.' in enst_id else ""

def annotate_fasta(fasta_file, mapping):
    out_file = os.path.splitext(fasta_file)[0] + "_annotated.faa"
    print(f"Annotating {fasta_file} → {out_file}")
    with open(out_file, "w") as out_f:
        for rec in SeqIO.parse(fasta_file, "fasta"):
            base_id = strip_version(rec.id)
            version = get_version(rec.id)
            gene_symbol = mapping.get(base_id)
            if gene_symbol:
                rec.id = f"{base_id}.{gene_symbol}.{version}" if version else f"{base_id}.{gene_symbol}"
            else:
                rec.id = rec.id
            rec.description = ""
            rec.name = ""
            SeqIO.write(rec, out_f, "fasta")
    print(f"Done writing {out_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python get_faas_gene_names.py transcript_to_symbol.tsv list_of_faa_files.txt")
        sys.exit(1)

    mapping_file = sys.argv[1]
    list_file = sys.argv[2]

    mapping = load_mapping(mapping_file)
    print(f"Loaded {len(mapping)} transcript→gene mappings")

    with open(list_file) as f:
        fasta_files = [line.strip() for line in f if line.strip()]

    for fasta in fasta_files:
        if not os.path.isfile(fasta):
            print(f"Warning: file {fasta} not found, skipping")
            continue
        annotate_fasta(fasta, mapping)

if __name__ == "__main__":
    main()
```

## Run Orthofinder
To confirm gene family membership, run Orthofinder on all of the protein fasta genomes, modified to include only the longest transcript (ENST) for each ENSG (`{species}_longest_transcript.faa`). 
==NOTE: I did not remove gaps from sequences==
The input argument for Orthofinder (`-f`) is simply the directory that contains all of our protein files

Run `orthofinder.sh`:
```
sbatch orthofinder.sh
```
==**NOTE: Orthofinder will NOT run if the output directory already exists!!!**==

#### `orthofinder.sh`
```bash!
#!/bin/bash
#SBATCH --output=ortho.out
#SBATCH --error=mv-ortho.err
#SBATCH --job-name=ortho
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --ntasks=20
#SBATCH --partition=cpu

module load conda/latest
conda activate compleasm

WORK_DIR=/scratch/workspace/nmcdonald1_access-ci_org-shared/nikki_work/cafe/6_orthofinder

orthofinder -a 4 -f $WORK_DIR/input_faas -o $WORK_DIR/results
```

## Format the orthogroup data for CAFE
The next step is to take our `Orthogroup.tsv` data from OrthoFinder and do some filtering and formatting in preparation for CAFE. No need to enrich the orthogroups (OG) at this stage.

We want to calculate: 
1) The number of genes/species in each OG
2) The number of genes/OG

Some filtering rules: 
1) Drop OGs with > 450 genes
2) Drop OGs with > 100 size variation "Families which varied in size by more than 100 members between the species with the highest and lowest count were also removed."
3) Merge OGs > 80% similarity, 
4) Drop all olfactory gene OGs

#### `3a_clean_pairs.py`
Because we did an all x all comparison of the orthogroups to measure similarity, we have duplicate Orthogroups (e.g., OG0 was compared to OG1, and OG1 was also compared to OG0). To remove these reduncancies, we use `3_clean_pairs.py`, which keeps track of all of the orthogroup pairs and writes a file called `clean_pairs.txt` free of duplicates. 
```python!
# This script takes the grep_similarity.txt summary, and eliminated duplicated comparisons. For example, if OG1 and OG2 are 80% similar, we do not also need OG2 vs. OG1 similarity to know that the two need to be merged. 

seen_pairs = set()

with open('grep_similarity.txt', 'r') as file:
    lines = file.readlines()

clean_lines = []
for i in range(0, len(lines), 3):
    orthogroup1_line = lines[i].strip()
    orthogroup2_line = lines[i + 1].strip()
    similarity_line = lines[i + 2].strip()

    orthogroup1 = orthogroup1_line.split(":")[1].strip()
    orthogroup2 = orthogroup2_line.split(":")[1].strip()

    pair = (orthogroup1, orthogroup2)
    reversed_pair = (orthogroup2, orthogroup1)

    if pair not in seen_pairs and reversed_pair not in seen_pairs:
        clean_lines.append(f"{orthogroup1_line}\n")
        clean_lines.append(f"{orthogroup2_line}\n")
        clean_lines.append(f"{similarity_line}\n")

        seen_pairs.add(pair)

with open('clean_pairs.txt', 'w') as file:
    file.writelines(clean_lines)
```
#### Time to merge
After removing duplicates, we have about 800 orthogroups that need to be joined. 

#### `grep_cleaner_pairs.sh`
```bash
grep -E "Orthogroup 1:|Orthogroup 2:" clean_pairs.txt > cleaner_pairs.txt
```

#### `calculate_minimum_mergers.py`
Ok now this bit is quite clever, because it takes `cleaner_pairs.txt` and creates a network with nodes and edges to determine the minimum number of mergers needed to join all of the OGs with similarity > 80. 
```python!
#Reads the clean_pairs.txt file and identifies the minimum number of mergers needed to achieve the desired task. To simplify this, we can create a graph where each orthogroup is a node, and each pair of orthogroups represents an edge. Then, we can find the connected components in the graph, which will give us the groups of orthogroups that need to be merged together.
## Usage: python calculate_min_mergers.py > min_mergers.txt
def read_pairs(file_path):
    pairs = []
    with open(file_path, 'r') as pairs_file:
        lines = pairs_file.readlines()

    i = 0
    while i + 1 < len(lines):
        orthogroup1_line = lines[i].strip()
        orthogroup2_line = lines[i + 1].strip()

        orthogroup1 = orthogroup1_line.split(":")[1].strip()
        orthogroup2 = orthogroup2_line.split(":")[1].strip()

        pairs.append((orthogroup1, orthogroup2))
        i += 2

    return pairs


def create_graph(pairs):
    graph = {}
    for orthogroup1, orthogroup2 in pairs:
        if orthogroup1 not in graph:
            graph[orthogroup1] = set()
        if orthogroup2 not in graph:
            graph[orthogroup2] = set()

        graph[orthogroup1].add(orthogroup2)
        graph[orthogroup2].add(orthogroup1)

    return graph


def find_connected_components(graph):
    components = []
    visited = set()

    def dfs(node, component):
        visited.add(node)
        component.add(node)

        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, component)

    for node in graph:
        if node not in visited:
            component = set()
            dfs(node, component)
            components.append(component)

    return components


def main():
    cleaner_pairs_file = 'cleaner_pairs.txt'
    min_mergers_file = 'min_mergers.txt'

    # Read pairs from the clean_pairs.txt file
    pairs = read_pairs(clean_pairs_file)

    # Create a graph from the pairs
    graph = create_graph(pairs)

    # Find connected components in the graph
    connected_components = find_connected_components(graph)

    # Open output file for writing
    with open(min_mergers_file, "w") as out:
        # Write the minimum number of mergers needed
        out.write(f"Minimum number of mergers needed: {len(connected_components)}\n\n")

        # Write each connected component group
        for i, component in enumerate(connected_components):
            out.write(f"Group {i + 1}: {', '.join(component)}\n")


if __name__ == "__main__":
    main()
```

#### `perform_min_mergers.py`
Now we need to actually perform the orthogroup merges, and create an output file with those merged orthogroups and groups that were not part of any mergers (non_merged): `cleaned_merged_orthogroups.tsv`

```python!
import csv
import re

# === Config ===
orthogroups_file = 'Orthogroups.tsv'
min_mergers_file = 'min_mergers.txt'
merged_output_file = 'merged_output.tsv'
non_merged_output_file = 'non_merged_output.tsv'
combined_output_file = 'cleaned_merged_orthogroups.tsv'
summary_log_file = 'summary_mergers.txt'

# === Helper: natural sort key ===
def natural_key(text):
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', text)]

# === Read header ===
with open(orthogroups_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    header = reader.fieldnames

# === Parse mergers ===
merged_rows = {}
merged_orthogroups_set = set()
summary_log = []

group_counter = 1

with open(min_mergers_file, 'r') as f:
    for line in f:
        if line.startswith("Group"):
            _, group_data = line.strip().split(":")
            orthogroups = [og.strip() for og in group_data.split(",")]
            group_name = f"og_{group_counter}"
            group_counter += 1
            merged_row = {col: "" for col in header}
            merged_row[header[0]] = group_name

            for og in orthogroups:
                merged_orthogroups_set.add(og)
                with open(orthogroups_file, 'r') as cf:
                    reader = csv.DictReader(cf, delimiter='\t')
                    for row in reader:
                        if row[header[0]] == og:
                            for col in header[1:]:
                                values = row[col].split(',') if row[col] else []
                                existing = merged_row[col].split(',') if merged_row[col] else []
                                merged = existing + values
                                # Clean and deduplicate
                                cleaned = [v.strip() for v in merged if v.strip()]
                                deduped = sorted(set(cleaned), key=str)
                                merged_row[col] = ', '.join(deduped)
                            break
            merged_rows[group_name] = merged_row
            summary_log.append(f"{group_name}: {', '.join(orthogroups)}")

# === Parse original Orthogroups.tsv and find non-merged ===
non_merged_rows = {}
with open(orthogroups_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        og_name = row[header[0]]
        if og_name not in merged_orthogroups_set:
            # Clean each cell
            for col in header[1:]:
                values = row[col].split(',') if row[col] else []
                cleaned = [v.strip() for v in values if v.strip()]
                row[col] = ', '.join(sorted(set(cleaned), key=str))
            non_merged_rows[og_name] = row

# === Write merged_output.tsv ===
with open(merged_output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
    writer.writeheader()
    for row in sorted(merged_rows.values(), key=lambda x: natural_key(x[header[0]])):
        writer.writerow(row)

# === Write non_merged_output.tsv ===
with open(non_merged_output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
    writer.writeheader()
    for row in sorted(non_merged_rows.values(), key=lambda x: natural_key(x[header[0]])):
        writer.writerow(row)

# === Write cleaned_merged_orthogroups.tsv ===
with open(combined_output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
    writer.writeheader()
    combined = list(merged_rows.values()) + list(non_merged_rows.values())
    for row in sorted(combined, key=lambda x: natural_key(x[header[0]])):
        writer.writerow(row)

# === Write summary_mergers.txt ===
with open(summary_log_file, 'w') as f:
    f.write("Merged Groups Summary:\n")
    for entry in summary_log:
        f.write(f"{entry}\n")

print("All outputs generated:")
print(f" - {merged_output_file}")
print(f" - {non_merged_output_file}")
print(f" - {combined_output_file}")
print(f" - {summary_log_file}")
```

#### `6a_filter_orthogroups.py`
Need to remove any orthogroups that are:
* Are very large: >450 unique genes total
* Are highly variable across species: >100 gene variation between any species
```
python!
import csv

def process_gene_name(gene):
    parts = gene.split('.')
    if len(parts) > 2:
        # remove first and last parts, join the middle parts back with '.'
        return '.'.join(parts[1:-1])
    else:
        # otherwise just take the first part
        return parts[0]

def count_genes_per_species(row, species_columns):
    species_gene_counts = {}
    all_genes = set()
    for species in species_columns:
        raw_genes = row[species].split(', ') if row[species].strip() else []
        processed_genes = [process_gene_name(g) for g in raw_genes]
        species_gene_counts[species] = len(raw_genes)            # total genes count (no uniqueness here)
        all_genes.update(processed_genes)                         # unique genes across species (processed)
    return species_gene_counts, len(all_genes)

def check_filters(species_gene_counts, total_unique_genes, max_genes_threshold=450, max_variation_threshold=100):
    max_genes = max(species_gene_counts.values())
    min_genes = min(species_gene_counts.values())
    variation = max_genes - min_genes
    reasons = []
    if max_genes > max_genes_threshold:
        reasons.append(f"species with >{max_genes_threshold} genes (max: {max_genes})")
    if variation > max_variation_threshold:
        reasons.append(f"variation > {max_variation_threshold} genes (variation: {variation})")
    if total_unique_genes > max_genes_threshold:
        reasons.append(f"total unique genes > {max_genes_threshold} (total: {total_unique_genes})")
    return reasons

def filter_orthogroups(input_tsv, output_tsv, log_file):
    with open(input_tsv, newline='') as infile, \
         open(output_tsv, 'w', newline='') as outfile, \
         open(log_file, 'w') as log:
        
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        
        species_columns = reader.fieldnames[1:]  # Assuming first column is orthogroup name
        
        removed_count = 0
        kept_count = 0
        
        for row in reader:
            species_gene_counts, total_unique_genes = count_genes_per_species(row, species_columns)
            
            reasons = check_filters(species_gene_counts, total_unique_genes)
            
            if reasons:
                removed_count += 1
                log.write(f"Removed {row[reader.fieldnames[0]]}: {', '.join(reasons)}\n")
            else:
                kept_count += 1
                writer.writerow(row)
        
        log.write(f"\nSummary:\nTotal orthogroups processed: {removed_count + kept_count}\n")
        log.write(f"Orthogroups removed: {removed_count}\n")
        log.write(f"Orthogroups kept: {kept_count}\n")

if __name__ == "__main__":
    input_tsv = 'cleaned_merged_orthogroups.tsv'  # Your input file
    output_tsv = 'filtered_merged_orthogroups.tsv'       # Output filtered orthogroups
    log_file = 'orthogroup_filter_log.txt'        # Log file
    
    filter_orthogroups(input_tsv, output_tsv, log_file)
```

#### `remove_olfactory_orthogroups.py`

```python!
import csv
import re
from collections import Counter

# Combined regex pattern for olfactory receptor genes
or_pattern = re.compile(r'^OR\d+(?:[A-Z]+\d*|\d+)$')

def extract_gene_symbols(row, fieldnames):
    """Extract unique gene symbols from all species columns."""
    genes = []
    for field in fieldnames[1:]:
        entries = row[field].split(', ')
        for entry in entries:
            parts = entry.split('.')
            if len(parts) > 2:
            # remove first and last parts, join the middle parts back with '.'
                gene = '.'.join(parts[1:-1])
                if gene:
                    genes.append(gene)
            else:
            # otherwise just take the first part
                gene = parts[0]
                if gene:
                    genes.append(gene)
    return genes

def filter_olfactory_rows(input_file, output_file, log_file):
    with open(input_file, 'r') as infile, \
         open(output_file, 'w', newline='') as outfile, \
         open(log_file, 'w', newline='') as logfile:

        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=reader.fieldnames)
        log_writer = csv.writer(logfile, delimiter='\t')
        
        writer.writeheader()
        log_writer.writerow(['Orthogroup', 'Olfactory_Percent', 'Unique_Genes', 'OR_Genes'])

        for row in reader:
            gene_symbols = extract_gene_symbols(row, reader.fieldnames)
            unique_genes = set(gene_symbols)
            or_genes = [g for g in unique_genes if or_pattern.match(g)]
            olfactory_percent = len(or_genes) / len(unique_genes) if unique_genes else 0

            if olfactory_percent > 0.5:
                log_writer.writerow([
                    row[reader.fieldnames[0]],
                    f"{olfactory_percent:.2f}",
                    "; ".join(sorted(unique_genes)),
                    "; ".join(sorted(or_genes))
                ])
            else:
                writer.writerow(row)

    print(f"Finished filtering olfactory orthogroups.")
    print(f"Filtered output → {output_file}")
    print(f"Log of removed rows → {log_file}")

if __name__ == "__main__":
    filter_olfactory_rows(
        input_file='filtered_merged_orthogroups.tsv',
        output_file='no_olfactory_filtered_merged_orthogroups.tsv',
        log_file='removed_olfactory_log.tsv'
    )
```

#### `7_convert_into_cafe_format.py`
```
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
```

Finish reformatting for CAFE

```bash!
awk 'BEGIN{OFS="\t"} NR==1{$0="Desc\t"$0} NR>1{$0="(null)\t"$0}1' cafe_input_file_w_olfactory.tsv > cafe_input_w_or.tsv
sed -i  's/_longest_transcript_annotated//g' cafe_input_w_or.tsv
awk 'BEGIN{OFS="\t"} NR==1{$0="Desc\t"$0} NR>1{$0="(null)\t"$0}1' cafe_input_file_no_olfactory.tsv > cafe_input_no_or.tsv
sed -i  's/_longest_transcript_annotated//g' cafe_input_no_or.tsv
```

## Step Run CAFE 
### Install CAFE via conda
```bash!
conda create -n cafe bioconda::cafe
```

### Generate ultrametric tree (required for CAFE)

#### Check that your tree looks approximately correct in terms of branch lengths (My), and ensure that the topology is the same as the original supermatrix_iqtree.contree (minus hg38)**

### 3. Estimate error rate with CAFE
```bash!
module load conda/latest
conda activate cafe

# Run CAFE
cafe5 -i cafe_input_w_or.tsv -t ultrametric_tree_clipped.nwk -e
# Copy output to work dir
cp results/Base_error_model.txt ../
```

### 4. Run CAFE analysis
CAFE recommends test multiple lambdas

#### See [CAFE Github](https://github.com/hahnlab/CAFE5) and [tutorial](github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md) for more information

#### `sbatch_cafe1.sh`
```bash!
#!/bin/bash

#SBATCH --job-name=cafe1
#SBATCH --output=cafe1.out
#SBATCH --error=cafe1.out
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH -p cpu

module load conda/latest
conda activate cafe

# Run CAFE with cafe_input_w_or.tsv
echo "Running CAFE with OR and lambda1"
cafe5 -i cafe_input_w_or.tsv -t ultrametric_tree_clipped.nwk -eBase_error_model.txt -o cafe1_run1 -y tree_lambda1.txt
echo "Running CAFE with OR and lambda2"
cafe5 -i cafe_input_w_or.tsv -t ultrametric_tree_clipped.nwk -eBase_error_model.txt -o cafe1_run2 -y tree_lambda2.txt
echo "Running CAFE with OR and lambda3"
cafe5 -i cafe_input_w_or.tsv -t ultrametric_tree_clipped.nwk -eBase_error_model.txt -o cafe1_run3 -y tree_lambda3.txt
```

### 4. Determine best model (best # lambdas)
Look at `Base_results.txt` for each run and determine if there is a substantial improvement in Model Base Final Likelihood (-lnL)

How to interpret:
See discussion posts: [one](https://github.com/hahnlab/CAFE5/discussions/164), [two](https://github.com/hahnlab/CAFE5/discussions/155)

### 4. Visualizing results with CAFEplotter and CAFE_fig
CafePlotter is very basic, it only produces trees with combined expansion/contraction results plotted as red/blue numbers with significance as ** or ***. 

### 1. Make CafePlotter environment

```
module load conda/latest
conda create -n cafeplotter pip
conda activate cafeplotter
pip install cafeplotter
```

### 2. Run CafePlotter
* `-i`: input directory is the cafe run folder
* `-o`: output directory for the plots
* `--ignore_branch_length`: flag to ignore branch length for plotting (default OFF)
```
cd 8_cafe
mkdir cafe1_run1_plots
cafeplotter -i cafe1_run1 -o cafe1_run1_plots/ --ignore_branch_length
```
