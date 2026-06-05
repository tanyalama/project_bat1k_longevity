---
title: HyPhy (aBSREL, MEME, RELAX)

---

# HyPhy (aBSREL, MEME, RELAX)

By Nikki Paulat & Tanya Lama & Bill Thomas & Graham Hughes & many more 
With #blessings from Sergei Kostakovsky Pond

### Objective: 
Run an exploratory round of selection tests on the full 23 bat & 7 nonbat species alignment for the bat1k longevity 

aBSREL (in the HyPhy package) performs branch-site model-based tests for episodic selection affecting a proportion of sites in the alignment along a proportion of branches in the tree. In other words, **is there evidence of selection *anywhere* in the alignment and *anywhere* on the tree**? 

### Required Software
* FASTconCAT
* Conda: iqtree, hyphy ([Github](https://github.com/veg/hyphy))
* R (packages: ape, phangorn)

## Table of Contents

[TOC]

## Generate Tree
==This pipeline assumes that you have already generated your gene alignment files, in this case from TOGA files run through the Hiller exon-by-exon pipeline==

### Copy and reformat alignment files
Make sure the alignment files have the extension .fasta or .fas
```bash!
cd /myworkingdirectory/masked_alignments/taper
for file in *.fa; do
  mv -- "$file" "${file%.fa}.fas"
done
```
Check for empty alignment files 
```
cp *.fas ../clean_fastas/
cd ../clean_fastas
find . -type f -empty
```
We also dropped alignments which include less than 5 species (see `alignments_insufficient_sp.txt`) (`rm <file>`)

Need to clean up the sequence headers for FASconCat to use -- only species abbreviations in the headers. Also need to change "REFERENCE" to "hg38".
```bash!
for file in *.fas; do sed -i '/^>/ s/[[:space:]].*//' "$file"; done
for file in *.fas; do sed -i 's/REFERENCE/hg38/g' "$file"; done
```

FASconCAT doesn't recognize "X" bases in nucleotide seqs, so need to change to Ns
```bash!
for file in *.fas; do sed -i '/^>/! s/X/N/g' "$file"; done
for file in *.fas; do sed -i '/^>/! s/x/n/g' "$file"; done
```

MACSE alignments denote frameshift mutations with "!" (see [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022594)), which aren't recognized by anything, and idk what to do with them, so swap with n

```bash!
for file in *.fas; do sed -i '/^>/! s/!/n/g' "$file"; done
```

### Concatenate gene alignments with FASconCat-G

Usage is: `perl FASconCAT-G_v1.05.pl [command options]` (Linux/Mac)
```!
Where command options:
  -p PHYLIP output (strict taxon names) #lets add this. We need it for iqtree
  -n NEXUS output #lets add this just in case we need it for iqtree  -- apparently this is the partition file. We need this too
  -l partitioned 
  -s START #include
```

These are the flags we have selected. You should end up with four files: `FcC_info.xls`, `FcC_supermatrix.fas`, `FcC_supermatrix.nex`,  `FcC_supermatrix.phy`
```
#perl FASconCAT-G_v1.05.pl -p -n -l -s
perl FASconCAT-G_v1.06.1.pl -n -p -l -s
```

Run fasconcat.sh, which needs to sbatched from within the your alignments folder. 
Will run in a couple of hours.

#### `fasconcat.sh`
```bash!
#! /bin/bash

#SBATCH --job-name=fasconcat
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --mem=200G #probably don't need this much, idk maybe 80G?
#SBATCH --time=04:00:00 #took 1.5 hours with 15.5k genes
#SBATCH -p cpu

#perl ./FASconCAT-G_v1.05.pl -n -p -l -s 
perl ./FASconCAT-G_v1.06.1.pl -n -p -l -s
```

### Review the outputs
The outputs we are most interested in are `supermatrix_partition.nex` and `supermatrix_partition.txt`. 
You can check that fasconcat ran to completion by seeing if # lines in supermatrix_partition.txt matches the number of input alignment files.
```
wc -l FcC_supermatrix_partition.txt
```

### Run iqtree: select best model of gene evolution and estimate phylogeny
We will create a conda environment to house iqtree and our other dependencies for this step. 
```bash!
module load conda/latest
conda create --name iqtree
#Activate the environment
conda activate iqtree
#Install any necessary modules
conda install -c bioconda iqtree
```

#### `slurm_iqtree.sh`
```bash!
#! /bin/bash
  
#SBATCH --job-name=iqtree-partitioned
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --mem=100G #15.5k genes needed min 47G
#SBATCH --time=24:00:00 #ran for 7 hours with 15.5k genes
#SBATCH -p cpu

iqtree -s *_supermatrix_partition.nex -p *_supermatrix_partition.txt --prefix supermatrix_partition_iqtree -B 1000 -T AUTO
```

### Review the outputs
The main output of interest is `supermatrix_partition_iqtree.contree`
But first let's check it out using figtree

#### Visualize your tree in figtree
Download compiled binaries for FigTree [here](https://github.com/rambaut/figtree/releases) and install accordingly. You will also need to install an updated version of java. FigTree is a desktop application. 
Download `supermatrix_partition_iqtree.contree` from the cluster and rename the suffix as **`supermatrix_partition_iqtree.tree`**
Open with FigTree
FigTree will notice immediately that you have an "added" field -- these are our bootstraps. Name 'bootstraps' accordingly
On the left, select Node Labels and view the aforementioned bootstraps on our tree: 
Save your rerooted tree as `supermatrix_partition_iqtree_reroot.tree`

### Drop bootstraps from `supermatrix_partition_iqtree_reroot.tree`

Save this version of your tree as **`supermatrix_partition_iqtree_reroot_no_bootstraps.tree`**

## Part 2: Create batches
Move all *.masked.fas files to a folder called fastas

Split all_fastas into batches of 1500 (max job count = 2000)
```bash!
sed -n -e '1,1500p' all_fastas.txt >batch1.txt
sed -n -e '1501,3000p' all_fastas.txt >batch2.txt
sed -n -e '3001,4500p' all_fastas.txt >batch3.txt
sed -n -e '4501,6000p' all_fastas.txt >batch4.txt
sed -n -e '6001,7500p' all_fastas.txt >batch5.txt
sed -n -e '7501,9000p' all_fastas.txt >batch6.txt
sed -n -e '9001,10500p' all_fastas.txt >batch7.txt
sed -n -e '10501,12000p' all_fastas.txt >batch8.txt
sed -n -e '12001,13500p' all_fastas.txt >batch9.txt
sed -n -e '13501,15000p' all_fastas.txt >batch10.txt
sed -n -e '15001,15484p' all_fastas.txt >batch11.txt
#last one is remainder, has 484 genes
```

## Part 3: Clip trees
Not all 15.5k genes are present for each species (e.g., a gene may be missing for myoMyo or artJam). So, for each gene we need to trim the tree and create an input for HyPhy. HyPhy aBSREL will then be run individually for each gene. 
### We have a folder for each batch of fastas. Copy your tree, `modclip.py`, and `cliptrees.R`

### Execute `modclip.py` and `cliptrees.R` from within the folder
Execute `modclip.py`
```bash!
module load conda/latest r-rocker-ml-verse/4.4.0+apptainer
WORK_DIR=/myworkingdirectory/hyphy
for i in {1..11}; do cd $WORK_DIR/batch${i}; python ./modclip.py; Rscript cliptrees.R; cd $WORK_DIR; done
```

#### `modclip.py` (Python conversion of Graham's `modclip.pl`)
```python!
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
```

#### `cliptrees.R`
```R!
tab<- read.table("details.txt")
       library("phangorn")
       library("ape")
       i<- 1
    while(i<nrow(tab)+1){
                          main<-read.tree("supermatrix_partition_iqtree_reroot_no_bootstraps.tree")
    a<- strsplit(as.matrix(tab[i,]),"\\|")
        genes<-a[[1]][1]
    tfl<-"tfl"
          j<-2
    while(j<length(a[[1]])+1){
                               main<-drop.tip(main,a[[1]][j])
    j<-j+1
           }
           rfoutput=paste(genes,tfl,sep=".")
           write.tree(file=rfoutput,main)
           i<-i+1
    }
```

#### `modclip.py` should have created a file called `details.txt` which looks like this:
```!
ENST00000323110.CYP11B2.cleanLb_hmm_manual.fasta|canFam4|bosTau9|equCab3|HLtadBra2 #missing four species
```

#### `cliptrees.R` should have created a .tfl file (reduced Newick tree) for each fasta!

### Add foreground tags to tree branches
Note: if your hyphy job doesn't include the `--branches <tag>` argument, tree tags are ignored, at least for aBSREL.

Here, the tag is `{FG}`, but it can be whatever makes sense.

```bash!
for i in {1..11}; do cd batch${i}; for tree in *.fas.tfl; do sed -e 's/myspecies/myspecies{FG}/' -e 's/myspecies/myspecies{FG}/' ${tree} > "${tree%.tfl}".labeled.tfl; done; cd ../; done
```
If you accidentally run this twice on a batch, here's the fix:
```bash
sed -E -i 's/(\{FG\})\1+/\1/g' *.labeled.tfl
```
Alternatively, HyPhy has ways to do this as well:
* Programatically: https://github.com/veg/hyphy-analyses/tree/master/LabelTrees
* Manually (for small trees or if you only have one tree for all genes): https://phylotree.hyphy.org/


### Create a gene list for each batch within the /batch folders
This is not necessary unless you want a redundant file to check in the folder.

```bash!
for i in {1..11}; do cp batch${i}.txt batch${i}/batch${i}_genelist.txt; done
#OR
#for i in {1..8}; do cd batch${i}; ls ENST*.fas > batch${i}_genelist.txt; done
```
## Part 4: Selection Testing with HyPhy ABSREL
HyPhy ([Github](https://github.com/veg/hyphy), [tutorials](https://www.hyphy.org/tutorials/CL-prompt-tutorial/)) has various analyses you can run, here's the breakdown of what we did:

1. aBSREL exploratory mode
2. MEME for genes under selection from aBSREL
3. RELAX 

All analyses with `--srv Yes` (see [paper](https://doi.org/10.1093/molbev/msad150) and Github [info](https://github.com/veg/hyphy-analyses/tree/master/SRV), [convo](https://github.com/veg/hyphy/issues/1647) for why)

Info about multiple test corrections in HyPhy: [convo1](https://github.com/veg/hyphy/issues/188), [convo2](https://github.com/veg/hyphy/issues/1707) (aBSREL does within gene error correction), [convo3](https://github.com/veg/hyphy/issues/1062), 

### Create a conda environment for hyphy

### Run HyPhy aBSREL jobs as array batches
:::info
**NOTE: aBSREL accepts rooted or unrooted trees, as it unroots in internally**
:::

:::spoiler non-working auto_sub_batches
After making the array script for each batch, generate a list of the submission scripts, and run the batching script `slurm_batch_hyphy.sh`

```bash!
ls slurm_absrel_hyp* > absrel_hyp_batch_list.txt
sbatch slurm_batch_hyphy.sh absrel_hyp_batch_list.txt
```
#### `slurm_batch_hyphy.sh`
```bash!
#!/bin/bash
#SBATCH --job-name=batch_hyphy
#SBATCH --output=batch_hyphy_%j.out
#SBATCH --error=batch_hyphy_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -p cpu

# === Arguments ===
# $1 - INPUT_LIST: Path to the input file listing all input files, one per line.

# === Configuration ===
MAX_REMAINING=500
JOB_LIST="$1"
PREV_JOB_ID=""

if [[ -z "$JOB_LIST" || ! -f "$JOB_LIST" ]]; then
    echo "Usage: $0 <job_list_file>"
    exit 1
fi

while IFS= read -r JOB_SCRIPT; do
    if [[ ! -f "$JOB_SCRIPT" ]]; then
        echo "Skipping missing job script: $JOB_SCRIPT"
        continue
    fi

    # Wait for previous job's array to drop below threshold
    if [[ -n "$PREV_JOB_ID" ]]; then
        echo "Waiting for fewer than $MAX_REMAINING total remaining tasks from job $PREV_JOB_ID..."

        # Get total number of array tasks from original job submission
        TOTAL_TASKS=$(sacct -j "$PREV_JOB_ID" --format=JobID,JobName%20,State --noheader -X | grep -E "^${PREV_JOB_ID}_[0-9]+" | wc -l)

        if [[ "$TOTAL_TASKS" -eq 0 ]]; then
            # fallback: try to parse from scontrol
            TOTAL_TASKS=$(scontrol show job "$PREV_JOB_ID" | grep -oP 'NumTasks=\K[0-9]+')
        fi

        [[ -z "$TOTAL_TASKS" || "$TOTAL_TASKS" -eq 0 ]] && TOTAL_TASKS=99999  # fail-safe

        while :; do
            # Count completed tasks
            COMPLETED_TASKS=$(sacct -j "$PREV_JOB_ID" --format=JobID,State --noheader -X \
                | grep -E "^${PREV_JOB_ID}_[0-9]+" | grep -E 'COMPLETED|CANCELLED|FAILED|TIMEOUT|BOOT_FAIL|PREEMPTED' | wc -l)

            REMAINING=$((TOTAL_TASKS - COMPLETED_TASKS))
            [[ $REMAINING -lt 0 ]] && REMAINING=0

            echo "  Remaining array tasks: $REMAINING"

            if (( REMAINING < MAX_REMAINING )); then
                echo "Below threshold; submitting next job."
                break
            else
                echo "Waiting 3 minutes..."
                sleep 180
            fi
        done
    fi

    # Submit the job script (which already includes #SBATCH --array)
    JOB_OUTPUT=$(sbatch "$JOB_SCRIPT")
    echo "Submitted: $JOB_SCRIPT → $JOB_OUTPUT"

    # Extract new job ID
    PREV_JOB_ID=$(echo "$JOB_OUTPUT" | awk '{print $4}')
    sleep 10
done < "$JOB_LIST"

```
:::

Make the template script and then replace `<NUM>` with batch numbers to generate scripts for each batch. Remember to manually modify the last batch's array length, since it will be the remainder and not a full batch.

```bash!
for i in {1..11}; do sed "s/<NUM>/$i/" template_slurm_absrel_array.sh > slurm_absrel_array${i}.sh; done
```

Submit each array job individually when you have <500 jobs already in queue.
```
sbatch slurm_absrel_array1.sh
```

#### `template_slurm_absrel_array.sh`
```bash!
#!/bin/bash
#SBATCH --job-name=absrel_array
#SBATCH --output=absrel_logs/absrel_%A_%a.out   # Log file per task (create 'logs' dir first)
#SBATCH --error=absrel_logs/absrel_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --array=0-1500         # Adjust based on number of lines in batch${i}.txt

# Load any necessary modules
module load conda/latest gcc/9.4.0
conda activate hyphy

# Read the input file corresponding to the SLURM array task ID
BATCH="batch<NUM>"
WORK_DIR=/myworkingdirectory/hyphy
LIST_FILE=$WORK_DIR/${BATCH}.txt
INPUT_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $LIST_FILE)

cd $WORK_DIR/$BATCH
#mkdir -p results/
echo "Processing file: $INPUT_FILE"

NAME="${INPUT_FILE%.exon*}"

hyphy absrel CPU=1 --pvalue 0.05 --alignment $INPUT_FILE --tree $INPUT_FILE.labeled.tfl --branches FG --srv Yes >> $WORK_DIR/absrel_hyp/$BATCH/${NAME}_absrel_output.txt
```

### Restarting disrupted hyphy runs
A batch was accidentally canceled. I want to restart it, without rerunning hyphy on all the genes for which there are completed results 

```bash!
ls results/* > completed
cat completed | grep -Po '(?<=(results/)).*(?=_absrel)' >> alldone
tail -n +37 batch1_genelist.txt > batch1_genelist_v2.txt #drop the first 40 lines of the batch1 genelist
```

If the run failed because of either alignment premature stop codons ("ASSERTION ERROR", need to fix alignment then rerun) or "internal errors" (???), grab the relevant files and make a redo folder

```bash!
mkdir redo_batch

suffix=".exon.aligned_masked.fas"  # Define your suffix
#suffix=".exon.aligned_masked.fas.labeled1.tfl"  # Define your suffix

while read -r partial_name; do
  target_file="${partial_name}${suffix}"
  for dir in batch{1..11}; do
    file=$(find "$dir" -type f -name "$target_file" -print -quit)
    if [[ -n "$file" ]]; then
      cp "$file" relax_redo/
      echo "Copied: $file"
      break  # Stop searching once found
    fi
  done
done < relax_redo_genes1.txt
```

### Get list of alignments with positive results
```bash!
# Get list of aBSREL result files with positive results
cd absrel_hyp/
#find batch{1..8} -name '*.txt' -exec grep -El 'found \*\*[1-9][0-9]*\*\*' {} + | xargs -n1 basename > absrel_positive_results.txt
find batch{1..11} -name '*.txt' -exec grep -El 'found \*\*[1-9][0-9]*\*\*' {} + | xargs -n1 basename > absrel_positive_results.txt

# Get list of alignment files that match up to those result files
WORK_DIR=/myworkingdirectory
> matched_files.txt
while IFS= read -r filename; do key=$(echo "$filename" | cut -d'_' -f1); match=$(find $WORK_DIR/masked_alignments/clean_fastas -type f -name "${key}.*.fas"); if [ -n "$match" ]; then echo "$match" >> matched_files.txt; fi; done < absrel_positive_results.txt
```

### Make a folder meme/fastas for genes with positive results
```bash!
WORK_DIR=/myworkingdirectory/hyphy
mkdir -p $WORK_DIR/meme
cd meme
cp $WORK_DIR/absrel_hyp/matched_files.txt .
mkdir -p fastas
# Copy over the gene alignment.fas files
for file in $(cat matched_files.txt); do cp $file fastas/; done
cd fastas/
# Make new list of alignments that are in the fastas/ folder
ls *.fas > ../alignment_list.txt
cd ..
# Copy over the relevant .tfl files for the .fas files in fastas/
for filename in $(cat alignment_list.txt); do for dir in batch{1..11}; do filepath=$(find $WORK_DIR/$dir -type f -name ${filename}.tfl -print -quit); if [ -n "$filepath" ]; then cp "$filepath" fastas/; break; fi; done; done
```

## Run HyPhy MEME for genes with positive aBSREL results 
```
sbatch slurm_meme_array.sh
```
#### `slurm_meme_array.sh`
```bash!
#!/bin/bash
#SBATCH --job-name=meme_array
#SBATCH --output=logs/meme_%A_%a.out   # Log file per task (create 'logs' dir first)
#SBATCH --error=logs/meme_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --array=0-1166         # Adjust based on number of lines in alignment_list.txt

# Load any necessary modules
module load conda/latest gcc/9.4.0
conda activate hyphy

# Read the input file corresponding to the SLURM array task ID
WORK_DIR=/myworkingdirectory/hyphy/meme
LIST_FILE=$WORK_DIR/alignment_list.txt
INPUT_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $LIST_FILE)

cd $WORK_DIR/fastas
#mkdir -p results/
echo "Processing file: $INPUT_FILE"

NAME="${INPUT_FILE%.exon*}"

hyphy meme CPU=1 --pvalue 0.05 --alignment $INPUT_FILE --tree $INPUT_FILE.tfl --srv Yes >> $WORK_DIR/results/${NAME}_MEME_output.txt

## --output /work/###_MEME_output.json
```

## Run RELAX on all genes for two hypotheses
#### **Test intensification or relaxation of selection with HyPhy RELAX**
### Add labels to clipped trees
This will be for two separate test runs, so you need to make 2 different labeled trees per gene.

### Set test and reference branches with {FG}

### Submit batches of RELAX jobs (1-11) for your hypothesis
```
sbatch slurm_relax_array1.sh
```
#### `slurm_relax_array1.sh`
```bash!
#!/bin/bash
#SBATCH --job-name=relax_array
#SBATCH --output=logs/relax1_%A_%a.out   # Log file per task (create 'logs' dir first)
#SBATCH --error=logs/relax1_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --array=0-1500         # Adjust based on number of lines in alignment_list.txt

# Load any necessary modules
module load conda/latest gcc/9.4.0
conda activate hyphy

# Read the input file corresponding to the SLURM array task ID
BATCH="batch1"
WORK_DIR=/myworkingdirectory/hyphy
LIST_FILE=$WORK_DIR/${BATCH}.txt
INPUT_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $LIST_FILE)

mkdir -p $WORK_DIR/relax/$BATCH

cd $WORK_DIR/$BATCH
echo "Processing file: $INPUT_FILE"

NAME="${INPUT_FILE%.exon*}"

hyphy relax CPU=1 --pvalue 0.05 --alignment $INPUT_FILE --tree $INPUT_FILE.labeled1.tfl --srv Yes --test FG1 --reference FG2 >> $WORK_DIR/relax/$BATCH/${NAME}_RELAX_output.txt

## --output /work/###_RELAX_output.json
```

### Redo any failed runs
If the run failed because of either alignment premature stop codons (need to fix alignment then rerun) or "internal errors" (???), grab the relevant files and make a redo folder

## Part 7: Parse RELAX results

### Pull the ENST.geneId and pvalue for each gene. 

### BH adjust pvalues 
Execute `padj.py` relax_results_parsed.txt`
You will need a conda env that has the `statsmodels` module installed!
```
module load conda/latest
conda activate stats_env
python padj_antPal2.py
```

#### `padj.py`
```python!
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
```

See `bh_adjusted_pvalues.txt` for output 

### Parse selection intensity and K parameter value
```
bash ./add_selection_antPal2.sh
```
#### `add_selection_antPal2.sh`
```bash!
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
```

See `relax_results_bh_adjusted_with_selection.txt` for output.
You should now have an tab-delimited file with geneID, geneSymbol, original pvalue, BH adjusted pvalue, selection intensity for p<0.05 genes, and K value (branch-dependent selection intensity parameter). 

### Create a dataframe with three columns for each results set: 
GeneID     pFDR    K
ABCD       0.01    2.10

```bash!
# Extract desired columns and rename header
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "GeneID","pFDR","K"} NR>1{print $2,$4,$6}' relax_results_bh_adjusted_with_selection.txt > relax_selection_summary.txt

# Extract desired columns and rename header
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "GeneID","pFDR","K"} NR>1{print $2,$4,$6}' relax_results_bh_adjusted_with_selection.txt > relax_selection_summary.txt
```

### Plot a selection intensity volcano in R  
```bash!
module load r-rocker-ml-verse/4.4.0+apptainer
Rscript relax_volcano_plot_antPal2.R
```

#### `relax_volcano_plot.R`
```R!
## Load libraries
library(ggplot2)
library(dplyr)
library(viridis)      # for viridis color palette
library(scales)        # for log10 axis formatting if needed

## Load data
data <- read.table("relax_selection_summary.txt", header = TRUE, sep = "\t")
#data <- read.table("relax_selection_summary.txt", header = TRUE, sep = "\t")

## Log2 transform K, but only if K > 0
data <- data %>%
  filter(K > 0) %>%  # Avoid taking log of zero or negative
  mutate(
    log2K = log2(K),
    Selection = case_when(
      pFDR < 0.05 & K > 1 ~ "Intensified",
      pFDR < 0.05 & K < 1 ~ "Relaxed",
      TRUE ~ "Not significant"
    ),
    Selection = factor(Selection, levels = c("Relaxed", "Not significant", "Intensified"))
  )

## Pick your genes of interest
intensified_genes <- c("FOXP3", "SLX4", "XPA", "BCCIP", "WDR45", "ARID1B", "TRIM65", "SMCHD1")
relaxed_genes <- c("PIK3CD", "CLSPN", "UMOD", "ERCC1", "EME1", "COA8", "HUWE1", "IL1RL1", "TNIP3", "ZMYM6", "KAT6A")

## Create your LabelGroup
data$LabelGroup <- case_when(
  data$GeneID %in% intensified_genes ~ "Intensified",
  data$GeneID %in% relaxed_genes ~ "Relaxed",
  TRUE ~ NA_character_
)

data$Label <- ifelse(!is.na(data$LabelGroup), data$GeneID, NA)
  
## Pick your color palette 
show_col(viridis_pal()(20))
show(viridis_pal()(20))


## Final figure with labels
ggplot(data, aes(x = log2K, y = pFDR)) +
    # Gray hollow background points (non-highlighted)
    geom_point(
        data = subset(data, is.na(LabelGroup)), 
        shape = 21, fill = "gray70", color = "gray70", 
        size = 2, stroke = 0.5, alpha = 0.3
    ) +
    
    # Colored solid points (highlighted genes only)
    geom_point(
        data = subset(data, !is.na(LabelGroup)), 
        aes(color = LabelGroup),
        size = 3
    ) +
    
    # Black labels for highlighted genes
    geom_text_repel(
        data = subset(data, !is.na(LabelGroup)),
        aes(label = Label),
        size = 3,
        max.overlaps = Inf,
        segment.color = "black",
        color = "black"
    ) +
    
    # Custom viridis-inspired colors
    scale_color_manual(values = c(
        "Relaxed" = "#32648EFF",        # viridis blue
        "Intensified" = "#B8DE29FF"     # viridis bright green
    )) +
    
    # Reference lines
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    
    # Layout and theme
    theme_minimal(base_size = 14) +
    coord_cartesian(ylim = c(1, -0.10)) +
    scale_y_reverse(limits = c(1, -0.10)) 
  ```
  
We also included the entrez summary in Supplementary Table 12 for each gene, and a note regarding whether each gene is associated with aging hallmarks or a priori lists (e.g., antagonistic pleiotropies of aging).

You can now move foward with investigations of invidividual genes

## Part 9: alphaFold & alphaMissense
### 1. Visualization with pymissense
Already installed in /bin
Usage is: 
```
python missense.py UNIPROTid output --maxacid --tsv
```

Where `--maxacid` is the number of aa you want displayed and `--tsv` is the alphaMissense predictions. You can use your own or the pre-computed predictions for all possible human aa substitutions and missense variants ([hosted here](https://console.cloud.google.com/storage/browser/dm_alphamissense;tab=objects?inv=1&invt=Ab2Nkg&prefix=&forceOnObjectsSortingFiltering=false)).
Example:
`python /work/pi_tlama_smith_edu/bin/pymissense/missense/missense.py Q9Y3I1 out --maxacid 409 --tsv Q9Y3I1_sites.tsv` 

### alphaMissense conda environment
At a minimum you will need:
```
conda install conda-forge::requests
conda install conda-forge::matplotlib-base
```

Honestly if you just run this command: `python /work/pi_tlama_smith_edu/bin/pymissense/missense/missense.py Q9Y3I1 out --maxacid 409 --tsv Q9Y3I1_sites.tsv` 

## Appendix and FAQ
