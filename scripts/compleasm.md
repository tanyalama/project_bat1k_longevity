---
title: 'hyphy meme'
disqus: hackmd
---

Calculating BUSCO scores with compleasm
===
By Blair Bentley

### Objective: 
Generating "BUSCO" annotation completeness scores using the compleasm software (https://github.com/huangnengCSU/compleasm) for the bat1k longevity paper. Note that the below uses the longest transcripts for each protein, and that these scores reflect annotation completeness rather than genome completeness. Therefore protein fastas were the input file rather than nucleotide fastas for these analyses. 

### Acknowledgements
Thanks to Bill Thomas for suggesting we try compleasm!

## Table of Contents

[TOC]

### Prepare inputs
Formatting the protein file for use with OrthoFinder primary_transcript.py script:
```
# Copy the protein fastas from Tanya's directory
cp /project/tlama_umass_edu/projects/project_bat1k_longevity/data/protein_fastas/ ./refs

# Make a file with all species prefixes:
ls ./refs/* > all_spp.txt
sed -i -e 's/.prot.fasta/g/' all_spp.txt

# Extract the query sequence for each species:
for q in {1..30}; do SPP=$(sed -n ${q}p ../../all_spp.txt); echo $SPP; grep -A 1 "QUERY" ${SPP}.prot.fasta > ${SPP}.query.prot.fasta; done

# Remove '-' and 'X' from files
for q in {1..30}; do SPP=$(sed -n ${q}p ../../../all_spp.txt); echo $SPP; sed -i -e 's/X//g' ${SPP}.query.prot.fasta
for q in {1..30}; do SPP=$(sed -n ${q}p ../../../all_spp.txt); echo $SPP; sed -i -e 's/-//g' ${SPP}.query.prot.fasta

# Remove blank lines that cause formatting issues:
for q in {1..30}; do SPP=$(sed -n ${q}p ../../../all_spp.txt); echo $SPP; sed -i -e '/^$/d' ${SPP}.query.prot.fasta; done

# Reformat the headers to include 'gene=' for the OrthoFinder script:
for q in {1..30}; do SPP=$(sed -n ${q}p ../../../all_spp.txt); echo $SPP; sed -E 's/>(ENST[0-9]+)\.([A-Z0-9_]+)\.([0-9]+)/>\1.\3; gene=\2/' ${SPP}.query.prot.fasta > ${SPP}.reform.query.prot.fasta; done
```
### Extracting the longest transcript per protein using the OrthoFinder tool:
```for q in {1..30}; do SPP=$(sed -n ${q}p ../../../all_spp.txt); echo $SPP; python ~/OrthoFinder/tools/primary_transcript.py ${SPP}.reform.query.prot.fasta; done
```
### Running Compleasm on the longest protein transcripts for each species:
```# Ensure the below script is run with the compleasm conda environment activated:
# Also download mammal database first (i.e. compleasm download mammalia)

#!/bin/bash
#SBATCH -J compleasm_protein_longest_allSpp
#SBATCH -o ./logs/all_spp/%x_%a_%A.log
#SBATCH -e ./logs/all_spp/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 02:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=64G  # Requested Memory
#SBATCH --array=1-30

SPP=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./all_spp.txt)
echo $SPP

compleasm protein -p /scratch3/workspace/bbentley_smith_edu-busco/refs/all_spp/cleaned/primary_transcripts/${SPP}.reform.query.prot.fasta -l mammalia -o compleasm/all_spp/${SPP} -t 16
```
### Download the log files to Drive:
```
scp -i ~/Documents/privkey.key bbentley_smith_edu@unity.rc.umass.edu:/scratch3/workspace/bbentley_smith_edu-busco/logs/all_spp/"*.log" .
```
### Use R to extract the data and write to file:
```
file_list<-list.files(path = "../bat1k_compleasm/", pattern = "log")

all_spp<-list()
for(q in 1:length(file_list)){
  df<-read.table(file=paste0("../bat1k_compleasm/",file_list[q]), sep = "\t")
  df2<-as.data.frame(t(df[c(1,5:8),]))
  df3<-as.data.frame(t(gsub("*.:","",gsub("%.*","",df2))))
  colnames(df3)<-c("SPP","S","D","F","M")
  all_spp[[q]]<-df3
}
all_spp_df<-do.call(rbind, all_spp)
all_spp_df[, 2:5] <- lapply(all_spp_df[, 2:5], as.numeric)
all_spp_df$C<-all_spp_df$S + all_spp_df$D
all_spp_df <- all_spp_df[, c(1, 6, 2:5)]

write.table(file="../bat1k_compleasm/all_species_compleasm_stats.txt", all_spp_df,
            quote = F, row.names = F, sep = "\t")
```
### Reformat the dataframe to plot:
```
# Format for plot:
library(ggpubr)
library(viridis)

lplot_df<-list()
for(s in 1:30){
  df<-all_spp_df[all_spp_df$SPP == levels(factor(all_spp_df$SPP))[s],]
  df2<-as.data.frame(t(df))
  df2$SPP<-df2[1,1]
  df2$SCORE<-rownames(df2)
  df3<-df2[c(3:6),]
  colnames(df3)[1]<-"NUM"
  lplot_df[[s]]<-df3
}

plot_df<-do.call(rbind, lplot_df)
plot_df$NUM<-as.numeric(plot_df$NUM)
plot_df$SCORE <- factor(plot_df$SCORE, levels = c("S", "D", "F", "M"))
plot_df_or<-plot_df %>%
  mutate(SPP = factor(SPP, levels = species_order)) %>%  # Set factor levels for ordering
  arrange(SPP)  # Order the dataframe by the new factor levels
plot_df_or<-plot_df_or[plot_df_or$SPP != "HLrhiSin2",]

bp1<-ggbarplot(plot_df_or, 
          x = "SPP", 
          y = "NUM", 
          fill = "SCORE",  # Fill by SCORE
          color = "black",  # Outline color of the bars
          position = position_stack(),  
          palette = c("#31688EFF","#35B779FF","#FDE725FF","#440154FF"),  # Using viridis with 4 categories
          lab.size = 10,  
          xlab = "Species",  
          ylab = "Score",
          orientation = "horiz")
bp1

ggsave(filename = "../bat1k_compleasm/all_species_compleasm.pdf", plot = bp1, width = 12, height = 14)
```