---
title: 'hyphy absrel'
disqus: hackmd
---

Selection tests with hyphy absrel
===
By Tanya Lama

### Objective: 
Run an exploratory round of selection tests on the full 23 bat & 7 nonbat species alignment for the bat1k longevity 

ABSREL (in the HYPHY package) performs branch-site model-based tests for episodic selection affecting a proportion of sites in the alignment along a proportion of branches in the tree. In other words, **is there evidence of selection *anywhere* in the alignment and *anywhere* on the tree**? 

### Acknowledgements
Thanks to Graham Hughes and Liliana Davalos for their expertise in creating this pipeline. 

## Table of Contents

[TOC]

# Address rhiSin hipLar issue
This is detailed elsewhere, but essentially Rhinolophus sinicus needed to be dropped from the alignment. We are awaiting a new alignment without Rhinolophus sinicus. For this pilot run of the pipeline, we will use 7,082 genes which are identical between rhiSin and hipLar (based on string comparison). See hackmd for details on how these 7,082 genes were identified. 

## Pull ~7k genes
We pulled all ~7k of our genes using the script pull_identical_genes.sh
```
#! /bin/bash
  
#
#SBATCH --job-name=identical_genes
#SBATCH --output=out_identical_genes.txt
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH -p gpu

#individuals_array=(ENST00000229135.IFNG   ENST00000297439.DEFB1   ENST00000367739.IFNGR1  ENST00000683941.IFNAR2)

for i in $(cat ./identical_genes.txt) #"${individuals_array[@]}"
do
cp -r ~/project_bat1k_longevity/data/alignments/bat_alignment_2021.12.07/rhiSin_hipLar_comparison/${i} ./${i}

done
```
We then confirmed that the list of pulled genes is the same length as the list of identical genes.
```
ls data/ > genes_pulled.txt 
wc-l genes_pulled.txt
wc -l identical_genes.txt
```
Both lists contain 7,082 genes so we are good to proceed.

## Drop rhiSin from with sed -e
This is a simple command that deletes the line including ">HLrhiSin2" and also deletes 1 line immediately following ">HLrhiSin2"
Let's try it on a gene called GDI1. We are only interested in the .cleanLb_hmm_manual.fasta for each gene.
`sed -e '/>HLrhiSin2/,+1d' ENST00000447750.GDI1.cleanLb_hmm_manual.fasta > ENST00000447750.GDI1.final.fasta`

This needs to be performed iteratively across all the genes we have selected. We will use the drop_rhisin.sh script
```
f#! /bin/bash

#
#SBATCH --job-name=drop_rhiSin
#SBATCH --output=drop_rhiSin.txt
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH -p gpu

for i in $(cat ./identical_genes.txt) 
do
sed -e '/>HLrhiSin2/,+1d' ./data/${i}/${i}.cleanLb_hmm_manual.fasta > ./7082_fastas/${i}.cleanLb_hmm_manual.fasta
done
```
Our fasta files have been saved in a folder called 7082_fastas. 
Let's check the length again to be sure we retained all 7,082 genes
```
ls 7082_fastas > final_fastas_genes.txt
wc -l final_fastas_genes.txt
```
Looks good. Now we have a folder full of fasta files that we can concatenate into a supermatrix with fasconcat

# Part 1: Tree
## fasconcat-g
FASconCAT is a tool for concatenating all of our gene alignments into a supermatrix. 
The perl script can be downloaded [here](https://github.com/PatrickKueck/FASconCAT-G). Note that this is a perl script that needs to be **in the same directory as your  fasta folder**.
FASconCAT will stitch all of the alignment files together in the same order, assuming all the fasta headers are the same across each file. It will also provide an output file detailing the coordinates of where each gene is in the supermatrix. 

I recommend testing this first with just a handful of genes in a separate test folder. Remember to run this on an interactive node

Usage is: 
perl FASconCAT-G_v1.05.pl [command options] 
	  (Linux/Mac)
      
      Where command options:
	  -p PHYLIP output (strict taxon names) #lets add this. We need it for iqtree
	  -n NEXUS output #lets add this just in case we need it for iqtree  -- apparently this is the partition file. We need this too
      -l partitioned 
	  -s START #include

These are the flags we have selected. You should end up with four files: FcC_info.xls  FcC_supermatrix.fas  FcC_supermatrix.nex  FcC_supermatrix.phy
```
perl FASconCAT-G_v1.05.pl -p -n -l -s
```
A couple of fastas were flagged by fasconcat as completely empty.
fasconcat/ENST00000390583.IGHD3-10.cleanLb_hmm_manual.fasta
fasconcat/ENST00000434970.TRDD2.cleanLb_hmm_manual.fasta
fasconcat/ENST00000454691.IGHD6-6.cleanLb_hmm_manual.fasta
fasconcat/ENST00000592904.ZNF266.cleanLb_hmm_manual.fasta
fasconcat/ENST00000632684.TRBD1.cleanLb_hmm_manual.fasta

We also dropped those which include 3 or fewer species:
fasconcat/ENST00000249007.RFPL3.cleanLb_hmm_manual.fasta
fasconcat/ENST00000342868.BTNL3.cleanLb_hmm_manual.fasta
fasconcat/ENST00000439466.CT47A5.cleanLb_hmm_manual.fasta
fasconcat/ENST00000592904.ZNF266.cleanLb_hmm_manual.fasta

Concatenating all 7,000 genes wasn't possible in an interactive session, so I wrote fasconcat.sh which will run for 24 hours on gpu-long. fasconcat.sh needs to sbatched from within the /7082_fastas folder. 
```
#! /bin/bash
  
#
#SBATCH --job-name=fasconcat
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=24:00:00 #likely excessive. Try 8h
#SBATCH -p gpu-long

perl ./FASconCAT-G_v1.05.pl -n -p -l -s 
```
### Review the outputs
The outputs we are most interested in are supermatrix_partition.nex and supermatrix_partition.txt. These are the files we will use for iqtree in the next step. 
The supermatrix is in NEXUS format.
The partition information is in raxml format: 
```
DNA, ENST00000000233.ARF5.cleanLb_hmm_manual.fasta = 1-540
DNA, ENST00000000412.M6PR.cleanLb_hmm_manual.fasta = 541-1374
DNA, ENST00000001146.CYP26B1.cleanLb_hmm_manual.fasta = 1375-2910
```
## iqtree: select best model of gene evolution and estimate phylogeny
We will create a conda environment to house iqtree and our other dependencies for this step. See [here](https://hackmd.io/@tlama/wgspipeline) for specific instructions on building conda environments, if they are new to you.
`/home/tlama_umass_edu/anaconda3/bin/conda create --name iqtree`
Environment location: `/home/tlama_umass_edu/anaconda3/envs/iqtree` 
#Activate the environment
`source activate iqtree` #you need to be in an interactive session for this
#Load any necessary modules
`conda install -c bioconda iqtree`

After much trial and error this is the final command we have landed on with iqtree: 
```
#! /bin/bash
  

#
#SBATCH --job-name=iqtree-partitioned
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=2
#SBATCH --time=12:00:00 #ran for 9 hours with 7k genes
#SBATCH -p gpu-long

iqtree -s *_supermatrix_partition.nex -p *_supermatrix_partition.txt --prefix supermatrix_partition_iqtree -B 1000 -T AUTO -nt AUTO

#-m GTR+I+G was an option we've used previously. I suggest leaving it as is and letting iqtree and ModelFinder step through all 299 models of sequence evolution, then pick the best model based on AIC. For expediency we may try -m TESTMODEL in the future to limit the number of potential models evaluated by ModelFinder.
#-nt is also optional. In our last successful run we only used -T AUTO (not both)

#step1 infers a concatentation-based species tree with 1k ultrafast bootstraps and an edge-linked partition model
```
Essentially what we have done here is provide a supermatrix and a partition file which says what range each gene is. iqtree will determine what model of sequence evolution is the best fit for each gene, then use maximum likelihood (model given the data) to optimize a tree from the supermatrix, then run 1k bootstraps to resolve the branchlength and topology of the best candidate tree and arrive at a consensus. 

### Review the outputs
The main output of interest is supermatrix_partition_iqtree.contree
This is the tree that has been selected based on the best model of sequence evolution for each partition, optimized by ML, and resolved by 1k ultrafast bootstraps. 
We will use this tree as our input for hyphy. But first let's check it out using figtree

## Visualize your tree in figtree
Download compiled binaries for Mac figtree [here](https://github.com/rambaut/figtree/releases) and install accordingly. You will also need to install an updated version of java. figtree is a desktop application. 
Download supermatrix_partition_iqtree.contree from the cluster and rename the suffix as **supermatrix_partition_iqtree.tree**
Open with figtree
figtree will notice immediately that you have an "added" field -- these are our boostraps. Name them boostraps accordingly
On the left, select Node Labels and view the aforementioned bootstraps on our tree: 
![](https://i.imgur.com/ZSqDE2a.png)
Note that we have one node that is relatively unresolved (Laurasiatheria). 

We can also use figtree to look at the branch lengths
![](https://i.imgur.com/U149KnL.png)

## Drop bootstraps from supermatrix_partition_iqtree.tree
We need to drop the bootstraps from our tree, because hyphy can't handle the added field. It will run, but will assume the boostrap is a branch length. 
We can remove this by manually deleting the 100 and 47 scores from the raw file, and leaving the ":" in place. You can check the output in figtree once more, being sure that the branch lengths and topology are the same. 

Save this vesion of your tree as **supermatrix_partition_iqtree_no_bootstraps.treefile**

# Part 2: Create batches
Move all *.cleanLb_hmm_manual.fasta files to a folder called fastas
```
ls /fastas > all_fastas
wc -l all_fastas #We have 10893 genes to run through hyphy
``` 

Split all_fastas into batches of 2000
```
sed -n -e '1,2000p' all_fastas.txt >batch1.txt
sed -n -e '2001,4000p' all_fastas.txt >batch2.txt
sed -n -e '4001,6000p' all_fastas.txt >batch3.txt
sed -n -e '6001,8000p' all_fastas.txt >batch4.txt
sed -n -e '8001,10893p' all_fastas.txt >batch5.txt #this one is a little longer, has 2893 genes
```
Move fastas accordingly to their folder in /work (repeat for all five batches)
`for i in $(cat ./batch1.txt); do cp fastas/${i} /work/talama_stonybrook_edu/project_bat1k_longevity/data/hyphy/batch1/; done`

# Part 3: Clip trees
Not all 17,000 genes are present for each species (e.g., a gene may be missing for myoNig or artJam). So, for each gene we need to trim the tree and create an input for hyphy. Hyphy absrel will then be run individually for each gene. 
## We have a folder for each batch of fastas. Copy your tree, modclip.pl, and cliptrees.R
## Execute perl modclip.pl from within the folder
```
@species=(); 
open(IN,"03.16.22_tree.tre");  ##opens the unrooted tree
while(<IN>){
    while($_=~s/([A-Za-z]+[0-9]+)\://){
    	push @species,$1;	    ##creates species in the supermatrix tree
    }
}
foreach $x(@species){
    #print $x."\n"; 
}
open(OUT,">details.txt"); 
    @array=(<*fasta>);                         ##for each alignment see what species are in the alignment. if there are species missing, make note of them in details.txt
    foreach $x(@array){
    $y=""; 
    if($x=~m/([\S]+)/){
    $gene=$1;
    }
    print OUT $gene;
    open(IN2, "$x"); 
    @keep=();
    while(<IN2>){
    if($_=~m/\>([\S]+)/){
    push @keep,$1;
    }
    }
    foreach $y(@species){
    if($y~~@keep){
    }
    else{
    print OUT "|" .$y;
    }
    }
    print OUT "\n";
}

`Rscript cliptrees.R`;

##note: Smartmatch is experimental at line 27 (S experimental::smartmatch) This warning is emitted if you use the smartmatch (~~) operator. This is currently an experimental feature, and its details are subject to change in future releases of Perl. 
```
## This should have created a file called details.txt which looks like this:
```
ENST00000323110.CYP11B2.cleanLb_hmm_manual.fasta|canFam4|bosTau9|equCab3|HLtadBra2 #missing four species
```
## Execute Rscript cliptrees.R from within the folder (modclip.pl now runs this on its own)
Before running this script, I needed to start an R session and install a couple of packages (ape, phangorn) and some dependencies. Unfortunately Unity uses R 3.6.3, which is rather outdated. I had to download an archived version of phangorn (see /project/tlama_umass_edu/bin) and install from source. 

‘quadprog’, ‘igraph’, ‘fastmatch’ are not available for package ‘phangorn’ install.packages("/project/tlama_umass_edu/bin/igraph_2.0.2.tar.gz", repos = NULL, type = "source")
install.packages("/project/tlama_umass_edu/bin/phangorn_2.5.5.tar.gz", repos = NULL, type = "source")

`tab<- read.table("details.txt")
       library("phangorn")
       library("ape")
       i<- 1
    while(i<nrow(tab)+1){
                          main<-read.tree("03.16.22_tree.tre")
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
    }`

## cliptrees.R should have created a .tfl file for each fasta!

## create a gene list for each batch 
ls ENST00000*.cleanLb_hmm_manual.fasta > batch1_genelist.txt
# Part 4: Selection
## Install anaconda3 if needed and export the $PATH to your bashrc 
```
export PATH="/home/talama_stonybrook_edu/bin/anaconda3/bin/:$PATH"
```
## Create a conda environment for hyphy
We use conda environments to load all of our dependences. See [here](https://hackmd.io/@tlama/wgspipeline) for specific instructions on building conda environments, if they are new to you.
This tutorial assumes the availability of HyPhy 2.2.x, either through conda or a local installation. Please download and install the program following the instructions at the [HyPhy project page](https://github.com/veg/hyphy) if you need to.
`conda create --name hyphy hyphy`
Environment location: `/home/tlama_umass_edu/anaconda3/envs/hyphy` 
#Activate the environment
`source activate hyphy`
#Load any necessary modules
```
module load gcc
```
## Submitting hyphy jobs
```
#!/bin/sh -l 
#SBATCH -job-name-reverselongbats
#SBATCH -t 312:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node 22
module load gcc

./hyphy absrel CPU=44 --pvalue 0.05 --alignment /ENST#####.cleanLb_hmm_manual.fasta --tree /ENST#####.cleanLb_hmm_manual.fasta.tfl >> /work/ENST#####.cleanLb_hmm_manual_absrel_output.txt

## --output /work/###_absrel_output.json
```

## Restarting disrupted hyphy runs
A batch was accidentally canceled. I want to restart it, without rerunning hyphy on all the genes for which there are completed results 
ls results/* > completed
cat completed | grep -Po '(?<=(results/)).*(?=_absrel)' >> alldone
```
tail -n +37 batch1_genelist.txt > batch1_genelist_v2.txt #drop the first 40 lines of the batch1 genelist
```
## Run hyphy with command line arguments
I recommend running hyphy -i (interactive mode) before using the command line arguments. A list of the available <additional_method_specific_arguments> can be seen by running hyphy <method_name> --help

Usage is: 
```hyphy <method_name> --alignment <path_to_alignment_file> <additional_method_specific_arguments>```

After much trial and error, this is what we have landed on as the best configuration for absrel 
```
#! /bin/bash

#
#SBATCH --job-name=7082-absrel-all
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=2
#SBATCH --time=336:00:00
#SBATCH -p gpu-long

hyphy absrel --alignment ./supermatrix_with_partitions.nex --tree ./supermatrix_iqtree_no_bootstraps.treefile --branches All --output 7082_genes_all_branches.ABSREL.json #Graham ran this on 44 cores on 2 nodes for 2 weeks

#hyphy <method_name> --alignment <path_to_alignment_file> <additional_method_specific_arguments>
#hyphy absrel --alignment ./FcC_supermatrix.phy --tree ./FcC_supermatrix.phy.treefile
```
## Alignment and tree
in fasta, nexus or phylib file format
/home/tlama_umass_edu/project_bat1k_longevity/data/alignments/bat_alignment_2021.12.07/out_longevityBats_v3/ENST00000000233.ARF5/ENST00000000233.ARF5.cleanLb_hmm.fasta

This is going to be a large quantity of data, so please run this section in an interactive session. 
```
bsub -q interactive -W 120 -Is bash ## on UMass
srun --partition short-40core --cpus-per-task 4 --mem-per-cpu 1000 --time 60:00 --pty bash ##on seawulf
```
Load dependencies
```
module load awscli/1.16.144`
```
Navigate to mPhyDis1 in the genomeark bucket
```
aws s3 --no-sign-request ls s3://genomeark/species/Phyllostomus_discolor/mPhyDis1/
```
Download the transcriptomic data
```
aws s3 --no-sign-request cp s3://genomeark/species/Phyllostomus_discolor/mPhyDis1/transcriptomic_data/ ./ --recursive
```
Download the curated mPhyDis1 assembly
```
aws s3 --no-sign-request cp s3://genomeark/species/Phyllostomus_discolor/mPhyDis1/assembly_curated/mPhyDis1.pri.cur.20200504.fasta.gz ./
```


Transfer all to Unity /
scp -r ./phyDis/ talama@login.seawulf.stonybrook.edu:/gpfs/scratch/talama/
scp -r ./phyDis/ tlama_umass_edu@unity.rc.umass.edu:/home/tlama_umass_edu/scratch/project_shrew_genome/data/
# Appendix and FAQ

## Alternative tree-building methods
ASTRAL estimates a tree for each locus and then combines the gene trees into a species tree. An interesting alternative is to use the gene tree species tree method described in ASTRAL.
To run ASTRAL you first need to generate unrooted gene trees as the input (makes sense) but I don’t know how to do that either. I might as well keep it simple for now.
it's unbelievably easy, split the aln into its loci and run raxml with these few taxa you'll get each tree in <10 s
5:22
by default most methods generate unrooted trees BT
:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `scripts`

# SCRATCH
Graham method of making tree and alignment input 

Can you tell me what method/tool you use to estimate the tree? OMA, SUMAC are a couple I’ve seen. I’m just trying to learn a bit about each step and haven’t done that before

1. concatenate the genes into a supermatris using FASconcat
2. FASconcat is a perl script that you just place in the same directory as your alignment files, and it will stitch them all together in the same order assuming all the fasta headers are the same across the files.It also gives you an output file that you can use to get the coordinates of where each gene is in the supermatrix
3. Use iqtree to make an ML tree. This required a partition file saying what range each gene is, and what model of sequence evolution was best fit for that gene. iqtree finds the best fit model for each gene first, before we make the ML tree using the supermatrix. 


## ASTRAL
We installed ASTRAL into the iqtree conda environment
conda install -c phylofisher astral

We also require java, so we loaded it as a module
module load java/11.0.2

Test it is running using
java -jar /home/tlama_umass_edu/anaconda3/envs/iqtree/ASTRAL-5.7.3/astral.5.7.3.jar

-i
a file containing input gene trees in newick format


## Running iqtree on ten genes
1.perl FASconCAT-G_v1.05.pl -s -n -p #just need p really
2.iqtree -s FcC_supermatrix.phy #we finally got the treefile
3.Take the gene treefile and use it as an input for astral 
java -jar /home/tlama_umass_edu/anaconda3/envs/iqtree/ASTRAL-5.7.3/astral.5.7.3.jar -i FcC_supermatrix.phy.treefile -o astral/FcC_supermatrix.phy.species.treefile

the stuff that worked is at: /home/tlama_umass_edu/project_bat1k_longevity/data/alignments/bat_alignment_2021.12.07/ten_genes/ten_genes/ten_genes


1. iqtree -p /home/tlama_umass_edu/project_bat1k_longevity/data/alignments/identical_genes/final_fastas/ --prefix concat2 -B 1000 -T AUTO #fail segmentation error

1. iqtree -s FcC_supermatrix.phy --prefix concat -B 1000 -T AUTO  #this worked!!!
3. iqtree -S /home/tlama_umass_edu/project_bat1k_longevity/data/alignments/identical_genes/final_fastas/ --prefix loci2 -T AUTO #step2 infer the locus trees
4. iqtree -t concat.treefile --gcf loci.treefile -p /home/tlama_umass_edu/project_bat1k_longevity/data/alignments/identical_genes/final_fastas/ --scf 100 --prefix concord -T 10 #step3 compute concordance factors

this worked!
iqtree -s FcC_supermatrix.phy --prefix concat -B 1000 -T AUTO
this worked!
iqtree -t concat.treefile --gcf loci.treefile -p data/ --scf 100 --prefix concord -T 2

this was all accomplished using the 10 genes dataset. Can we do the same with the complete dataset? Start with fasconcat (scripted) to get *.phy and go from there.



1. re-run fasconcat and iqtree on 10-genes with rhiSin dropped
2. use FcC_supermatrix.phy.treefile as the input for absrel

Execute the hyphy_absrel.sh script. The interactive version of hyphy is great for troubleshooting. Options can be reviewed using hyphy absrel --help. 
```
#hyphy <method_name> --alignment <path_to_alignment_file> <additional_method_specific_arguments>
hyphy absrel --alignment ./FcC_supermatrix.phy --tree ./FcC_supermatrix.phy.treefile
```

Final changes
1. We needed to add the -l flag to fasconcatg-g to generate the gene partition file. This flag generates two types -- a partition.txt file in raxml format and a partition.nex file in nexus format. 
2. We needed to add the -p flag to iqtree to identify the gene partition file and allow modelfinder to find a unique model of sequence evolution for each gene. 
iqtree -s FcC_supermatrix.phy -p FcC_supermatrix_partition.txt -m GTR+I+G

iqtree -s FcC_supermatrix_partition.nex -p FcC_supermatrix_partition.txt -m GTR+I+G

3. Substition model -m GTR+I+G 
Note that using RAxML-style partition.txt file, all partitions will use the **same rate heterogeneity model** given in -m option. If you want to specify different models for different genes you can do that with the partition.nex file.

**What is it doing?** Here, IQ-TREE partitions the alignment FcC_supermatrix.phy into n sub-alignments containing sites (columns) 1-100 and 101-384, etc. IQ-TREE then applies the subtitution model GTR+I+G to each partition, respectively. Substitution model parameters and trees with branch lengths can be found in the result file FcC_supermatrix.nex.iqtree.

hyphy absrel --alignment ./FcC_supermatrix.fas --tree ./FcC_supermatrix.phy.treefile --branches Internal --output internal.ABSREL.json

#looks like a nexus file is necessary for performing a partitioned analysis with hyphy

## Looking for recombination breakpoints with GARD
Motivation: Phylogenetic and evolutionary inference can be severely misled if recombination is not accounted for, hence screening for it should be an essential component of nearly every comparative study. The evolution of recombinant sequences can not be properly explained by a single phylogenetic tree, but several phylogenies may be used to correctly model the evolution of non-recombinant fragments.

![](https://i.imgur.com/Z8igmoC.png)
![](https://i.imgur.com/CqoVHwD.png)
If GARD indicates that you have recombination in your alignment and different partitions have different tree topologies, then you need to interpret aBSREL results carefully.

This is because aBSREL operates on branches, and if you have different topologies, some of the internal branches may exist in one tree but not the other. The only branches that are guaranteed to exist in all these topologies are leaf branches. For example, the blue branches highlighted in the each trees below do not have a match in the other tree.

My recommendation would be to run aBSREL on individual partitions separately, and without selecting branches to test in one partition based on the results in the other (this would induce a testing bias, calling HARK : hypothesizing after results are known). You may then identify branches that appear to be under significant positive selection in both partitions, assuming they exist in both partitions: this may be limited to terminal branches only, depending on what the topologies look like.

hyphy absrel --alignment ./FcC_supermatrix_partition.nex --tree ./FcC_supermatrix_partition.treefile --output gard.ABSREL.json

-m TESTONLY
Perform standard model selection like jModelTest (for DNA) and ProtTest (for protein). Moreover, IQ-TREE also works for codon, binary and morphogical data. If a partition file is specified, IQ-TREE will find the best-fit model for each partition.

-m is a powerful option to specify substitution models, state frequency and rate heterogeneity type. The general syntax is:
-m MODEL+FreqType+RateType

from login 
ssh-keygen ##this needs to be generated WITHOUT a passphrase
/home/tlama_umass_edu/.ssh/id_rsa
key fingerprint SHA256:+1sykLNwcZEIjtzE803GxYd4NwWXDb+wzeFz5Lxm0Sc
cat /home/tlama_umass_edu/.ssh/id_rsa.pub
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABgQCtPDsOXs02o/wgjsFCH0879lr4FSl5RhodS3XYXwqkwhiFX9EbR2hCFRI5c8N6uR8RVrXNakVFIq2do30eBisQh+HEtMRPOiuuuofYvAUYHjq+8SAMHLaKn1qTjMBvF/fAa0dqgBE4Nbh9fM74lyuw7fxq/S10xVYnHzwwrt7mhDshdi5vcHVZrC0MO7qrhj3AmBVjFfNF/eIpHjpgqZe610iU+vRXNqGp9IeIhk0nQ+akdOsUDpT92UP2kxzLzWIS1R3vNlzj1ZmT6i6IrpD68Y3Xkkn2OrGeHJ6SHEjqxRvZslSCBorr1ZJktKG+x8grI5sqT5L2RdpGYJravR5pXSd6b5VW3G/81DMzUe9goepHHrz/giNSeTqw528NQMevUEjI+fVL2oMlGu3myDIBVwi4EWfRpgkUMjtw6wLRW6VDCj5oO52fwfIVjgB6hhEFJ2IomGtVMKNYJi0OVwuhqooPKPM7f4ZSI5WT6mo+Nu136vYd0dSiBctQ1+Whq7s= tlama_umass_edu@login #add to unity

/bin/bash
ssh login

from an interactive session on cpu
/home/tlama_umass_edu/scratch/project_shrew_genome/scripts/modified_doBlastzChain.pl `pwd`/DEF -verbose=2 -noDbNameCheck -workhorse=login -bigClusterHub=login -skipDownload -dbHost=login -smallClusterHub=login -trackHub -fileServer=login -syntenicNet > do.log

## Errors
hgsql: command not found --ignore?
> [/home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz/doPartition.bash: line 25: twoBitToFa: command not found
/home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz/doPartition.bash: line 26: faToTwoBit: command not found]

This is a little weird, when I installed all of the scripts and tools (see bill's instructions above) it created an unnamed folder in /bin and saved them there. And now scripts are having a little bit of trouble finding them (/home/tlama_umass_edu/bin//). So I've gone through the script and just put in hard paths which is not a perfect solution, but I'm not sure how to direct the scripts to the tools otherwise. This has resolved the above errors about twoBitToFa and faToTwoBit.

Note that you need to clear the /tmp and /trackData folders every time you test this script.
> [Command failed:
> + /scratch/tlama_umass_edu/project_shrew_genome/scripts/partitionSequence.pl 175000000 0 /home/tlama_umass_edu/project_bat1k_longevity/data/annotations/hg38/hg38.2bit /home/tlama_umass_edu/project_bat1k_longevity/data/annotations/hg38/hg38.chrom.sizes -xdir xdir.sh -rawDir ../psl 18 -lstDir tParts
lstDir tParts must be empty, but seems to have files  (part000.lst ...)]


ssh -x -o 'StrictHostKeyChecking = no' -o 'BatchMode = yes' login nice /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz/doPartition.bash

-x option: The executable is run on the remote machine and displayed (drawn) on the local machine. What ssh -X remote does is start up a proxy X11 server on the remote machine.
-o option

the problem was in setting path to /home/tlama_umass_edu/bin/
should be /home/tlama_umass_edu/bin
#delete and restart?? or just change path?

echo 'export PATH=/home/tlama_umass_edu/bin:/home/tlama_umass_edu/scratch/project_shrew_genome/scripts:$PATH' >> $HOME/.bashrc

## What I did
run to failure
edit doPartition.bash (faTwoBit)
clear /trackData/tparts and /trackData/qparts
run ssh -x -o 'StrictHostKeyChecking = no' -o 'BatchMode = yes' login nice /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz/doPartition.bash

/home/tlama_umass_edu/scratch/project_shrew_genome/scripts/modified_doBlastzChain.pl `pwd`/DEF -verbose=2 -noDbNameCheck -workhorse=login -bigClusterHub=login -skipDownload -dbHost=login -smallClusterHub=login -trackHub -fileServer=login -syntenicNet -continue blastz > do.log

/home/tlama_umass_edu/scratch/project_shrew_genome/scripts/tanya_modified_doBlastzChain.pl `pwd`/DEF -verbose=2 -noDbNameCheck -workhorse=login -bigClusterHub=login -skipDownload -dbHost=login -smallClusterHub=login -trackHub -fileServer=login -syntenicNet > do.log


gensub2

copy modified_do_blast and edit it
line 585 add the path to faToTwoBit
586
594
595
648-651

cd /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz

gensub2 hg38.lst myoNig_masked.lst gsub jobList
cp /home/tlama_umass_edu/scratch/project_shrew_genome/scripts/diana_scripts/scriptSetup.sh .
sh scriptSetup.sh

cp /home/tlama_umass_edu/scratch/project_shrew_genome/scripts/diana_scripts/scriptSetup.sh .
sh scriptSetup.sh
scriptSetup.sh: 8: [[: not found
sbatch: error: Batch job submission failed: Invalid job array specification

HgStepManager: executing step 'blastz' Sun Feb 13 16:02:55 2022.#we got through the partition step successfully
-continue blastz 

/home/tlama_umass_edu/scratch/project_shrew_genome/scripts/tanya_modified_doBlastzChain.pl `pwd`/DEF -verbose=2 -noDbNameCheck -workhorse=login -bigClusterHub=login -skipDownload -dbHost=login -smallClusterHub=login -trackHub -fileServer=login -syntenicNet -continue blastz > do.log

Command that failed: ssh -x -o 'StrictHostKeyChecking = no' -o 'BatchMode = yes' login nice /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz/doClusterRun.csh

COUNT=`wc -l jobList`
ARRAY=$(echo $COUNT | awk '{print $1}')
cp /home/tlama_umass_edu/scratch/project_shrew_genome/scripts/diana_scripts/doArrayJob.sh . #copy the slurm script for 1-1990 jobs
cp /home/tlama_umass_edu/scratch/project_shrew_genome/scripts/diana_scripts/doArrayJob_2nd.sh . #copy the slurm script for jobs 2000-COUNT

sed -i "s/CHANGE/$ARRAY/g" doArrayJob_2nd.sh #Set the number of total jobs to be submited

if '2448' -le '1990'; then
        sed -i "s/1990/$ARRAY/g" doArrayJob.sh
        sbatch doArrayJob.sh ;
else
        echo "submiting two array jobs" 
        sbatch doArrayJob.sh
        sleep 10m
        X=`squeue -u tlama_umass_edu | wc -l`
        while [ $X -gt 1000 ]; do sleep 1m; X=`squeue -u tlama_umass_edu | wc -l`; done
        sbatch --dependency=singleton --job-name=Array doArrayJob_2nd.sh
fi


NAMESFILE=jobList

#Creates a job for each line in the jobList file
CHRANGE=$(awk $NAMESFILE)

## submitting array jobs
This actually worked: 
```
#!/bin/bash
#SBATCH --output=batch_%A_%a.out
#SBATCH --error=batch_%A.%a.err
#SBATCH --job-name=Array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --array=1-10

#The jobList generated by gsub script
NAMESFILE=/home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz/jobList

#Creates a job for each line in the jobList file
CHRANGE=$(awk'{print $NAMESFILE}')

#This will run the line for each jobList, that will call blastz-run-ucsc
${CHRANGE}
```

## jobList permissions
run chmod +x on jobList before running doArray.sh
copy /scripts to /bin/scripts

Command failed:
export PATH=/home/tlama_umass_edu/bin:/home/tlama_umass_edu/scratch/project_shrew_genome/scripts:/bin:/usr/bin; 
/bin/cp /dev/null /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/tmp/blastz.fsgEJx/local.lav
temp files in /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/tmp/blastz.3WeMfv -- will attempt cleanup before die-ing.
/bin/sh: 1: set: Illegal option -o pipefail

blastz-run-ucsc target query DEF out

/scratch/tlama_umass_edu/project_shrew_genome/scripts/blastz-run-ucsc /home/tlama_umass_edu/project_bat1k_longevity/data/annotations/hg38/hg38.2bit:chr1:0-175000000 /home/tlama_umass_edu/project_bat1k_longevity/data/annotations/myoNig/myoNig_masked.2bit:myoNig_3:0-50010000 ../DEF ../psl/hg38.2bit:chr1:0-175000000/hg38.2bit:chr1:0-175000000_myoNig_masked.2bit:myoNig_3:0-50010000.psl

1. changed paths in DEF
2. changed parameters to match Bill's README
3. changed myoNig files (chrom.sizes) and "2bit" to match Bill's example DEF file

We still have issues with the output 
Command failed:
export PATH=/home/tlama_umass_edu/bin:/home/tlama_umass_edu/bin/scripts:/bin:/usr/bin; 
/bin/cp /dev/null /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/tmp/blastz.QWlFpv/local.lav
temp files in /home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/tmp/blastz.GVMj2e -- will attempt cleanup before die-ing.
/bin/sh: 1: set: Illegal option -o pipefail


/usr/bin

PATH=/home/tlama_umass_edu/bin:/usr/bin;/bin/cp /dev/null

## doArray.sh
#!/bin/bash
#SBATCH --output=batch_%A_%a.out
#SBATCH --error=batch_%A.%a.err
#SBATCH --job-name=Array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --array=1-2

#The jobList generated by gsub script
NAMESFILE=/home/tlama_umass_edu/scratch/project_shrew_genome/data/genomes/myoNig/trackData/run.blastz/jobList

#Creates a job for each line in the jobList file
CHRANGE=$NAMESFILE

#This will run the line for each jobList, that will call blastz-run-ucsc
${CHRANGE}


## scriptsetup.sh
COUNT=`wc -l jobList`
ARRAY=$(echo $COUNT | awk '{print $1}')
cp /home/tlama_umass_edu/scratch/project_shrew_genome/scripts/diana_scripts/doArrayJob.sh . #copy the slurm script for 1-1990 jobs
cp /home/tlama_umass_edu/scratch/project_shrew_genome/scripts/diana_scripts/doArrayJob_2nd.sh . #copy the slurm script for jobs 2000-COUNT

sed -i "s/CHANGE/$ARRAY/g" doArrayJob_2nd.sh #Set the number of total jobs to be submited

if $ARRAY -le '1990'; then
        sed -i "s/1990/$ARRAY/g" doArrayJob.sh
        sbatch doArrayJob.sh ;
else
        echo "submiting two array jobs" 
        sbatch doArrayJob.sh
        sleep 10m
        X=`squeue -u tlama_umass_edu | wc -l`
        while [ $X -gt 1000 ]; do sleep 1m; X=`squeue -u tlama_umass_edu | wc -l`; done
        sbatch --dependency=singleton --job-name=Array doArrayJob_2nd.sh
fi

echo Wait until all submitted jobs are complete

X=`squeue -u tlama_umass_edu | wc -l`

while [ $X -gt 1 ]
do
        echo There are $X jobs running on `date`
        sleep 20m
        X=`squeue -u tlama_umass_edu | wc -l`
done

echo Check if all psl files exists, as a control point

awk '{print $NF}' jobList | while read file; do if [[ ! -f ${file} ]]; then echo WARNING: ${file} does not exist, check error manually >> failed_jobs.txt; fi; done

rm *err *out

## Graham April 7th
- supermatrix is done
- patterns didn't match -- horse, human, mouse kind of screwy
- also forced the topology from the jebb et al. paper
- Graham will send the treefile -- from the 17,970 genes supermatrix (final alignment data)
- for i in j, run hyphy absrel
