#!/bin/bash

#move this to results folder for branch shift then run it
#example inputs bash hypothal.output.converter -f BSThetaTestLRTsadult_hypothal_BS2ALV.res -s 16 -g genelist -o result
while getopts f:s:g:o: flag
do
    case "${flag}" in
        f) filename=${OPTARG};;
        s) species_num=${OPTARG};;
        g) genelist=${OPTARG};;
        o) output=${OPTARG};;
    esac
done
var=1
header=$(($species_num+$var))

#how to convert the output
# take all of the spaces of LRT file and turn them into linebreaks, then save as a new file .LTR
tr ' ' '\n' < "../analysis/results_branchshift/"$filename > "../analysis/results_branchshift/"$filename".LRTs"
# use awk to convert it out of scientific notation
awk '{printf "%8.6f\n", $1}' "../analysis/results_branchshift/"$filename".LRTs" > "../analysis/results_branchshift/"$filename".LRTs.nosci"
#paste it next to a file with just the ENST gene names
paste "../data/10_genelist/"$genelist "../analysis/results_branchshift/"$filename".LRTs.nosci" > "../analysis/results_branchshift/"$output
