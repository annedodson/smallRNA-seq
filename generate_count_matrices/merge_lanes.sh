#!/bin/bash/

# This script merges fastq files from samples that were run over multiple Illumina lanes.
# USAGE: sh merge_lanes.sh <full pathname of the main directory>

# Assign a variable to the input directory
input_dir=${1}/results/trim3_trim5

# Create output directory if it doesn't already exist
mkdir -p ${1}/results/trim3_trim5_merge

# Create a list of output sample names
## Use metadata.txt to creat a list of the input sample names
## Remove the lane designation of the sample name (_L*) as well as anything that follows (e.g. .fastq)
## Remove duplicate names with sort
samples=$(awk 'FNR > 1 {print $1}' ${1}/meta/metadata.txt | awk -F '_L' '{print $1}' | sort -u)

for sample in $samples
do
	# Assign variable to the output filename
	output=${1}/results/trim3_trim5_merge/${sample}_trim3_trim5_merge.fastq
	# Assign variable to the input filenames
	input=$(ls -d -1 $input_dir/*.fastq | grep $sample)
	# Merge all input files into one output file
	cat $input > $output
done
