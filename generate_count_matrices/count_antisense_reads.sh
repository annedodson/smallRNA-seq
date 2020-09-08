#!/bin/bash/

# This script takes .bam files and counts the number of reads that are antisense to each gene using featureCounts.
# USAGE: sh count_reads.sh <full pathname of bam file> <number of cores> <full pathname of the main directory>

# Assign variable to the gene annotation file  
## Make sure it matches the version of the genome used for mapping
gtf=/n/groups/kennedy/anne/WBcel235/Annotation/Genes/genes.gtf

# Assign variables to the outputs
base=$(basename $1 _trim3_trim5_merge.bam)
output_dir=${3}/results/counts_antisense
log_dir=${3}/logs/counts_antisense

# Create output directories if they don't already exist
mkdir -p $output_dir $log_dir

# Add featureCounts to $PATH
export PATH=/n/app/bcbio/tools/bin:$PATH

# Count antisense reads
## -T option specificies the number of cores
## -s 2 option counts only the antisense reads
featureCounts -T $2 -s 2 -a $gtf -o ${output_dir}/${base}_featureCounts.txt $1

# Move the summary file to the logs directory
mv ${output_dir}/${base}_featureCounts.txt.summary $log_dir

# Make a simplified count matrix for future analysis in R 
## First, cut out all the columns except for the gene IDs and counts
## Then, delete the first line
## Then, strip the pathname in the column 2 header to just the sample name
## Lastly, save as an .Rmatrix.txt file
cut -f1,7,8,9,10,11,12 ${output_dir}/${base}_featureCounts.txt | sed '1d' | awk 'NR==1{$2=a}1' a=$base > ${output_dir}/${base}_featureCounts.Rmatrix.txt
