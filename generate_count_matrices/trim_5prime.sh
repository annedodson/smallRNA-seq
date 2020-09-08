#!/bin/bash/

# This script takes a fastq file and trims off the 5' barcode.
# USAGE: sh trim_5prime.sh <full pathname of fastq file> <full pathname of the main directory>

# Load modules
module load gcc/6.2.0
module load python/2.7.12
module load cutadapt/1.14  # Need to load python before loading cutadapt

# Make the output directories if they don't already exist
mkdir -p ${2}/results/trim3_trim5
mkdir -p ${2}/logs/trim5

# Assign variables to the output filename
base=$(basename $1 .fastq)
output=${2}/results/trim3_trim5/${base}_trim5.fastq

# Trim the 5' barcode using cutadapt
## The -u option specifies the number of bases to remove from the beginning of the read
## Note, this command trims off the first 4 bases, no matter what the sequence

cutadapt -u 4 -o $output $1 > ${2}/logs/trim5/${base}.txt