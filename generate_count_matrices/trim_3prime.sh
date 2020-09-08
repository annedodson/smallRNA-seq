#!/bin/bash/

# This script takes a fastq file and trims off the 3' adapter.
# USAGE: sh trim_3prime.sh <full pathname of fastq file> <adaptor sequence> <full pathname of the main directory>

# Load modules
module load gcc/6.2.0
module load python/2.7.12
module load cutadapt/1.14  # Need to load python before loading cutadapt

# Make the output directories if they don't already exist
mkdir -p ${3}/results/trim3
mkdir -p ${3}/logs/trim3

# Assign variables to the output filename
base=$(basename $1 .fastq)
output=${3}/results/trim3/${base}_trim3.fastq

# Trim the 3' adapter using cutadapt
## The -a option finds the given adapter sequence (whole or partial) and removes that sequence + any sequence that follows
## The -m option removes reads that are shorter than 14bp following adapter removal
## --discard-untrimmed option removes reads in which no adapter was found
cutadapt -a $2 -m 14 --discard-untrimmed -o $output $1 > ${3}/logs/trim3/${base}.txt