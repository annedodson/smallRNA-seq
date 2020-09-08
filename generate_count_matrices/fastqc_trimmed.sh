#!/bin/bash/

# This script takes a fastq file and runs quality control with fastqc.
# USAGE: sh fastqc_trimmed.sh <full pathname of fastq file> <full pathname of the main directory>

# Load modules
module load gcc/6.2.0
module load fastqc/0.11.5

# Make the output directory if it doesn't already exist
mkdir -p ${2}/results/fastqc_trimmed

# Run quality control
fastqc -o ${2}/results/fastqc_trimmed $1