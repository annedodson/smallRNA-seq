#!/bin/bash/

# This script takes fastq files and maps them using bowtie.
# USAGE: sh map_reads.sh <full pathname of fastq file> <number of cores> <full pathname of the main directory>

# Assign variable to the reference genome
genome=/n/groups/kennedy/anne/WBcel235/Sequence/BowtieIndex/genome

# Assign variables to the outputs
base=$(basename $1 .fastq)
output_dir=${3}/results/alignment
log_dir=${3}/logs/alignment

# Create output directories if they don't already exist
mkdir -p $output_dir $log_dir

# Load modules
module load gcc/6.2.0
module load bowtie/1.2.2
module load samtools/1.9
module load python/2.7.12
module load deeptools/3.0.2

# Map to the genome using Bowtie
## -q option specifies that the reads file is in fastq format
## -n option specifies the number of allowed mismatches
## -p option specifies the number of cores requested
bowtie -q -n 0 -p $2 $genome $1 ${output_dir}/${base}.sam -S 2> ${log_dir}/${base}.txt

# Use samtools to generate a .bam file and sort by position
samtools sort ${output_dir}/${base}.sam > ${output_dir}/${base}.bam

## Use samtools to generate a bam index file (.bai) for visualization with IGV
samtools index ${output_dir}/${base}.bam ${output_dir}/${base}.bai

## Use deeptools to generate a normalized bedgraph file for each strand (for plotting with Sushi in R)
### -bs option specifices the bin size in bases
### normalizes using counts per million (CPM)
#### First, generate bedgraph for the forward strand
#### --samFlagExclude 16 filters out reads that map to the reverse strand
bamCoverage -b ${output_dir}/${base}.bam -p $2 -o ${output_dir}/${base}.norm.forward.bedGraph -of bedgraph -bs 1 --normalizeUsing CPM --samFlagExclude 16
#### Second, generate bedgraph for the reverse strand
#### --samFlagInclude 16 selects only the reads the map to the reverse strand
bamCoverage -b ${output_dir}/${base}.bam -p $2 -o ${output_dir}/${base}.norm.reverse.bedGraph -of bedgraph -bs 1 --normalizeUsing CPM --samFlagInclude 16

