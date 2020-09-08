#!/bin/bash/

# This script takes a fastq file and selects reads based on their 5' barcode.
# USAGE: sh select_5prime_barcode.sh <full pathname of fastq file> <full pathname of the main directory>

# This script assumes the following about the structure of the metadata file:
## 1. The metadata file is named "metadata.txt"
## 2. The output filename must be located in column 1
## 3. There must be a column containing the input filenames (e.g. LIB037835_SRN00139199_S1_L001_R1.fastq)

# Assign variables to the 5' barcodes
barcode1=AGCG #AF-PP-333
barcode2=CGTC #AF-PP-334

# Assign variables to the input and output directories
input_path=$(dirname $1)/
input_file=$(basename $1)
output_path=${2}/results/sort_5prime/
meta=${2}/meta/metadata.txt
log_file=${2}/logs/sort_5prime.log

# Make the output directory if it doesn't already exist
mkdir -p $output_path

# Use metadata to assign a new, meaningful name to each fastq file
## Searches within metadata for the line containing the specified BPF filename
## Then assigns the first field of that line to a variable
new_name=$(grep $input_file $meta | awk '{print $1}')

# Create new fastq file containing only reads that begin with any one of the specified barcodes.
## -B 1 prints 1 line of leading context before the matching line
## -A 2 prints 2 lines of trailing context after the matching line
## -e option specifies each pattern to search for, acts as an OR operator when used more than once.
## ^ option restricts the search to the beginning of each line
## sed command deletes the "--" lines that separate each read
grep -B 1 -A 2 -e ^$barcode1 -e ^$barcode2 $input_path$input_file | sed '/^--/d' > $output_path$new_name

# Report the number of reads that were sorted
## First, count the number of lines in each file and divide by 4 to get the number of reads
total_reads=$(wc -l $input_path$input_file | awk '{print $1/4}')
barcoded_reads=$(wc -l $output_path$new_name | awk '{print $1/4}')
## Then, determine the number of reads that didn't begin with either barcode
unassigned=$(($total_reads - $barcoded_reads))
## Finally, print the above information and save it to a log file
echo -e "$barcoded_reads reads were assigned to $new_name.\n$unassigned reads from $input_file were unassigned.\n" >> $log_file

