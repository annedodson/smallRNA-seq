#!/bin/bash/

# This script creates a summary of the cutadapt results for all samples.
# USAGE: sh trim_3prime_stats.sh <full pathname of the main directory>

# Create a new file that will contain the stats
## Assign a variable to the new file
summary=${1}/logs/trim3/summary.txt
## Fill in column 1 with the sample names using metadata.txt
awk '{print $1}' ${1}/meta/metadata.txt | sed 's/.fastq//' > $summary
## Assign a variable to the text in the header
sample_header=$(awk 'FNR == 1 {print $1}' $summary)
## Add 4 new columns (separated by tabs) and fill in the header
sed -i "s/$sample_header/$sample_header\ttotal_reads_processed\treads_with_adapters\treads_that_were_too_short\treads_written_(passing_filters)/" $summary
## Get sample names
samples=$(awk 'FNR > 1 {print $1}' $summary)

## For each sample, find the main stats and incorporate them into the summary file
for sample in $samples
do
	# Assign variable to the sample name's cutadapt log file
	log_file=${1}/logs/trim3/${sample}.txt
	# Find the relevant stats within the log file and assign them to variables
	total_reads=$(grep 'Total reads processed:' $log_file | awk '{print $NF}')
	adapter_reads=$(grep 'Reads with adapters:' $log_file | awk '{print $(NF-1), $NF}')
	too_short=$(grep 'Reads that were too short:' $log_file | awk '{print $(NF-1), $NF}')
	written_reads=$(grep 'Reads written (passing filters):' $log_file | awk '{print $(NF-1), $NF}')
	# Incorporate the stats into the appropriate row in summary file
	sed -i "s/$sample/$sample\t$total_reads\t$adapter_reads\t$too_short\t$written_reads/" $summary
done
