#! /bin/bash

# This is the main workflow for sorting, trimming, mapping, and counting small RNA-seq reads

# Assign a variable to the pathname of the project (this is the main folder). Change this to fit your own path.
main_dir=/n/groups/kennedy/anne/smallRNAseq/hangover_official

# Assign variables to the subdirectories
logs=${main_dir}/logs
meta=${main_dir}/meta
raw_data=${main_dir}/raw_data
results=${main_dir}/results
scripts=${main_dir}/scripts

# Make the above directories if they don't already exist
mkdir -p $logs $meta $raw_data $results $scripts

# Before continuing, make sure that:
## 1. The raw fastq files are in the "raw_data" directory
## 2. The scripts called below are in the "scripts" directory
## 3. metadata.txt is in the "meta" directory (see select_5prime_barcode.sh for more instructions on metadata)


# Select reads based on their 5' barcode (AGCG or CGTC) and assign meaningful names to the resulting fastq files
## Store the resulting .fastq files in ${results}/sort_5prime
for fq in ${raw_data}/*.fastq
do
	sbatch -p short -t 0-12:00 --job-name sort-5prime --wrap="sh ${scripts}/select_5prime_barcode.sh $fq $main_dir"
done

# Trim off the 3' adapter using cutadapt
## First, assign a variable to the 3' adapter sequence
adapter_3prime=CTGTAGGCACCATCAATAGATCGGAAGAGCAC
## Then, trim off the adapter (+ anything that follows) and discard reads that didn't contain any adapter sequence
## Store the resulting .fastq files in ${results}/trim3
## Store the log reports in ${logs}/trim3
for fq in ${results}/sort_5prime/*.fastq
do
	sbatch -p short -t 0-12:00 --job-name trim-3prime --wrap="sh ${scripts}/trim_3prime.sh $fq $adapter_3prime $main_dir"
done
## Create a file (summary.txt) that summarizes the main stats of the 3' trimming and store it in ${logs}/trim3
sbatch -p short -t 0-1:00 --wrap="sh ${scripts}/trim_3prime_stats.sh $main_dir"

# Trim off the 5' barcodes (AGCG and CGTC) from the 3'-trimmed files using cutadapt
## Store the resulting .fastq files in ${results}/trim3_trim5
## Store the log reports in ${logs}/trim5
for fq in ${results}/trim3/*.fastq
do
	sbatch -p short -t 0-12:00 --job-name trim-5prime --wrap="sh ${scripts}/trim_5prime.sh $fq $main_dir"
done

# Run quality control on the trimmed reads with fastQC and store in ${results}/fastqc_trimmed
for fq in ${results}/trim3_trim5/*.fastq
do
	sbatch -p short -t 0-2:00 --job-name fastqc --wrap="sh ${scripts}/fastqc_trimmed.sh $fq $main_dir"
done

# Merge the fastq files from different lanes and store them in ${results}/trim3_trim5_merge
## Note, merge_lanes.sh assumes that the lane number indicator immediately precedes the .fastq extension in the sample file name (see metadata.txt example)
sbatch -p short -t 0-4:00 --job-name merge --wrap="sh ${scripts}/merge_lanes.sh $main_dir"

# Map to the genome (WBcel235) using bowtie
## Note, edit map_reads.sh to change the reference genome or to alter the alignment settings
## First, specify the number of cores to use
cores=4
## Then, align reads to the genome
## Use samtools to generate .bam and .bai files
## Store .sam, .bam, and .bai files in ${results}/alignment
## Store the log reports in ${logs}/alignment
for fq in ${results}/trim3_trim5_merge/*.fastq
do
	sbatch -p short -c $cores -t 0-12:00 --job-name align --wrap="sh ${scripts}/map_reads.sh $fq $cores $main_dir"
done
## Create a file (summary.txt) that summarizes the main stats of the mapping and store it in ${logs}/alignment
sbatch -p short -t 0-1:00 --wrap="sh ${scripts}/alignment_stats.sh ${logs}/alignment"

# Count the number of reads that are antisense to each gene using featureCounts
## Note, edit count_antisense_reads.sh to change the gene annotation file or to alter the featureCounts settings
## Store .txt files in ${results}/counts_antisense
## Store the log reports in ${logs}/counts_antisense
## Create a simplified count file for R with the extension .Rmatrix.txt and store in ${results}/counts_antisense
for bam in ${results}/alignment/*.bam
do
	sbatch -p short -c $cores -t 0-12:00 --job-name count --wrap="sh ${scripts}/count_antisense_reads.sh $bam $cores $main_dir"
done
## Create a file (summary.txt) that summarizes the main stats of the counting and store it in ${logs}/counts_antisense
sbatch -p short -t 0-1:00 --wrap="sh ${scripts}/counts_stats.sh ${logs}/counts_antisense"


