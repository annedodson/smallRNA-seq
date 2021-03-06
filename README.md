# *C. elegans* small RNA-seq analysis

This is a pipeline to analyze next-generation sequencing of small RNAs in *C. elegans*. The pipeline can be broken down into two major parts:

1. **Generate count matrices.** Trims and maps reads to the *C. elegans* genome, then generates count matrices of the number of reads mapping antisense to each gene. This first part of the pipeline is designed to run in a high-performance computing cluster based on Linux and Slurm.
	- <code>generate_count_matrices.sh</code> is the main file for Part 1 and calls all the other Part 1 scripts. On line 6 of <code>generate_count_matrices.sh</code>, specify the full pathname of your project directory:
		```
		# Assign a variable to the pathname of the project (this is the main directory). Change this to fit your own path.
		main_dir=<project pathname>
		```
	- Then execute lines 8-16 of <code>generate_count_matrices.sh</code> to generate the following directory structure:
		```
		project name
		  ├── logs
		  ├── meta
		  ├── raw_data
		  ├── results
		  └── scripts
		```
	- Before continuing on with the rest of <code>generate_count_matrices.sh</code>, make sure that:
		- You've copied your raw, demultiplexed fastq files into the <code>raw_data</code> directory.
		- All Part 1 scripts are in the <code>scripts</code> directory.
		- Your metadata file <code>metadata.txt</code> is in the <code>meta</code> directory. Column 1 of <code>metadata.txt</code> must contain the desired output filename, and there must also be a column containing the input filename. See <code>metadata.txt</code> in this repository for an example.

	- Note, this pipeline assumes the reads contain a 4-nucleotide-long barcode at the 5' end. If your reads do not contain a 5' barcode and instead begin immediately with the insert, make the following two changes:

		- Change line 36 in <code>select_5prime_barcode.sh</code> from:
			```bash
			grep -B 1 -A 2 -e ^$barcode1 -e ^$barcode2 $input_path$input_file | sed '/^--/d' > $output_path$new_name
			```
			to:
			```bash
			cp $input_path$input_file $output_path$new_name
			```
			With this change, running <code>select_5prime_barcode.sh</code> will simply assign new, meaningful names to the fastq files using <code>metadata.txt</code> and place them in a new directory in <code>results</code> called <code>sort_5prime</code>.

		- Change line 23 in <code>trim_5prime.sh</code> from:
			```bash
			cutadapt -u 4 -o $output $1 > ${2}/logs/trim5/${base}.txt
			```
			to:
			```bash
			cp $1 $output
			```
			With this change, running <code>trim_5prime.sh</code> will simply copy the fastq files into a new directory in <code>results</code> called <code>trim3_trim5</code> and add "_trim5" to the end of each filename.

2. **Differential analysis and visualization.** Uses the count matrices generated in Part 1 to perform a simple **wild type vs. mutant** analysis to identify genes that are differentially targeted by small RNAs. This part of the pipeline is designed to run as an RStudio project (<code>DA_and_visualization.Rproj</code>).
	- <code>main_script.R</code> is the main file for Part 2.

	- Before beginning, make sure the count matrices are in the <code>data</code> directory and that the metadata file (.csv format) is in the <code>meta</code> directory. Examples of these files can be found in the <code>example_files</code> directory.

	- Part 2 outputs include the following:
		- A table of normalized counts (median of ratios method)
		- A list of differentially targeted genes and their corresponding log2 fold changes and adjusted p-values
		- A biplot of the top two principal components determined by principal component analysis
		- A volcano plot of log2 fold change vs. significance, with labels for the top 10 significant genes
		- The option to plot normalized counts for any given gene (specified by WormBase Gene ID)

## Software requirements

| Software                    | Version      | Used in                                         |
| --------------------------- | ------------ | ----------------------------------------------- |
| <code>gcc</code>            | 6.2.0        | Part 1: Generate count matrices                 |
| <code>python</code>         | 2.7.12       | Part 1: Generate count matrices                 |
| <code>cutadapt</code>       | 1.14         | Part 1: Generate count matrices                 |
| <code>fastqc</code>         | 0.11.5       | Part 1: Generate count matrices                 |
| <code>bowtie</code>         | 1.2.2        | Part 1: Generate count matrices                 |
| <code>samtools</code>       | 1.9          | Part 1: Generate count matrices                 |
| <code>deeptools</code>      | 3.0.2        | Part 1: Generate count matrices                 |
| <code>featureCounts</code>  | 2.0.0        | Part 1: Generate count matrices                 |
| <code>R</code>              | 3.5.1        | Part 2: Differential analysis and visualization |
| <code>DESeq2</code>         | 1.22.2       | Part 2: Differential analysis and visualization |
| <code>tidyverse</code>      | 1.2.1        | Part 2: Differential analysis and visualization |
| <code>ggrepel</code>        | 0.8.1        | Part 2: Differential analysis and visualization |