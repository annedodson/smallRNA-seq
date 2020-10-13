# Gene-level analysis of small RNA targets (wild type vs. mutant)

# Starting with raw count matrices, this script normalizes counts using DESeq2, tests for differential abundance using DESeq2, and generates plots for quality control and data visualization


####################################################
# LOAD TOOLS
library(DESeq2)
library(tidyverse)
library(ggrepel)


####################################################
# SOURCE FILES
source("functions/data_import.R")
source("functions/data_analysis.R")
source("functions/plots.R")


####################################################
# DATA IMPORT

# Read in the metadata file (.csv) using "import_meta" function from "functions/data_import.R"
## The sample names in column 1 should match the basename of the data file, but without the "_featureCounts.Rmatrix.txt" (e.g. if the data file is named "sampleA_featureCounts.Rmatrix.txt", then the sample name in the metadata should be "sampleA")
metadata <- import_meta("meta/metadata_for_DESeq.csv")

# Read in a key for the gene IDs using "import_gene_names" function from "functions/data_import.R"
gene_names <- import_gene_names("meta/gene_names.txt")

# Read in the count matrices and merge them into one dataframe using "import_counts" function from "functions/data_import.R"
raw_counts <- import_counts(metadata$Sample)


####################################################
# DATA ANALYSIS

# Normalize counts with DESeq2
## First, make sure the sample names in "metadata" match and are in the same order as in "raw_counts"
(all(colnames(raw_counts) == metadata$Sample))
## Next, create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ Strain)
## Calculate and view the normalization factor (median of ratios method)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
## Retrieve and save the normalized counts matrix from "dds"
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Run DESeq analysis
dds <- DESeq(dds)

# Pairwise comparison between wild type (N2) and mutant (JH3225)
## Define contrast
contrast <- c("Strain", "JH3225", "N2")
## Specify the cutoffs
significance <- 0.05
log_fc <- 1
## Extract the results table using a significance 
res_table_unshrunken <- results(dds, contrast = contrast, alpha = significance, lfcThreshold = log_fc)
## Shrink the log2 fold changes
res_table <- lfcShrink(dds, contrast = contrast, res = res_table_unshrunken)
## Summarize results
summary(res_table)
## Convert the results table into a tibble and add a column with meaningful gene names using "new_res_table" function from "functions/data_analysis.R"
res_table_tb <- new_res_table(res_table)
## Create a list of just the significant genes using "sig_genes" function from "functions/data_analysis.R"
sig_tb <- sig_genes(res_table_tb, padj_cutoff = significance, lfc_cutoff = log_fc)
## Save
write.table(sig_tb, file="data/sig_genes.txt", sep="\t", quote=F, col.names=NA)


####################################################
# PLOTS

# Quality control - PCA
## log2-transform the counts and execute PCA
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup = "Genotype", returnData=TRUE)
## Plot
### First, reorder the factor levels of pcaData$Genotype to control sample order in the legend
pcaData$Genotype <- factor(pcaData$Genotype, levels = c("wild_type", "mutant"))
### Plot using "PCA_plot" function from "functions/plots.R" and save
tiff("figures/PCAplot.tiff", width = 5, height = 2.5, units = 'in', res = 300)
PCA_plot(pcaData, category="Genotype")
dev.off()

# Create a volcano plot that labels the top 10 genes using "volcano" function from "functions/plots.R"
tiff("figures/volcano.tiff", width = 6, height = 4, units = 'in', res = 300, compression = "lzw")
volcano(res_table_tb, padj_cutoff = significance, lfc_cutoff = log_fc)
dev.off()

# Plot levels for a single gene
## Specify the wormbase gene ID
gene <- "WBGene00023421"
## Specify the order of the samples on the x-axis
sample_order <- c("wild_type", "mutant")
## Plot using "gene_counts" function from "functions/plots.R" and save
tiff(sprintf("figures/%s.tiff", gene), width = 3, height = 4, units = 'in', res = 300)
gene_counts(dds, gene, sample_order)
dev.off()


####################################################
# Save session info to a file
writeLines(capture.output(sessionInfo()), "data/sessionInfo.txt")
