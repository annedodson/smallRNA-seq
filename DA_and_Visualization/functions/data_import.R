# Functions for importing files:
# 1. import_meta
# 2. import_gene_names
# 3. import_counts


####################################################
# 1. import_meta (reads in the .csv metadata file and converts it to a tibble)
import_meta <- function(metadata_file) {
  metadata <- read.csv(file = metadata_file, header=T, row.names=1)
  metadata <- metadata %>%
    rownames_to_column(var="Sample") %>%
    as_tibble()
  return(metadata)
}


####################################################
# 2. import_gene_names (reads in list of gene IDs and their corresponding gene names)
import_gene_names <- function(gene_list) {
  gene_names <- read.delim(gene_list, header=F, colClasses=c("character", "character"))
  ## Rename the columns
  colnames(gene_names)[1] <- "gene"
  colnames(gene_names)[2] <- "gene_name"
  ## Remove spaces from gene_name strings
  gene_names$gene_name <- gsub('\\s+', '', gene_names$gene_name)
  return(gene_names)
}


####################################################
# 3. import_counts (reads in count matrices and merges them into one table called "raw_counts")
import_counts <- function(sample_list) {
  ## Create a for loop to load each file for every sample
  for (sample in sample_list) {
    ### Assign a variable to the filename
    file <- sprintf("data/%s_featureCounts.Rmatrix.txt", sample)
    cat("Reading in", file, "\n")
    ### Read in data file
    temp <- read.table(file, header=T, row.names=1)
    ### If this is the first sample, use its file to initialize the "raw_counts" table
    if (sample == sample_list[1]) {
      raw_counts <- temp 
    }
    ### Otherwise, just append the current counts column to the "raw_counts" table
    ### But first, make sure the gene IDs of "temp" are all ordered the same as the gene IDs of "raw_counts"
    else if (all(rownames(temp) == rownames(raw_counts))) {
      raw_counts <- cbind(raw_counts, temp[ ,1, drop=FALSE]) 
    }
    ### If something went wrong with the first two options, print an error and stop the for loop
    else {
      cat("Could not complete the loading of the data files")
      break
    }
  }
  rm(temp)
  return(raw_counts)
}
